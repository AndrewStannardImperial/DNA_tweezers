import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import expit
from dna_calculator import DNACalculator, DNACalculatorTH as DNA_th
from joblib import Parallel, delayed
from tqdm.auto import tqdm

BASE_DIR = Path(__file__).resolve().parent

#--------------------------------
#------- Plotting Configs -------
#--------------------------------

mpl.rcParams.update({
    'text.usetex':       False,
    'mathtext.fontset':  'dejavusans',        # CM for all math
    'font.family':       'helvetica',     # DejaVu Serif or your system serif
    'font.size':         15,
    'axes.titlesize':    22,
    'axes.labelsize':    15,
    'xtick.labelsize':   16,
    'ytick.labelsize':   16,
    'axes.linewidth':    1.2,
    'xtick.major.width': 1.5,
    'ytick.major.width': 1.5,
    'xtick.major.size': 6, 
    'xtick.minor.size': 3, 
    'ytick.major.size': 6, 
})

#-------------------------
#------- Constants -------
#-------------------------

qe = 1.609e-19                          # fundamental charge
kb = 1.38e-23                           # boltzmann constant
T=298                                   # temperature
eps0 = 8.85e-12                         # permittivity of free space
lb = 1e10*qe**2 / (4*np.pi*eps0 * kb*T) # Bjerrum length in A

Lam  = 5            # nonlocal correlation length
epsb = 80           # dielectric constant of water
epss = 5            # short range dielectric constant
epsc = 2            # dielectric constant of helix

g_b = epsb/epss     # ratio of bulk/short-range epsilon 
g_c = epsc/epss     # ratio of core/short-range epsilon

cb = 0.100          # electrolyte concentration n moldm-3
Na = 6.023*1e23     # avogadro constant in mol-1
nb = cb*Na/1e27     # unit conversions
kappa = np.sqrt(2*nb*qe**2/(1e-10*eps0*epsb*kb*T)) # debye length in A-1

phis = 0.8 * np.pi  # width of minor groove
H    = 34           # helical pitch
hr   = 3.4          # rise per base-pair

f3   = 0.0          # strands

a = 10              # radius of first helix
b = 6               # radius of second helix


N_phi = 500
N_psi = 25

phi_grid = np.linspace(-np.pi, np.pi, N_phi)

Nmax, Mmax = 5, 5
n_idx = np.arange(-Nmax, Nmax+1)
m_idx = np.arange(-Mmax, Mmax+1)

cos_m_grid = np.cos(np.outer(phi_grid, m_idx))  # (N_phi, Mmax)
sin_m_grid = np.sin(np.outer(phi_grid, m_idx))  # (N_phi, Mmax)


l_c = 100

L_bp = 32
L_DNA = L_bp*hr

#------------------------------
#------- Fitting R_F(c) -------
#------------------------------

conc_exp = np.array([0, 1, 2, 5, 10, 20, 50, 100, 200])

rel_FR_Mg = np.array([1, 0.995, 0.994, 0.991, 0.986, 0.983, 0.976, 0.969, 0.963])
rel_FR_Ca = np.array([1, 0.993, 0.990, 0.984, 0.979, 0.972, 0.966, 0.956, 0.950])
rel_FR_Mg_uncert = np.array([0, 0.002, 0.002, 0.003, 0.002, 0.002, 0.002, 0.002, 0.002])
rel_FR_Ca_uncert = np.array([0, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.003])


FR_Mg = np.array([5.71, 5.68, 5.67, 5.66, 5.63, 5.61, 5.57, 5.53, 5.50])
FR_Ca = np.array([5.71, 5.67, 5.65, 5.62, 5.59, 5.55, 5.52, 5.46, 5.42])
FR_Mg_uncert = np.array([0.4, 0.4, 0.4, 0.4, 0.39, 0.39, 0.39, 0.39, 0.39])
FR_Ca_uncert = np.array([0.4, 0.4, 0.4, 0.39, 0.39, 0.39, 0.39, 0.38, 0.38])

def R_c(c, a, k):
    return 1-a*np.log(1+c/k)

params_Mg = curve_fit(R_c, conc_exp, rel_FR_Mg)[0]
params_Ca = curve_fit(R_c, conc_exp, rel_FR_Ca)[0]

L_nm = 12.2

def psi_s_from_RF(RF_nm):
    x = RF_nm / (2.0 * L_nm)
    # clip for numerical safety
    x = np.clip(x, 0.0, 0.999999999)
    return 2.0 * np.arcsin(x)

def psi_s_mean_Mg(c_mM):
    return psi_s_from_RF(FR_Mg[0]*R_c(c_mM, params_Mg[0], params_Mg[1]))

def psi_s_mean_Ca(c_mM):
    return psi_s_from_RF(FR_Ca[0]*R_c(c_mM, params_Ca[0], params_Ca[1]))


#-------------------------------------------
#------- Setting up curve calculator -------
#-------------------------------------------


def compute_curves_nominal(R, f2, K, ion,
                           concentrations,
                           a, b, L_DNA, l_c,
                           N_psi, N_phi,
                           n_idx, m_idx,
                           phi_grid, cos_m_grid, sin_m_grid,
                           Lam, g_b, kappa, epsb, g_c, epsc, lb, phis, H, hr, f3,
                           DNACalculator, DNA_th):
    """
    Returns (deltaG_hom, deltaG_nonhom, deltaG_recog) arrays of shape (len(concentrations),)
    for a single (R,f2,K) using nominal psi_s_mean(c).
    """
    f1 = 1.0 - f2

    def create_calculator(TH):
        calc = DNACalculator(Lam, g_b, kappa, epsb, g_c, epsc, lb, phis, H, hr, f1, f2, f3, TH)
        orig_p = calc.p
        calc.p = lambda k, TH=None: orig_p(int(abs(k)), TH)
        return calc

    omega_vals = np.empty_like(concentrations, dtype=float)

    for i, c in enumerate(concentrations):
        TH = K * c / (1.0 + K * c)
        if ion == "Mg":
            psi_s = psi_s_mean_Mg(c)
        elif ion == "Ca":
            psi_s = psi_s_mean_Ca(c)

        psi_grid = np.linspace(psi_s, np.pi, N_psi)

        calc = create_calculator(TH)
        integrand = np.empty(N_psi, dtype=float)

        for j, psi in enumerate(psi_grid):
            w_mat = np.array([[calc.w_nm(R, psi, a, b, n, m) for m in m_idx] for n in n_idx])

            C2_all = cos_m_grid.dot(w_mat.T)
            S2_all = sin_m_grid.dot(w_mat.T)

            Umin = np.inf
            for phi1 in phi_grid:
                u1 = np.cos(n_idx * phi1)
                v1 = np.sin(n_idx * phi1)
                E_vals = C2_all.dot(u1) + S2_all.dot(v1)
                Umin = min(Umin, E_vals.min())

            integrand[j] = np.exp(-Umin / 2.0) * np.sin(psi) / (1.0 - np.cos(psi_s))

        omega_vals[i] = np.trapezoid(integrand, psi_grid)

    calc_th = DNA_th(Lam, g_b, kappa, epsb, g_c, epsc, lb, phis, H, hr, f1, f2, f3)
    TH_list = K * concentrations / (1.0 + K * concentrations)

    u_hom = calc_th.u_hom_min(R, a, b, L_DNA, TH_list, method='nonlocal')
    u_nonhom = calc_th.u_nonhom_min(R, a, b, L_DNA, l_c, TH_list, method='nonlocal')

    deltaG_hom = np.array(u_hom) + np.log(omega_vals)
    deltaG_nonhom = np.array(u_nonhom) + np.log(omega_vals)
    deltaG_recog = deltaG_hom - deltaG_nonhom

    return deltaG_hom, deltaG_nonhom, deltaG_recog


def sample_thetas(theta_0, Sigma_0, n_samples, seed=0,
                           max_draws=1_000_000):
    """
    Sample theta ~ N(theta_0, Sigma_0) and reject non-physical draws.

    Physical domain enforced:
      - R > 0
      - 0 < f2 < 1
      - K > 0
    where theta = [R, f2, K] (adjust order if yours differs).

    Returns
    -------
    thetas : (n_samples, 3) array
    """
    rng = np.random.default_rng(seed)
    out = np.empty((n_samples, len(theta_0)), dtype=float)

    kept = 0
    draws = 0

    while kept < n_samples:
        if draws >= max_draws:
            raise RuntimeError(
                f"Reached max_draws={max_draws} with kept={kept}. "
                "Acceptance rate likely too low; consider log/logit sampling."
            )

        th = rng.multivariate_normal(mean=theta_0, cov=Sigma_0)
        draws += 1

        R, f2, K = th  # <-- adjust if your ordering differs
        if (R > 0.0) and (K > 0.0) and (0 < f2 < 1.0):
            out[kept] = th
            kept += 1

    return out

# def sample_thetas(theta0, Sigma_theta, n_samples, seed=0, eps=1e-20):
#     """
#     Approximate transform sampling:
#       z = [log R, logit f2, log K]
#     using Sigma_z ≈ J Sigma_theta J^T evaluated at theta0.
#     """
#     R0, f20, K0 = map(float, theta0)

#     # guard against exactly 0/1 in logit
#     f20 = float(np.clip(f20, eps, 1.0 - eps))

#     z0 = np.array([
#         np.log(R0),
#         np.log(f20 / (1.0 - f20)),
#         np.log(K0),
#     ])

#     # Jacobian of z wrt theta at theta0
#     J = np.diag([
#         1.0 / R0,
#         1.0 / (f20 * (1.0 - f20)),
#         1.0 / K0,
#     ])

#     Sigma_z = J @ Sigma_theta @ J.T

#     rng = np.random.default_rng(seed)
#     z = rng.multivariate_normal(mean=z0, cov=Sigma_z, size=n_samples)

#     R = np.exp(z[:, 0])
#     f2 = expit(z[:, 1])
#     f2 = np.clip(f2, eps, 1.0 - eps)
#     K = np.exp(z[:, 2])

#     return np.column_stack([R, f2, K])

# def envelope(Y, interval=0.90):
#     """
#     Y: (n_samples, n_points)
#     interval: e.g. 0.90, 0.95, 0.68
#     Returns (lo, median, hi)
#     """
#     alpha = 1.0 - interval
#     lo = np.quantile(Y, alpha/2, axis=0)
#     md = np.quantile(Y, 0.50, axis=0)
#     hi = np.quantile(Y, 1.0 - alpha/2, axis=0)
#     return lo, md, hi

def _one_sample(idx, theta, ion, concentrations):
    R_s, f2_s, K_s = theta

    y_h, y_nh, y_r = compute_curves_nominal(
        R_s, f2_s, K_s, ion,
        concentrations,
        a, b, L_DNA, l_c,
        N_psi, N_phi,
        n_idx, m_idx,
        phi_grid, cos_m_grid, sin_m_grid,
        Lam, g_b, kappa, epsb, g_c, epsc, lb, phis, H, hr, f3,
        DNACalculator, DNA_th
    )
    return idx, y_h, y_nh, y_r

def run_mc_and_save(
    ion: str,
    thetas: np.ndarray,
    concentrations: np.ndarray,
    out_dir: str,
    n_jobs: int,
    backend: str = "loky",
    batch_size: int = 1,
    flush_every: int = 25,
):
    """
    Runs MC in parallel and saves Yh/Ynh/Yr incrementally as memmaps.

    Produces:
      out_dir/{ion}_Yh.dat, {ion}_Ynh.dat, {ion}_Yr.dat
      out_dir/{ion}_thetas.npy
      out_dir/{ion}_meta.npz
    """
    os.makedirs(out_dir, exist_ok=True)

    N = thetas.shape[0]
    M = len(concentrations)

    # Save sampled thetas for reproducibility
    np.save(os.path.join(out_dir, f"{ion}_thetas.npy"), thetas)

    # Memmaps for incremental writing
    Yh_path  = os.path.join(out_dir, f"{ion}_Yh.dat")
    Ynh_path = os.path.join(out_dir, f"{ion}_Ynh.dat")
    Yr_path  = os.path.join(out_dir, f"{ion}_Yr.dat")

    Yh  = np.memmap(Yh_path,  dtype="float64", mode="w+", shape=(N, M))
    Ynh = np.memmap(Ynh_path, dtype="float64", mode="w+", shape=(N, M))
    Yr  = np.memmap(Yr_path,  dtype="float64", mode="w+", shape=(N, M))

    # Mark incomplete rows as NaN initially (so you can resume/diagnose)
    Yh[:]  = np.nan
    Ynh[:] = np.nan
    Yr[:]  = np.nan
    Yh.flush(); Ynh.flush(); Yr.flush()

    # Run jobs
    iterator = (delayed(_one_sample)(i, thetas[i], ion, concentrations) for i in range(N))

    completed = 0
    for idx, y_h, y_nh, y_r in tqdm(
        Parallel(n_jobs=n_jobs, backend=backend, batch_size=batch_size, return_as="generator")(iterator),
        total=N,
        desc=f"MC {ion}"
    ):
        Yh[idx, :]  = y_h
        Ynh[idx, :] = y_nh
        Yr[idx, :]  = y_r
        completed += 1

        if (completed % flush_every) == 0:
            Yh.flush(); Ynh.flush(); Yr.flush()

    # final flush
    Yh.flush(); Ynh.flush(); Yr.flush()

    # Write metadata (handy later)
    np.savez(
        os.path.join(out_dir, f"{ion}_meta.npz"),
        ion=ion,
        N=N,
        M=M,
        concentrations=concentrations,
        Yh_path=Yh_path,
        Ynh_path=Ynh_path,
        Yr_path=Yr_path,
    )

    # Return memmaps (still backed by disk)
    return Yh, Ynh, Yr
#----------------------------------------
#------- Sampling Mg2+ data first -------
#----------------------------------------

N_samples = 1500
seed = 0

def load_fit_results(path=BASE_DIR / "fit_results.npz"):
    if not path.exists():
        raise FileNotFoundError(
            f"Missing fitted parameters: {path}\n"
            "Run `python fit_data.py` before running `python MC_sampling.py`."
        )

    fit = np.load(path)
    required = ["phat_Mg", "cov_Mg", "phat_Ca", "cov_Ca"]
    missing = [key for key in required if key not in fit.files]
    if missing:
        raise KeyError(
            f"{path} is missing required arrays: {missing}. "
            "Regenerate it with `python fit_data.py`."
        )

    return fit["phat_Mg"], fit["cov_Mg"], fit["phat_Ca"], fit["cov_Ca"]


theta_0_Mg, Sigma_0_Mg, theta_0_Ca, Sigma_0_Ca = load_fit_results()
print("Loaded fitted parameters from:", BASE_DIR / "fit_results.npz")
print("theta_0_Mg:", theta_0_Mg)
print("theta_0_Ca:", theta_0_Ca)


thetas_Mg = sample_thetas(theta_0_Mg, Sigma_0_Mg, N_samples, seed=seed)
thetas_Ca = sample_thetas(theta_0_Ca, Sigma_0_Ca, N_samples, seed=1)

conc_low = np.logspace(0, 1, 51)
conc_high = np.linspace(10, 250, 50)
concentrations = np.concatenate((conc_low[:-1], conc_high))

def get_n_jobs():
    default = min(60, max(1, (os.cpu_count() or 2) - 1))
    value = os.environ.get("MC_N_JOBS", default)
    try:
        n_jobs = int(value)
    except ValueError as exc:
        raise ValueError("MC_N_JOBS must be an integer.") from exc
    if n_jobs < 1:
        raise ValueError("MC_N_JOBS must be at least 1.")
    return n_jobs


n_jobs = get_n_jobs()
print(f"Running Monte Carlo with n_jobs={n_jobs} parallel worker(s)")
out_dir = BASE_DIR / "mc_outputs"

# --- Mg ---

Yh_Mg, Ynh_Mg, Yr_Mg = run_mc_and_save(
    ion="Mg",
    thetas=thetas_Mg,
    concentrations=concentrations,
    out_dir=out_dir,
    n_jobs=n_jobs,
    batch_size=8,
    flush_every=100,
)

# --- Ca ---

Yh_Ca, Ynh_Ca, Yr_Ca = run_mc_and_save(
    ion="Ca",
    thetas=thetas_Ca,
    concentrations=concentrations,
    out_dir=out_dir,
    n_jobs=n_jobs,
    batch_size=8,
    flush_every=100,
)

def band_about_best(Y, y_best, interval=0.68):
    """
    Y: (n_samples, n_points) sampled curves
    y_best: (n_points,) best-fit curve
    interval: e.g. 0.68, 0.90, 0.95

    Returns (lo, center, hi) where center == y_best and
    hi/lo are symmetric around y_best using |deviation| quantiles.
    """
    dev = Y - y_best[None, :]            # (S, M)
    rad = np.quantile(np.abs(dev), interval, axis=0)  # (M,)
    lo = y_best - rad
    hi = y_best + rad
    return lo, y_best, hi

def compute_best_fit_curves(theta_0_Mg, theta_0_Ca, concentrations):
    """
    Compute best-fit curves using the nominal parameter vectors theta_0_*.
    Returns dict with keys:
      best["Mg"]["h"], best["Mg"]["nh"], best["Mg"]["r"]
      best["Ca"]["h"], best["Ca"]["nh"], best["Ca"]["r"]
    """
    best_h_Mg, best_nh_Mg, best_r_Mg = compute_curves_nominal(
        theta_0_Mg[0], theta_0_Mg[1], theta_0_Mg[2], "Mg",
        concentrations,
        a, b, L_DNA, l_c,
        N_psi, N_phi,
        n_idx, m_idx,
        phi_grid, cos_m_grid, sin_m_grid,
        Lam, g_b, kappa, epsb, g_c, epsc, lb, phis, H, hr, f3,
        DNACalculator, DNA_th,
    )
    best_h_Ca, best_nh_Ca, best_r_Ca = compute_curves_nominal(
        theta_0_Ca[0], theta_0_Ca[1], theta_0_Ca[2], "Ca",
        concentrations,
        a, b, L_DNA, l_c,
        N_psi, N_phi,
        n_idx, m_idx,
        phi_grid, cos_m_grid, sin_m_grid,
        Lam, g_b, kappa, epsb, g_c, epsc, lb, phis, H, hr, f3,
        DNACalculator, DNA_th,
    )
    return {
        "Mg": {"h": best_h_Mg, "nh": best_nh_Mg, "r": best_r_Mg},
        "Ca": {"h": best_h_Ca, "nh": best_nh_Ca, "r": best_r_Ca},
    }

best = compute_best_fit_curves(theta_0_Mg, theta_0_Ca, concentrations)

#----------------------------------------
#------- Saving and plotting data -------
#----------------------------------------

def interval_tag(interval: float) -> str:
    # 0.68 -> "68", 0.95 -> "95"
    return f"{int(round(interval * 100)):02d}"

# def save_bands_for_interval(out_dir, concentrations, interval,
#                             h_lo_Mg, h_md_Mg, h_hi_Mg,
#                             nh_lo_Mg, nh_md_Mg, nh_hi_Mg,
#                             r_lo_Mg, r_md_Mg, r_hi_Mg,
#                             h_lo_Ca, h_md_Ca, h_hi_Ca,
#                             nh_lo_Ca, nh_md_Ca, nh_hi_Ca,
#                             r_lo_Ca, r_md_Ca, r_hi_Ca):
#     os.makedirs(out_dir, exist_ok=True)
#     tag = interval_tag(interval)

#     # --- Mg ---
#     df_Mg = pd.DataFrame({
#         "c [mM]": concentrations,
#         f"MC {tag}% Homo lo": h_lo_Mg,
#         f"MC {tag}% Homo med": h_md_Mg,
#         f"MC {tag}% Homo hi": h_hi_Mg,
#         f"MC {tag}% Hetero lo": nh_lo_Mg,
#         f"MC {tag}% Hetero med": nh_md_Mg,
#         f"MC {tag}% Hetero hi": nh_hi_Mg,
#     })
#     df_Mg.to_csv(os.path.join(out_dir, f"Mg2+_theory_MC_{tag}.csv"), index=False)

#     # --- Ca ---
#     df_Ca = pd.DataFrame({
#         "c [mM]": concentrations,
#         f"MC {tag}% Homo lo": h_lo_Ca,
#         f"MC {tag}% Homo med": h_md_Ca,
#         f"MC {tag}% Homo hi": h_hi_Ca,
#         f"MC {tag}% Hetero lo": nh_lo_Ca,
#         f"MC {tag}% Hetero med": nh_md_Ca,
#         f"MC {tag}% Hetero hi": nh_hi_Ca,
#     })
#     df_Ca.to_csv(os.path.join(out_dir, f"Ca2+_theory_MC_{tag}.csv"), index=False)

#     # --- Recognition (both ions in one file) ---
#     df_recog = pd.DataFrame({
#         "c [mM]": concentrations,
#         f"MC {tag}% recog lo (Mg)": r_lo_Mg,
#         f"MC {tag}% recog med (Mg)": r_md_Mg,
#         f"MC {tag}% recog hi (Mg)": r_hi_Mg,
#         f"MC {tag}% recog lo (Ca)": r_lo_Ca,
#         f"MC {tag}% recog med (Ca)": r_md_Ca,
#         f"MC {tag}% recog hi (Ca)": r_hi_Ca,
#     })
#     df_recog.to_csv(os.path.join(out_dir, f"Recog_Energy_MC_{tag}.csv"), index=False)

# def save_bands_for_interval(out_dir, concentrations, interval,
#                             # Mg bands
#                             h_lo_Mg, h_md_Mg, h_hi_Mg,
#                             nh_lo_Mg, nh_md_Mg, nh_hi_Mg,
#                             r_lo_Mg, r_md_Mg, r_hi_Mg,
#                             # Ca bands
#                             h_lo_Ca, h_md_Ca, h_hi_Ca,
#                             nh_lo_Ca, nh_md_Ca, nh_hi_Ca,
#                             r_lo_Ca, r_md_Ca, r_hi_Ca,
#                             # NEW: best-fit curves
#                             best_h_Mg, best_nh_Mg, best_r_Mg,
#                             best_h_Ca, best_nh_Ca, best_r_Ca):
#     """
#     Writes 3 CSVs (Mg, Ca, Recog) and includes BOTH:
#       - MC lo/med/hi bands
#       - best-fit curves (single line)
#     """
#     os.makedirs(out_dir, exist_ok=True)
#     tag = interval_tag(interval)

#     # --- Mg ---
#     df_Mg = pd.DataFrame({
#         "c [mM]": concentrations,

#         f"MC {tag}% Homo lo":  h_lo_Mg,
#         f"MC {tag}% Homo med": h_md_Mg,
#         f"MC {tag}% Homo hi":  h_hi_Mg,
#         "BestFit Homo":        best_h_Mg,

#         f"MC {tag}% Hetero lo":  nh_lo_Mg,
#         f"MC {tag}% Hetero med": nh_md_Mg,
#         f"MC {tag}% Hetero hi":  nh_hi_Mg,
#         "BestFit Hetero":        best_nh_Mg,
#     })
#     df_Mg.to_csv(os.path.join(out_dir, f"Mg2+_theory_MC_{tag}.csv"), index=False)

#     # --- Ca ---
#     df_Ca = pd.DataFrame({
#         "c [mM]": concentrations,

#         f"MC {tag}% Homo lo":  h_lo_Ca,
#         f"MC {tag}% Homo med": h_md_Ca,
#         f"MC {tag}% Homo hi":  h_hi_Ca,
#         "BestFit Homo":        best_h_Ca,

#         f"MC {tag}% Hetero lo":  nh_lo_Ca,
#         f"MC {tag}% Hetero med": nh_md_Ca,
#         f"MC {tag}% Hetero hi":  nh_hi_Ca,
#         "BestFit Hetero":        best_nh_Ca,
#     })
#     df_Ca.to_csv(os.path.join(out_dir, f"Ca2+_theory_MC_{tag}.csv"), index=False)

#     # --- Recognition (both ions in one file) ---
#     df_recog = pd.DataFrame({
#         "c [mM]": concentrations,

#         f"MC {tag}% recog lo (Mg)":  r_lo_Mg,
#         f"MC {tag}% recog med (Mg)": r_md_Mg,
#         f"MC {tag}% recog hi (Mg)":  r_hi_Mg,
#         "BestFit recog (Mg)":        best_r_Mg,

#         f"MC {tag}% recog lo (Ca)":  r_lo_Ca,
#         f"MC {tag}% recog med (Ca)": r_md_Ca,
#         f"MC {tag}% recog hi (Ca)":  r_hi_Ca,
#         "BestFit recog (Ca)":        best_r_Ca,
#     })
#     df_recog.to_csv(os.path.join(out_dir, f"Recog_Energy_MC_{tag}.csv"), index=False)

def save_bands_for_interval(out_dir, concentrations, interval,
                            h_lo_Mg, h_ctr_Mg, h_hi_Mg,
                            nh_lo_Mg, nh_ctr_Mg, nh_hi_Mg,
                            r_lo_Mg, r_ctr_Mg, r_hi_Mg,
                            h_lo_Ca, h_ctr_Ca, h_hi_Ca,
                            nh_lo_Ca, nh_ctr_Ca, nh_hi_Ca,
                            r_lo_Ca, r_ctr_Ca, r_hi_Ca,
                            best_h_Mg, best_nh_Mg, best_r_Mg,
                            best_h_Ca, best_nh_Ca, best_r_Ca):
    os.makedirs(out_dir, exist_ok=True)
    tag = f"{int(round(interval*100)):02d}"

    df_Mg = pd.DataFrame({
        "c [mM]": concentrations,
        "best Homo": best_h_Mg,
        "best Hetero": best_nh_Mg,
        f"MC {tag}% Homo lo": h_lo_Mg,
        f"MC {tag}% Homo hi": h_hi_Mg,
        f"MC {tag}% Hetero lo": nh_lo_Mg,
        f"MC {tag}% Hetero hi": nh_hi_Mg,
    })
    df_Mg.to_csv(os.path.join(out_dir, f"Mg2+_theory_MC_{tag}.csv"), index=False)

    df_Ca = pd.DataFrame({
        "c [mM]": concentrations,
        "best Homo": best_h_Ca,
        "best Hetero": best_nh_Ca,
        f"MC {tag}% Homo lo": h_lo_Ca,
        f"MC {tag}% Homo hi": h_hi_Ca,
        f"MC {tag}% Hetero lo": nh_lo_Ca,
        f"MC {tag}% Hetero hi": nh_hi_Ca,
    })
    df_Ca.to_csv(os.path.join(out_dir, f"Ca2+_theory_MC_{tag}.csv"), index=False)

    df_R = pd.DataFrame({
        "c [mM]": concentrations,
        "best recog (Mg)": best_r_Mg,
        "best recog (Ca)": best_r_Ca,
        f"MC {tag}% recog lo (Mg)": r_lo_Mg,
        f"MC {tag}% recog hi (Mg)": r_hi_Mg,
        f"MC {tag}% recog lo (Ca)": r_lo_Ca,
        f"MC {tag}% recog hi (Ca)": r_hi_Ca,
    })
    df_R.to_csv(os.path.join(out_dir, f"Recog_Energy_MC_{tag}.csv"), index=False)


intervals = [0.68, 0.75, 0.9, 0.95, 0.98]
out_dir = BASE_DIR / "mc_band_csvs"

exp_c = np.array([20, 50, 100, 200])

Mg_homo_E = np.array([2.3918, 1.1692, 0.4255, -0.0440])
Mg_homo_E_e = np.array([0.0873, 0.0354, 0.0289, 0.0456])
Mg_hetero_E = np.array([3.0343, 1.7963, 1.0130, 0.5215])
Mg_hetero_E_e = np.array([0.3311, 0.2417, 0.1602, 0.1467])

Ca_homo_E = np.array([1.68223306987079, 0.405881792142694, -0.248470487508777, -0.643747074316099])
Ca_homo_E_e = np.array([0.0779041455291794, 0.0995916367860504, 0.100332779864682, 0.100093609669645])
Ca_hetero_E = np.array([2.29305234124575, 1.04285238181179, 0.317028891171946, -0.0908624617070255])
Ca_hetero_E_e = np.array([0.220911931657532, 0.0742133745579114, 0.105404745198588, 0.163136249912870])


best_h_Mg, best_nh_Mg, best_r_Mg = best["Mg"]["h"], best["Mg"]["nh"], best["Mg"]["r"]
best_h_Ca, best_nh_Ca, best_r_Ca = best["Ca"]["h"], best["Ca"]["nh"], best["Ca"]["r"]

for interval in intervals:
    # h_lo_Mg, h_md_Mg, h_hi_Mg   = envelope(Yh_Mg,  interval=interval)
    # nh_lo_Mg, nh_md_Mg, nh_hi_Mg = envelope(Ynh_Mg, interval=interval)
    # r_lo_Mg, r_md_Mg, r_hi_Mg   = envelope(Yr_Mg,  interval=interval)

    # h_lo_Ca, h_md_Ca, h_hi_Ca   = envelope(Yh_Ca,  interval=interval)
    # nh_lo_Ca, nh_md_Ca, nh_hi_Ca = envelope(Ynh_Ca, interval=interval)
    # r_lo_Ca, r_md_Ca, r_hi_Ca   = envelope(Yr_Ca,  interval=interval)

    h_lo_Mg, h_ctr_Mg, h_hi_Mg     = band_about_best(Yh_Mg,  best_h_Mg,  interval)
    nh_lo_Mg, nh_ctr_Mg, nh_hi_Mg  = band_about_best(Ynh_Mg, best_nh_Mg, interval)
    r_lo_Mg, r_ctr_Mg, r_hi_Mg     = band_about_best(Yr_Mg,  best_r_Mg,  interval)

    h_lo_Ca, h_ctr_Ca, h_hi_Ca     = band_about_best(Yh_Ca,  best_h_Ca,  interval)
    nh_lo_Ca, nh_ctr_Ca, nh_hi_Ca  = band_about_best(Ynh_Ca, best_nh_Ca, interval)
    r_lo_Ca, r_ctr_Ca, r_hi_Ca     = band_about_best(Yr_Ca,  best_r_Ca,  interval)

    save_bands_for_interval(
        out_dir, concentrations, interval,
        h_lo_Mg, h_ctr_Mg, h_hi_Mg,
        nh_lo_Mg, nh_ctr_Mg, nh_hi_Mg,
        r_lo_Mg, r_ctr_Mg, r_hi_Mg,
        h_lo_Ca, h_ctr_Ca, h_hi_Ca,
        nh_lo_Ca, nh_ctr_Ca, nh_hi_Ca,
        r_lo_Ca, r_ctr_Ca, r_hi_Ca,
        best_h_Mg, best_nh_Mg, best_r_Mg,
        best_h_Ca, best_nh_Ca, best_r_Ca
        #best["Mg"]["h"], best["Mg"]["nh"], best["Mg"]["r"],
        #best["Ca"]["h"], best["Ca"]["nh"], best["Ca"]["r"],
    )
    
    fig = plt.figure(figsize=(8,4))
    gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1.2, 1.2], wspace=0.05)

    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    ax2 = fig.add_subplot(gs[0, 3])
    ax = [ax0, ax1, ax2]

    # --- Mg best-fit (dashed), CI shading from lo/hi ---
    ax[0].plot(concentrations, best["Mg"]["h"],  ls='--', lw=1.5, label=r'$\mathrm{Mg}^{2+}$ homo.',   color='#07baec')
    ax[0].fill_between(concentrations, h_lo_Mg,  h_hi_Mg,  alpha=0.3, color='#07baec')

    ax[0].plot(concentrations, best["Mg"]["nh"], ls='--', lw=1.5, label=r'$\mathrm{Mg}^{2+}$ hetero.', color='#05a876')
    ax[0].fill_between(concentrations, nh_lo_Mg, nh_hi_Mg, alpha=0.3, color='#05a876')

    ax[0].errorbar(exp_c, Mg_homo_E,   yerr=Mg_homo_E_e,   marker='^', markersize=8, c='#07baec', ls='', capsize=5,
                   elinewidth=1.5, capthick=2.5, markeredgewidth=1, markeredgecolor='black')
    ax[0].errorbar(exp_c, Mg_hetero_E, yerr=Mg_hetero_E_e, marker='^', markersize=8, c='#05a876', ls='', capsize=5,
                   elinewidth=1.5, capthick=2.5, markeredgewidth=1, markeredgecolor='black')

    # --- Ca best-fit (dashed), CI shading from lo/hi ---
    ax[1].plot(concentrations, best["Ca"]["h"],  ls='--', lw=1.5, label=r'$\mathrm{Ca}^{2+}$ homo.',   color='#f59621')
    ax[1].fill_between(concentrations, h_lo_Ca,  h_hi_Ca,  alpha=0.3, color='#f59621')

    ax[1].plot(concentrations, best["Ca"]["nh"], ls='--', lw=1.5, label=r'$\mathrm{Ca}^{2+}$ hetero.', color='#eadf3d')
    ax[1].fill_between(concentrations, nh_lo_Ca, nh_hi_Ca, alpha=0.3, color='#eadf3d')

    ax[1].errorbar(exp_c, Ca_homo_E,   yerr=Ca_homo_E_e,   marker='^', markersize=8, c='#f59621', ls='', capsize=5,
                   elinewidth=1.5, capthick=2.5, markeredgewidth=1, markeredgecolor='black')
    ax[1].errorbar(exp_c, Ca_hetero_E, yerr=Ca_hetero_E_e, marker='^', markersize=8, c='#eadf3d', ls='', capsize=5,
                   elinewidth=1.5, capthick=2.5, markeredgewidth=1, markeredgecolor='black')

    # --- Recognition best-fit + CI shading ---
    ax[2].plot(concentrations, best["Mg"]["r"], ls='--', lw=1.5, color='#07baec')
    ax[2].fill_between(concentrations, r_lo_Mg, r_hi_Mg, alpha=0.3, color='#07baec')

    ax[2].plot(concentrations, best["Ca"]["r"], ls='--', lw=1.5, color='#f59621')
    ax[2].fill_between(concentrations, r_lo_Ca, r_hi_Ca, alpha=0.3, color='#f59621')

    ax[2].errorbar(exp_c, Mg_homo_E - Mg_hetero_E,
                   yerr=np.sqrt(Mg_homo_E_e**2 + Mg_hetero_E_e**2),
                   marker='^', markersize=8, c='#07baec', ls='', capsize=5,
                   elinewidth=1.5, capthick=2.5, markeredgewidth=1, markeredgecolor='black')
    ax[2].errorbar(exp_c, Ca_homo_E - Ca_hetero_E,
                   yerr=np.sqrt(Ca_homo_E_e**2 + Ca_hetero_E_e**2),
                   marker='^', markersize=8, c='#f59621', ls='', capsize=5,
                   elinewidth=1.5, capthick=2.5, markeredgewidth=1, markeredgecolor='black')

    # --- formatting (unchanged from your script) ---
    ax[0].set_ylim(9.5, -1)
    ax[1].set_ylim(9.5, -1)
    ax[2].set_ylim(0.5, -1.5)
    ax[0].spines[['top', 'right']].set_visible(False)
    ax[1].spines[['top', 'right', 'left']].set_visible(False)
    ax[2].spines[['top', 'right']].set_visible(False)
    ax[1].tick_params(axis='y', which='both', left=False, labelleft=False)

    ax[0].set_xlabel(r"[Mg$^{2+}$] (mM)")
    ax[1].set_xlabel(r"[Ca$^{2+}$] (mM)")
    ax[2].set_xlabel(r"[M$^{2+}$] (mM)")
    ax[0].set_ylabel(r"Coalignment free energy, $\Delta G~(k_BT)$")
    ax[2].set_ylabel(r"Recognition energy, $\Delta\Delta G~(k_BT)$")

    for axt in ax:
        axt.set_xscale('log')
        axt.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
        axt.set_xticks([1, 10, 100])
        axt.tick_params(which='both', direction='in', top=False, right=False)
        axt.set_xlim(1, 300)

    fig.legend(fontsize='x-small', fancybox=False, edgecolor='white', bbox_to_anchor=(0.58, 0.45))

    tag = int(round(interval * 100))
    plt.savefig(BASE_DIR / f'uncert_band_{tag}CI.pdf', format='pdf', dpi=300, bbox_inches='tight')
