import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from dna_calculator import DNACalculator, DNACalculatorTH as DNA_th
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


def project_theta(theta):
    R, f2, K = theta
    R = max(R, 1e-9)
    f2 = float(np.clip(f2, 0.0, 1.0))
    K = max(K, 1e-12)
    return np.array([R, f2, K])



theta_0_Mg = np.array([27.91284911,  0.99134506,  0.06163151])
Sigma_0_Mg = np.array([[7.67070122e-02, 1.43708700e-02, 3.96326919e-04],
                    [1.43708700e-02, 2.88677370e-03, 1.40687928e-05],
                    [3.96326919e-04, 1.40687928e-05, 3.00636906e-05]])

theta_0_Ca = np.array([27.65372212,  0.9999797,   0.08336363])
Sigma_0_Ca = np.array([[2.38163902e-03,  4.84967809e-04, -4.13614841e-05],
                        [ 4.84967809e-04,  1.16605885e-04, -2.41717000e-05],
                        [-4.13614841e-05, -2.41717000e-05,  2.08282753e-05]])

conc_low = np.logspace(0, 1, 51)
conc_high = np.linspace(10, 250, 50)
concentrations = np.concatenate((conc_low[:-1], conc_high))

#----------------------------------------
#------- Plotting data -------
#----------------------------------------

def load_bands_for_interval(out_dir: str, interval: float):
    """
    Loads the CSVs written by save_bands_for_interval() and returns arrays:
      concentrations,
      Mg: h_lo, h_med, h_hi, nh_lo, nh_med, nh_hi,
      Ca: h_lo, h_med, h_hi, nh_lo, nh_med, nh_hi,
      Recog: r_lo_Mg, r_med_Mg, r_hi_Mg, r_lo_Ca, r_med_Ca, r_hi_Ca
    """
    tag = interval_tag(interval)

    out_dir = Path(out_dir)
    if not out_dir.is_absolute():
        out_dir = BASE_DIR / out_dir

    mg_path = out_dir / f"Mg2+_theory_MC_{tag}.csv"
    ca_path = out_dir / f"Ca2+_theory_MC_{tag}.csv"
    r_path = out_dir / f"Recog_Energy_MC_{tag}.csv"

    print(f"Loading {tag}% CI band CSVs from {out_dir}")

    missing = [path for path in (mg_path, ca_path, r_path) if not path.exists()]
    if missing:
        missing_names = ", ".join(path.name for path in missing)
        raise FileNotFoundError(
            f"Missing Monte Carlo band CSVs for {tag}% CI: {missing_names}. "
            "Run `python MC_sampling.py` before `python load_and_plot.py`."
        )

    df_Mg = pd.read_csv(mg_path)
    df_Ca = pd.read_csv(ca_path)
    df_R  = pd.read_csv(r_path)
    print(f"Loaded: {mg_path.name}, {ca_path.name}, {r_path.name}")

    c = df_Mg["c [mM]"].to_numpy()
    if not np.allclose(c, df_Ca["c [mM]"].to_numpy()):
        raise ValueError("Mg/Ca concentration grids do not match.")
    if not np.allclose(c, df_R["c [mM]"].to_numpy()):
        raise ValueError("Recog concentration grid does not match Mg/Ca.")
    print(f"Using {len(c)} concentration points from {c[0]:.3g} to {c[-1]:.3g} mM")

    def center_or_best(df, med_col, best_col):
        if med_col in df.columns:
            return df[med_col].to_numpy()
        if best_col in df.columns:
            return df[best_col].to_numpy()
        raise KeyError(f"Missing both '{med_col}' and '{best_col}' in CSV.")

    # Mg
    h_lo_Mg = df_Mg[f"MC {tag}% Homo lo"].to_numpy()
    h_med_Mg = center_or_best(df_Mg, f"MC {tag}% Homo med", "best Homo")
    h_hi_Mg = df_Mg[f"MC {tag}% Homo hi"].to_numpy()
    nh_lo_Mg = df_Mg[f"MC {tag}% Hetero lo"].to_numpy()
    nh_med_Mg = center_or_best(df_Mg, f"MC {tag}% Hetero med", "best Hetero")
    nh_hi_Mg = df_Mg[f"MC {tag}% Hetero hi"].to_numpy()

    # Ca
    h_lo_Ca = df_Ca[f"MC {tag}% Homo lo"].to_numpy()
    h_med_Ca = center_or_best(df_Ca, f"MC {tag}% Homo med", "best Homo")
    h_hi_Ca = df_Ca[f"MC {tag}% Homo hi"].to_numpy()
    nh_lo_Ca = df_Ca[f"MC {tag}% Hetero lo"].to_numpy()
    nh_med_Ca = center_or_best(df_Ca, f"MC {tag}% Hetero med", "best Hetero")
    nh_hi_Ca = df_Ca[f"MC {tag}% Hetero hi"].to_numpy()

    # recognition
    r_lo_Mg = df_R[f"MC {tag}% recog lo (Mg)"].to_numpy()
    r_med_Mg = center_or_best(df_R, f"MC {tag}% recog med (Mg)", "best recog (Mg)")
    r_hi_Mg = df_R[f"MC {tag}% recog hi (Mg)"].to_numpy()
    r_lo_Ca = df_R[f"MC {tag}% recog lo (Ca)"].to_numpy()
    r_med_Ca = center_or_best(df_R, f"MC {tag}% recog med (Ca)", "best recog (Ca)")
    r_hi_Ca = df_R[f"MC {tag}% recog hi (Ca)"].to_numpy()

    bands = {
        "Mg": {"h": (h_lo_Mg, h_med_Mg, h_hi_Mg), "nh": (nh_lo_Mg, nh_med_Mg, nh_hi_Mg)},
        "Ca": {"h": (h_lo_Ca, h_med_Ca, h_hi_Ca), "nh": (nh_lo_Ca, nh_med_Ca, nh_hi_Ca)},
        "recog": {"Mg": (r_lo_Mg, r_med_Mg, r_hi_Mg), "Ca": (r_lo_Ca, r_med_Ca, r_hi_Ca)},
    }
    return c, bands

# ----------------------------
# Your best-fit curve function
# ----------------------------

# Assumes you already have compute_curves_nominal(...) defined exactly as in your main script.
def compute_best_fit_curves(theta0_Mg, theta0_Ca, concentrations):
    """
    Returns best-fit curves (not MC median) for each ion:
      best["Mg"]["h"], best["Mg"]["nh"], best["Mg"]["r"]
      best["Ca"]["h"], best["Ca"]["nh"], best["Ca"]["r"]
    """
    best_h_Mg, best_nh_Mg, best_r_Mg = compute_curves_nominal(
        theta0_Mg[0], theta0_Mg[1], theta0_Mg[2], "Mg",
        concentrations,
        a, b, L_DNA, l_c,
        N_psi, N_phi,
        n_idx, m_idx,
        phi_grid, cos_m_grid, sin_m_grid,
        Lam, g_b, kappa, epsb, g_c, epsc, lb, phis, H, hr, f3,
        DNACalculator, DNA_th
    )
    best_h_Ca, best_nh_Ca, best_r_Ca = compute_curves_nominal(
        theta0_Ca[0], theta0_Ca[1], theta0_Ca[2], "Ca",
        concentrations,
        a, b, L_DNA, l_c,
        N_psi, N_phi,
        n_idx, m_idx,
        phi_grid, cos_m_grid, sin_m_grid,
        Lam, g_b, kappa, epsb, g_c, epsc, lb, phis, H, hr, f3,
        DNACalculator, DNA_th
    )

    return {
        "Mg": {"h": best_h_Mg, "nh": best_nh_Mg, "r": best_r_Mg},
        "Ca": {"h": best_h_Ca, "nh": best_nh_Ca, "r": best_r_Ca},
    }

def interval_tag(interval: float) -> str:
    # 0.68 -> "68", 0.95 -> "95"
    return f"{int(round(interval * 100)):02d}"

def save_bands_for_interval(out_dir, concentrations, interval,
                            h_lo_Mg, h_md_Mg, h_hi_Mg,
                            nh_lo_Mg, nh_md_Mg, nh_hi_Mg,
                            r_lo_Mg, r_md_Mg, r_hi_Mg,
                            h_lo_Ca, h_md_Ca, h_hi_Ca,
                            nh_lo_Ca, nh_md_Ca, nh_hi_Ca,
                            r_lo_Ca, r_md_Ca, r_hi_Ca):
    os.makedirs(out_dir, exist_ok=True)
    tag = interval_tag(interval)

    # --- Mg ---
    df_Mg = pd.DataFrame({
        "c [mM]": concentrations,
        f"MC {tag}% Homo lo": h_lo_Mg,
        f"MC {tag}% Homo med": h_md_Mg,
        f"MC {tag}% Homo hi": h_hi_Mg,
        f"MC {tag}% Hetero lo": nh_lo_Mg,
        f"MC {tag}% Hetero med": nh_md_Mg,
        f"MC {tag}% Hetero hi": nh_hi_Mg,
    })
    df_Mg.to_csv(os.path.join(out_dir, f"Mg2+_theory_MC_{tag}.csv"), index=False)

    # --- Ca ---
    df_Ca = pd.DataFrame({
        "c [mM]": concentrations,
        f"MC {tag}% Homo lo": h_lo_Ca,
        f"MC {tag}% Homo med": h_md_Ca,
        f"MC {tag}% Homo hi": h_hi_Ca,
        f"MC {tag}% Hetero lo": nh_lo_Ca,
        f"MC {tag}% Hetero med": nh_md_Ca,
        f"MC {tag}% Hetero hi": nh_hi_Ca,
    })
    df_Ca.to_csv(os.path.join(out_dir, f"Ca2+_theory_MC_{tag}.csv"), index=False)

    # --- Recognition (both ions in one file) ---
    df_recog = pd.DataFrame({
        "c [mM]": concentrations,
        f"MC {tag}% recog lo (Mg)": r_lo_Mg,
        f"MC {tag}% recog med (Mg)": r_md_Mg,
        f"MC {tag}% recog hi (Mg)": r_hi_Mg,
        f"MC {tag}% recog lo (Ca)": r_lo_Ca,
        f"MC {tag}% recog med (Ca)": r_md_Ca,
        f"MC {tag}% recog hi (Ca)": r_hi_Ca,
    })
    df_recog.to_csv(os.path.join(out_dir, f"Recog_Energy_MC_{tag}.csv"), index=False)

exp_c = np.array([20, 50, 100, 200])

Mg_homo_E = np.array([2.3918, 1.1692, 0.4255, -0.0440])
Mg_homo_E_e = np.array([0.0873, 0.0354, 0.0289, 0.0456])
Mg_hetero_E = np.array([3.0343, 1.7963, 1.0130, 0.5215])
Mg_hetero_E_e = np.array([0.3311, 0.2417, 0.1602, 0.1467])

Ca_homo_E = np.array([1.68223306987079, 0.405881792142694, -0.248470487508777, -0.643747074316099])
Ca_homo_E_e = np.array([0.0779041455291794, 0.0995916367860504, 0.100332779864682, 0.100093609669645])
Ca_hetero_E = np.array([2.29305234124575, 1.04285238181179, 0.317028891171946, -0.0908624617070255])
Ca_hetero_E_e = np.array([0.220911931657532, 0.0742133745579114, 0.105404745198588, 0.163136249912870])

def plot_ci_with_best_fit(out_dir: str,
                          interval: float,
                          theta0_Mg: np.ndarray,
                          theta0_Ca: np.ndarray,
                          exp_c: np.ndarray,
                          Mg_homo_E: np.ndarray, Mg_homo_E_e: np.ndarray,
                          Mg_hetero_E: np.ndarray, Mg_hetero_E_e: np.ndarray,
                          Ca_homo_E: np.ndarray, Ca_homo_E_e: np.ndarray,
                          Ca_hetero_E: np.ndarray, Ca_hetero_E_e: np.ndarray,
                          savepath: str | None = None):
    """
    Loads CI bands from CSV, computes best-fit curves, plots:
      - best-fit as dashed line
      - CI shading from loaded bands (lo/hi)
      - experimental points
    """
    concentrations, bands = load_bands_for_interval(out_dir, interval)
    print("Computing best-fit curves for plotting")
    best = compute_best_fit_curves(theta0_Mg, theta0_Ca, concentrations)
    print("Finished best-fit curve calculation")

    # Unpack bands (lo, med, hi). We only use lo/hi for shading.
    h_lo_Mg, _, h_hi_Mg   = bands["Mg"]["h"]
    nh_lo_Mg, _, nh_hi_Mg = bands["Mg"]["nh"]
    h_lo_Ca, _, h_hi_Ca   = bands["Ca"]["h"]
    nh_lo_Ca, _, nh_hi_Ca = bands["Ca"]["nh"]
    r_lo_Mg, _, r_hi_Mg   = bands["recog"]["Mg"]
    r_lo_Ca, _, r_hi_Ca   = bands["recog"]["Ca"]

    fig = plt.figure(figsize=(8, 4))
    gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1.2, 1.2], wspace=0.05)
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    ax2 = fig.add_subplot(gs[0, 3])
    ax = [ax0, ax1, ax2]

    # --- Mg (panel 0) ---
    ax[0].plot(concentrations, best["Mg"]["h"],  ls="--", lw=1.5, label=r'$\mathrm{Mg}^{2+}$ homo.', color="#07baec")
    ax[0].fill_between(concentrations, h_lo_Mg, h_hi_Mg, alpha=0.3, color="#07baec")
    ax[0].plot(concentrations, best["Mg"]["nh"], ls="--", lw=1.5, label=r'$\mathrm{Mg}^{2+}$ hetero.', color="#05a876")
    ax[0].fill_between(concentrations, nh_lo_Mg, nh_hi_Mg, alpha=0.3, color="#05a876")
    ax[0].errorbar(exp_c, Mg_homo_E,   yerr=Mg_homo_E_e,   marker="^", markersize=8, c="#07baec",
                   ls="", capsize=5, elinewidth=1.5, capthick=2.5, markeredgewidth=1, markeredgecolor="black")
    ax[0].errorbar(exp_c, Mg_hetero_E, yerr=Mg_hetero_E_e, marker="^", markersize=8, c="#05a876",
                   ls="", capsize=5, elinewidth=1.5, capthick=2.5, markeredgewidth=1, markeredgecolor="black")

    # --- Ca (panel 1) ---
    ax[1].plot(concentrations, best["Ca"]["h"],  ls="--", lw=1.5, label=r'$\mathrm{Ca}^{2+}$ homo.', color="#f59621")
    ax[1].fill_between(concentrations, h_lo_Ca, h_hi_Ca, alpha=0.3, color="#f59621")
    ax[1].plot(concentrations, best["Ca"]["nh"], ls="--", lw=1.5, label=r'$\mathrm{Ca}^{2+}$ hetero.', color="#eadf3d")
    ax[1].fill_between(concentrations, nh_lo_Ca, nh_hi_Ca, alpha=0.3, color="#eadf3d")
    ax[1].errorbar(exp_c, Ca_homo_E,   yerr=Ca_homo_E_e,   marker="^", markersize=8, c="#f59621",
                   ls="", capsize=5, elinewidth=1.5, capthick=2.5, markeredgewidth=1, markeredgecolor="black")
    ax[1].errorbar(exp_c, Ca_hetero_E, yerr=Ca_hetero_E_e, marker="^", markersize=8, c="#eadf3d",
                   ls="", capsize=5, elinewidth=1.5, capthick=2.5, markeredgewidth=1, markeredgecolor="black")

    # --- Recognition (panel 2) ---
    ax[2].plot(concentrations, best["Mg"]["r"], ls="--", lw=1.5, color="#07baec")
    ax[2].fill_between(concentrations, r_lo_Mg, r_hi_Mg, alpha=0.3, color="#07baec")
    ax[2].plot(concentrations, best["Ca"]["r"], ls="--", lw=1.5, color="#f59621")
    ax[2].fill_between(concentrations, r_lo_Ca, r_hi_Ca, alpha=0.3, color="#f59621")

    ax[2].errorbar(exp_c, Mg_homo_E - Mg_hetero_E,
                   yerr=np.sqrt(Mg_hetero_E_e**2 + Mg_homo_E_e**2),
                   marker="^", markersize=8, c="#07baec", ls="",
                   capsize=5, elinewidth=1.5, capthick=2.5, markeredgewidth=1, markeredgecolor="black")
    ax[2].errorbar(exp_c, Ca_homo_E - Ca_hetero_E,
                   yerr=np.sqrt(Ca_hetero_E_e**2 + Ca_homo_E_e**2),
                   marker="^", markersize=8, c="#f59621", ls="",
                   capsize=5, elinewidth=1.5, capthick=2.5, markeredgewidth=1, markeredgecolor="black")

    # --- Formatting (matches your existing style) ---
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
        axt.set_xscale("log")
        axt.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
        axt.set_xticks([1, 10, 100])
        axt.tick_params(which="both", direction="in", top=False, right=False)
        axt.set_xlim(1, 300)

    fig.legend(fontsize="x-small", fancybox=False, edgecolor="white", bbox_to_anchor=(0.58, 0.45))

    if savepath is not None:
        savepath = Path(savepath)
        if not savepath.is_absolute():
            savepath = BASE_DIR / savepath
        fig.savefig(savepath, format="pdf", dpi=300, bbox_inches="tight")
        print(f"Saved figure: {savepath}")

    return fig, ax

out_dir = "mc_band_csvs"
interval = 0.68
print(f"Starting load_and_plot.py for {int(round(interval * 100))}% CI")

theta_0_Mg = np.array([27.91284911, 0.99134506, 0.06163151])
theta_0_Ca = np.array([27.65372212, 0.9999797,  0.08336363])

tag = int(round(interval*100))
plot_ci_with_best_fit(
    out_dir=out_dir,
    interval=interval,
    theta0_Mg=theta_0_Mg,
    theta0_Ca=theta_0_Ca,
    exp_c=exp_c,
    Mg_homo_E=Mg_homo_E, Mg_homo_E_e=Mg_homo_E_e,
    Mg_hetero_E=Mg_hetero_E, Mg_hetero_E_e=Mg_hetero_E_e,
    Ca_homo_E=Ca_homo_E, Ca_homo_E_e=Ca_homo_E_e,
    Ca_hetero_E=Ca_hetero_E, Ca_hetero_E_e=Ca_hetero_E_e,
    savepath=f"uncert_band_bestfit_{tag}CI.pdf",
)
