import argparse
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy.optimize import differential_evolution, minimize
from scipy.optimize import curve_fit
from tqdm.auto import tqdm

from dna_calculator import DNACalculator, DNACalculatorTH as DNA_th


BASE_DIR = Path(__file__).resolve().parent


# -------------------------
# Physical constants
# -------------------------

qe = 1.609e-19
kb = 1.38e-23
T = 298
eps0 = 8.85e-12
lb = 1e10 * qe**2 / (4 * np.pi * eps0 * kb * T)

Lam = 5
epsb = 80
epss = 5
epsc = 2

g_b = epsb / epss
g_c = epsc / epss

cb = 0.100
Na = 6.023e23
nb = cb * Na / 1e27
kappa = np.sqrt(2 * nb * qe**2 / (1e-10 * eps0 * epsb * kb * T))

phis = 0.8 * np.pi
H = 34
hr = 3.4
f3 = 0.0

PARAM_BOUNDS = np.array([
    [26.0, 30.0],
    [0.0, 1.0],
    [0.01, 0.2],
])

a = 10.0
b = 6.0
l_c = 100.0
L_bp = 32
L_DNA = L_bp * hr


# -------------------------
# Experimental data
# -------------------------

conc_exp = np.array([0, 1, 2, 5, 10, 20, 50, 100, 200])
rel_FR_Mg = np.array([1, 0.995, 0.994, 0.991, 0.986, 0.983, 0.976, 0.969, 0.963])
rel_FR_Ca = np.array([1, 0.993, 0.990, 0.984, 0.979, 0.972, 0.966, 0.956, 0.950])
FR_Mg = np.array([5.71, 5.68, 5.67, 5.66, 5.63, 5.61, 5.57, 5.53, 5.50])
FR_Ca = np.array([5.71, 5.67, 5.65, 5.62, 5.59, 5.55, 5.52, 5.46, 5.42])

exp_c = np.array([20, 50, 100, 200])

Mg_homo_E = np.array([2.3918, 1.1692, 0.4255, -0.0440])
Mg_homo_E_e = np.array([0.0873, 0.0354, 0.0289, 0.0456])
Mg_hetero_E = np.array([3.0343, 1.7963, 1.0130, 0.5215])
Mg_hetero_E_e = np.array([0.3311, 0.2417, 0.1602, 0.1467])

Ca_homo_E = np.array([1.68223306987079, 0.405881792142694, -0.248470487508777, -0.643747074316099])
Ca_homo_E_e = np.array([0.0779041455291794, 0.0995916367860504, 0.100332779864682, 0.100093609669645])
Ca_hetero_E = np.array([2.29305234124575, 1.04285238181179, 0.317028891171946, -0.0908624617070255])
Ca_hetero_E_e = np.array([0.220911931657532, 0.0742133745579114, 0.105404745198588, 0.163136249912870])

WEIGHTS = {
    "Mg": {
        "homo": np.array([0.15, 1.2, 1.0, 1.2]),
        "hetero": np.array([0.15, 1.2, 1.0, 1.2]),
    },
    "Ca": {
        "homo": np.array([0.15, 1.2, 15.0, 15.0]),
        "hetero": np.array([0.15, 1.2, 2.0, 2.0]),
    },
}


def R_c(c, a_fit, k_fit):
    return 1 - a_fit * np.log(1 + c / k_fit)


params_Mg = curve_fit(R_c, conc_exp, rel_FR_Mg)[0]
params_Ca = curve_fit(R_c, conc_exp, rel_FR_Ca)[0]

L_nm = 12.2


def psi_s_from_RF(RF_nm):
    x = RF_nm / (2.0 * L_nm)
    x = np.clip(x, 0.0, 0.999999999)
    return 2.0 * np.arcsin(x)


def psi_s_mean_Mg(c_mM):
    return psi_s_from_RF(FR_Mg[0] * R_c(c_mM, params_Mg[0], params_Mg[1]))


def psi_s_mean_Ca(c_mM):
    return psi_s_from_RF(FR_Ca[0] * R_c(c_mM, params_Ca[0], params_Ca[1]))


# -------------------------
# Model and fit helpers
# -------------------------

@dataclass
class FitConfig:
    a: float = a
    b: float = b
    L_DNA: float = L_DNA
    l_c: float = l_c
    psi_s: float = 0.478
    Nmax: int = 10
    Mmax: int = 10
    concentrations: np.ndarray = field(default_factory=lambda: np.linspace(20, 200, 50))
    phi_grid: np.ndarray = field(default_factory=lambda: np.linspace(-np.pi, np.pi, 500))
    n_idx: np.ndarray = field(default_factory=lambda: np.arange(-5, 6))
    m_idx: np.ndarray = field(default_factory=lambda: np.arange(-5, 6))
    cos_m_grid: np.ndarray = None
    sin_m_grid: np.ndarray = None
    phi1_grid: np.ndarray = None
    U1c: np.ndarray = None
    U1s: np.ndarray = None

    def __post_init__(self):
        if self.cos_m_grid is None:
            self.cos_m_grid = np.cos(np.outer(self.phi_grid, self.m_idx))
        if self.sin_m_grid is None:
            self.sin_m_grid = np.sin(np.outer(self.phi_grid, self.m_idx))
        if self.phi1_grid is None:
            self.phi1_grid = self.phi_grid
        if self.U1c is None:
            self.U1c = np.cos(np.outer(self.phi1_grid, self.n_idx))
        if self.U1s is None:
            self.U1s = np.sin(np.outer(self.phi1_grid, self.n_idx))


def make_coarse_cfg():
    phi_c = np.linspace(-np.pi, np.pi, 150)
    n_c = np.arange(-5, 6)
    m_c = np.arange(-5, 6)
    return FitConfig(
        Nmax=5,
        Mmax=5,
        concentrations=np.linspace(20, 200, 6),
        phi_grid=phi_c,
        n_idx=n_c,
        m_idx=m_c,
        cos_m_grid=np.cos(np.outer(phi_c, m_c)),
        sin_m_grid=np.sin(np.outer(phi_c, m_c)),
        phi1_grid=phi_c,
        U1c=np.cos(np.outer(phi_c, n_c)),
        U1s=np.sin(np.outer(phi_c, n_c)),
    )


CFG_COARSE = make_coarse_cfg()


def psi_quadrature(n_nodes, psi_s):
    x, w = np.polynomial.legendre.leggauss(n_nodes)
    lo, hi = psi_s, np.pi
    psi = 0.5 * (x + 1) * (hi - lo) + lo
    weights = 0.5 * (hi - lo) * w
    return psi, weights


def create_calculator(TH, f2):
    calc = DNACalculator(Lam, g_b, kappa, epsb, g_c, epsc, lb, phis, H, hr, 1 - f2, f2, f3, TH)
    orig_p = calc.p
    calc.p = lambda k, TH=None: orig_p(int(abs(k)), TH)
    return calc


def emin_for_psi(w_mat, cfg):
    A = cfg.U1c @ w_mat
    B = cfg.U1s @ w_mat
    E = A @ cfg.cos_m_grid.T + B @ cfg.sin_m_grid.T
    return E.min()


def run_model(R, f2, K, cfg, n_psi_nodes, psi_s_fn):
    TH_list = K * cfg.concentrations / (1 + K * cfg.concentrations)
    omega_vals = np.empty_like(cfg.concentrations, dtype=float)

    for i, c in enumerate(cfg.concentrations):
        TH = TH_list[i]
        calc = create_calculator(TH, f2)

        psi_s = float(psi_s_fn(c))
        psi_nodes, psi_wts = psi_quadrature(n_psi_nodes, psi_s)
        denom = 1.0 - np.cos(psi_s)

        integrand = np.empty_like(psi_nodes)
        for j, psi in enumerate(psi_nodes):
            w_mat = np.array([[calc.w_nm(R, psi, cfg.a, cfg.b, n, m) for m in cfg.m_idx] for n in cfg.n_idx])
            Umin = emin_for_psi(w_mat, cfg)
            integrand[j] = np.exp(-Umin / 2.0) * np.sin(psi) / denom

        omega_vals[i] = np.sum(integrand * psi_wts)

    calc_th = DNA_th(Lam, g_b, kappa, epsb, g_c, epsc, lb, phis, H, hr, 1 - f2, f2, f3)
    u_hom = calc_th.u_hom_min(R, cfg.a, cfg.b, cfg.L_DNA, TH_list, method="nonlocal")
    u_non = calc_th.u_nonhom_min(R, cfg.a, cfg.b, cfg.L_DNA, cfg.l_c, TH_list, method="nonlocal")

    pred_hom = np.asarray(u_hom) + np.log(omega_vals)
    pred_heter = np.asarray(u_non) + np.log(omega_vals)
    return pred_hom, pred_heter


def ion_data(ion):
    if ion == "Mg":
        return Mg_homo_E, Mg_homo_E_e, Mg_hetero_E, Mg_hetero_E_e, psi_s_mean_Mg
    if ion == "Ca":
        return Ca_homo_E, Ca_homo_E_e, Ca_hetero_E, Ca_hetero_E_e, psi_s_mean_Ca
    raise ValueError("ion must be 'Mg' or 'Ca'")


def residual_vector(params, ion, cfg, n_psi_nodes):
    R, f2, K = params
    if np.any(params < PARAM_BOUNDS[:, 0]) or np.any(params > PARAM_BOUNDS[:, 1]):
        return np.full(8, np.nan)

    homo_E, homo_E_e, hetero_E, hetero_E_e, psi_s_fn = ion_data(ion)
    pred_hom, pred_heter = run_model(R, f2, K, cfg, n_psi_nodes, psi_s_fn)

    model_hom = np.interp(exp_c, cfg.concentrations, pred_hom)
    model_heter = np.interp(exp_c, cfg.concentrations, pred_heter)

    weights = WEIGHTS[ion]
    r_h = np.sqrt(weights["homo"]) * (model_hom - homo_E) / homo_E_e
    r_x = np.sqrt(weights["hetero"]) * (model_heter - hetero_E) / hetero_E_e
    return np.concatenate([r_h, r_x])


def chi2(params, ion, cfg, n_psi_nodes):
    r = residual_vector(params, ion, cfg, n_psi_nodes)
    if np.any(~np.isfinite(r)):
        return np.inf
    return float(r @ r)


class ObjWithBar:
    def __init__(self, func, total=None, desc="objective evals"):
        self.func = func
        self.bar = tqdm(total=total, desc=desc)
        self.count = 0

    def __call__(self, x, *args, **kwargs):
        self.count += 1
        self.bar.update(1)
        return self.func(x, *args, **kwargs)

    def close(self):
        self.bar.close()


def estimate_de_evals(dim, popsize, maxiter):
    return popsize * dim + popsize * maxiter


def fit_parameters(ion, seed, popsize, maxiter, cfg, n_psi_nodes, polish=False):
    bounds = [tuple(row) for row in PARAM_BOUNDS]
    est = estimate_de_evals(len(bounds), popsize, maxiter) + 200
    obj = ObjWithBar(
        lambda x: chi2(x, ion, cfg, n_psi_nodes),
        total=est,
        desc=f"{ion} objective evals",
    )
    try:
        de_res = differential_evolution(
            obj,
            bounds,
            seed=seed,
            popsize=popsize,
            maxiter=maxiter,
            polish=False,
            workers=1,
            updating="deferred",
        )
    finally:
        obj.close()

    if not polish:
        return de_res

    polish_res = minimize(
        lambda x: chi2(x, ion, cfg, n_psi_nodes),
        de_res.x,
        method="L-BFGS-B",
        bounds=bounds,
        options={"maxiter": 100},
    )
    if polish_res.fun < de_res.fun:
        return polish_res
    return de_res


def jacobian_bounded(fun, p, bounds=PARAM_BOUNDS, rel_step=1e-6, abs_step=1e-8):
    p = np.asarray(p, dtype=float)
    r0 = fun(p)
    if np.any(~np.isfinite(r0)):
        raise ValueError("Residuals are not finite at p.")

    J = np.zeros((r0.size, p.size), dtype=float)
    for j in range(p.size):
        h = max(abs_step, rel_step * abs(p[j]))
        lo, hi = bounds[j]

        for _ in range(12):
            dp = np.zeros_like(p)
            dp[j] = h
            can_step_down = p[j] - h >= lo
            can_step_up = p[j] + h <= hi

            if can_step_down and can_step_up:
                rp = fun(p + dp)
                rm = fun(p - dp)
                if np.all(np.isfinite(rp)) and np.all(np.isfinite(rm)):
                    J[:, j] = (rp - rm) / (2 * h)
                    break

            if can_step_up:
                rp = fun(p + dp)
                if np.all(np.isfinite(rp)):
                    J[:, j] = (rp - r0) / h
                    break

            if can_step_down:
                rm = fun(p - dp)
                if np.all(np.isfinite(rm)):
                    J[:, j] = (r0 - rm) / h
                    break

            h *= 0.25
        else:
            raise ValueError(f"Non-finite residuals when stepping parameter {j}.")
    return J


def fit_covariance(p_hat, ion, cfg, n_psi_nodes):
    fun = lambda p: residual_vector(p, ion, cfg, n_psi_nodes)
    J = jacobian_bounded(fun, p_hat)
    return np.linalg.pinv(J.T @ J)


def main():
    parser = argparse.ArgumentParser(description="Fit Mg/Ca scissors data and save phat/cov.")
    parser.add_argument("--out", default=str(BASE_DIR / "fit_results.npz"))
    parser.add_argument("--popsize", type=int, default=8)
    parser.add_argument("--maxiter", type=int, default=60)
    parser.add_argument("--seed-mg", type=int, default=42)
    parser.add_argument("--seed-ca", type=int, default=1)
    parser.add_argument(
        "--raw-de",
        action="store_true",
        help="Save the raw differential-evolution endpoint instead of the polished final fit.",
    )
    parser.add_argument(
        "--eval-mg",
        type=float,
        nargs=3,
        metavar=("R", "F2", "K"),
        help="Print the Mg coarse objective at a supplied (R, f2, K) and exit.",
    )
    parser.add_argument(
        "--eval-ca",
        type=float,
        nargs=3,
        metavar=("R", "F2", "K"),
        help="Print the Ca coarse objective at a supplied (R, f2, K) and exit.",
    )
    args = parser.parse_args()

    cfg = CFG_COARSE

    if args.eval_mg is not None:
        p = np.array(args.eval_mg, dtype=float)
        print("chi2_Mg:", chi2(p, "Mg", cfg, n_psi_nodes=6))
        print("residuals_Mg:", residual_vector(p, "Mg", cfg, n_psi_nodes=6))
        return

    if args.eval_ca is not None:
        p = np.array(args.eval_ca, dtype=float)
        print("chi2_Ca:", chi2(p, "Ca", cfg, n_psi_nodes=4))
        print("residuals_Ca:", residual_vector(p, "Ca", cfg, n_psi_nodes=4))
        return

    # Matches the nominal notebook settings: Mg objective uses 6 psi nodes, Ca uses 4.
    fit_Mg = fit_parameters(
        "Mg",
        args.seed_mg,
        args.popsize,
        args.maxiter,
        cfg,
        n_psi_nodes=6,
        polish=not args.raw_de,
    )
    phat_Mg = fit_Mg.x
    cov_Mg = fit_covariance(phat_Mg, "Mg", cfg, n_psi_nodes=6)

    print("phat_Mg:", phat_Mg)
    print("cov_Mg:\n", cov_Mg)

    fit_Ca = fit_parameters(
        "Ca",
        args.seed_ca,
        args.popsize,
        args.maxiter,
        cfg,
        n_psi_nodes=4,
        polish=not args.raw_de,
    )
    phat_Ca = fit_Ca.x
    print("phat_Ca:", phat_Ca)
    cov_Ca = fit_covariance(phat_Ca, "Ca", cfg, n_psi_nodes=6)

    print("cov_Ca:\n", cov_Ca)
    
    out_path = Path(args.out)
    if not out_path.is_absolute():
        out_path = BASE_DIR / out_path
    out_path.parent.mkdir(parents=True, exist_ok=True)

    np.savez(
        out_path,
        phat_Mg=phat_Mg,
        cov_Mg=cov_Mg,
        phat_Ca=phat_Ca,
        cov_Ca=cov_Ca,
        chi2_Mg=fit_Mg.fun,
        chi2_Ca=fit_Ca.fun,
    )
    print("Saved:", out_path)
    

if __name__ == "__main__":
    main()
