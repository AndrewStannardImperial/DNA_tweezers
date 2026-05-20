import numpy as np
from scipy.special import kn, iv, ivp, kvp

__all__ = ["DNACalculator", "DNACalculatorTH"]

class DNACalculator:
    """
    Calculator for classical and nonlocal Bessel-based physical functions,
    surface charge, helical expansion harmonics, and minimal interaction metrics.
    """
    def __init__(self,
                 Lam, g_b, kappa, epsb, g_c, epsc,
                 lb, phis, H, hr, f1, f2, f3, TH=None):
        # Physical parameters
        self.Lam = Lam
        self.g_b = g_b
        self.kappa = kappa
        self.epsb = epsb
        self.g_c = g_c
        self.epsc = epsc
        # Helical / surface charge parameters
        self.lb = lb
        self.phis = phis
        self.H = H
        self.hr = hr
        self.f1 = f1
        self.f2 = f2
        self.f3 = f3
        self.TH = TH
        # Derived wavenumber
        self.g = 2 * np.pi / H

    def _resolve_TH(self, TH):
        if TH is not None:
            return TH
        if self.TH is not None:
            return self.TH
        raise ValueError("TH must be provided either when constructing DNACalculator or as a method argument.")

    # Surface charge distribution
    def p(self, n, TH=None):
        TH = self._resolve_TH(TH)
        term1 = np.cos(n * self.phis / 2)
        term2 = self.f1 + self.f2 * ((-1)**n) + self.f3 * np.cos(n * self.phis / 2)
        return term1 - TH * term2

    # Nonlocal helper functions
    def Q1(self, q):
        t1 = 2 * q**2
        t2 = self.g_b * self.kappa**2
        discr = np.sqrt((1 - 2*self.kappa*self.Lam + self.g_b*self.kappa**2*self.Lam**2)
                        * (1 + 2*self.kappa*self.Lam + self.g_b*self.kappa**2*self.Lam**2))
        t3 = (1/self.Lam**2) * (1 + discr)
        return np.sqrt((t1 + t2 + t3) / 2)

    def Q2(self, q):
        t1 = 2 * q**2
        t2 = self.g_b * self.kappa**2
        discr = np.sqrt((1 - 2*self.kappa*self.Lam + self.g_b*self.kappa**2*self.Lam**2)
                        * (1 + 2*self.kappa*self.Lam + self.g_b*self.kappa**2*self.Lam**2))
        t3 = (1/self.Lam**2) * (1 - discr)
        return np.sqrt((t1 + t2 + t3) / 2)

    def g1(self):
        top = self.g_b * self.Lam**2 * self.Q1(0)**2 - 1
        bot = np.sqrt((1 - 2*self.kappa*self.Lam + self.g_b*self.kappa**2*self.Lam**2)
                      * (1 + 2*self.kappa*self.Lam + self.g_b*self.kappa**2*self.Lam**2))
        return top / bot

    def g2(self):
        top = 1 - self.g_b * self.Lam**2 * self.Q2(0)**2
        bot = np.sqrt((1 - 2*self.kappa*self.Lam + self.g_b*self.kappa**2*self.Lam**2)
                      * (1 + 2*self.kappa*self.Lam + self.g_b*self.kappa**2*self.Lam**2))
        return top / bot

    # Core Bessel-based functions
    def Omega(self, q, x, y, z, n, m, method='nonlocal'):
        if method == 'classical':
            kq = np.sqrt(self.kappa**2 + q**2)
            bess = kn(n-m, kq*x) * iv(n, kq*y) * iv(m, kq*z)
            return (4*np.pi * y * z / self.epsb) * ((-1)**m) * bess
        elif method == 'nonlocal':
            Q1, Q2 = self.Q1(q), self.Q2(q)
            b1 = self.g1() * kn(n-m, Q1*x) * iv(n, Q1*y) * iv(m, Q1*z)
            b2 = self.g2() * kn(n-m, Q2*x) * iv(n, Q2*y) * iv(m, Q2*z)
            return (4*np.pi * y * z / self.epsb) * ((-1)**m) * (b1 + b2)
        else:
            raise ValueError("method must be 'classical' or 'nonlocal'")

    def Aa(self, x, y, q, m, method='nonlocal'):
        if method == 'classical':
            kq = np.sqrt(self.kappa**2 + q**2)
            bess = kn(m, kq*x) * iv(m, kq*y)
            return 4*np.pi * bess / self.epsb
        elif method == 'nonlocal':
            Q1, Q2 = self.Q1(q), self.Q2(q)
            b1 = self.g1() * kn(m, Q1*x) * iv(m, Q1*y)
            b2 = self.g2() * kn(m, Q2*x) * iv(m, Q2*y)
            return 4*np.pi * (b1 + b2) / self.epsb
        else:
            raise ValueError("method must be 'classical' or 'nonlocal'")

    def Aa_p(self, x, y, q, m, method='nonlocal'):
        if method == 'classical':
            kq = np.sqrt(self.kappa**2 + q**2)
            term = kq * kn(m, kq*x) * ivp(m, kq*y)
            return 4*np.pi * term / self.epsb
        elif method == 'nonlocal':
            Q1, Q2 = self.Q1(q), self.Q2(q)
            t1 = self.g1()*Q1 * kn(m, Q1*x) * ivp(m, Q1*y)
            t2 = self.g2()*Q2 * kn(m, Q2*x) * ivp(m, Q2*y)
            return 4*np.pi * (t1 + t2) / self.epsb
        else:
            raise ValueError("method must be 'classical' or 'nonlocal'")

    def Ab(self, x, y, q, m, method='nonlocal'):
        if method == 'classical':
            kq = np.sqrt(self.kappa**2 + q**2)
            bess = kn(m, kq*y) * iv(m, kq*x)
            return 4*np.pi * bess / self.epsb
        elif method == 'nonlocal':
            Q1, Q2 = self.Q1(q), self.Q2(q)
            b1 = self.g1() * kn(m, Q1*y) * iv(m, Q1*x)
            b2 = self.g2() * kn(m, Q2*y) * iv(m, Q2*x)
            return 4*np.pi * (b1 + b2) / self.epsb
        else:
            raise ValueError("method must be 'classical' or 'nonlocal'")

    def Ab_p(self, x, y, q, m, method='nonlocal'):
        if method == 'classical':
            kq = np.sqrt(self.kappa**2 + q**2)
            term = kq * kvp(m, kq*y) * iv(m, kq*x)
            return 4*np.pi * term / self.epsb
        elif method == 'nonlocal':
            Q1, Q2 = self.Q1(q), self.Q2(q)
            t1 = self.g1()*Q1 * kvp(m, Q1*y) * iv(m, Q1*x)
            t2 = self.g2()*Q2 * kvp(m, Q2*y) * iv(m, Q2*x)
            return 4*np.pi * (t1 + t2) / self.epsb
        else:
            raise ValueError("method must be 'classical' or 'nonlocal'")

    def Xi(self, q, m, a, b, method='nonlocal'):
        prefac = -a / b
        bess_frac = ivp(m, abs(q)*b) / iv(m, abs(q)*b)
        if method == 'classical':
            g_cst = self.epsc / self.epsb
            top = self.Aa_p(a, b, q, m, method) - g_cst * abs(q) * bess_frac * self.Aa(a, b, q, m, method)
            bot = self.Ab_p(b, b, q, m, method) - g_cst * abs(q) * bess_frac * self.Ab(b, b, q, m, method)
        elif method == 'nonlocal':
            top = self.Aa_p(a, b, q, m, method) - self.g_c * abs(q) * bess_frac * self.Aa(a, b, q, m, method)
            bot = self.Ab_p(b, b, q, m, method) - self.g_c * abs(q) * bess_frac * self.Ab(b, b, q, m, method)
        else:
            raise ValueError("method must be 'classical' or 'nonlocal'")
        return prefac * (top / bot)
    
    def W(self, q, R, a, b, n, m, method='nonlocal'):
        term1 = self.Omega(q, R, a, a, n, m, method)
        term2 = self.Xi(q, m, a, b, method) * self.Omega(q, R, a, b, n, m, method)
        term3 = self.Xi(-q, -n, a, b, method) * self.Omega(q, R, b, a, n, m, method)
        term4 = self.Xi(-q, -n, a, b, method) * self.Xi(q, m, a, b, method) * self.Omega(q, R, b, b, n, m, method)
        return term1 + term2 + term3 + term4

    # Helical expansion harmonics
    def a0(self, R, a, b, TH=None, method='nonlocal'):
        TH = self._resolve_TH(TH)
        prefac = 2*self.lb*(1-TH)**2 / (np.pi*a**2*self.hr**2)
        term = self.W(0, R, a, b, 0, 0, method=method)
        return prefac * term

    def an(self, R, a, b, n, TH=None, method='nonlocal'):
        TH = self._resolve_TH(TH)
        prefac = 4*self.lb / (np.pi*a**2*self.hr**2)
        surf = self.p(n, TH)**2
        term = self.W(n*self.g, R, a, b, n, n, method=method)
        return prefac * surf * term

    def u(self, R, dphi, a, b, nmax, TH=None, method='nonlocal'):
        TH = self._resolve_TH(TH)
        energy = self.a0(R, a, b, TH, method)
        for n in range(1, nmax+1):
            energy += self.an(R, a, b, n, TH, method) * np.cos(n*dphi)
        return energy

    # Minimal xi and interaction metrics
    def xi_min_nh(self, R, a, b, L, l_c, TH=None, method='nonlocal'):
        TH = self._resolve_TH(TH)
        Lbar = L / l_c
        a1 = self.an(R, a, b, 1, TH, method)
        a2 = self.an(R, a, b, 2, TH, method)
        cos_xi = -a1 * (1 - np.exp(-Lbar/2))
        cos_xi /= a2 * (1 - np.exp(-2*Lbar))
        cos_xi = np.clip(cos_xi, -1.0, 1.0)
        return np.arccos(cos_xi)

    def xi_min_h(self, R, a, b, TH=None, method='nonlocal'):
        TH = self._resolve_TH(TH)
        a1 = self.an(R, a, b, 1, TH, method)
        a2 = self.an(R, a, b, 2, TH, method)
        cos_xi = -a1 / (4 * a2)
        cos_xi = np.clip(cos_xi, -1.0, 1.0)
        return np.arccos(cos_xi)

    def nu1_min(self, R, a, b, L, l_c, TH=None, method='nonlocal'):
        TH = self._resolve_TH(TH)
        Lbar = L / l_c
        xi_s = self.xi_min_nh(R, a, b, L, l_c, TH, method)
        return (2/Lbar) * (1 - np.exp(-Lbar/2)) * np.cos(xi_s)

    def nu2_min(self, R, a, b, L, l_c, TH=None, method='nonlocal'):
        TH = self._resolve_TH(TH)
        Lbar = L / l_c
        xi_s = self.xi_min_nh(R, a, b, L, l_c, TH, method)
        return (1/(2*Lbar)) * (1 - np.exp(-2*Lbar)) * np.cos(2*xi_s)

    def u_nonhom_min(self, R, a, b, L, l_c, TH=None, method='nonlocal'):
        TH = self._resolve_TH(TH)
        u0 = self.a0(R, a, b, TH, method)
        u1 = self.an(R, a, b, 1, TH, method) * self.nu1_min(R, a, b, L, l_c, TH, method)
        u2 = self.an(R, a, b, 2, TH, method) * self.nu2_min(R, a, b, L, l_c, TH, method)
        return L * (u0 + u1 + u2)

    def u_hom_min(self, R, a, b, L, TH=None, method='nonlocal'):
        TH = self._resolve_TH(TH)
        u0 = self.a0(R, a, b, TH, method)
        xi_s = self.xi_min_h(R, a, b, TH, method)
        u1 = self.an(R, a, b, 1, TH, method) * np.cos(xi_s)
        u2 = self.an(R, a, b, 2, TH, method) * np.cos(2*xi_s)
        return L * (u0 + u1 + u2)

    def u_recog_min(self, R, a, b, L, l_c, TH=None, method='nonlocal'):
        TH = self._resolve_TH(TH)
        return self.u_nonhom_min(R, a, b, L, l_c, TH, method) - self.u_hom_min(R, a, b, L, TH, method)
    

    ### SKEWED MOLECULES ###

    def u_nm(self, psi, n, m):
        top = n*self.g - m*self.g*np.cos(psi)
        bot = np.sin(psi)
        return top/bot

    def I_nm_1(self, R, psi, x, y, n, m):
        Q1m = self.Q1(m*self.g)
        Q1n = self.Q1(n*self.g)
        unm = self.u_nm(psi, n, m)
        umn = self.u_nm(psi, m, n)
        prefac = np.pi*(-1.0)**(n+m)*(1-self.g_b*(self.Lam**2)*(self.Q1(0)**2))
        t1 = np.exp(-R*np.sqrt(unm**2+Q1m**2))
        t2 = iv(n, Q1n*x)*iv(m, Q1m*y) / (np.sqrt(unm**2+Q1m**2))
        with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
            t3 = ((np.sqrt(umn**2 + Q1n**2)+umn)/(np.sqrt(umn**2 + Q1n**2)-umn))**(n/2)
            t4 = ((np.sqrt(unm**2 + Q1m**2)+unm)/(np.sqrt(unm**2 + Q1m**2)-unm))**(m/2)
            return prefac*t1*t2*t3*t4
    
    def I_nm_2(self, R, psi, x, y, n, m):
        Q2m = self.Q2(m*self.g)
        Q2n = self.Q2(n*self.g)
        unm = self.u_nm(psi, n, m)
        umn = self.u_nm(psi, m, n)
        prefac = np.pi*(-1.0)**(n+m)*(1-self.g_b*(self.Lam**2)*(self.Q2(0)**2))
        t1 = np.exp(-R*np.sqrt(unm**2+Q2m**2))
        t2 = iv(n, Q2n*x)*iv(m, Q2m*y) / (np.sqrt(unm**2+Q2m**2))
        with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
            t3 = ((np.sqrt(umn**2 + Q2n**2)+umn)/(np.sqrt(umn**2 + Q2n**2)-umn))**(n/2)
            t4 = ((np.sqrt(unm**2 + Q2m**2)+unm)/(np.sqrt(unm**2 + Q2m**2)-unm))**(m/2)
            return prefac*t1*t2*t3*t4
    
    def I_nm(self, R, psi, x, y, n, m):
        Q1m = self.Q1(m*self.g)
        Q2m = self.Q2(m*self.g)
        prefac = (4*np.pi)/(self.epsb*(self.Lam**2))
        t1 = (1/(Q2m**2 - Q1m**2))
        t2 = self.I_nm_1(R, psi, x, y, n, m) - self.I_nm_2(R, psi, x, y, n, m)
        return prefac*t1*t2
    
    def w_nm(self, R, psi, a, b, n, m, TH=None):
        TH = self._resolve_TH(TH)
        prefac = 2*self.lb / (np.pi*(a**2)*(self.hr**2) * np.abs(np.sin(psi)))
        t1 = self.p(n, TH)*self.p(-m, TH)
        t2 = a**2 * self.I_nm(R, psi, a, a, n, m)
        t3 = a*b*(self.Xi(m*self.g, m, a, b)*self.I_nm(R, psi, a, b, n, m) + self.Xi(-n*self.g, -n, a, b)*self.I_nm(R, psi, b, a, n, m))
        t4 = b**2 * self.Xi(-n*self.g, -n, a, b)*self.Xi(m*self.g, m, a, b)*self.I_nm(R, psi, b, b, n, m)
        return prefac*t1*(t2+t3+t4)
    
    def U_psi(self, R, psi, phi1, phi2, a, b, Nmax = 10, Mmax = 10):
        total_U = 0
        for n in range(Nmax):
            for m in range(Mmax):
                c = np.cos(n*phi1 - m*phi2)
                wnm = self.w_nm(R, psi, a, b, n, m)
                t_nm = c*wnm
                total_U+=t_nm
        return total_U
    
    
DNACalculatorTH = DNACalculator
