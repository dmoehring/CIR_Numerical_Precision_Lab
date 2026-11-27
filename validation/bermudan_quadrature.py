#!/usr/bin/env python3
"""
High-precision (mpmath, 30 digits) Bermudan reference for the standard contract, written independently of the Java
class lab.cir.analytic.CirQuadratureBermudan. Shared with it is only the mathematics:

    C_j(r) = P(e_j, e_{j+1}; r) * E^{e_{j+1}}[ max(C_{j+1}, E_{j+1})(r_{e_{j+1}}) | r_{e_j} = r ]

with the forward-measure law of r (scaled noncentral chi-square, Brigo & Mercurio 2006, Sect. 3.2.3).
Different from the Java implementation:

* density: Bessel form  f(y) = 1/2 e^{-(y+d)/2} (y/d)^{nu/2} I_nu(sqrt(d y)),  nu = v/2 - 1  (mpmath.besseli)
  instead of the Poisson mixture;
* quadrature: mpmath tanh-sinh (mp.quad) between break points instead of adaptive Gauss-Kronrod;
* interpolation: Chebyshev points of the second kind (Lobatto) on fewer, wider panels, root finding with
  mp.findroot (Anderson) instead of bisection;
* zero bonds: closed form A(tau) exp(-B(tau) r) in 30-digit arithmetic.

The script also checks the forward-measure density (mass 1, tower property P(0,T) E^T[P(T,S)] = P(0,S)) and prices
the four European control cases.

Usage (from the project root):
    python3 -m venv /tmp/venv && /tmp/venv/bin/pip install mpmath
    /tmp/venv/bin/python validation/bermudan_quadrature.py > validation/bermudan_quadrature.out
Runtime: a few minutes.
"""
import time

import mpmath as mp

mp.mp.dps = 30

THETA, MU, SIGMA, R0 = mp.mpf('0.50'), mp.mpf('0.04'), mp.mpf('0.10'), mp.mpf('0.03')
CASHFLOWS = [(mp.mpf(t), mp.mpf('4.00') + (mp.mpf('100.00') if t == 5 else 0)) for t in (1, 2, 3, 4, 5)]
STRIKE = mp.mpf('100.00')
EXERCISE = [mp.mpf(1), mp.mpf(2), mp.mpf(3), mp.mpf(4)]

R_MAX = mp.mpf('0.6')
PANELS = 24
NODES = 15  # Chebyshev-Lobatto points per panel

H = mp.sqrt(THETA ** 2 + 2 * SIGMA ** 2)
DOF = 4 * THETA * MU / SIGMA ** 2


def zero_bond(tau, r):
    if tau == 0:
        return mp.mpf(1)
    em1 = mp.expm1(H * tau)
    den = 2 * H + (THETA + H) * em1
    a = (2 * H * mp.exp((THETA + H) * tau / 2) / den) ** (2 * THETA * MU / SIGMA ** 2)
    return a * mp.exp(-2 * em1 / den * r)


def bond(e, r):
    return mp.fsum(c * zero_bond(t - e, r) for t, c in CASHFLOWS if t > e)


def exercise_value(e, r):
    return max(bond(e, r) - STRIKE, mp.mpf(0))


class Transition:
    """Law of r_{t+tau} given r_t under the (t+tau)-forward measure: c r ~ ncx2(DOF, k r)."""

    def __init__(self, tau):
        rho = 2 * H / (SIGMA ** 2 * mp.expm1(H * tau))
        psi = (THETA + H) / SIGMA ** 2
        self.c = 2 * (rho + psi)
        self.k = 2 * rho ** 2 * mp.exp(H * tau) / (rho + psi)
        self.tau = tau

    def density(self, x, r):
        y = self.c * x
        d = self.k * r
        if y <= 0:
            return mp.mpf(0)
        nu = DOF / 2 - 1
        if d == 0:
            return self.c * y ** nu * mp.exp(-y / 2) / (2 ** (DOF / 2) * mp.gamma(DOF / 2))
        return self.c * mp.mpf(1) / 2 * mp.exp(-(y + d) / 2) * (y / d) ** (nu / 2) * mp.besseli(nu, mp.sqrt(d * y))

    def support(self, r):
        d = self.k * r
        mean = DOF + d
        sd = mp.sqrt(2 * (DOF + 2 * d))
        return max(mp.mpf(0), (mean - 20 * sd) / self.c), (mean + 25 * sd + 120) / self.c


class PiecewiseChebyshev:
    def __init__(self, f):
        self.width = R_MAX / PANELS
        self.t = [mp.cos(mp.pi * k / (NODES - 1)) for k in range(NODES)]
        self.w = [(-1) ** k * (mp.mpf(1) / 2 if k in (0, NODES - 1) else 1) for k in range(NODES)]
        self.values = []
        for p in range(PANELS):
            mid = (p + mp.mpf(1) / 2) * self.width
            self.values.append([f(mid + self.width / 2 * tk) for tk in self.t])

    def __call__(self, x):
        x = min(max(x, mp.mpf(0)), R_MAX)
        p = min(PANELS - 1, int(x / self.width))
        u = (x - (p + mp.mpf(1) / 2) * self.width) / (self.width / 2)
        num = den = mp.mpf(0)
        for tk, wk, yk in zip(self.t, self.w, self.values[p]):
            if u == tk:
                return yk
            q = wk / (u - tk)
            num += q * yk
            den += q
        return num / den

    def breaks(self):
        return [p * self.width for p in range(PANELS + 1)]


def strike_kink(e):
    if bond(e, 0) <= STRIKE:
        return None
    return mp.findroot(lambda r: bond(e, r) - STRIKE, (mp.mpf(0), mp.mpf(1)), solver='anderson')


def continuation_kinks(e, cont):
    g = lambda r: cont(r) - exercise_value(e, r)
    xs = [R_MAX * i / (PANELS * 60) for i in range(PANELS * 60 + 1)]
    roots = []
    prev_x, prev_g = xs[0], g(xs[0])
    for x in xs[1:]:
        gx = g(x)
        if (gx > 0) != (prev_g > 0):
            roots.append(mp.findroot(g, (prev_x, x), solver='anderson'))
        prev_x, prev_g = x, gx
    return roots


def value_function(e, cont):
    kinks = [k for k in [strike_kink(e)] if k is not None]
    breaks = []
    if cont is None:
        f = lambda x: exercise_value(e, x)
        ck = []
    else:
        ck = continuation_kinks(e, cont)
        kinks += ck
        breaks = cont.breaks()
        f = lambda x: max(cont(x), exercise_value(e, x))
    return f, sorted(kinks), breaks, ck


def expect(vf, tr, r):
    f, kinks, breaks, _ = vf
    lo, hi = tr.support(r)
    pts = sorted(set([lo, hi] + [x for x in kinks + breaks if lo < x < hi]))
    return mp.quad(lambda x: f(x) * tr.density(x, r), pts)


def bermudan(exercise):
    m = len(exercise)
    vf = value_function(exercise[-1], None)
    boundaries = {m - 1: vf[1]}
    for j in range(m - 2, -1, -1):
        tr = Transition(exercise[j + 1] - exercise[j])
        cur = vf
        t0 = time.time()
        cont = PiecewiseChebyshev(lambda r: zero_bond(tr.tau, r) * expect(cur, tr, r))
        vf = value_function(exercise[j], cont)
        boundaries[j] = vf[3]
        print(f"#   date e{j + 1}: continuation interpolant built in {time.time() - t0:.0f} s", flush=True)
    tr = Transition(exercise[0])
    return zero_bond(exercise[0], R0) * expect(vf, tr, R0), boundaries


def main():
    print(f"# mpmath {mp.__version__}, {mp.mp.dps} digits; panels {PANELS} x {NODES} Lobatto nodes on [0, {R_MAX}]")
    tr = Transition(mp.mpf(1))
    mass = mp.quad(lambda x: tr.density(x, R0), [0, 0.05, 0.1, 0.3, 1])
    tower = zero_bond(1, R0) * mp.quad(lambda x: zero_bond(2, x) * tr.density(x, R0), [0, 0.05, 0.1, 0.3, 1])
    print(f"# density check: mass - 1 = {mp.nstr(mass - 1, 3)}, tower P(0,1)E[P(1,3)] - P(0,3) = "
          f"{mp.nstr(tower - zero_bond(3, R0), 3)}")
    for j, e in enumerate(EXERCISE):
        v, _ = bermudan([e])
        print(f"european e{j + 1} {mp.nstr(v, 25)}")
    t0 = time.time()
    v, boundaries = bermudan(EXERCISE)
    print(f"bermudan {mp.nstr(v, 25)}")
    for j in sorted(boundaries):
        print(f"boundary e{j + 1} " + " ".join(mp.nstr(b, 20) for b in boundaries[j]))
    print(f"# Bermudan runtime {time.time() - t0:.0f} s")


if __name__ == '__main__':
    main()
