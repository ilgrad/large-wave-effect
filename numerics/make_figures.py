# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2", "matplotlib>=3", "scipy>=1.11"]
# ///
"""Generate the figures for paper/main.tex (vector PDF, headless).

Fig 1 (ring): A_N vs ln N. The ceiling U_N=(1/2N)sum csc(pi r/N)~(1/pi)ln N is exact; for primes
and 2^m, A_N=U_N (Thms); for composite N the exact subtorus value sits a little below (points
taken from verify_exact_AN.py). Establishes A_N=Theta(ln N) for all N.

Fig 2 (Schrodinger): B_N=sup_t ||e^{-itL_N}||_{inf->inf} vs sqrt(N), with a linear fit (R^2=1.0):
the ell^inf amplification grows ballistically as Theta(sqrt N).
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
from scipy.special import airy, gamma, jv

OUT = __file__.rsplit("/", 2)[0] + "/paper/figures"


def _save(fig: Figure, name: str) -> None:
    fig.tight_layout()
    # dpi controls only rasterized artists (the pcolormesh density plots); vector
    # elements (axes, text, lines) stay vector. Keeps fig_branched et al. small.
    fig.savefig(f"{OUT}/{name}.pdf", dpi=200)
    plt.close(fig)


def U(n: int) -> float:
    r = np.arange(1, n)
    return float(np.sum(1.0 / np.sin(np.pi * r / n)) / (2 * n))


def schrodinger_B(n: int) -> float:
    lam = 4.0 * np.sin(np.pi * np.arange(n) / n) ** 2
    ts = np.linspace(0.01, 6.0 * n, 8000)
    return max(float(np.sum(np.abs(np.fft.ifft(np.exp(-1j * t * lam))))) for t in ts)


def fig_ring() -> None:
    ns = np.arange(6, 1025)
    u = np.array([U(int(n)) for n in ns])
    primes2m = {5: 0.5506, 7: 0.6585, 8: 0.7012, 11: 0.8029, 13: 0.8562, 16: 0.9224}
    composite = {6: 0.5594, 9: 0.6717, 12: 0.7546, 15: 0.7542, 18: 0.8747, 24: 0.9618}

    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    ax.plot(np.log(ns), u, color="0.4", lw=1.4, label=r"ceiling $U_N$")
    ax.plot(np.log(ns), np.log(ns) / np.pi + (np.euler_gamma + np.log(2 / np.pi)) / np.pi,
            "k--", lw=1.0, label=r"$\frac{1}{\pi}\ln N + 0.040$")
    ax.scatter([np.log(n) for n in primes2m], list(primes2m.values()),
               marker="o", color="C0", zorder=5, label=r"prime / $2^m$  ($A_N=U_N$)")
    ax.scatter([np.log(n) for n in composite], list(composite.values()),
               marker="x", color="C3", zorder=5, label=r"composite (exact $A_N$, subtorus)")
    ax.set_xlabel(r"$\ln N$")
    ax.set_ylabel(r"$A_N$")
    ax.set_title(r"Ring large wave: $A_N=\Theta(\ln N)$ for all $N$")
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(alpha=0.3)
    _save(fig, "fig_ring")


def fig_schrodinger() -> None:
    ns = np.array([16, 32, 64, 128, 256, 512])
    b = np.array([schrodinger_B(int(n)) for n in ns])
    a, c = np.polyfit(np.sqrt(ns), b, 1)
    xs = np.linspace(0, np.sqrt(ns[-1]) * 1.02, 100)

    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    ax.plot(xs, a * xs + c, "k--", lw=1.0, label=rf"fit $B_N\approx{a:.3f}\sqrt{{N}}{c:+.2f}$")
    ax.scatter(np.sqrt(ns), b, marker="s", color="C2", zorder=5, label=r"$B_N$ (operator norm)")
    for n, bb in zip(ns, b, strict=True):
        ax.annotate(f"N={n}", (np.sqrt(n), bb), textcoords="offset points", xytext=(4, -9), fontsize=7)
    ax.set_xlabel(r"$\sqrt{N}$")
    ax.set_ylabel(r"$B_N$")
    ax.set_title(r"Schrodinger ring: $B_N=\Theta(\sqrt{N})$ (ballistic)")
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(alpha=0.3)
    _save(fig, "fig_schrodinger")


def fig_bessel() -> None:
    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    for t, c in [(20, "C0"), (80, "C1"), (320, "C3")]:
        n = np.arange(int(2 * t * 1.45))
        ax.plot(n, np.abs(jv(n, 2 * t)), color=c, lw=1.0, label=rf"$t={t}$ (front $n\approx{2 * t}$)")
        ax.axvline(2 * t, color=c, ls=":", lw=0.8)
    ax.set_xlabel(r"site $n$")
    ax.set_ylabel(r"$|K_\infty(n,t)|=|J_n(2t)|$")
    ax.set_title(r"Bessel kernel $|J_n(2t)|$: ballistic front, $(2t)^{-1/3}$ decay")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    _save(fig, "fig_bessel")


def fig_dispersion() -> None:
    th = np.linspace(0, np.pi, 400)
    exact = 2 * np.sin(th / 2)
    cont = th
    taylor = th - th**3 / 24  # 2 sin(theta/2) = theta - theta^3/24 + ...
    a2 = 0.25 - 1 / np.pi**2
    pade = th / np.sqrt(1 + a2 * th**2)
    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    ax.plot(th, exact, "k", lw=1.8, label=r"exact lattice $2\sin(\theta/2)$")
    ax.plot(th, pade, "C2--", lw=1.2, label=r"two-point Padé $\theta/\sqrt{1+a^2\theta^2}$")
    ax.plot(th, taylor, "C1-.", lw=1.0, label=r"Taylor $\theta-\theta^3/24$")
    ax.plot(th, cont, "C0:", lw=1.0, label=r"continuum $\theta$")
    ax.scatter([np.pi], [2.0], color="k", zorder=5)
    ax.set_xlabel(r"$\theta$")
    ax.set_ylabel(r"$\omega(\theta)$")
    ax.set_title("Dispersion: Padé matches both band edges")
    ax.set_xlim(0, np.pi)
    ax.set_ylim(0, np.pi)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    _save(fig, "fig_dispersion")


def fig_zeta() -> None:
    s = np.linspace(0.02, 0.49, 240)
    i_s = 4.0 ** (-s) * gamma(0.5 - s) / (np.sqrt(np.pi) * gamma(1.0 - s))
    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    ax.plot(s, i_s, "k", lw=1.6, label=r"density $I(s)=4^{-s}\Gamma(\frac{1}{2}-s)/(\sqrt{\pi}\,\Gamma(1-s))$")
    for nn, c in [(32, "C0"), (256, "C3")]:
        r = np.arange(1, nn)
        lam = 4 * np.sin(np.pi * r / nn) ** 2
        zn = np.array([np.sum(lam ** (-si)) / nn for si in s])
        ax.plot(s, zn, c + "--", lw=0.9, label=rf"$(1/N)\zeta_{{L_N}}(s),\ N={nn}$")
    ax.axvline(0.5, ls=":", color="gray")
    ax.annotate(r"$s=\frac{1}{2}$: $U_N=\zeta(\frac{1}{2})/N\sim\frac{1}{\pi}\ln N$",
                xy=(0.5, 3.2), xytext=(0.27, 3.4), fontsize=8,
                arrowprops={"arrowstyle": "->", "color": "gray"})
    ax.set_xlabel(r"$s$")
    ax.set_ylabel(r"$(1/N)\,\zeta_{L_N}(s)$")
    ax.set_title(r"Spectral zeta density: pole at $s=\frac{1}{2}$ drives the log ceiling")
    ax.set_xlim(0, 0.55)
    ax.set_ylim(0, 4)
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(alpha=0.3)
    _save(fig, "fig_zeta")


def _totient(n: int) -> int:
    r, m, p = n, n, 2
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            r -= r // p
        p += 1
    if m > 1:
        r -= r // m
    return r


def _green(n: int, t: float) -> np.ndarray:
    om = 2.0 * np.sin(np.pi * np.arange(n) / n)
    nz = om > 1e-12
    c = np.zeros(n)
    c[nz] = np.sin(om[nz] * t) / om[nz]
    return np.real(np.fft.ifft(c))


def fig_dos() -> None:
    n = 4000
    lam = 4.0 * np.sin(np.pi * np.arange(1, n) / n) ** 2
    x = np.linspace(0.02, 3.98, 500)
    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    ax.hist(lam, bins=70, density=True, color="C0", alpha=0.55, label=r"empirical $\lambda_r$")
    ax.plot(x, 1.0 / (np.pi * np.sqrt(x * (4 - x))), "k", lw=1.6,
            label=r"arcsine $1/(\pi\sqrt{\lambda(4-\lambda)})$")
    ax.set_xlabel(r"$\lambda$")
    ax.set_ylabel("density of states")
    ax.set_title(r"Spectral density of states: arcsine law ($N=4000$)")
    ax.set_ylim(0, 1.6)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    _save(fig, "fig_dos")


def fig_value_dist() -> None:
    n = 61
    om = 2.0 * np.sin(np.pi * np.arange(1, n) / n)
    ts = np.linspace(0.0, 150000.0, 300000)
    vals = (np.sin(np.outer(ts, om)) * (1.0 / n / om)).sum(axis=1)
    sd = float(vals.std())
    x = np.linspace(vals.min(), vals.max(), 400)
    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    ax.hist(vals, bins=90, density=True, color="C3", alpha=0.55, label=r"values of $u_0(t)$")
    ax.plot(x, np.exp(-x * x / (2 * sd * sd)) / (sd * np.sqrt(2 * np.pi)), "k--", lw=1.2,
            label="Gaussian (same variance)")
    ax.set_xlabel(r"$u_0(t)$")
    ax.set_ylabel("probability density")
    ax.set_title(r"Bohr-Jessen value law: var $=1/12$, platykurtic (kurt $\approx2.4$)")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    _save(fig, "fig_value_dist")


def fig_heatmap() -> None:
    n = 64
    ts = np.linspace(0.0, 6.0 * n, 500)
    field = np.array([np.abs(_green(n, float(t))) for t in ts])  # |u_j(t)|, shape (T, N)
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    im = ax.pcolormesh(np.arange(n), ts, field, shading="auto", cmap="magma", rasterized=True)
    fig.colorbar(im, ax=ax, label=r"$|u_j(t)|$")
    ax.set_xlabel(r"site $j$")
    ax.set_ylabel(r"time $t$")
    ax.set_title(r"Space-time $|u_j(t)|$ on the ring ($N=64$): impulse spreading")
    _save(fig, "fig_heatmap")


def fig_qrank() -> None:
    ns = np.arange(3, 49)
    rank = np.array([_totient(2 * int(x)) // 2 for x in ns])
    half = ns // 2

    def cls(x: int) -> str:
        if all(x % d for d in range(2, int(x**0.5) + 1)):
            return "prime"
        return "2^m" if (x & (x - 1)) == 0 else "composite"

    kinds = np.array([cls(int(x)) for x in ns])
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.plot(ns, half, "k:", lw=1.0, label=r"full rank $\lfloor N/2\rfloor$")
    for k, c, mk in [("prime", "C0", "o"), ("2^m", "C2", "s"), ("composite", "C3", "x")]:
        sel = kinds == k
        ax.scatter(ns[sel], rank[sel], c=c, marker=mk, s=28, label=k, zorder=5)
    ax.set_xlabel(r"$N$")
    ax.set_ylabel(r"$\mathrm{rank}_\mathbb{Q}=\dim\mathbb{T}_N=\varphi(2N)/2$")
    ax.set_title(r"Subtorus dimension $\varphi(2N)/2$; full rank iff $N$ prime or $2^m$")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    _save(fig, "fig_qrank")


def _peregrine(x: np.ndarray, t: np.ndarray) -> np.ndarray:
    xx, tt = np.meshgrid(x, t)
    return np.exp(2j * tt) * (1.0 - 4.0 * (1.0 + 4j * tt) / (1.0 + 4.0 * xx**2 + 16.0 * tt**2))


def fig_peregrine() -> None:
    x = np.linspace(-4, 4, 400)
    t = np.linspace(-2, 2, 400)
    amp = np.abs(_peregrine(x, t))
    fig, (axl, axr) = plt.subplots(1, 2, figsize=(8.6, 3.7), width_ratios=[1.25, 1])
    im = axl.pcolormesh(x, t, amp, shading="auto", cmap="inferno", vmin=0, vmax=3, rasterized=True)
    fig.colorbar(im, ax=axl, label=r"$|\psi(x,t)|$")
    axl.set_xlabel(r"$x$")
    axl.set_ylabel(r"$t$")
    axl.set_title(r"Peregrine soliton: rogue peak from nowhere")
    axr.plot(x, np.abs(_peregrine(x, np.array([0.0]))[0]), "C3", lw=1.8, label=r"$t=0$ (peak $=3$)")
    axr.axhline(1.0, color="0.5", ls=":", lw=1.0, label="background $=1$")
    axr.axhline(3.0, color="0.7", ls="--", lw=0.8)
    axr.set_xlabel(r"$x$")
    axr.set_ylabel(r"$|\psi(x,0)|$")
    axr.set_title(r"Amplification factor $3$")
    axr.legend(fontsize=8)
    axr.grid(alpha=0.3)
    _save(fig, "fig_peregrine")


def fig_mi() -> None:
    kappa = np.linspace(1e-3, np.pi - 1e-3, 500)
    lam_k = 4.0 * np.sin(kappa / 2) ** 2  # gamma=A=1: sigma^2 = lam_k(2 - lam_k)
    sigma = np.sqrt(np.clip(lam_k * (2.0 - lam_k), 0.0, None))
    cont = kappa * np.sqrt(np.clip(2.0 - kappa**2, 0.0, None))  # continuum k sqrt(2A^2 - k^2)
    sim = {20: 0.3436, 40: 0.6465, 60: 0.8766}  # measured slopes, verify_rogue_dnls.py (N=512)
    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    ax.plot(kappa, sigma, "k", lw=1.8, label=r"lattice $\sqrt{\lambda_\kappa(2-\lambda_\kappa)}$")
    ax.plot(kappa, cont, "C0:", lw=1.1, label=r"continuum $\kappa\sqrt{2-\kappa^2}$")
    ax.scatter([2 * np.pi * k / 512 for k in sim], list(sim.values()),
               marker="s", color="C3", zorder=5, label="DNLS simulation")
    ax.axvline(np.pi / 2, color="0.6", ls="--", lw=0.8)
    ax.annotate(r"band edge $\kappa=\pi/2$", xy=(np.pi / 2, 0.1), xytext=(np.pi / 2 + 0.1, 0.25),
                fontsize=8, color="0.4")
    ax.scatter([np.pi / 3], [1.0], color="k", zorder=6)
    ax.annotate(r"$\sigma_{\max}=\gamma A^2$", xy=(np.pi / 3, 1.0), xytext=(np.pi / 3 + 0.15, 0.85),
                fontsize=8)
    ax.set_xlabel(r"$\kappa$")
    ax.set_ylabel(r"growth rate $\sigma(\kappa)$")
    ax.set_title("Modulational (Benjamin-Feir) instability gain")
    ax.set_xlim(0, np.pi)
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(alpha=0.3)
    _save(fig, "fig_mi")


def _dnls_step(z: np.ndarray, dt: float, gamma: float, lam: np.ndarray) -> np.ndarray:
    z = z * np.exp(-1j * gamma * np.abs(z) ** 2 * dt)
    return np.fft.ifft(np.exp(1j * lam * dt) * np.fft.fft(z))


def fig_rogue_dnls() -> None:
    rng = np.random.default_rng(0)
    n = 256
    lam = 4.0 * np.sin(np.pi * np.arange(n) / n) ** 2
    z = (1.0 + 0.05 * rng.standard_normal(n)) * np.exp(1j * 0.05 * rng.standard_normal(n))
    dt, gamma = 0.004, 1.0
    snaps, times, amps = [], [], []
    for s in range(6000):
        z = _dnls_step(z, dt, gamma, lam)
        a = np.abs(z)
        amps.append(a.copy())
        if s % 30 == 0:
            snaps.append(a.copy())
            times.append(s * dt)
    field = np.array(snaps)
    amps = np.concatenate(amps)
    rms = float(np.sqrt(np.mean(amps**2)))

    fig, (axl, axr) = plt.subplots(1, 2, figsize=(8.8, 3.8), width_ratios=[1.15, 1])
    im = axl.pcolormesh(np.arange(n), times, field, shading="auto", cmap="magma", rasterized=True)
    fig.colorbar(im, ax=axl, label=r"$|z_j(t)|$")
    axl.set_xlabel(r"site $j$")
    axl.set_ylabel(r"time $t$")
    axl.set_title(r"DNLS on the ring ($N=256$): rogue spikes")

    xs = np.linspace(0, amps.max() / rms, 200)
    emp = np.array([np.mean(amps / rms > x) for x in xs])  # exceedance P(a/rms > x)
    rayleigh = np.exp(-(xs**2))  # Gaussian-field envelope: <a^2>=rms^2 -> P=exp(-x^2)
    axr.semilogy(xs, emp, "C3", lw=1.8, label="DNLS amplitudes")
    axr.semilogy(xs, rayleigh, "k--", lw=1.1, label=r"Rayleigh $e^{-x^2}$ (linear sea)")
    axr.axvline(2.2, color="0.6", ls=":", lw=1.0)
    axr.annotate("rogue\nthreshold", xy=(2.2, 1e-4), xytext=(1.4, 1e-5), fontsize=7, color="0.4")
    axr.set_xlabel(r"$|z|/\mathrm{rms}$")
    axr.set_ylabel(r"exceedance $P(|z|/\mathrm{rms}>x)$")
    axr.set_title("Heavy tail above the linear prediction")
    axr.set_ylim(1e-6, 1.5)
    axr.legend(fontsize=8)
    axr.grid(alpha=0.3, which="both")
    _save(fig, "fig_rogue_dnls")


def _ceiling_d(n: int, d: int) -> float:
    s = 2.0 * np.sin(np.pi * np.arange(n) / n)
    if d == 1:
        w2 = s**2
    elif d == 2:
        w2 = s[:, None] ** 2 + s[None, :] ** 2
    else:
        w2 = s[:, None, None] ** 2 + s[None, :, None] ** 2 + s[None, None, :] ** 2
    w = np.sqrt(w2)
    inv = np.zeros_like(w)
    inv[w > 1e-12] = 1.0 / w[w > 1e-12]
    return float(inv.sum() / n**d)


def fig_dimension() -> None:
    ns1 = np.array([8, 16, 32, 64, 128, 256, 512, 1024])
    ns = np.array([8, 16, 32, 64, 128])
    u1 = np.array([_ceiling_d(int(n), 1) for n in ns1])
    u2 = np.array([_ceiling_d(int(n), 2) for n in ns])
    u3 = np.array([_ceiling_d(int(n), 3) for n in ns])
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.plot(np.log(ns1), u1, "o-", color="C3", lw=1.3, label=r"$d=1$ (grows $\sim\frac{1}{\pi}\ln N$)")
    ax.plot(np.log(ns1), np.log(ns1) / np.pi + (np.euler_gamma + np.log(2 / np.pi)) / np.pi,
            "k--", lw=0.9, label=r"$\frac{1}{\pi}\ln N+C_0$")
    ax.plot(np.log(ns), u2, "s-", color="C0", lw=1.3, label=r"$d=2$ (bounded)")
    ax.plot(np.log(ns), u3, "^-", color="C2", lw=1.3, label=r"$d=3$ (bounded)")
    ax.set_xlabel(r"$\ln N$")
    ax.set_ylabel(r"ceiling $U_N^{(d)}$")
    ax.set_title(r"The large wave is one-dimensional: $U_N^{(d)}=\Theta(\ln N)$ only for $d=1$")
    ax.legend(fontsize=8, loc="center left")
    ax.grid(alpha=0.3)
    _save(fig, "fig_dimension")


def _fput_force(u: np.ndarray, beta: float) -> np.ndarray:
    dr = np.roll(u, -1) - u
    return (dr - np.roll(dr, 1)) + beta * (dr**3 - np.roll(dr**3, 1))


def _fput_varforce(u: np.ndarray, du: np.ndarray, beta: float) -> np.ndarray:
    dr = np.roll(u, -1) - u
    ddr = np.roll(du, -1) - du
    g = dr**2 * ddr
    return (ddr - np.roll(ddr, 1)) + 3.0 * beta * (g - np.roll(g, 1))


def fig_fput() -> None:
    n, beta, amp, dt = 32, 0.25, 4.0, 0.05
    o2 = 4.0 * np.sin(np.pi * np.arange(n) / n) ** 2
    j = np.arange(n)
    u = amp * np.sin(2 * np.pi * j / n)
    p = np.zeros(n)
    f = _fput_force(u, beta)
    ts, em = [], []
    for s in range(400000):
        p += 0.5 * dt * f
        u += dt * p
        f = _fput_force(u, beta)
        p += 0.5 * dt * f
        if s % 200 == 0:
            uh = np.fft.fft(u) / np.sqrt(n)
            ph = np.fft.fft(p) / np.sqrt(n)
            me = 0.5 * (np.abs(ph) ** 2 + o2 * np.abs(uh) ** 2)
            ts.append(s * dt)
            em.append([me[1] + me[n - 1], me[2] + me[n - 2], me[3] + me[n - 3]])
    ts, em = np.array(ts), np.array(em)
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    for k, c, lab in [(0, "C3", "mode 1"), (1, "C0", "mode 2"), (2, "C2", "mode 3")]:
        ax.plot(ts, em[:, k] / em[0, 0], color=c, lw=0.9, label=lab)
    ax.set_xlabel(r"time $t$")
    ax.set_ylabel(r"mode energy $E_r/E_1(0)$")
    ax.set_title(r"FPUT-$\beta$ recurrence ($N=32$): energy returns to mode 1")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    _save(fig, "fig_fput")


def fig_selftrapping() -> None:
    n, dt, steps = 128, 0.01, 40000
    lam = 4.0 * np.sin(np.pi * np.arange(n) / n) ** 2
    gammas = np.array([1, 2, 3, 4, 5, 6, 8, 12], dtype=float)
    pf = []
    for g in gammas:
        psi = np.zeros(n, dtype=complex)
        psi[0] = 1.0
        for _ in range(steps):
            psi = _dnls_step(psi, dt, float(g), lam)
        a2 = np.abs(psi) ** 2
        pf.append(a2.sum() ** 2 / np.sum(a2**2))
    pf = np.array(pf)
    psi = np.zeros(n, dtype=complex)
    psi[0] = 1.0
    snaps, tt = [], []
    for s in range(steps):
        psi = _dnls_step(psi, dt, 8.0, lam)
        if s % 400 == 0:
            snaps.append(np.abs(psi).copy())
            tt.append(s * dt)
    field = np.roll(np.array(snaps), n // 2, axis=1)  # center the breather for display

    fig, (axl, axr) = plt.subplots(1, 2, figsize=(8.8, 3.8), width_ratios=[1, 1.1])
    axl.plot(gammas, pf, "o-", color="C0")
    axl.axvline(4.0, ls="--", color="0.6", lw=0.8)
    axl.annotate(r"band top $\gamma\approx4$", (4.2, pf.max() * 0.6), fontsize=8, color="0.4")
    axl.set_xlabel(r"nonlinearity $\gamma$")
    axl.set_ylabel(r"participation ratio $P$")
    axl.set_title("Self-trapping transition")
    axl.grid(alpha=0.3)
    im = axr.pcolormesh(np.arange(n) - n // 2, tt, field, shading="auto", cmap="viridis", rasterized=True)
    fig.colorbar(im, ax=axr, label=r"$|\psi_j|$")
    axr.set_xlabel(r"site $j$")
    axr.set_ylabel(r"time $t$")
    axr.set_xlim(-20, 20)
    axr.set_title(r"Discrete breather ($\gamma=8$)")
    _save(fig, "fig_selftrapping")


def _lyap(n: int, beta: float, amp: float, dt: float, t_steps: int) -> tuple[float, float]:
    rng = np.random.default_rng(0)
    j = np.arange(n)
    u = amp * np.sin(2 * np.pi * j / n)
    p = np.zeros(n)
    dr = np.roll(u, -1) - u
    h0 = float(np.sum(p**2) / 2 + np.sum(dr**2) / 2 + beta * np.sum(dr**4) / 4)
    d0 = 1e-7
    du = rng.standard_normal(n)
    dp = rng.standard_normal(n)
    nrm = np.sqrt(np.sum(du**2) + np.sum(dp**2))
    du *= d0 / nrm
    dp *= d0 / nrm
    f, vf = _fput_force(u, beta), _fput_varforce(u, du, beta)
    re = max(1, round(1.0 / dt))
    ls, nc = 0.0, 0
    for s in range(1, t_steps + 1):
        p += 0.5 * dt * f
        dp += 0.5 * dt * vf
        u += dt * p
        du += dt * dp
        f, vf = _fput_force(u, beta), _fput_varforce(u, du, beta)
        p += 0.5 * dt * f
        dp += 0.5 * dt * vf
        if s % re == 0:
            d = np.sqrt(np.sum(du**2) + np.sum(dp**2))
            ls += np.log(d / d0)
            du *= d0 / d
            dp *= d0 / d
            nc += 1
    return ls / (nc * re * dt), h0


def fig_lyapunov() -> None:
    n, beta, dt = 32, 0.25, 0.02
    energies, lams = [], []
    for amp in (1, 2, 3, 4, 6, 8, 10):
        lam, h0 = _lyap(n, beta, float(amp), dt, 120000)
        lams.append(lam)
        energies.append(h0)
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.plot(energies, lams, "o-", color="C3")
    ax.axhline(0.0, color="0.6", lw=0.8)
    ax.set_xlabel(r"energy $H$")
    ax.set_ylabel(r"largest Lyapunov exponent $\lambda_{\max}$")
    ax.set_title(r"Route to chaos: $\lambda_{\max}$ turns on with energy")
    ax.grid(alpha=0.3)
    _save(fig, "fig_lyapunov")


def _dnls2d_step(psi: np.ndarray, dt: float, gamma: float, lam: np.ndarray) -> np.ndarray:
    psi = psi * np.exp(-1j * gamma * np.abs(psi) ** 2 * dt)
    return np.fft.ifft2(np.exp(1j * lam * dt) * np.fft.fft2(psi))


def fig_breather2d() -> None:
    n, dt, steps, gamma = 48, 0.01, 15000, 16.0
    d = 4.0 * np.sin(np.pi * np.arange(n) / n) ** 2
    lam = d[:, None] + d[None, :]
    psi = np.zeros((n, n), dtype=complex)
    psi[0, 0] = 1.0
    for _ in range(steps):
        psi = _dnls2d_step(psi, dt, gamma, lam)
    amp = np.roll(np.roll(np.abs(psi), n // 2, 0), n // 2, 1)  # center the breather
    w, c = 12, n // 2
    sub = amp[c - w:c + w + 1, c - w:c + w + 1]
    xx, yy = np.meshgrid(np.arange(-w, w + 1), np.arange(-w, w + 1))
    fig = plt.figure(figsize=(6.4, 4.7))
    ax = fig.add_subplot(projection="3d")
    ax.plot_surface(xx, yy, sub, cmap="viridis", linewidth=0, antialiased=True)
    ax.set_xlabel(r"site $x$")
    ax.set_ylabel(r"site $y$")
    ax.set_zlabel(r"$|\psi_{x,y}|$")
    ax.set_title(r"2D discrete breather (DNLS, $\gamma=16$, $N=48$)")
    _save(fig, "fig_breather2d")


def fig_phase3d() -> None:
    n, beta, dt = 32, 0.25, 0.05
    o2 = 4.0 * np.sin(np.pi * np.arange(n) / n) ** 2

    def run(amp: float, steps: int) -> np.ndarray:
        j = np.arange(n)
        u = amp * np.sin(2 * np.pi * j / n)
        p = np.zeros(n)
        f = _fput_force(u, beta)
        em = []
        for s in range(steps):
            p += 0.5 * dt * f
            u += dt * p
            f = _fput_force(u, beta)
            p += 0.5 * dt * f
            if s % 20 == 0:
                uh = np.fft.fft(u) / np.sqrt(n)
                ph = np.fft.fft(p) / np.sqrt(n)
                me = 0.5 * (np.abs(ph) ** 2 + o2 * np.abs(uh) ** 2)
                # mode-1 seed + cubic coupling excites only ODD modes (parity); use 1,3,5
                em.append([me[1] + me[n - 1], me[3] + me[n - 3], me[5] + me[n - 5]])
        return np.array(em)

    e_lo = run(4.0, 200000)   # regular (quasi-periodic torus)
    e_hi = run(10.0, 200000)  # chaotic (space-filling)
    fig = plt.figure(figsize=(8.8, 4.2))
    for idx, (em, title, c) in enumerate(
        [(e_lo, "regular (low energy)", "C0"), (e_hi, "chaotic (high energy)", "C3")]
    ):
        ax = fig.add_subplot(1, 2, idx + 1, projection="3d")
        ax.plot(em[:, 0], em[:, 1], em[:, 2], color=c, lw=0.3)
        ax.set_xlabel(r"$E_1$")
        ax.set_ylabel(r"$E_3$")
        ax.set_zlabel(r"$E_5$")
        ax.set_title(title)
    fig.suptitle("FPUT phase portrait in mode-energy space: torus vs chaos")
    _save(fig, "fig_phase3d")


def fig_dispersion3d() -> None:
    th = np.linspace(0, 2 * np.pi, 130)
    tx, ty = np.meshgrid(th, th)
    omega = 2.0 * np.sqrt(np.sin(tx / 2) ** 2 + np.sin(ty / 2) ** 2)
    fig = plt.figure(figsize=(6.6, 5.0))
    ax = fig.add_subplot(projection="3d")
    ax.plot_surface(tx, ty, omega, cmap="plasma", linewidth=0, antialiased=True, rcount=90, ccount=90)
    ax.set_xlabel(r"$\theta_x$")
    ax.set_ylabel(r"$\theta_y$")
    ax.set_zlabel(r"$\omega$")
    ax.set_title(r"Dispersion surface $\omega(\theta_x,\theta_y)$ on the 2-torus")
    _save(fig, "fig_dispersion3d")


def fig_caustic() -> None:
    fig, (axl, axr) = plt.subplots(1, 2, figsize=(8.8, 3.9))
    for t, c in [(50, "C0"), (200, "C1"), (800, "C3")]:
        nn = np.arange(int(2 * t) - 55, int(2 * t) + 15)
        s = 2.0 ** (1 / 3) * (nn - 2 * t) / (2 * t) ** (1 / 3)
        axl.plot(s, t ** (1 / 3) * jv(nn, 2 * t), color=c, lw=1.0, label=rf"$t={t}$")
    ss = np.linspace(-6, 3, 300)
    axl.plot(ss, airy(ss)[0], "k--", lw=1.3, label=r"$\mathrm{Ai}(s)$")
    axl.set_xlabel(r"$s=2^{1/3}(n-2t)/(2t)^{1/3}$")
    axl.set_ylabel(r"$t^{1/3}J_n(2t)$")
    axl.set_title("Discrete front collapses onto the Airy caustic")
    axl.legend(fontsize=8)
    axl.grid(alpha=0.3)
    ts = np.array([30, 60, 120, 240, 480, 960], dtype=float)
    peaks = np.array([np.abs(jv(np.arange(int(2 * t) - 55, int(2 * t) + 15), 2 * t)).max() for t in ts])
    axr.loglog(ts, peaks, "o", color="C3", label="front peak")
    axr.loglog(ts, 0.5357 * ts ** (-1 / 3), "k--", lw=1.0, label=r"$0.536\,t^{-1/3}$")
    axr.set_xlabel(r"$t$")
    axr.set_ylabel(r"$\max_n|J_n(2t)|$")
    axr.set_title(r"Fold-caustic scaling $\sim t^{-1/3}$")
    axr.legend(fontsize=8)
    axr.grid(alpha=0.3, which="both")
    _save(fig, "fig_caustic")


def fig_branched() -> None:
    rng = np.random.default_rng(0)
    n, vrms, corr, dz, zsteps = 1024, 0.12, 15.0, 0.05, 3000
    k = 2 * np.pi * np.fft.fftfreq(n)
    v = np.real(np.fft.ifft(np.fft.fft(rng.standard_normal(n)) * np.exp(-0.5 * (k * corr) ** 2)))
    v = vrms * v / np.std(v)
    kin = np.exp(-1j * (2 * np.pi * np.fft.fftfreq(n)) ** 2 * dz / 2)
    half = np.exp(-1j * v * dz / 2)
    psi = np.ones(n, dtype=complex)
    fields, zs, sci = [], [], []
    for s in range(zsteps):
        psi = half * psi
        psi = np.fft.ifft(kin * np.fft.fft(psi))
        psi = half * psi
        inten = np.abs(psi) ** 2
        sci.append(inten.max() / inten.mean())
        if s % 6 == 0:
            fields.append(inten.copy())
            zs.append(s * dz)
    field, zs, sci = np.array(fields), np.array(zs), np.array(sci)
    fig, (axl, axr) = plt.subplots(1, 2, figsize=(9.0, 4.0), width_ratios=[1.3, 1])
    im = axl.pcolormesh(np.arange(n), zs, np.log10(field + 1e-3), shading="auto",
                        cmap="inferno", vmin=-1, vmax=1.2, rasterized=True)
    fig.colorbar(im, ax=axl, label=r"$\log_{10}|\psi|^2$")
    axl.set_xlabel(r"transverse $x$")
    axl.set_ylabel(r"propagation $z$")
    axl.set_title("Branched flow: caustic hot spots")
    axr.plot(np.arange(zsteps) * dz, sci, color="C3", lw=1.0)
    axr.set_xlabel(r"propagation $z$")
    axr.set_ylabel(r"$\max_x|\psi|^2/\mathrm{mean}$")
    axr.set_title("Scintillation: focal hot-spot onset")
    axr.grid(alpha=0.3)
    _save(fig, "fig_branched")


def fig_cusp() -> None:
    th = np.linspace(0.06, np.pi - 0.06, 60)
    hitx, hity = np.cos(th), np.sin(th)
    nx, ny = np.cos(th), np.sin(th)
    dx, dy = 1 - 2 * nx * nx, -2 * nx * ny
    cx = 0.75 * np.cos(th) - 0.25 * np.cos(3 * th)
    cy = 0.75 * np.sin(th) - 0.25 * np.sin(3 * th)
    fig = plt.figure(figsize=(9.0, 4.2))
    axl = fig.add_subplot(1, 2, 1)
    circ = np.linspace(0, 2 * np.pi, 200)
    axl.plot(np.cos(circ), np.sin(circ), "0.6", lw=1.0)
    for i in range(0, len(th), 2):
        t = np.array([0.0, 2.2])
        axl.plot(hitx[i] + dx[i] * t, hity[i] + dy[i] * t, "C0", lw=0.3, alpha=0.6)
    axl.plot(cx, cy, "C3", lw=2.0, label="nephroid caustic")
    axl.plot(-cx, -cy, "C3", lw=2.0)
    axl.set_aspect("equal")
    axl.set_xlim(-1.6, 1.6)
    axl.set_ylim(-1.6, 1.6)
    axl.set_title("Teacup caustic (reflected rays $\\to$ nephroid)")
    axl.legend(fontsize=8, loc="lower right")
    xs = np.linspace(-4.5, 2.5, 90)
    ys = np.linspace(-4.0, 4.0, 110)
    s = np.arange(-7.0, 7.0, 0.004)
    win = np.exp(-(s / 6.4) ** 12)
    pe = np.empty((len(ys), len(xs)))
    for i, xv in enumerate(xs):
        integ = np.exp(1j * (s**4 + xv * s**2 + np.outer(ys, s))) * win
        pe[:, i] = np.abs(integ.sum(axis=1) * (s[1] - s[0])) ** 2
    axr = fig.add_subplot(1, 2, 2, projection="3d")
    xx, yy = np.meshgrid(xs, ys)
    axr.plot_surface(xx, yy, pe, cmap="inferno", linewidth=0, antialiased=True, rcount=90, ccount=80)
    axr.set_xlabel(r"$X$")
    axr.set_ylabel(r"$Y$")
    axr.set_zlabel(r"$|\mathrm{Pe}|^2$")
    axr.set_title("Pearcey cusp diffraction")
    _save(fig, "fig_cusp")


def fig_raycaustic() -> None:
    k, sigma, dz, zsteps = 0.04, 12.0, 0.1, 950
    x0 = np.linspace(-30, 30, 41)
    x = x0.copy()
    v = np.zeros_like(x0)
    acc = -k * x * np.exp(-x**2 / (2 * sigma**2))
    traj = [x.copy()]
    for _ in range(zsteps):
        v += 0.5 * dz * acc
        x = x + dz * v
        acc = -k * x * np.exp(-x**2 / (2 * sigma**2))
        v += 0.5 * dz * acc
        traj.append(x.copy())
    traj = np.array(traj)
    z = np.arange(zsteps + 1) * dz
    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    for j in range(traj.shape[1]):
        ax.plot(z, traj[:, j], "C0", lw=0.4, alpha=0.7)
    ax.set_xlabel(r"propagation $z$")
    ax.set_ylabel(r"transverse $x$")
    ax.set_title(r"Ray tracing in a current jet: caustic from a ray bundle")
    ax.grid(alpha=0.3)
    _save(fig, "fig_raycaustic")


def fig_caustic_mi() -> None:
    n = 4096
    xg = np.linspace(-80, 80, n, endpoint=False)
    k2 = (2 * np.pi * np.fft.fftfreq(n, d=(xg[1] - xg[0]))) ** 2
    amp, w, b, dt, steps = 0.42, 14.0, 0.28, 0.01, 2500

    def step(psi, gamma):
        psi = psi * np.exp(1j * gamma * np.abs(psi) ** 2 * dt / 2)
        psi = np.fft.ifft(np.exp(-1j * k2 * dt / 2) * np.fft.fft(psi))
        return psi * np.exp(1j * gamma * np.abs(psi) ** 2 * dt / 2)

    def run(psi0, gamma):
        psi = psi0.astype(complex).copy()
        peaks = [float(np.abs(psi).max())]
        best = np.abs(psi).copy()
        for _ in range(steps):
            psi = step(psi, gamma)
            peaks.append(float(np.abs(psi).max()))
            if np.abs(psi).max() > best.max():
                best = np.abs(psi).copy()
        return np.array(peaks), best

    lens = amp * np.exp(-xg**2 / (2 * w**2)) * np.exp(-1j * b * xg**2 / 2)
    plane = amp * np.ones(n, dtype=complex)
    p_lin, prof_lin = run(lens, 0.0)
    p_nl, prof_nl = run(lens, 1.0)
    p_pw, _ = run(plane, 1.0)
    zz = np.arange(steps + 1) * dt
    fig, (axl, axr) = plt.subplots(1, 2, figsize=(9.0, 4.0))
    axl.plot(zz, p_lin, "C0", lw=1.2, label="lens, linear (caustic)")
    axl.plot(zz, p_nl, "C3", lw=1.2, label="lens + NLS (self-focusing)")
    axl.plot(zz, p_pw, "0.6", lw=1.2, label="plane wave + NLS")
    axl.set_xlabel(r"propagation $z$")
    axl.set_ylabel(r"$\max_x|\psi|$")
    axl.set_title("Caustic seeds self-focusing")
    axl.legend(fontsize=8)
    axl.grid(alpha=0.3)
    axr.plot(xg, prof_lin, "C0", lw=1.0, label="caustic (linear)")
    axr.plot(xg, prof_nl, "C3", lw=1.0, label="caustic + self-focusing")
    axr.set_xlim(-12, 12)
    axr.set_xlabel(r"$x$")
    axr.set_ylabel(r"$|\psi|$ at focus")
    axr.set_title("Focal profile: nonlinearity sharpens the peak")
    axr.legend(fontsize=8)
    axr.grid(alpha=0.3)
    _save(fig, "fig_caustic_mi")


def _laplace(leader: np.ndarray) -> np.ndarray:
    ny, nx = leader.shape
    phi = np.zeros((ny, nx))
    phi[-1, :] = 1.0
    free = ~leader.copy()
    free[-1, :] = False
    free[0, :] = False
    idx = -np.ones((ny, nx), dtype=int)
    coords = np.argwhere(free)
    for m, (i, j) in enumerate(coords):
        idx[i, j] = m
    a = lil_matrix((len(coords), len(coords)))
    b = np.zeros(len(coords))
    for m, (i, j) in enumerate(coords):
        a[m, m] = -4.0
        for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            ii, jj = i + di, (j + dj) % nx
            if 0 <= ii < ny:
                if free[ii, jj]:
                    a[m, idx[ii, jj]] += 1.0
                else:
                    b[m] -= phi[ii, jj]
    sol = spsolve(a.tocsr(), b)
    for m, (i, j) in enumerate(coords):
        phi[i, j] = sol[m]
    return phi


def fig_lightning() -> None:
    rng = np.random.default_rng(0)
    ny, nx, steps = 60, 71, 700
    leader = np.zeros((ny, nx), dtype=bool)
    leader[0, nx // 2] = True
    for _ in range(steps):
        phi = _laplace(leader)
        cand = np.zeros((ny, nx), dtype=bool)
        for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            cand |= np.roll(leader, (di, dj), (0, 1))
        adj = cand & ~leader
        adj[0, :] = adj[-1, :] = False
        ci, cj = np.where(adj)
        if len(ci) == 0:
            break
        wgt = np.maximum(phi[ci, cj], 0.0)
        wgt = wgt if wgt.sum() > 0 else np.ones_like(wgt)
        p = rng.choice(len(ci), p=wgt / wgt.sum())
        leader[ci[p], cj[p]] = True
        if ci[p] == ny - 2:
            break
    field = _laplace(leader)
    field[leader] = np.nan
    fig, ax = plt.subplots(figsize=(5.6, 4.8))
    im = ax.pcolormesh(field, shading="auto", cmap="cividis", rasterized=True)
    fig.colorbar(im, ax=ax, label=r"potential $\phi$ (field $=\nabla\phi$)")
    yy, xx = np.where(leader)
    ax.scatter(xx + 0.5, yy + 0.5, c="white", s=2)
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"cloud $\to$ ground")
    ax.set_title("Lightning as Laplacian growth (discrete Laplacian)")
    ax.invert_yaxis()
    _save(fig, "fig_lightning")


def _u_response(n: int, t: np.ndarray) -> np.ndarray:
    r = np.arange(1, n)
    omega = 2.0 * np.sin(np.pi * r / n)
    coeff = 1.0 / (2 * n * np.sin(np.pi * r / n))
    return (np.sin(np.outer(t, omega)) * coeff).sum(axis=1)


def fig_reachability() -> None:
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(9.6, 4.0))

    primes = [5, 7, 11, 13, 17, 19]
    teps = []
    grid = np.arange(0.0, 20000.0, 0.05)
    for n in primes:
        rmax = np.maximum.accumulate(_u_response(n, grid))
        target = 0.8 * U(n)
        hit = np.argmax(rmax >= target)
        teps.append(grid[hit] if rmax[hit] >= target else np.nan)
    axL.semilogy(primes, teps, "o-", color="C0")
    axL.set_xlabel(r"$N$ (prime)")
    axL.set_ylabel(r"$T_{0.2}(N)$  (time to reach $0.8\,U_N$)")
    axL.set_title("Recurrence time grows exponentially in $N$")
    axL.grid(alpha=0.3, which="both")

    horizons = [(300.0, "C3", r"$T=300$"), (2000.0, "C1", r"$T=2000$"), (20000.0, "C2", r"$T=2\cdot10^4$")]
    ns = [5, 11, 17, 23, 29]
    for tmax, c, lab in horizons:
        g = np.arange(0.0, tmax, 0.05)
        eff = [float(np.maximum.accumulate(_u_response(n, g))[-1] / U(n)) for n in ns]
        axR.plot(ns, eff, "o-", color=c, label=lab)
    axR.axhline(1.0, ls=":", color="k", lw=0.8)
    axR.set_xlabel(r"$N$")
    axR.set_ylabel(r"finite-time efficiency $A_N(T)/U_N$")
    axR.set_title("At fixed horizon, efficiency falls with $N$")
    axR.legend(fontsize=8)
    axR.grid(alpha=0.3)
    _save(fig, "fig_reachability")


def fig_spectral_dim() -> None:
    def cyc(n: int) -> np.ndarray:
        return 2.0 - 2.0 * np.cos(2.0 * np.pi * np.arange(n) / n)

    def prod(*sizes: int) -> np.ndarray:
        e = np.zeros(1)
        for s in sizes:
            e = (e[:, None] + cyc(s)[None, :]).ravel()
        return e

    def ug(e: np.ndarray) -> float:
        pos = e[e > 1e-9]
        return float(np.sum(pos ** -0.5) / e.size)

    sizes = [50, 100, 200, 400, 800, 1600]
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    families = [
        ("ring $C_N$ ($d_s{=}1$)", "C0", "o", lambda n: cyc(n)),
        (r"ladder $C_N\times C_2$", "C1", "s", lambda n: prod(n, 2)),
        (r"tube $C_N\times C_3$", "C2", "^", lambda n: prod(n, 3)),
        (r"square torus ($d_s{=}2$)", "C3", "x", lambda n: prod(n, n)),
    ]
    lnx = np.log(sizes)
    for lab, c, mk, mk_eigs in families:
        vals = [ug(mk_eigs(n)) for n in sizes]
        ax.plot(lnx, vals, mk + "-", color=c, label=lab, ms=5)
    ax.set_xlabel(r"$\ln N$")
    ax.set_ylabel(r"ceiling $U_G=|V|^{-1}\sum\lambda^{-1/2}$")
    ax.set_title(r"Large wave at spectral dimension $d_s=1$ (incl.\ quasi-1D)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    _save(fig, "fig_spectral_dim")


def fig_filimonov() -> None:
    """Filimonov's discrete-Schrodinger large wave on the Dirichlet segment (alpha_j=1)."""
    def d_jN(n: int) -> np.ndarray:
        k = np.arange(1, n)
        a = (2.0 / n) * np.sin(np.outer(k, np.arange(1, n)) * np.pi / n).sum(axis=1)
        j = np.arange(1, n)
        return np.abs(a[None, :] * np.sin(np.outer(j, k) * np.pi / n)).sum(axis=1)

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(9.6, 4.0))

    # left: d_j^N vs ln j, the C = 4/pi^2 law
    n = 4096
    d = d_jN(n)
    j = np.arange(2, 400)
    c = 4.0 / np.pi**2
    intc = np.polyfit(np.log(j), d[j - 1], 1)[1]
    axL.plot(np.log(j), d[j - 1], color="C0", lw=1.3, label=r"$d_j^N$ ($N=2^{12}$)")
    axL.plot(np.log(j), c * np.log(j) + intc, "k--", lw=1.0,
             label=rf"$\frac{{4}}{{\pi^2}}\ln j + {intc:.2f}$")
    axL.set_xlabel(r"$\ln j$")
    axL.set_ylabel(r"site ceiling $d_j^N=\sum_k|a_k\sin\frac{\pi k j}{N}|$")
    axL.set_title(r"Theorem 2: $d_j^N\sim\frac{4}{\pi^2}\ln j$")
    axL.legend(fontsize=8, loc="upper left")
    axL.grid(alpha=0.3)

    # right: spatial profile |z_j(t*)| growing with j, under the d_j^N envelope (the large wave)
    n = 512
    k = np.arange(1, n)
    a = (2.0 / n) * np.sin(np.outer(k, np.arange(1, n)) * np.pi / n).sum(axis=1)
    w = np.cos(np.pi * k / n)
    dprof = d_jN(n)
    js = np.arange(1, n)
    sinjk = np.sin(np.outer(js, k) * np.pi / n)
    ts = np.linspace(0, 6000.0, 6000)
    best_t, best_amp = 0.0, 0.0
    for t in ts:                                   # pick a time with a large far-site amplitude
        amp = np.abs(sinjk[n // 2 - 1] @ (a * np.exp(-2j * w * t)))
        if amp > best_amp:
            best_amp, best_t = amp, t
    prof = np.abs(sinjk @ (a * np.exp(-2j * w * best_t)))
    axR.plot(js, prof, color="C3", lw=1.0, label=rf"$|z_j(t^*)|$, $t^*={best_t:.0f}$")
    axR.plot(js, dprof, color="0.5", lw=1.2, ls="--", label=r"ceiling $d_j^N$ (sup over $t$)")
    axR.set_xlabel(r"site $j$")
    axR.set_ylabel(r"$|z_j(t)|$")
    axR.set_title(r"Large wave: amplitude grows with site $j$ ($N=512$)")
    axR.legend(fontsize=8, loc="upper left")
    axR.grid(alpha=0.3)
    _save(fig, "fig_filimonov")


def fig_largewave() -> None:
    """Time-domain picture: the rare large-wave spike of G_N(0,t) standing above the thin-tailed bulk."""
    n = 13
    r = np.arange(1, n)
    om = 2.0 * np.sin(np.pi * r / n)
    coef = 1.0 / (n * om)
    Uc = float(np.sum(1.0 / np.sin(np.pi * r / n)) / (2 * n))
    rms = np.sqrt((n * n - 1) / (12.0 * n * n))
    ts = np.arange(0.0, 150000.0, 0.1)
    g = (np.sin(np.outer(ts, om)) * coef).sum(axis=1)
    ipk = int(np.argmax(g))
    tpk = ts[ipk]

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(9.6, 4.0))
    # left: a long stretch with the rare spike; sub-sample for plotting
    sl = slice(0, 200000)                                   # first 20000 time units
    axL.plot(ts[sl], g[sl], color="0.55", lw=0.4)
    axL.axhline(Uc, color="C3", ls="--", lw=1.0, label=rf"$U_{{{n}}}={Uc:.3f}$")
    axL.axhline(rms, color="C0", ls=":", lw=0.8, label=r"$\pm\,$RMS$=\pm0.29$")
    axL.axhline(-rms, color="C0", ls=":", lw=0.8)
    axL.set_xlabel(r"$t$")
    axL.set_ylabel(r"$G_N(0,t)$")
    axL.set_title(r"Thin-tailed bulk with rare spikes ($N=13$)")
    axL.legend(fontsize=8, loc="lower right")
    axL.grid(alpha=0.3)
    # right: zoom on the largest near-alignment peak -- the large wave forming
    w = (ts > tpk - 12) & (ts < tpk + 12)
    axR.plot(ts[w], g[w], color="C3", lw=1.2)
    axR.axhline(Uc, color="C3", ls="--", lw=1.0, label=rf"ceiling $U_{{{n}}}$")
    axR.axhline(rms, color="C0", ls=":", lw=0.8)
    axR.axhline(-rms, color="C0", ls=":", lw=0.8)
    axR.scatter([tpk], [g[ipk]], color="k", zorder=5, s=18)
    axR.annotate(rf"$\approx{g[ipk] / Uc:.2f}\,U_{{{n}}}$", (tpk, g[ipk]), textcoords="offset points",
                 xytext=(6, -4), fontsize=8)
    axR.set_xlabel(r"$t$")
    axR.set_ylabel(r"$G_N(0,t)$")
    axR.set_title(r"Zoom: the large wave at a near-alignment time")
    axR.legend(fontsize=8, loc="lower right")
    axR.grid(alpha=0.3)
    _save(fig, "fig_largewave")


def main() -> int:
    import os

    os.makedirs(OUT, exist_ok=True)
    for f in (fig_ring, fig_schrodinger, fig_bessel, fig_dispersion, fig_zeta,
              fig_dos, fig_value_dist, fig_heatmap, fig_qrank,
              fig_peregrine, fig_mi, fig_rogue_dnls, fig_dimension,
              fig_fput, fig_selftrapping, fig_lyapunov,
              fig_breather2d, fig_phase3d, fig_dispersion3d,
              fig_caustic, fig_branched,
              fig_cusp, fig_raycaustic, fig_caustic_mi, fig_lightning,
              fig_reachability, fig_spectral_dim, fig_filimonov, fig_largewave):
        f()
    print("wrote 29 figures to", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
