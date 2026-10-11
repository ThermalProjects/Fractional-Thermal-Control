"""Figures and numbers for Chapter 4 -- A Thermal Model of the Plant."""
import json
import numpy as np
from numba import njit
from scipy.optimize import least_squares
from scipy.special import gamma as Gamma
from style import plt, save, W, K, G1, G2, C, simlabel
from style import HERE
import vplant as v

TS = v.TS
levels = [0, 25, 50, 75, 100, 75, 50, 25, 0]
pwm = np.concatenate([np.zeros(200), np.repeat(levels[1:], 2000), np.zeros(1000)]).astype(float)
rep1 = v.open_loop(pwm, Tamb=23.7, noise=0.03, seed=1)
rep2 = v.open_loop(pwm, Tamb=23.9, noise=0.03, seed=2)
t = rep1["t"]


@njit(cache=False)
def ol3(u, Tamb, R1, R2, Rc, Cc, Cm, Ch, phi, td, kP):
    N = u.shape[0]
    d = max(1, int(round(td / 0.1)))
    Tc = Tamb; Tm = Tamb; Th = Tamb
    oc = np.empty(N); oh = np.empty(N)
    for k in range(N):
        j = k - d
        uu = u[j] if j >= 0 else u[0]
        Pe = kP * uu
        dTc = (Tm - Tc) / (R1 * Cc) - phi * Pe / Cc
        dTm = (Tc - Tm) / (R1 * Cm) + (Th - Tm) / (R2 * Cm)
        dTh = (Tm - Th) / (R2 * Ch) + (1 - phi) * Pe / Ch - (Th - Tamb) / (Rc * Ch)
        Tc += 0.1 * dTc; Tm += 0.1 * dTm; Th += 0.1 * dTh
        oc[k] = Tc; oh[k] = Th
    return oc, oh


NAMES = ["R1", "R2", "Rconv", "Cc", "Cm", "Ch", "phi", "td"]
TRUE = np.array([v.P[n if n != "Cm" else "Cm"] for n in ["R1", "R2", "Rconv", "Cc", "Cm", "Ch", "phi", "td"]])


def model(x, Tamb):
    return ol3(pwm, Tamb, *x, 0.12)


def res(x, data, Tamb):
    oc, oh = model(x, Tamb)
    return np.concatenate([(oc - data["Tc_meas"])[::5], (oh - data["Th_meas"])[::5]])


lb = [0.1, 0.1, 0.1, 0.5, 10, 2, 0.05, 0.1]
ub = [50, 50, 50, 200, 5000, 1000, 0.95, 3.0]
fits = []
rng = np.random.default_rng(7)
for s in range(6):
    x0 = TRUE * np.exp(rng.normal(0, 0.5, TRUE.size))
    x0[6] = np.clip(x0[6], 0.1, 0.8)
    x0 = np.clip(x0, np.array(lb) * 1.01, np.array(ub) * 0.99)
    sol = least_squares(res, x0, args=(rep1, 23.7), bounds=(lb, ub), x_scale=TRUE, max_nfev=400)
    oc, oh = model(sol.x, 23.7)
    rm1 = np.sqrt(np.mean((oc - rep1["Tc_meas"]) ** 2))
    oc2, oh2 = model(sol.x, 23.9)
    rm2 = np.sqrt(np.mean((oc2 - rep2["Tc_meas"]) ** 2))
    fits.append(dict(x=sol.x.tolist(), rmse_fit=rm1, rmse_val=rm2,
                     rmse_fit_h=np.sqrt(np.mean((oh - rep1["Th_meas"]) ** 2))))
    print(s, np.round(sol.x, 3), round(rm1, 4), round(rm2, 4))
json.dump(fits, open(HERE + "/ch04_fits.json", "w"), indent=1)
best = min(fits, key=lambda f: f["rmse_fit"])
xb = np.array(best["x"])

# one-node model fitted to the first step (for contrast)
step = (t >= 20) & (t < 220)
y = rep1["Tc_meas"][step]
tt = t[step] - 20


def one(p, tt):
    K1, tau, L = p
    return 23.7 + K1 * 25 * (1 - np.exp(-np.clip(tt - L, 0, None) / tau))


s1 = least_squares(lambda p: one(p, tt) - y, [-0.15, 20, 1], bounds=([-1, 0.5, 0], [0, 500, 5]))
K1, tau1, L1 = s1.x
print("one node:", s1.x)

# ---------------------------------------------------------------- tau reading
fig, ax = plt.subplots(figsize=(W, 2.1))
ax.plot(tt, y, ".", ms=1.5, color=G2, label="cold face, first step (0 → 25)")
ax.plot(tt, one(s1.x, tt), color=C[1], ls="--", label=fr"one-node fit, $\tau$ = {tau1:.0f} s")
Tinf = 23.7 + K1 * 25
T63 = 23.7 + 0.632 * (Tinf - 23.7)
ax.axhline(T63, color=G1, lw=0.6, ls=":")
ax.axvline(L1 + tau1, color=G1, lw=0.6, ls=":")
ax.annotate("63 %", (L1 + tau1, T63), (L1 + tau1 + 12, T63 + 0.6), fontsize=7)
ax.set_xlabel("time since the step (s)")
ax.set_ylabel(r"$T_{\mathrm{sup}}$ (°C)")
ax.legend(loc="upper right")
simlabel(ax)
save(fig, "fig04_tau")

# ---------------------------------------------------------------- Peltier balance (reference module, from its datasheet)
# ETH-241-10-13-E-H1: Imax 3.5 A, Vmax 29.8 V, dTmax 73 K at Th = 298 K
Imax, Vmax, dTmax, Thd = 3.5, 29.8, 73.0, 298.15
S = Vmax / Thd                                   # V/K
Rm = Vmax / Imax * (Thd - dTmax) / Thd           # ohm
Km = Imax * Vmax * (Thd - dTmax) / (2 * Thd * dTmax)   # W/K
print("module from datasheet: S=%.4f V/K  R=%.2f ohm  K=%.3f W/K" % (S, Rm, Km))
I = np.linspace(0, Imax, 300)
Tcold = 288.15
fig, ax = plt.subplots(figsize=(W, 2.2))
for dT, ls in [(0, "-"), (20, "--"), (40, ":")]:
    Qc = S * I * Tcold - 0.5 * I ** 2 * Rm - Km * dT
    ax.plot(I, Qc, ls=ls, color=K, label=fr"$\Delta T$ = {dT} K")
ax.plot(I, S * I * Tcold, color=C[0], lw=0.7, ls="-.", label="Peltier term only")
ax.axvspan(0, 12 / Rm, color="0.9")
ax.text(0.05, 72, "12 V supply", fontsize=6.5, color=G1)
ax.set_xlabel("current I (A)")
ax.set_ylabel(r"heat pumped from cold face $Q_c$ (W)")
ax.set_ylim(-25, 105)
ax.legend(loc="upper left", ncol=2, bbox_to_anchor=(0.18, 1.0))
save(fig, "fig04_balance")

# ---------------------------------------------------------------- one-node global fit
@njit(cache=False)
def ol1(u, Tamb, tau, Kg, td):
    N = u.shape[0]; d = max(1, int(round(td / 0.1))); T = Tamb
    o = np.empty(N)
    for k in range(N):
        j = k - d
        uu = u[j] if j >= 0 else u[0]
        T += 0.1 * ((Tamb - T) + Kg * uu) / tau
        o[k] = T
    return o
g1 = least_squares(lambda p: (ol1(pwm, 23.7, *p) - rep1["Tc_meas"])[::5], [40, -0.15, 1.0], bounds=([1, -1, 0.1], [2000, 0, 5]))
one_glob = ol1(pwm, 23.7, *g1.x)
print("one-node global", g1.x, np.sqrt(np.mean((one_glob - rep1["Tc_meas"]) ** 2)))
RM1G = float(np.sqrt(np.mean((one_glob - rep1["Tc_meas"]) ** 2)))

# ---------------------------------------------------------------- fit and validation
oc1, oh1 = model(xb, 23.7)
oc2, oh2 = model(xb, 23.9)
fig, ax = plt.subplots(3, 1, figsize=(W, 3.6), sharex=True, gridspec_kw=dict(height_ratios=[2, 1, 2]))
ax[0].plot(t, rep1["Tc_meas"], color=G2, lw=0.6, label="run 1 (fitting)")
ax[0].plot(t, oc1, color=K, lw=0.9, ls="--", label="three-node fit")
ax[0].plot(t, rep1["Th_meas"], color=G2, lw=0.6)
ax[0].plot(t, oh1, color=K, lw=0.9, ls="--")
ax[0].set_ylabel("°C")
ax[0].legend(loc="upper right", ncol=2)
ax[1].plot(t, rep1["Tc_meas"] - one_glob, color=C[1], lw=0.5, label="one node")
ax[1].plot(t, rep1["Tc_meas"] - oc1, color=K, lw=0.4, label="three nodes")
ax[1].set_ylabel("residual (°C)")
ax[1].legend(loc="upper right", ncol=2)
ax[2].plot(t, rep2["Tc_meas"], color=G2, lw=0.6, label="run 2 (validation)")
ax[2].plot(t, oc2, color=C[1], lw=0.9, ls="--", label="prediction, no refit")
ax[2].set_ylabel(r"$T_{\mathrm{sup}}$ (°C)")
ax[2].set_xlabel("time (s)")
ax[2].legend(loc="lower right", ncol=2)
simlabel(ax[2])
save(fig, "fig04_fit")

# ---------------------------------------------------------------- memory kernels
tk = np.linspace(0.01, 60, 3000)
tau = 10.0
beta = 0.69


def mittag_leffler_relax(beta, tau, t, h=0.01):
    """x(t) solving D^beta x = -x / tau^beta, x(0)=1, by Grunwald-Letnikov (full memory)."""
    n = int(t[-1] / h) + 1
    g = v.gl_weights(beta, n)
    x = np.zeros(n)
    x[0] = 1.0
    a = h ** beta / tau ** beta
    for k in range(1, n):
        # sum_j g_j (x_{k-j} - x0) = -a x_k
        s = np.dot(g[1:k + 1], x[k - 1::-1] - 1.0)
        x[k] = (1.0 - s) / (1.0 + a)
    tt_ = np.arange(n) * h
    return np.interp(t, tt_, x)


ml = mittag_leffler_relax(beta, tau, tk, h=0.05)
fig, ax = plt.subplots(1, 2, figsize=(W, 1.9))
ax[0].plot(tk, np.exp(-tk / tau), color=K, label=r"exponential ($\beta$ = 1)")
ax[0].plot(tk, ml, color=C[1], ls="--", label=fr"power-law tail ($\beta$ = {beta})")
ax[0].set_xlabel("time (s)")
ax[0].set_ylabel("normalised response")
ax[0].legend(loc="upper right")
ax[1].loglog(tk, np.exp(-tk / tau), color=K)
ax[1].loglog(tk, ml, color=C[1], ls="--")
ax[1].set_ylim(1e-3, 1.2)
ax[1].set_xlabel("time (s)")
fig.tight_layout()
save(fig, "fig04_memory")

out = dict(one_global=g1.x.tolist(), rmse_one_global=RM1G, tau1=tau1, K1=K1, L1=L1, best=best,
           rmse_one=float(np.sqrt(np.mean((one(s1.x, tt) - y) ** 2))))
json.dump(out, open(HERE + "/ch04_numbers.json", "w"), indent=1, default=float)
print(out)
