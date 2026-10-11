"""Figures for Chapter 5 -- From PID to FOPID."""
import json
import numpy as np
from style import plt, save, W, K, G1, G2, C, simlabel
from style import HERE
import vplant as v

nums = {}
# ---------------------------------------------------------------- P only
fig, ax = plt.subplots(figsize=(W, 2.1))
for Kp, ls in [(10, ":"), (25, "--"), (50, "-")]:
    o = v.simulate(v.FOPID(Kp=Kp, Ki=0, Kd=0), N=3001, noise=0.03, seed=1)
    ax.plot(o["t"], o["Tc"], ls=ls, color=K, label=fr"$K_p$ = {Kp}")
    off = o["Tc"][-200:].mean() - 18.2
    nums[f"P{Kp}"] = dict(offset=float(off), pwm=float(o["pwm"][-200:].mean()))
    ax.annotate(f"{off:.1f} °C", (300, 18.2 + off), (306, 18.2 + off), fontsize=6.5, va="center")
ax.axhline(18.2, color=C[1], lw=0.8, ls="-.", label=r"$T_{\mathrm{ref}}$")
ax.set_xlim(0, 330)
ax.set_ylim(17.5, 24)
ax.set_xlabel("time (s)")
ax.set_ylabel(r"$T_{\mathrm{sup}}$ (°C)")
ax.legend(loc="upper right", ncol=4)
simlabel(ax)
save(fig, "fig05_ponly")

# ---------------------------------------------------------------- windup
fig, ax = plt.subplots(2, 1, figsize=(W, 2.8), sharex=True)
for aw, ls, lab in [("none", "--", "no anti-windup"), ("ramp", "-", "conditional anti-windup")]:
    c = v.FOPID(Kp=50, Ki=1.6, Kd=0, antiwindup=aw)
    o = v.simulate(c, N=1201, Tref=12.0, noise=0.03, seed=1)
    m = v.metrics(o)
    nums[f"windup_{aw}"] = dict(undershoot=float(m["undershoot"]), settle=float(m["settle"]), IAE=float(m["IAE"]))
    ax[0].plot(o["t"], o["Tc"], ls=ls, color=K if aw == "ramp" else C[1], label=lab)
    ax[1].plot(o["t"], o["pwm"], ls=ls, color=K if aw == "ramp" else C[1])
ax[0].axhline(12.0, color=G2, lw=0.7, ls=":")
ax[0].set_ylabel(r"$T_{\mathrm{sup}}$ (°C)")
ax[0].legend(loc="upper right")
ax[1].set_ylabel("PWM")
ax[1].set_xlabel("time (s)")
simlabel(ax[1])
save(fig, "fig05_windup")

# ---------------------------------------------------------------- GL weights
M = 30
j = np.arange(M + 1)
lam, mu = v.TUNE_FOPID["lam"], v.TUNE_FOPID["mu"]
fig, ax = plt.subplots(1, 2, figsize=(W, 2.0))
wI_pid = v.gl_weights(0.0, M)
wI = v.gl_weights(1 - lam, M)
ax[0].stem(j + 0.15, wI_pid, linefmt="0.6", markerfmt="o", basefmt=" ", label=r"PID ($\lambda$ = 1)")
ax[0].stem(j - 0.15, wI, linefmt="k-", markerfmt="kD", basefmt=" ", label=fr"FOPID ($\lambda$ = {lam:.2f})")
ax[0].set_title(r"weights $g_j^{(1-\lambda)}$, integral path")
ax[0].set_xlabel(r"sample age $j$")
ax[0].set_xlim(-0.8, 12)
ax[0].legend(loc="upper right")
wD_pid = v.gl_weights(1.0, M)
wD = v.gl_weights(mu, M)
ax[1].stem(j + 0.15, wD_pid, linefmt="0.6", markerfmt="o", basefmt=" ", label=r"PID ($\mu$ = 1)")
ax[1].stem(j - 0.15, wD, linefmt="k-", markerfmt="kD", basefmt=" ", label=fr"FOPID ($\mu$ = {mu:.2f})")
ax[1].set_title(r"weights $g_j^{(\mu)}$, derivative path")
ax[1].set_xlabel(r"sample age $j$")
ax[1].set_xlim(-0.8, 12)
ax[1].legend(loc="lower right")
for a in ax:
    a.set_xticks(range(0, 13, 2))
    for line in a.get_lines():
        line.set_markersize(2.5)
fig.tight_layout()
save(fig, "fig05_glweights")
nums["wI"] = wI[:6].tolist()
nums["wD"] = wD[:6].tolist()

# ---------------------------------------------------------------- memory horizon
Jl = np.arange(1, 2001)
fig, ax = plt.subplots(figsize=(W, 2.0))
for a, ls, lab in [(1 - lam, "-", fr"integral path, order $1-\lambda$ = {1 - lam:.2f}"),
                   (mu, "--", fr"derivative path, order $\mu$ = {mu:.2f}")]:
    g = np.abs(v.gl_weights(a, 2000))[1:]
    ax.loglog(Jl * v.TS, g, ls=ls, color=K, label=lab)
ax.axvline(M * v.TS, color=C[1], lw=0.8)
ax.text(M * v.TS * 1.1, 2e-1, r"$M T_s$ = 3 s", fontsize=7, color=C[1])
ax.set_xlabel(r"age of the sample, $j\,T_s$ (s)")
ax.set_ylabel(r"$|g_j|$")
ax.legend(loc="lower left")
save(fig, "fig05_horizon")
gI = np.abs(v.gl_weights(1 - lam, 100000))
nums["tail_frac_I"] = float(gI[M + 1:].sum() / gI[1:].sum())
gD = np.abs(v.gl_weights(mu, 100000))
nums["tail_frac_D"] = float(gD[M + 1:].sum() / gD[1:].sum())

# ---------------------------------------------------------------- PID vs FOPID
fig, ax = plt.subplots(2, 1, figsize=(W, 2.8), sharex=True, gridspec_kw=dict(height_ratios=[1.6, 1]))
for name, x, ls, col in [("PID", v.TUNE_PID, "--", C[1]), ("FOPID", v.TUNE_FOPID, "-", K)]:
    o = v.simulate(v.FOPID(**x), N=1201, noise=0.03, seed=5)
    m = v.metrics(o)
    nums[name] = {k: float(val) for k, val in m.items()}
    ax[0].plot(o["t"], o["Tc"], ls=ls, color=col, label=name)
    ax[1].plot(o["t"], o["pwm"], ls=ls, color=col)
ax[0].axhline(18.2, color=G2, lw=0.7, ls=":")
ax[0].set_ylabel(r"$T_{\mathrm{sup}}$ (°C)")
ax[0].legend(loc="upper right")
ins = ax[0].inset_axes([0.33, 0.35, 0.4, 0.55])
for name, x, ls, col in [("PID", v.TUNE_PID, "--", C[1]), ("FOPID", v.TUNE_FOPID, "-", K)]:
    o = v.simulate(v.FOPID(**x), N=1201, noise=0.03, seed=5)
    ins.plot(o["t"], o["Tc"], ls=ls, color=col, lw=0.8)
ins.axhline(18.2, color=G2, lw=0.6, ls=":")
ins.set_xlim(5, 40)
ins.set_ylim(17.9, 18.8)
ins.tick_params(labelsize=6)
ax[1].set_ylabel("PWM")
ax[1].set_xlabel("time (s)")
simlabel(ax[1])
save(fig, "fig05_pidfopid")
json.dump(nums, open(HERE + "/ch05_numbers.json", "w"), indent=1)
print(json.dumps(nums, indent=1))
