"""Figures for Chapter 3 -- Building the Plant and Closing the First Loop."""
import numpy as np
from style import plt, save, W, K, G1, G2, G3, C, simlabel
import vplant as v

# ---------------------------------------------------------------- duty cycle
fig, ax = plt.subplots(3, 1, figsize=(W, 2.3), sharex=True)
t = np.linspace(0, 3, 3000)
for a, d in zip(ax, [0.25, 0.5, 0.75]):
    sq = 12.0 * ((t % 1.0) < d)
    a.plot(t, sq, color=K, lw=0.8)
    a.axhline(12 * d, color=C[1], ls="--", lw=0.9)
    a.set_ylim(-1, 14)
    a.set_yticks([0, 12])
    a.text(3.12, 12 * d, f"PWM = {int(round(255 * d))}\nmean {12 * d:.0f} V", fontsize=6.5, va="center")
    a.grid(False)
ax[1].set_ylabel("voltage across the cell (V)")
ax[-1].set_xlabel("time (PWM periods)")
fig.subplots_adjust(right=0.78, hspace=0.25)
save(fig, "fig03_duty")

# ---------------------------------------------------------------- staircase
levels = [0, 25, 50, 75, 100, 75, 50, 25, 0]
dwell = 200.0
N = int(len(levels) * dwell / v.TS)
pwm = np.concatenate([np.zeros(200), np.repeat(levels[1:], int(dwell / v.TS)), np.zeros(1000)]).astype(float)
reps = []
for r, (Ta, sd) in enumerate([(23.7, 1), (23.9, 2), (23.5, 3)]):
    reps.append(v.open_loop(pwm, Tamb=Ta, noise=0.03, seed=sd))
o = reps[0]
fig, ax = plt.subplots(2, 1, figsize=(W, 2.9), sharex=True, gridspec_kw=dict(height_ratios=[1, 2.2]))
ax[0].step(o["t"], pwm, where="post", color=K)
ax[0].set_ylabel("PWM")
ax[0].set_ylim(-5, 115)
ax[1].plot(o["t"], o["Th_meas"], color=C[1], lw=0.8, label=r"hot face $T_h$")
ax[1].plot(o["t"], o["Tc_meas"], color=C[0], lw=0.8, label=r"cold face $T_{\mathrm{sup}}$")
ax[1].axhline(23.7, color=G2, ls=":", lw=0.8)
ax[1].text(5, 24.6, r"$T_{\mathrm{amb}}$", fontsize=7, color=G1)
ax[1].set_ylabel("°C")
ax[1].set_xlabel("time (s)")
ax[1].legend(loc="upper right", ncol=2)
# annotate fast and slow parts of the first step
ax[1].annotate("most of the change within a minute,\nthen a slow drift", (40, o["Tc"][400]), (330, 27.0), fontsize=6.5,
               arrowprops=dict(arrowstyle="->", lw=0.5))
simlabel(ax[1])
save(fig, "fig03_staircase")

# ---------------------------------------------------------------- repeatability
Tc = np.vstack([r["Tc_meas"] for r in reps])
fig, ax = plt.subplots(2, 1, figsize=(W, 2.5), sharex=True, gridspec_kw=dict(height_ratios=[2, 1]))
for i, ls in enumerate(["-", "--", ":"]):
    ax[0].plot(o["t"], Tc[i], ls=ls, color=K, lw=0.7, label=f"run {i + 1}")
ax[0].set_ylabel(r"$T_{\mathrm{sup}}$ (°C)")
ax[0].legend(loc="lower left", ncol=3)
sig = Tc.std(axis=0, ddof=1)
ax[1].plot(o["t"], sig, color=K, lw=0.7)
ax[1].set_ylabel(r"$\sigma(t)$ (°C)")
ax[1].set_xlabel("time (s)")
simlabel(ax[1])
save(fig, "fig03_repeat")
print("mean sigma", sig.mean(), "max", sig.max())

# ---------------------------------------------------------------- hysteresis
h = 0.3
fig = plt.figure(figsize=(W, 2.4))
gs = fig.add_gridspec(2, 3, width_ratios=[1, 0.08, 2.4], height_ratios=[1.6, 1])
a0 = fig.add_subplot(gs[:, 0])
a0.plot([-1.2, -h], [255, 255], color=K)
a0.plot([h, 1.2], [0, 0], color=K)
a0.plot([-h, h], [0, 0], color=K, ls="--")
a0.plot([-h, h], [255, 255], color=K, ls="--")
a0.annotate("", (-h, 0), (-h, 255), arrowprops=dict(arrowstyle="->", lw=0.7))
a0.annotate("", (h, 255), (h, 0), arrowprops=dict(arrowstyle="->", lw=0.7))
a0.set_xticks([-h, 0, h])
a0.set_xticklabels([r"$-h$", "0", r"$+h$"])
a0.set_xlabel(r"$e$ (°C)")
a0.set_ylabel("PWM")
a0.set_ylim(-20, 280)
a0.set_title("relay with dead band", fontsize=7)
c = v.Hysteresis(h=h, on=255)
o = v.simulate(c, N=3001, noise=0.03, seed=4)
a1 = fig.add_subplot(gs[0, 2])
a1.plot(o["t"], o["Tc"], color=K, lw=0.8, label=r"$T_{\mathrm{sup}}$")
a1.axhline(18.2, color=C[1], ls="--", lw=0.8, label=r"$T_{\mathrm{ref}}$")
a1.axhspan(18.2 - h, 18.2 + h, color="0.88")
a1.set_ylabel("°C")
a1.set_ylim(14.5, 24)
a1.legend(loc="upper right", ncol=2)
a2 = fig.add_subplot(gs[1, 2], sharex=a1)
a2.plot(o["t"], o["pwm"], color=K, lw=0.7)
a2.set_ylabel("PWM")
a2.set_xlabel("time (s)")
plt.setp(a1.get_xticklabels(), visible=False)
simlabel(a2)
fig.subplots_adjust(hspace=0.12, wspace=0.05)
save(fig, "fig03_hysteresis")
m = o["t"] > 100
e = o["e"][m]
print("hyst mean|e|", np.mean(np.abs(e)), "pp", np.ptp(o["Tc"][m]), "duty", np.mean(o["pwm"][m]) / 2.55)
sw = np.where(np.diff((o["pwm"] > 0).astype(int)) > 0)[0]
print("period", np.diff(o["t"][sw])[2:].mean())
