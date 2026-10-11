"""Figures for Chapter 2 -- Instrumentation and Monitoring."""
import numpy as np
from style import plt, save, W, K, G1, G2, C, simlabel
import vplant as v

B, R0, T0, RF = 3435.0, 10e3, 298.15, 10e3


def r_ntc(Tc):
    return R0 * np.exp(B * (1 / (Tc + 273.15) - 1 / T0))


def counts(Tc, bits=12):
    """NTC on the high side, fixed resistor to ground: V rises with T."""
    return (2 ** bits - 1) * RF / (RF + r_ntc(Tc))


# ---------------------------------------------------------------- NTC curve
T = np.linspace(0, 50, 501)
fig, ax = plt.subplots(1, 2, figsize=(W, 1.9))
ax[0].plot(T, r_ntc(T) / 1e3, color=C[0])
ax[0].set_xlabel(r"temperature (°C)")
ax[0].set_ylabel(r"$R_{\mathrm{NTC}}$ (k$\Omega$)")
ax[0].plot([25], [10], "o", color=K, ms=3)
ax[0].annotate(r"10 k$\Omega$ at 25 °C", (25, 10), (28, 18), fontsize=7,
               arrowprops=dict(arrowstyle="-", lw=0.5))
V = 3.3 * RF / (RF + r_ntc(T))
ax[1].plot(T, V, color=C[0])
ax[1].plot(T, 3.3 / 2 + (T - 25) * (V[251] - V[249]) / 0.2, ls=":", color=G2, lw=0.8)
ax[1].set_xlabel(r"temperature (°C)")
ax[1].set_ylabel(r"$V_{\mathrm{out}}$ (V)")
ax[1].set_ylim(0, 3.3)
ax[1].text(2, 2.85, "tangent at 25 °C", fontsize=6.5, color=G1)
fig.tight_layout()
save(fig, "fig02_ntc_curve")

# ---------------------------------------------------------------- resolution
fig, ax = plt.subplots(figsize=(W, 1.8))
for bits, ls, lab in [(12, "-", "12-bit ADC (Nano 33 IoT)"), (10, "--", "10-bit ADC (Uno)")]:
    dNdT = np.gradient(counts(T, bits), T)
    ax.plot(T, 1 / dNdT, ls=ls, color=K, label=lab)
for Tp in (5, 25):
    for bits in (12, 10):
        d = 1 / np.interp(Tp, T, np.gradient(counts(T, bits), T))
        ax.plot(Tp, d, "o", color=G1, ms=2.5)
        ax.annotate(f"{d * 1000:.0f} mK", (Tp, d), (Tp + 1.2, d * 1.08), fontsize=6.5)
ax.set_xlabel(r"temperature (°C)")
ax.set_ylabel("°C per ADC count")
ax.set_ylim(0, 0.16)
ax.legend(loc="upper center")
save(fig, "fig02_resolution")

# ---------------------------------------------------------------- filter
rng = np.random.default_rng(3)
t = np.arange(0, 4.0, 0.1)
x = np.where(t >= 1.0, 21.0, 22.0) + rng.normal(0, 0.03, t.size)
x[25] += 1.2                      # isolated spike
f = v.MedEMA(x[0])
y = np.array([f(xi) for xi in x])
ema_only = np.zeros_like(x)
ema_only[0] = x[0]
for k in range(1, x.size):
    ema_only[k] = ema_only[k - 1] + 0.2 * (x[k] - ema_only[k - 1])
fig, ax = plt.subplots(figsize=(W, 1.9))
ax.plot(t, x, ".", color=G2, ms=3, label="raw samples")
ax.plot(t, ema_only, ls="--", color=C[1], label="EMA only")
ax.plot(t, y, color=K, label="median of 3 + EMA")
ax.axvline(1.0, color=G2, lw=0.6, ls=":")
ax.annotate("spike", (2.5, x[25]), (2.9, 22.5), fontsize=7, arrowprops=dict(arrowstyle="->", lw=0.5))
ax.set_xlabel("time (s)")
ax.set_ylabel(r"temperature (°C)")
ax.legend(loc="upper right", ncol=3, bbox_to_anchor=(1, 1.18))
save(fig, "fig02_filter")

# ---------------------------------------------------------------- dew point
Ta = np.linspace(10, 35, 200)
fig, ax = plt.subplots(figsize=(W, 2.1))
for rh, ls in [(30, ":"), (50, "-"), (70, "--"), (90, "-.")]:
    ax.plot(Ta, v.magnus_dew(Ta, rh), ls=ls, color=K, label=f"RH = {rh} %")
ax.plot(Ta, Ta, color=G2, lw=0.6)
ax.text(30.5, 31.6, r"$T_{\mathrm{dew}}=T_{\mathrm{amb}}$ (saturated)", fontsize=6.5, color=G1, ha="right")
ax.plot(25, v.magnus_dew(25, 50), "o", color=C[1], ms=3.5)
ax.annotate(r"25 °C, 50 % → 13.9 °C", (25, 13.87), (26.5, 3.5), fontsize=7,
            arrowprops=dict(arrowstyle="->", lw=0.5))
ax.set_xlabel(r"ambient temperature $T_{\mathrm{amb}}$ (°C)")
ax.set_ylabel(r"dew point $T_{\mathrm{dew}}$ (°C)")
ax.legend(loc="upper left")
save(fig, "fig02_dewpoint")

# ---------------------------------------------------------------- monitor view
o = v.simulate(v.FOPID(**v.TUNE_FOPID), N=1801, noise=0.03, seed=2)
fig, ax = plt.subplots(2, 1, figsize=(W, 2.6), sharex=True, gridspec_kw=dict(height_ratios=[1.4, 1]))
ax[0].plot(o["t"], o["Tsup"], color=K, label=r"$T_{\mathrm{sup}}$ (filtered)")
ax[0].plot(o["t"], o["Tref"], ls="--", color=C[1], label=r"$T_{\mathrm{ref}}$")
ax[0].set_ylabel("°C")
ax[0].legend(loc="upper right")
ax[1].plot(o["t"], o["e"], color=K)
ax[1].axhline(0, color=G2, lw=0.6)
ax[1].fill_between(o["t"], 0, o["e"], color="0.85")
ax[1].set_ylabel(r"$e=T_{\mathrm{ref}}-T_{\mathrm{sup}}$ (°C)")
ax[1].set_xlabel("time (s)")
iae = np.sum(np.abs(o["e"])) * 0.1
ax[1].text(100, -3.5, f"IAE = shaded area = {iae:.0f} °C·s", fontsize=7)
simlabel(ax[1])
save(fig, "fig02_monitor")
print("IAE monitor", iae)
print("counts/K at 25 (12 bit)", np.gradient(counts(T), T)[250])
for Tp in (5, 25):
    for bits in (12, 10):
        print(Tp, bits, 1 / np.interp(Tp, T, np.gradient(counts(T, bits), T)))
