"""Figures and numbers for Chapter 10 -- Dew-Point Control: Making Water from Air."""
import json
import numpy as np
from style import plt, save, W, K, G1, G2, C, simlabel
from style import HERE
import vplant as v

nums = {}
F = v.TUNE_FOPID
TS = v.TS
RV = 461.5                     # J/(kg K), water vapour
AREA = 0.040 * 0.040           # m^2, cold plate (assumed 40 mm x 40 mm)
HM = 0.006                     # m/s, mass-transfer coefficient, natural convection (indicative)


def rho_v(p_hpa, T):
    return p_hpa * 100 / (RV * (T + 273.15))


def cond_rate(Ts, Tamb, RH):
    """Condensation rate on the plate, g/h (negative = evaporation potential)."""
    pv = RH / 100 * v.p_sat(Tamb)
    Tfilm = 0.5 * (Ts + Tamb)
    return HM * AREA * (pv - v.p_sat(Ts)) * 100 / (RV * (Tfilm + 273.15)) * 3.6e6


# ---------------------------------------------------------------- vapour pressure and condensation zone
T = np.linspace(0, 35, 300)
fig, ax = plt.subplots(figsize=(W, 2.2))
ax.plot(T, v.p_sat(T), color=K, lw=1.2, label="saturation pressure $p_s(T)$")
Ta, RHa = 24.0, 65.0
pv = RHa / 100 * v.p_sat(Ta)
Td = float(v.magnus_dew(Ta, RHa))
ax.axhline(pv, color=C[0], lw=0.8, ls="--", label=fr"vapour pressure of the air, $T_{{amb}}$ = {Ta:.0f} °C, RH = {RHa:.0f} %")
ax.fill_between(T, v.p_sat(T), pv, where=T < Td, color="0.85", lw=0, label="surface condenses")
ax.plot([Ta], [pv], "o", color=C[0], ms=4)
ax.plot([Td], [pv], "s", color=K, ms=4)
ax.annotate(r"$T_{\mathrm{dew}}$", (Td, pv), xytext=(Td - 7, pv + 9), arrowprops=dict(arrowstyle="->", lw=0.6), fontsize=7)
ax.annotate("air", (Ta, pv), xytext=(Ta + 2, pv - 8), arrowprops=dict(arrowstyle="->", lw=0.6), fontsize=7)
Tsurf = Td - 2.0
ax.annotate("", xy=(Tsurf, v.p_sat(Tsurf)), xytext=(Tsurf, pv), arrowprops=dict(arrowstyle="<->", lw=0.8, color=C[1]))
ax.text(Tsurf - 0.5, (pv + v.p_sat(Tsurf)) / 2, r"$\Delta p$", ha="right", va="center", fontsize=7, color=C[1])
ax.set_xlabel("temperature (°C)")
ax.set_ylabel("pressure (hPa)")
ax.set_xlim(0, 35); ax.set_ylim(0, 50)
ax.legend(loc="upper left", fontsize=6)
save(fig, "fig10_vapour")
nums["example"] = dict(Tamb=Ta, RH=RHa, Tdew=Td, pv=float(pv))

# ---------------------------------------------------------------- uncertainty of T_dew
def dTd(Tamb, RH, dT, dRH):
    gT = (v.magnus_dew(Tamb + 1e-3, RH) - v.magnus_dew(Tamb - 1e-3, RH)) / 2e-3
    gR = (v.magnus_dew(Tamb, RH + 1e-3) - v.magnus_dew(Tamb, RH - 1e-3)) / 2e-3
    return np.sqrt((gT * dT) ** 2 + (gR * dRH) ** 2), gT, gR


RHs = np.linspace(30, 90, 121)
unc = {}
fig, ax = plt.subplots(figsize=(W, 1.9))
for dT_, dRH_, ls in [(0.5, 3.0, "-"), (0.2, 2.0, "--"), (0.1, 1.5, ":")]:
    s = np.array([dTd(24.0, r, dT_, dRH_)[0] for r in RHs])
    ax.plot(RHs, s, ls=ls, color=K, label=fr"$\pm${dT_} °C, $\pm${dRH_:.1f} % RH")
    unc[f"{dT_},{dRH_}"] = dict(at65=float(dTd(24.0, 65.0, dT_, dRH_)[0]), at40=float(dTd(24.0, 40.0, dT_, dRH_)[0]),
                                at85=float(dTd(24.0, 85.0, dT_, dRH_)[0]))
_, gT, gR = dTd(24.0, 65.0, 0, 0)
nums["sens"] = dict(dTd_dT=float(gT), dTd_dRH=float(gR))
ax.set_xlabel("relative humidity (%)")
ax.set_ylabel(r"uncertainty of $T_{\mathrm{dew}}$ (°C)")
ax.legend(title=r"sensor tolerance at $T_{amb}$ = 24 °C", fontsize=6, title_fontsize=6)
ax.set_ylim(0, None)
save(fig, "fig10_uncertainty")
nums["unc"] = unc

# ---------------------------------------------------------------- closed loop with a moving reference
def run_dew(Tamb, RH, b, N, mode="moving", tau_ref=20.0, fixed=None, rh_noise=0.3, t_noise=0.03,
            noise=0.03, seed=0, ctrl=None):
    """Closed loop with the ambient sensor read every 2 s, Magnus, EMA on T_dew and
    T_ref = T_dew - b (moving) or a fixed reference."""
    rng = np.random.default_rng(seed)
    c = ctrl if ctrl is not None else v.FOPID(**F)
    pl = v.Plant(Tamb[0]); fl = v.MedEMA(pl.Tc)
    Tdew_f = float(v.magnus_dew(Tamb[0], RH[0]))
    a_ref = 2.0 / tau_ref if tau_ref else 1.0
    out = {k: np.zeros(N) for k in ["Tc", "Tref", "Tdew_true", "Tdew_meas", "pwm"]}
    Tdew_m = Tdew_f
    for k in range(N):
        if k % 20 == 0:
            Tdew_m = float(v.magnus_dew(Tamb[k] + rng.normal(0, t_noise), RH[k] + rng.normal(0, rh_noise)))
            Tdew_f += a_ref * (Tdew_m - Tdew_f)
        Tref = (Tdew_f - b) if mode == "moving" else fixed
        Tf = fl(pl.Tc + rng.normal(0, noise))
        u = c.step(Tref - Tf, Tref)
        pl.step(u, Tamb[k])
        out["Tc"][k] = pl.Tc; out["Tref"][k] = Tref; out["pwm"][k] = u
        out["Tdew_true"][k] = v.magnus_dew(Tamb[k], RH[k]); out["Tdew_meas"][k] = Tdew_m
    return out


# design curve: steady condensation, drive and margin as a function of b_dew
N = 6001
bs = [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0]
Tamb_c = np.full(N, 24.0); RH_c = np.full(N, 65.0)
design = []
for b in bs:
    o = run_dew(Tamb_c, RH_c, b, N, seed=11)
    m = np.arange(N) * TS > 300
    rate = cond_rate(o["Tc"][m], 24.0, 65.0)
    below = np.mean(o["Tc"][m] < o["Tdew_true"][m])
    design.append(dict(b=b, rate=float(np.mean(rate)), below=float(below),
                       pwm=float(np.mean(o["pwm"][m])), sd_Tc=float(np.std(o["Tc"][m]))))
nums["design"] = design
sig_td = unc["0.5,3.0"]["at65"]
fig, ax = plt.subplots(1, 2, figsize=(W, 1.9))
ax[0].plot(bs, [d["rate"] for d in design], "o-", color=K, ms=3)
ax[0].axhline(0, color=G2, lw=0.5)
ax[0].axvspan(0, sig_td, color="0.88", lw=0)
ax[0].text(sig_td / 2, max(d["rate"] for d in design) * 0.85, "within\ntolerance", ha="center", fontsize=6)
ax[0].set_xlabel(r"$b_{\mathrm{dew}}$ (°C)"); ax[0].set_ylabel("condensation (g/h)")
ax[1].plot(bs, [d["pwm"] / 2.55 for d in design], "s-", color=K, ms=3)
ax[1].set_xlabel(r"$b_{\mathrm{dew}}$ (°C)"); ax[1].set_ylabel("mean duty cycle (%)")
fig.tight_layout()
simlabel(ax[1])
save(fig, "fig10_design")

# rate uncertainty: if the true T_dew differs from the measured by +-sigma
nums["sig_td"] = float(sig_td)
nums["rate_vs_bias"] = {f"{b}": [float(cond_rate(17.03 - b + d, 24.0, 65.0)) for d in (-sig_td, 0, sig_td)] for b in (0.5, 1.0, 1.5, 2.0)}

# ---------------------------------------------------------------- filtering the moving reference
filt = {}
N2 = 6001
for tau in (0.0, 20.0, 60.0):
    o = run_dew(Tamb_c[:N2], RH_c[:N2], 1.5, N2, tau_ref=tau if tau else None, seed=12)
    m = np.arange(N2) * TS > 300
    filt[str(tau)] = dict(sd_ref=float(np.std(o["Tref"][m])), sd_pwm=float(np.std(o["pwm"][m])),
                          sd_Tc=float(np.std(o["Tc"][m])))
nums["ref_filter"] = filt

# ---------------------------------------------------------------- a simulated day: fixed vs moving reference
H = 8.0
N3 = int(H * 3600 / TS) + 1
t3 = np.arange(N3) * TS
ph = t3 / (H * 3600)
Tamb_d = 22.0 + 4.0 * np.sin(np.pi * ph) ** 1.5           # warms up towards the afternoon
RH_d = 80.0 - 25.0 * np.sin(np.pi * ph * 0.9) ** 1.2     # air dries
Td_d = v.magnus_dew(Tamb_d, RH_d)
b = 1.5
runs = {}
fixed_start = float(Td_d[0] - b)
fixed_cons = float(Td_d.min() - b)
for name, kw in [("moving", dict(mode="moving")), ("fixed (start)", dict(mode="fixed", fixed=fixed_start)),
                 ("fixed (lowest)", dict(mode="fixed", fixed=fixed_cons))]:
    o = run_dew(Tamb_d, RH_d, b, N3, seed=13, **kw)
    rate = np.maximum(cond_rate(o["Tc"], Tamb_d, RH_d), 0) * 0 + cond_rate(o["Tc"], Tamb_d, RH_d)
    mass = float(np.sum(rate) * TS / 3600)
    runs[name] = dict(mass=mass, below=float(np.mean(o["Tc"] < Td_d)), duty=float(np.mean(o["pwm"]) / 2.55),
                      gpkwh=None, o=o)
for name in runs:
    o = runs[name].pop("o")
    runs[name]["o_ds"] = dict(t=t3[::600] / 3600, Tc=o["Tc"][::600], Tref=o["Tref"][::600], pwm=o["pwm"][::600],
                              rate=cond_rate(o["Tc"], Tamb_d, RH_d)[::600])
nums["day"] = {k: {kk: vv for kk, vv in r.items() if kk != "o_ds"} for k, r in runs.items()}
o = run_dew(Tamb_d, RH_d, 2.4, N3, seed=13)
nums["day"]["moving b=2.4"] = dict(mass=float(np.sum(cond_rate(o["Tc"], Tamb_d, RH_d)) * TS / 3600),
                                   duty=float(np.mean(o["pwm"]) / 2.55))
nums["day"]["fixed_start"] = fixed_start; nums["day"]["fixed_lowest"] = fixed_cons
nums["day"]["Td_range"] = [float(Td_d.min()), float(Td_d.max())]
nums["day"]["Tamb_range"] = [float(Tamb_d.min()), float(Tamb_d.max())]
nums["day"]["RH_range"] = [float(RH_d.min()), float(RH_d.max())]

fig, ax = plt.subplots(3, 1, figsize=(W, 3.4), sharex=True)
th = t3[::600] / 3600
ax[0].plot(t3 / 3600, Tamb_d, color=G2, lw=0.9, label=r"$T_{\mathrm{amb}}$")
ax[0].plot(t3 / 3600, Td_d, color=K, lw=1.0, label=r"$T_{\mathrm{dew}}$")
a0 = ax[0].twinx(); a0.grid(False); a0.plot(t3 / 3600, RH_d, color=C[0], lw=0.8, ls="--"); a0.set_ylabel("RH (%)", color=C[0])
ax[0].set_ylabel("°C"); ax[0].legend(loc="center left", fontsize=6)
sty = {"moving": dict(color=K, ls="-"), "fixed (start)": dict(color=C[1], ls="--"), "fixed (lowest)": dict(color=G2, ls=":")}
for name, r in runs.items():
    ax[1].plot(r["o_ds"]["t"], r["o_ds"]["Tc"], lw=0.9, label=name, **sty[name])
    ax[2].plot(r["o_ds"]["t"], r["o_ds"]["rate"], lw=0.9, **sty[name])
ax[1].plot(t3 / 3600, Td_d, color=K, lw=0.5, alpha=0.5)
ax[1].set_ylabel(r"$T_{\mathrm{sup}}$ (°C)"); ax[1].legend(loc="upper right", fontsize=6, ncol=3)
ax[1].set_ylim(13.5, 20.5)
ax[2].axhline(0, color=G2, lw=0.5)
ax[2].set_ylabel("condensation (g/h)"); ax[2].set_xlabel("time (h)")
ax[2].set_ylim(-0.06, 0.15)
simlabel(ax[2])
save(fig, "fig10_day")

json.dump(nums, open(HERE + "/ch10_numbers.json", "w"), indent=1, default=float)
print(json.dumps({k: nums[k] for k in ("design", "rate_vs_bias", "day")}, default=float))
