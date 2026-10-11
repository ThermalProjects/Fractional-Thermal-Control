"""Figures and numbers for Chapter 11 -- The Digital Twin and the Final Experiment.

One simulated two-hour experiment with everything on: dew-point reference with
b_dew = 1.5 C and a filtered dew point, FOPID-AS, two-probe observer with an
innovation alarm, and a digital twin (the virtual plant fed with the measured
ambient temperature and the applied PWM) with a residual band and a
recalibration.  From minute 60 the heat sink degrades slowly (R_conv rises
linearly to x3 at minute 120), as when a fan slows down or dust builds up.
"""
import json
import numpy as np
from scipy.optimize import minimize_scalar
from style import plt, save, W, K, G1, G2, C, simlabel
from style import HERE
import vplant as v

nums = {}
TS = v.TS
RV, AREA, HM = 461.5, 0.040 * 0.040, 0.006


def cond_rate(Ts, Tamb, RH):
    pv = RH / 100 * v.p_sat(Tamb)
    return HM * AREA * (pv - v.p_sat(Ts)) * 100 / (RV * (0.5 * (Ts + Tamb) + 273.15)) * 3.6e6


Hh = 2.0
N = int(Hh * 3600 / TS) + 1
t = np.arange(N) * TS
tm = t / 60
rng = np.random.default_rng(21)
Tamb = 23.0 + 2.0 * np.clip(t / (Hh * 3600), 0, 1) + 0.15 * np.sin(2 * np.pi * t / 1500)
RH = 70.0 - 10.0 * np.clip(t / (Hh * 3600), 0, 1) + 1.0 * np.sin(2 * np.pi * t / 2100)
rconv = np.where(tm < 60, 1.0, 1.0 + 2.0 * (tm - 60) / 60)
TRUE = dict(kP=v.P["kP"] * 1.01, R2=v.P["R2"] * 1.03)      # the plant is never exactly the model
BDEW, TAU_R = 1.5, 20.0

# controller, observer, twin
ctrl = v.Supervisor(v.FOPID(**v.TUNE_FOPID), near=v.SUP_NEAR, far=v.SUP_FAR, **v.SUP_OPTS)
plant = v.Plant(Tamb[0], TRUE)
filt = v.MedEMA(plant.Tc); filth = v.MedEMA(plant.Th)
LB = v.observer_gain([[1, 0, 0], [0, 0, 1]], 5.0)
obs = v.Observer([Tamb[0]] * 3, LB)
twin_p = dict(v.P)
twin = v.Plant(Tamb[0], twin_p)

rec = {k: np.zeros(N) for k in ["Tc", "Tm", "Th", "ysup", "yh", "Tref", "Tdew", "pwm", "Tm_hat", "inn_c", "inn_h",
                                 "tw_c", "tw_h", "w", "Q", "mass", "below"]}
Tdew_f = float(v.magnus_dew(Tamb[0], RH[0]))
a_ref = 2.0 / TAU_R
Q = 0.0; mass = 0.0; nbelow = 0
recal_times = []
band = None
outside = 0
win_res = []
pwm_log, tamb_log, yh_log, yc_log = [], [], [], []

for k in range(N):
    if k % 20 == 0:
        Tdm = float(v.magnus_dew(Tamb[k] + rng.normal(0, 0.03), RH[k] + rng.normal(0, 0.3)))
        Tdew_f += a_ref * (Tdm - Tdew_f)
    Tref = Tdew_f - BDEW if t[k] > 30 else Tdew_f - BDEW
    yc = plant.Tc + rng.normal(0, 0.03); yh = plant.Th + rng.normal(0, 0.03)
    ys = filt(yc); yhf = filth(yh)
    e = Tref - ys
    u = ctrl.step(e, Tref)
    xh, inn = obs.step(np.array([ys, yhf]), u, Tamb[k])
    # twin: same inputs, same moment, no correction
    tc_, tm_, th_ = twin.Tc, twin.Tm, twin.Th
    rec["tw_c"][k], rec["tw_h"][k] = tc_, th_
    twin.step(u, Tamb[k])
    # condensation bookkeeping (uses the measured, filtered dew point)
    Q += max(0.0, Tdew_f - ys) * TS / 3600
    nbelow += ys < Tdew_f
    mass += cond_rate(plant.Tc, Tamb[k], RH[k]) * TS / 3600
    rec["Tc"][k], rec["Tm"][k], rec["Th"][k] = plant.Tc, plant.Tm, plant.Th
    rec["ysup"][k], rec["yh"][k] = ys, yhf
    rec["Tref"][k], rec["Tdew"][k], rec["pwm"][k] = Tref, Tdew_f, u
    rec["Tm_hat"][k] = xh[1]; rec["inn_c"][k], rec["inn_h"][k] = inn
    rec["w"][k] = ctrl.f
    rec["Q"][k], rec["mass"][k], rec["below"][k] = Q, mass, nbelow / (k + 1)
    pwm_log.append(u); tamb_log.append(Tamb[k]); yh_log.append(yhf); yc_log.append(ys)
    plant.step(u, Tamb[k], 0.0, rconv[k])
    # twin residual supervision (hot face), band from the first 20 healthy minutes
    res = yhf - th_
    if band is None and t[k] >= 20 * 60:
        r0 = rec["yh"][3000:k] - rec["tw_h"][3000:k]
        band = (float(np.mean(r0)), float(4 * np.std(r0) + 0.3))
    if band is not None:
        outside = outside + 1 if abs(res - band[0]) > band[1] else 0
        if outside >= 600 and (not recal_times or t[k] - recal_times[-1] > 900):     # 60 s outside: recalibrate
            j0 = max(0, k - 9000)                                       # last 15 min
            P_ = np.array(pwm_log[j0 - 10:k + 1]); A_ = np.array(tamb_log[j0 - 10:k + 1])
            Yh = np.array(yh_log[j0:k + 1]); Yc = np.array(yc_log[j0:k + 1])

            def sse(mult):
                pp = dict(twin_p); pp["Rconv"] = v.P["Rconv"] * mult
                pl = v.Plant(A_[0], pp)
                pl.Tc, pl.Tm, pl.Th = rec["ysup"][j0 - 10], rec["Tm_hat"][j0 - 10], rec["yh"][j0 - 10]
                out = []
                for i in range(len(P_)):
                    out.append(pl.Th); pl.step(P_[i], A_[i])
                out = np.array(out[10:])
                return float(np.mean((out - Yh) ** 2))
            r = minimize_scalar(sse, bounds=(0.8, 3.0), method="bounded", options=dict(xatol=1e-3))
            twin_p["Rconv"] = v.P["Rconv"] * r.x
            twin = v.Plant(Tamb[k], twin_p)
            twin.Tc, twin.Tm, twin.Th = ys, xh[1], yhf               # restart from measured faces + observed interior
            recal_times.append(float(t[k]))
            nums.setdefault("recal", []).append(dict(t_min=float(tm[k]), Rconv_mult=float(r.x),
                                                     true_mult=float(rconv[k])))
            outside = 0

res_h = rec["yh"] - rec["tw_h"]
res_c = rec["ysup"] - rec["tw_c"]
m_ss = tm > 5
nums["band"] = band
nums["summary"] = dict(
    mean_abs_err_to_ref=float(np.mean(np.abs(rec["ysup"][m_ss] - rec["Tref"][m_ss]))),
    max_abs_err_after_5min=float(np.max(np.abs(rec["ysup"][m_ss] - rec["Tref"][m_ss]))),
    time_below=float(rec["below"][-1]), Q_degCh=float(Q), mass_g=float(mass),
    duty=float(np.mean(rec["pwm"]) / 2.55), duty_first_h=float(np.mean(rec["pwm"][tm < 60]) / 2.55),
    duty_last_10=float(np.mean(rec["pwm"][tm > 110]) / 2.55),
    Th_first=float(np.mean(rec["Th"][(tm > 50) & (tm < 60)])), Th_last=float(np.mean(rec["Th"][tm > 110])),
    Tm_rmse=float(np.sqrt(np.mean((rec["Tm_hat"][m_ss] - rec["Tm"][m_ss]) ** 2))),
    res_h_rmse_first_h=float(np.sqrt(np.mean(res_h[(tm > 5) & (tm < 60)] ** 2))),
    res_c_rmse_first_h=float(np.sqrt(np.mean(res_c[(tm > 5) & (tm < 60)] ** 2))),
    sup_inside=float(np.mean(rec["w"][m_ss] == 0)),
    Tdew_range=[float(rec["Tdew"].min()), float(rec["Tdew"].max())])
# observer innovation alarm (B30 rule) on the hot face
win = 100
ma = lambda x: np.convolve(x, np.ones(win) / win, mode="full")[:len(x)]
ih = ma(rec["inn_h"]); healthy = (tm > 5) & (tm < 60)
thr = 5 * np.std(ih[healthy]) + abs(np.mean(ih[healthy]))
al = np.abs(ih) > thr
runs5 = np.convolve(al.astype(float), np.ones(50), mode="full")[:N] >= 50
nums["obs_alarm"] = dict(thr=float(thr), first_min=float(tm[np.argmax(runs5 & (tm > 60))]) if np.any(runs5 & (tm > 60)) else None)
# first time residual leaves band
if band:
    ob = np.abs(res_h - band[0]) > band[1]
    nums["res_first_out_min"] = float(tm[np.argmax(ob & (tm > 20))]) if np.any(ob & (tm > 20)) else None

# ---------------------------------------------------------------- figure: twin vs plant
fig, ax = plt.subplots(2, 1, figsize=(W, 2.7), sharex=True, gridspec_kw=dict(height_ratios=[1.3, 1]))
ax[0].plot(tm, rec["yh"], color=K, lw=0.8, label=r"$T_h$ measured")
ax[0].plot(tm, rec["tw_h"], color=C[0], lw=0.9, ls="--", label=r"$T_h$ of the twin")
ax[0].plot(tm, rec["ysup"], color=G2, lw=0.8, label=r"$T_{\mathrm{sup}}$ measured")
ax[0].plot(tm, rec["tw_c"], color=C[1], lw=0.8, ls=":", label=r"$T_{\mathrm{sup}}$ of the twin")
ax[0].set_ylabel("°C"); ax[0].legend(loc="center left", fontsize=6, ncol=2)
ax[1].plot(tm, res_h, color=C[0], lw=0.8, label="hot-face residual")
ax[1].plot(tm, res_c, color=K, lw=0.8, label="cold-face residual")
if band:
    ax[1].axhspan(band[0] - band[1], band[0] + band[1], color="0.9", lw=0)
for r in recal_times:
    for a in ax:
        a.axvline(r / 60, color=C[1], lw=0.7, ls="--")
ax[1].axvline(60, color=G2, lw=0.6, ls=":")
ax[1].text(60.5, ax[1].get_ylim()[1] * 0.8 if ax[1].get_ylim()[1] > 0 else 1, "heat sink starts to degrade", fontsize=6, color="0.3")
ax[1].set_ylabel("residual (°C)"); ax[1].set_xlabel("time (min)")
ax[1].legend(loc="upper left", fontsize=6)
simlabel(ax[1])
save(fig, "fig11_twin")

# ---------------------------------------------------------------- figure: final experiment panel
fig, ax = plt.subplots(5, 1, figsize=(W, 5.6), sharex=True)
ax[0].plot(tm, Tamb, color=G2, lw=0.8, label=r"$T_{\mathrm{amb}}$")
ax[0].plot(tm, rec["Tdew"], color=K, lw=0.9, label=r"$T_{\mathrm{dew}}$ (filtered)")
ax[0].set_ylabel("°C")
a0 = ax[0].twinx(); a0.grid(False); a0.plot(tm, RH, color=C[0], lw=0.7, ls="--"); a0.set_ylabel("RH (%)", color=C[0])
ax[0].legend(loc="center right", fontsize=6)
ax[1].plot(tm, rec["Tdew"], color=K, lw=0.6, alpha=0.6, label=r"$T_{\mathrm{dew}}$")
ax[1].plot(tm, rec["Tref"], color=C[1], lw=0.7, label=r"$T_{\mathrm{ref}} = T_{\mathrm{dew}} - 1.5$")
ax[1].plot(tm, rec["ysup"], color=C[0], lw=0.7, label=r"$T_{\mathrm{sup}}$")
ax[1].set_ylabel("°C"); ax[1].set_ylim(12, 19); ax[1].legend(loc="lower right", fontsize=6, ncol=3)
ax[2].plot(tm, rec["pwm"], color=K, lw=0.5)
ax[2].set_ylabel("PWM")
ax[2].set_ylim(0, 120)
a2 = ax[2].twinx(); a2.grid(False); a2.fill_between(tm, 0, runs5.astype(float), step="pre", color=C[1], alpha=0.25, lw=0); a2.set_ylabel("observer alarm", color=C[1]); a2.set_ylim(0, 1.05); a2.set_yticks([0, 1])
ax[3].plot(tm, rec["yh"], color=K, lw=0.8, label=r"$T_h$ measured")
ax[3].plot(tm, rec["Tm_hat"], color=C[0], lw=0.8, ls="--", label=r"$\hat T_m$ (observer)")
ax[3].set_ylabel("°C"); ax[3].legend(loc="upper left", fontsize=6, ncol=2)
ax[4].plot(tm, rec["Q"], color=K, lw=0.9, label=r"$Q_{\mathrm{dew}}$ (°C·h)")
ax[4].set_ylabel(r"$Q_{\mathrm{dew}}$ (°C·h)")
a4 = ax[4].twinx(); a4.grid(False); a4.plot(tm, rec["mass"] * 1000, color=C[0], lw=0.8, ls="--"); a4.set_ylabel("mass (mg)", color=C[0])
ax[4].set_xlabel("time (min)")
for r in recal_times:
    ax[3].axvline(r / 60, color=C[1], lw=0.6, ls="--")
simlabel(ax[4])
save(fig, "fig11_final")

json.dump(nums, open(HERE + "/ch11_numbers.json", "w"), indent=1, default=float)
print(json.dumps(nums, indent=1, default=float))
