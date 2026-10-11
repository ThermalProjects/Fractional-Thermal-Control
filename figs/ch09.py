"""Figures and numbers for Chapter 9 -- Seeing What Is Not Measured (observer)."""
import json
import numpy as np
from scipy.signal import place_poles
from style import plt, save, W, K, G1, G2, C, simlabel
from style import HERE
import vplant as v

nums = {}
F = v.TUNE_FOPID
TS = v.TS


def model(p):
    A = np.array([[-1 / (p["R1"] * p["Cc"]), 1 / (p["R1"] * p["Cc"]), 0],
                  [1 / (p["R1"] * p["Cm"]), -1 / (p["R1"] * p["Cm"]) - 1 / (p["R2"] * p["Cm"]), 1 / (p["R2"] * p["Cm"])],
                  [0, 1 / (p["R2"] * p["Ch"]), -1 / (p["R2"] * p["Ch"]) - 1 / (p["Rconv"] * p["Ch"])]])
    B = np.array([-p["phi"] * p["kP"] / p["Cc"], 0, (1 - p["phi"]) * p["kP"] / p["Ch"]])
    E = np.array([0, 0, 1 / (p["Rconv"] * p["Ch"])])
    return np.eye(3) + TS * A, TS * B, TS * E, A


P0 = dict(v.P)
Ad, Bd, Ed, Ac = model(P0)
eig = np.linalg.eigvals(Ac)
nums["plant_eig"] = sorted([float(-1 / e.real) for e in eig])      # time constants (s)
delay = max(1, int(round(P0["td"] / TS)))

# ---------------------------------------------------------------- observability per sensor
def gram(Cm, T=300.0):
    n = int(T / TS)
    Wo = np.zeros((3, 3)); Ak = np.eye(3)
    for _ in range(n):
        Wo += Ak.T @ Cm.T @ Cm @ Ak * TS
        Ak = Ad @ Ak
    return Wo


obs = {}
for name, Cm in [("cold face", np.array([[1, 0, 0]])), ("hot face", np.array([[0, 0, 1]])),
                 ("both faces", np.array([[1, 0, 0], [0, 0, 1]]))]:
    O = np.vstack([Cm @ np.linalg.matrix_power(Ad, i) for i in range(3)])
    s = np.linalg.svd(gram(Cm), compute_uv=False)
    # how well each state is "seen": diagonal of the inverse Gramian = error variance per unit noise
    Winv = np.linalg.inv(gram(Cm))
    obs[name] = dict(rank=int(np.linalg.matrix_rank(O)), cond=float(s[0] / s[-1]),
                     smin=float(s[-1]), sd=[float(np.sqrt(Winv[i, i])) for i in range(3)])
nums["observability"] = obs

# ---------------------------------------------------------------- observer gain
def observer_gain(Cm, tau):
    """Luenberger gain placing the three observer poles at continuous time
    constants tau, 1.25 tau, 1.5 tau (s)."""
    zp = np.exp(-TS / (tau * np.array([1.0, 1.25, 1.5])))
    return place_poles(Ad.T, Cm.T, zp).gain_matrix.T


class Observer:
    def __init__(self, x0, L, p=P0, Cm=np.array([[1.0, 0, 0]])):
        self.Ad, self.Bd, self.Ed, _ = model(p)
        self.x = np.array(x0, float); self.L = L; self.C = Cm
        self.ubuf = [0.0] * (delay + 1)

    def step(self, y, u, Tamb):
        """Predict-correct: y measured now, u applied now, Tamb measured.
        Returns the estimate of the current state and the innovation."""
        inn = np.atleast_1d(y) - self.C @ self.x
        xk = self.x + self.L @ inn                       # corrected estimate of x[k]
        self.ubuf.append(float(u)); ud = self.ubuf.pop(0)
        self.x = self.Ad @ xk + self.Bd * ud + self.Ed * Tamb
        return xk, inn


CC = np.array([[1.0, 0, 0]])
CB = np.array([[1.0, 0, 0], [0, 0, 1.0]])
N = 9001
t = np.arange(N) * TS
Tref = np.where(t < 300, 18.2, np.where(t < 600, 15.0, 20.0))
TRUE = dict(R2=P0["R2"] * 1.10, Ch=P0["Ch"] * 0.85, kP=P0["kP"] * 1.05)
XWARM = (25.0, 29.0, 30.0)          # restarted while the module interior is still warm


def closed_loop(p_true=TRUE, x_init=XWARM, rconv=None, noise=0.03, seed=3, Tamb=23.7):
    rng = np.random.default_rng(seed)
    c = v.FOPID(**F)
    pl = v.Plant(Tamb, p_true)
    pl.Tc, pl.Tm, pl.Th = x_init
    fl = v.MedEMA(pl.Tc)
    rc = np.ones(N) if rconv is None else rconv
    out = {k: np.zeros(N) for k in ["Tc", "Tm", "Th", "yc", "yh", "pwm"]}
    for k in range(N):
        yc = pl.Tc + rng.normal(0, noise); yh = pl.Th + rng.normal(0, noise)
        u = c.step(Tref[k] - fl(yc), Tref[k])
        out["pwm"][k] = u; out["yc"][k] = yc; out["yh"][k] = yh
        out["Tc"][k], out["Tm"][k], out["Th"][k] = pl.Tc, pl.Tm, pl.Th
        pl.step(u, Tamb, 0.0, rc[k])
    return out


x0 = np.array([23.7, 23.7, 23.7])   # observer knows nothing: starts at ambient


def run_obs(o, L, Cm=CB, meas=("yc", "yh"), x0=x0):
    ob = Observer(x0, L, P0, Cm)
    X = np.zeros((N, 3)); I = np.zeros((N, Cm.shape[0]))
    for k in range(N):
        y = np.array([o[mm][k] for mm in meas])
        X[k], I[k] = ob.step(y, o["pwm"][k], 23.7)
    return X, I


m = t > 300
rmse = lambda a, b, mm=m: float(np.sqrt(np.mean((a[mm] - b[mm]) ** 2)))


def tconv(err, thr=0.5):
    bad = np.where(np.abs(err) > thr)[0]
    return float((bad[-1] + 1) * TS) if len(bad) else 0.0


# ---------------------------------------------------------------- copy vs observer (both faces measured)
TAU = 5.0
LB = observer_gain(CB, TAU)
nums["tau_obs"] = TAU
nums["L_both"] = LB.tolist()
o = closed_loop()
Xol, _ = run_obs(o, np.zeros((3, 2)))
Xob, Iob = run_obs(o, LB)
oe = closed_loop(p_true={})
Xobe, _ = run_obs(oe, LB)
nums["est"] = dict(copy=dict(Tm=rmse(Xol[:, 1], o["Tm"]), Tm_conv=tconv(Xol[:, 1] - o["Tm"]),
                             Tm_err_end=float(Xol[-1, 1] - o["Tm"][-1])),
                   observer=dict(Tm=rmse(Xob[:, 1], o["Tm"]), Tm_conv=tconv(Xob[:, 1] - o["Tm"]),
                                 Tm_bias=float(np.mean(Xob[m, 1] - o["Tm"][m]))),
                   observer_exact=dict(Tm=rmse(Xobe[:, 1], oe["Tm"]), Tm_conv=tconv(Xobe[:, 1] - oe["Tm"])))
fig, ax = plt.subplots(2, 1, figsize=(W, 2.8), sharex=True, gridspec_kw=dict(height_ratios=[1, 1.4]))
ax[0].plot(t, o["yc"], color=K, lw=0.6, label=r"$T_c$ measured")
ax[0].plot(t, o["yh"], color=G2, lw=0.6, label=r"$T_h$ measured")
ax[0].set_ylabel("°C"); ax[0].legend(loc="center right", ncol=2, fontsize=6)
ax[1].plot(t, o["Tm"], color=K, lw=1.1, label=r"true $T_m$ (no probe)")
ax[1].plot(t, Xol[:, 1], color=G2, lw=1.0, ls=":", label="open-loop copy")
ax[1].plot(t, Xob[:, 1], color=C[0], lw=1.0, ls="--", label="observer")
ax[1].set_ylabel(r"$T_m$ (°C)"); ax[1].set_xlabel("time (s)")
ax[1].legend(loc="lower right", ncol=3, fontsize=6, bbox_to_anchor=(1.0, 0.1))
ax[1].set_ylim(22.5, 31)
simlabel(ax[1])
save(fig, "fig09_estimate")

# ---------------------------------------------------------------- choosing L: speed vs noise and model error
taus = [2, 5, 10, 20, 40, 80]
trade = []
on = closed_loop(p_true={}, seed=5)                 # exact model, noise only
for tau in taus:
    Lt = observer_gain(CB, tau)
    Xe, _ = run_obs(on, Lt)
    Xm, _ = run_obs(o, Lt)
    trade.append(dict(tau=tau, conv=tconv(Xe[:, 1] - on["Tm"]),
                      noise_sd=float(np.std(np.diff(Xe[m, 1])) / np.sqrt(2)),
                      rmse_mismatch=rmse(Xm[:, 1], o["Tm"]),
                      Lmax=float(np.abs(Lt).max())))
nums["L_trade"] = trade
fig, ax = plt.subplots(1, 2, figsize=(W, 1.8))
ax[0].plot(taus, [r["conv"] for r in trade], "o-", color=K, ms=3)
ax[0].set_xscale("log"); ax[0].set_xlabel(r"observer time constant $\tau_o$ (s)"); ax[0].set_ylabel(r"$\hat T_m$ convergence (s)")
ax[1].plot(taus, [r["noise_sd"] for r in trade], "s-", color=K, ms=3, label="noise in $\\hat T_m$")
ax[1].plot(taus, [r["rmse_mismatch"] for r in trade], "^--", color=C[0], ms=3, label="error, model 5--15 % off")
ax[1].set_xscale("log"); ax[1].set_yscale("log"); ax[1].set_xlabel(r"observer time constant $\tau_o$ (s)"); ax[1].set_ylabel("°C")
ax[1].legend(fontsize=6, loc="center right")
for a in ax:
    a.axvline(TAU, color=C[1], lw=0.7, ls="--")
fig.tight_layout()
simlabel(ax[1])
save(fig, "fig09_gain")

# ---------------------------------------------------------------- observability per sensor (bar chart)
sig = 0.03
fig, ax = plt.subplots(figsize=(W, 1.8))
labels = [r"$T_c$", r"$T_m$", r"$T_h$"]
xx = np.arange(3)
for i, (name, col, hatch) in enumerate([("cold face", K, ""), ("hot face", G2, "///"), ("both faces", C[0], "")]):
    sd = np.array(obs[name]["sd"]) * sig * np.sqrt(TS)
    ax.bar(xx + (i - 1) * 0.26, sd, width=0.25, color=col, hatch=hatch, edgecolor="0.2", lw=0.3, label=f"probe on {name}" if name != "both faces" else "both probes")
ax.set_yscale("log"); ax.set_xticks(xx); ax.set_xticklabels(labels)
ax.set_ylabel("best attainable s.d. (°C)")
ax.legend(fontsize=6, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.18))
save(fig, "fig09_observability")
nums["obs_sd_scaled"] = {k: [float(x * sig * np.sqrt(TS)) for x in obs[k]["sd"]] for k in obs}

# ---------------------------------------------------------------- hide a sensor: estimate T_h from T_c only
hide = {}
fig, ax = plt.subplots(2, 1, figsize=(W, 3.0), sharex=True)
for j, (lab, oo) in enumerate([("model = plant", oe), ("model 5--15 % off", o)]):
    for tau, ls, col in [(80, "--", C[0]), (20, ":", C[1])]:
        Lc = observer_gain(CC, tau)
        Xh, _ = run_obs(oo, Lc, Cm=CC, meas=("yc",))
        hide[f"{lab}, tau={tau}"] = dict(rmse=rmse(Xh[:, 2], oo["yh"]), bias=float(np.mean(Xh[m, 2] - oo["yh"][m])),
                                          peak=float(np.max(np.abs(Xh[:, 2] - oo["Th"]))), L=Lc.ravel().tolist())
        ax[j].plot(t, Xh[:, 2], color=col, ls=ls, lw=0.9, label=fr"$\hat T_h$, $\tau_o$ = {tau} s")
    ax[j].plot(t, oo["yh"], color=K, lw=0.7, label=r"$T_h$ logged (hidden from the observer)")
    ax[j].set_ylim(20, 50); ax[j].set_ylabel(r"$T_h$ (°C)")
    ax[j].set_title(lab, fontsize=7)
ax[1].legend(loc="lower right", fontsize=6, ncol=3)
ax[1].set_xlabel("time (s)")
simlabel(ax[0])
fig.tight_layout()
save(fig, "fig09_hide")
nums["hide"] = hide

# ---------------------------------------------------------------- when the plant changes: fan failure at 450 s
rc = np.where(t >= 450, 3.0, 1.0)
of = closed_loop(rconv=rc, seed=4)
Xf, If = run_obs(of, LB)
win = 100
ma = lambda x: np.convolve(x, np.ones(win) / win, mode="full")[:len(x)]   # causal 10 s mean
rc_, rh_ = ma(If[:, 0]), ma(If[:, 1])
healthy = (t > 200) & (t < 440)
thr_c = 5 * np.std(rc_[healthy]) + abs(np.mean(rc_[healthy]))
thr_h = 5 * np.std(rh_[healthy]) + abs(np.mean(rh_[healthy]))
det_h = t[np.where((t > 450) & (np.abs(rh_) > thr_h))[0]]
det_c = t[np.where((t > 450) & (np.abs(rc_) > thr_c))[0]]
Xn, In = run_obs(o, LB)
fa = int(np.sum((t > 200) & ((np.abs(ma(In[:, 1])) > thr_h) | (np.abs(ma(In[:, 0])) > thr_c))))
Lc80 = observer_gain(CC, 20.0)
_, Ifc = run_obs(of, Lc80, Cm=CC, meas=("yc",))
rcc = ma(Ifc[:, 0]); thr_cc = 5 * np.std(rcc[healthy]) + abs(np.mean(rcc[healthy]))
det_cc = t[np.where((t > 450) & (np.abs(rcc) > thr_cc))[0]]
nums["fault_cold_only"] = dict(thr=float(thr_cc), detect=float(det_cc[0] - 450) if len(det_cc) else None,
                               max_after=float(np.max(np.abs(rcc[t > 450]))))
nums["fault"] = dict(thr_c=float(thr_c), thr_h=float(thr_h),
                     detect_hot=float(det_h[0] - 450) if len(det_h) else None,
                     detect_cold=float(det_c[0] - 450) if len(det_c) else None,
                     false_alarm_samples=fa,
                     pwm_before=float(np.mean(of["pwm"][(t > 380) & (t < 450)])),
                     pwm_after=float(np.mean(of["pwm"][(t > 560) & (t < 600)])),
                     Th_before=float(np.mean(of["Th"][(t > 380) & (t < 450)])),
                     Th_after=float(np.mean(of["Th"][(t > 560) & (t < 600)])))
fig, ax = plt.subplots(3, 1, figsize=(W, 3.2), sharex=True)
ax[0].plot(t, of["yc"], color=K, lw=0.6); ax[0].plot(t, Tref, color=C[1], lw=0.6)
ax[0].set_ylabel(r"$T_{\mathrm{sup}}$ (°C)")
a0 = ax[0].twinx(); a0.grid(False); a0.plot(t, of["pwm"], color=G2, lw=0.5); a0.set_ylabel("PWM", color=G2)
ax[1].plot(t, of["yh"], color=K, lw=0.7, label=r"$T_h$ measured")
ax[1].plot(t, Xf[:, 2], color=C[0], lw=0.8, ls="--", label=r"$T_h$ predicted by the observer model")
ax[1].set_ylabel(r"$T_h$ (°C)"); ax[1].legend(loc="upper left", fontsize=6)
ax[2].plot(t, rh_, color=C[0], lw=0.8, ls="--", label="hot-face innovation")
ax[2].plot(t, rc_, color=K, lw=0.8, label="cold-face innovation")
for th, col in [(thr_h, C[0]), (thr_c, K)]:
    ax[2].axhline(th, color=col, lw=0.5, ls=":"); ax[2].axhline(-th, color=col, lw=0.5, ls=":")
ax[2].set_ylabel("innovation, 10 s mean (°C)"); ax[2].set_xlabel("time (s)")
ax[2].legend(loc="upper left", fontsize=6)
ax[2].set_ylim(-0.4, 0.8)
for a in ax:
    a.axvline(450, color=C[1], lw=0.6, ls="--")
ax[0].set_xlim(200, 900)
simlabel(ax[2])
save(fig, "fig09_fault")

json.dump(nums, open(HERE + "/ch09_numbers.json", "w"), indent=1, default=float)
print(json.dumps({k: nums[k] for k in nums if k not in ("observability",)}, indent=1, default=float))

# persistence check for the alarm rule
def runs(mask):
    r, c = [], 0
    for x in mask:
        if x: c += 1
        elif c: r.append(c); c = 0
    if c: r.append(c)
    return r
fa_mask = (t > 200) & ((np.abs(ma(In[:, 1])) > thr_h) | (np.abs(ma(In[:, 0])) > thr_c))
det_mask = (t > 450) & (np.abs(rh_) > thr_h)
nums["alarm_runs"] = dict(false_runs_s=[r * TS for r in runs(fa_mask)], false_times=t[fa_mask][:5].tolist(),
                          det_first_run_s=runs(det_mask)[0] * TS if runs(det_mask) else 0)
json.dump(nums, open(HERE + "/ch09_numbers.json", "w"), indent=1, default=float)
print(nums["alarm_runs"])
