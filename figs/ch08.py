"""Figures and numbers for Chapter 8 -- Adaptive Supervision of the FOPID."""
import json
import numpy as np
from style import plt, save, W, K, G1, G2, C, simlabel
from style import HERE
import vplant as v

nums = {}
F, Pd = v.TUNE_FOPID, v.TUNE_PID
NEAR, FAR, OPT = v.SUP_NEAR, v.SUP_FAR, v.SUP_OPTS


def make_as(**over):
    o = dict(OPT); o.update(over)
    return v.Supervisor(v.FOPID(**F), near=NEAR, far=FAR, **o)


# ---------------------------------------------------------------- scheduling functions
s = make_as()
ea = np.linspace(0, 4, 801)
up = [];
s.inside = False
for a in ea:                       # error magnitude rising
    up.append(dict(s.schedule(a), f=s.f))
dn = []
for a in ea[::-1]:                 # error magnitude falling
    dn.append(dict(s.schedule(a), f=s.f))
dn = dn[::-1]
fig, ax = plt.subplots(1, 3, figsize=(W, 1.75))
for i, (key, lab) in enumerate([("f", r"far weight $w_{\mathrm{far}}$"), ("Kp", r"$K_p$"), ("lam", r"$\lambda$")]):
    ax[i].plot(ea, [d[key] for d in up], color=K, lw=1.0, label="|e| rising")
    ax[i].plot(ea, [d[key] for d in dn], color=C[1], lw=0.9, ls="--", label="|e| falling")
    ax[i].set_xlabel(r"$|e|$ (°C)")
    ax[i].set_title(lab, fontsize=7)
    ax[i].axvspan(OPT["e_low"], OPT["e_low"] + OPT["hyst"], color="0.85", lw=0)
    for x in (OPT["e_low"], OPT["e_high"]):
        ax[i].axvline(x, color=G2, lw=0.5, ls=":")
ax[0].legend(loc="lower right", fontsize=6)
fig.tight_layout()
save(fig, "fig08_schedule")

# ---------------------------------------------------------------- scenario: setpoint change + load
N = 6001
t = np.arange(N) * v.TS
Tref = np.where(t < 200, 18.2, 15.0)
Qd = np.where((t >= 400) & (t < 460), 2.0, 0.0)


def run(ctrl, p=None, seed=1):
    return v.simulate(ctrl, N=N, Tref=Tref, Qd=Qd, noise=0.03, seed=seed, p=p)


def seg_metrics(o):
    e = o["Tref"] - o["Tc"]
    r = {}
    for name, a, b in [("S1 start", 0, 200), ("S2 step", 200, 400), ("S3 load", 400, 600)]:
        m = (t >= a) & (t < b)
        iae = float(np.sum(np.abs(e[m])) * v.TS)
        under = float(max(0.0, np.max(e[m])))
        peak = float(np.max(np.abs(e[m] if a >= 400 else e[m][t[m] > a + 5])))
        r[name] = dict(IAE=iae, undershoot=under, peak=peak,
                       duty=float(np.mean(o["pwm"][m]) / 2.55),
                       sdpwm=float(np.std(o["pwm"][m & (t > a + 120)])) if a < 400 else float("nan"))
    return r


res = {}
outs = {}
for name, mk in [("PID", lambda: v.FOPID(**Pd)), ("FOPID", lambda: v.FOPID(**F)), ("FOPID-AS", make_as)]:
    c = mk()
    o = run(c)
    outs[name] = o
    res[name] = seg_metrics(o)
nums["scenario"] = res

fig, ax = plt.subplots(3, 1, figsize=(W, 3.6), sharex=False)
sty = {"PID": dict(color=G2, ls=":", lw=1.0), "FOPID": dict(color=K, ls="--", lw=0.9),
       "FOPID-AS": dict(color=C[0], ls="-", lw=1.0)}
for name, o in outs.items():
    ax[0].plot(t, o["Tc"], label=name, **sty[name])
    ax[1].plot(t, o["Tc"], **sty[name])
    ax[2].plot(t, o["Tc"], **sty[name])
for a in ax:
    a.plot(t, Tref, color=C[1], lw=0.6, ls="-", alpha=0.7)
    a.set_ylabel(r"$T_{\mathrm{sup}}$ (°C)")
ax[0].set_xlim(0, 600); ax[0].set_ylim(13.5, 24.5)
ax[0].axvspan(400, 460, color="0.9", lw=0)
ax[0].legend(loc="upper right", ncol=3)
ax[0].set_title("whole run: start, reference step at 200 s, 2 W load 400--460 s", fontsize=7)
ax[1].set_xlim(195, 260); ax[1].set_ylim(14.5, 18.6)
ax[1].set_title("reference step (zoom)", fontsize=7)
ax[2].set_xlim(395, 540); ax[2].set_ylim(14.1, 16.0)
ax[2].axvspan(400, 460, color="0.9", lw=0)
ax[2].set_title("load (zoom)", fontsize=7)
ax[2].set_xlabel("time (s)")
fig.tight_layout()
simlabel(ax[2])
save(fig, "fig08_compare")

# drive and scheduled parameters for the AS run
o = outs["FOPID-AS"]
c = make_as()
Kp_h, lam_h, f_h = [], [], []
pl = v.Plant(23.7); fl = v.MedEMA(pl.Tc); rng = np.random.default_rng(1)
for k in range(N):
    Tf = fl(pl.Tc + rng.normal(0, 0.03))
    e = Tref[k] - Tf
    p = c.step(e, Tref[k])
    pl.step(p, 23.7, Qd[k])
    Kp_h.append(c.c.Kp); lam_h.append(c.c.lam); f_h.append(c.f)
fig, ax = plt.subplots(3, 1, figsize=(W, 2.9), sharex=True)
ax[0].plot(t, outs["FOPID"]["pwm"], **sty["FOPID"], label="FOPID")
ax[0].plot(t, o["pwm"], **sty["FOPID-AS"], label="FOPID-AS")
ax[0].set_ylabel("PWM (counts)")
ax[0].legend(loc="upper right", ncol=2)
ax[1].plot(t, Kp_h, color=K, lw=0.8)
ax[1].set_ylabel(r"$K_p$")
ax[2].plot(t, lam_h, color=K, lw=0.9)
ax[2].set_ylabel(r"$\lambda$")
ax[2].set_xlabel("time (s)")
ax[2].set_xlim(0, 600)
simlabel(ax[2])
save(fig, "fig08_schedule_run")
nums["sched_fraction_far"] = float(np.mean(np.array(f_h) > 0))
nums["sched_fraction_inside"] = float(np.mean(np.array(f_h) == 0))

# ---------------------------------------------------------------- changing the orders online: keep vs reset
def order_jump(mode):
    c = v.FOPID(**F)
    pl = v.Plant(23.7); fl = v.MedEMA(pl.Tc)
    P = []
    for k in range(2401):
        Tf = fl(pl.Tc)
        e = 18.2 - Tf
        if k == 1500:
            if mode == "keep":
                c.set_orders(1.0, 1.24)
            elif mode == "reset":
                c.set_orders(1.0, 1.24); c.Iacc = 0.0; c.eh[:] = 0.0
        p = c.step(e, 18.2)
        pl.step(p, 23.7)
        P.append((pl.Tc, p, c.u))
    return np.array(P)


jk, jr = order_jump("keep"), order_jump("reset")
tt = np.arange(2401) * v.TS
nums["orders"] = dict(keep_du=float(np.max(np.abs(np.diff(jk[1495:1510, 2])))),
                      reset_du=float(np.max(np.abs(np.diff(jr[1495:1510, 2])))),
                      reset_dT=float(np.max(np.abs(jr[1500:, 0] - 18.2))),
                      keep_dT=float(np.max(np.abs(jk[1500:, 0] - 18.2))))
fig, ax = plt.subplots(1, 2, figsize=(W, 1.8))
ax[0].plot(tt, jk[:, 1], color=K, lw=0.9, label="keep integral and history")
ax[0].plot(tt, jr[:, 1], color=C[1], lw=0.9, ls="--", label="reset them")
ax[0].set_xlim(145, 175); ax[0].set_ylabel("PWM (counts)"); ax[0].set_xlabel("time (s)")
ax[0].axvline(150, color=G2, lw=0.5, ls=":")
ax[0].legend(loc="upper right", fontsize=6)
ax[1].plot(tt, jk[:, 0], color=K, lw=0.9)
ax[1].plot(tt, jr[:, 0], color=C[1], lw=0.9, ls="--")
ax[1].set_xlim(145, 240); ax[1].set_ylabel(r"$T_{\mathrm{sup}}$ (°C)"); ax[1].set_xlabel("time (s)")
ax[1].axvline(150, color=G2, lw=0.5, ls=":")
fig.tight_layout()
simlabel(ax[1])
save(fig, "fig08_orders")

# ---------------------------------------------------------------- frozen-parameter margins along the schedule
# Linearised discrete loop: controller (exact GL implementation, linear derivative), plant
# discretised (Euler, as in the virtual plant), transport delay, EMA filter.
p = v.P
A = np.array([[-1 / (p["R1"] * p["Cc"]), 1 / (p["R1"] * p["Cc"]), 0],
              [1 / (p["R1"] * p["Cm"]), -1 / (p["R1"] * p["Cm"]) - 1 / (p["R2"] * p["Cm"]), 1 / (p["R2"] * p["Cm"])],
              [0, 1 / (p["R2"] * p["Ch"]), -1 / (p["R2"] * p["Ch"]) - 1 / (p["Rconv"] * p["Ch"])]])
B = np.array([-p["phi"] * p["kP"] / p["Cc"], 0, (1 - p["phi"]) * p["kP"] / p["Ch"]])
Ad = np.eye(3) + v.TS * A
Bd = v.TS * B
delay = max(1, int(round(p["td"] / v.TS)))
w = np.logspace(-4, np.log10(np.pi / v.TS) - 1e-3, 6000)
z = np.exp(1j * w * v.TS)


def G(z):
    return np.array([(np.array([1, 0, 0]) @ np.linalg.solve(zi * np.eye(3) - Ad, Bd)) for zi in z]) * z ** (-delay)


Gz = G(z)
Fz = 0.2 / (1 - 0.8 / z)                      # EMA


def Cz(Kp, Ki, Kd, lam, mu, M=30):
    wv = v.gl_weights(1 - lam, M); wd = v.gl_weights(mu, M)
    zz = np.array([z ** (-j) for j in range(M + 1)])
    V = (wv @ zz) / v.TS ** (1 - lam)
    I = v.TS / (1 - 1 / z) * V
    D = (wd @ zz) / v.TS ** mu
    return Kp + Ki * I + Kd * D


def margins(L):
    mag = np.abs(L); ph = np.unwrap(np.angle(L))
    i = np.where((mag[:-1] >= 1) & (mag[1:] < 1))[0]
    pm = np.nan; wc = np.nan
    if len(i):
        i = i[-1]; wc = w[i]; pm = 180 + np.degrees(ph[i]) % 360 - 360 if False else np.degrees(np.angle(-L[i])) * -1
        pm = 180 - abs(np.degrees(np.angle(L[i])))
    # gain margin: |L| where phase crosses -180
    ang = np.degrees(np.angle(L))
    j = np.where((np.sign(np.imag(L[:-1])) != np.sign(np.imag(L[1:]))) & (np.real(L[:-1]) < 0))[0]
    gm = (20 * np.log10(1 / np.abs(L[j[0]]))) if len(j) else np.inf
    return pm, gm, wc


rows = []
ff = np.linspace(0, 1, 21)
for A_ in (1.0, OPT["A"]):
    for f in ff:
        q = {k: NEAR[k] + f * (FAR[k] - NEAR[k]) for k in NEAR}
        L = -Cz(q["Kp"] * A_, F["Ki"], F["Kd"], q["lam"], q["mu"]) * Gz * Fz   # PWM = -C(e), G < 0
        pm, gm, wc = margins(L)
        rows.append(dict(A=A_, f=float(f), PM=float(pm), GM=float(gm), wc=float(wc)))
Lf = -Cz(**F) * Gz * Fz
Lp = -Cz(**Pd) * Gz * Fz
nums["margins_fixed"] = dict(FOPID=margins(Lf), PID=margins(Lp))
nums["margins_sched"] = rows
fig, ax = plt.subplots(1, 2, figsize=(W, 1.8))
for A_, ls, col in [(1.0, "-", K), (OPT["A"], "--", C[0])]:
    r = [x for x in rows if x["A"] == A_]
    ax[0].plot([x["f"] for x in r], [x["PM"] for x in r], ls=ls, color=col, label=f"A = {A_}")
    ax[1].plot([x["f"] for x in r], [x["GM"] for x in r], ls=ls, color=col)
ax[0].set_xlabel(r"far weight $w_{\mathrm{far}}$"); ax[0].set_ylabel("phase margin (°)")
ax[1].set_xlabel(r"far weight $w_{\mathrm{far}}$"); ax[1].set_ylabel("gain margin (dB)")
ax[0].legend(loc="best", fontsize=6)
ax[0].set_ylim(0, None); ax[1].set_ylim(0, None)
fig.tight_layout()
save(fig, "fig08_margins")

# ---------------------------------------------------------------- robustness: plant variations
rob = {}
for lab, pv in [("nominal", None), ("kP -20%", dict(kP=0.096)), ("kP +20%", dict(kP=0.144)),
                ("Rconv x1.5", dict(Rconv=v.P["Rconv"] * 1.5)), ("delay x2", dict(td=2 * v.P["td"]))]:
    rob[lab] = {}
    for name, mk in [("PID", lambda: v.FOPID(**Pd)), ("FOPID", lambda: v.FOPID(**F)), ("FOPID-AS", make_as)]:
        o = run(mk(), p=pv)
        e = o["Tref"] - o["Tc"]
        rob[lab][name] = dict(IAE=float(np.sum(np.abs(e)) * v.TS), undershoot=float(max(0, e[t > 5].max())))
nums["robust"] = rob

json.dump(nums, open(HERE + "/ch08_numbers.json", "w"), indent=1, default=float)
print(json.dumps({k: nums[k] for k in nums if k != "margins_sched"}, indent=1, default=float))
for r in rows[::4]:
    print(r)
