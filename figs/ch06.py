"""Figures and numbers for Chapter 6 -- Tuning and a Fair Comparison."""
import json, sys, os
import numpy as np
from style import plt, save, W, K, G1, G2, C, simlabel
from style import HERE
import vplant as v
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "vplant"))
import pso_v2_port as pp

nums = {}
cache = HERE + "/ch06_pso.json"
if os.path.exists(cache):
    R = json.load(open(cache))
else:
    R = {}
    lb5, ub5 = [5, 0.1, 0.0, 0.7, 0.7], [50, 3.0, 5.0, 1.3, 1.3]
    for tag, (l, u) in {"std": (lb5, ub5), "wide": ([5, 0.1, 0.0, 0.7, 0.7], [100, 6.0, 5.0, 1.3, 1.3])}.items():
        xF, fF, hF = pp.pso(pp.cost, l, u, seed=42)
        xP, fP, hP = pp.pso(pp.cost, l[:3], u[:3], seed=42)
        R[tag] = dict(FOPID=dict(x=list(map(float, xF)), f=float(fF), h=list(map(float, hF))),
                      PID=dict(x=list(map(float, xP)), f=float(fP), h=list(map(float, hP))))
        print(tag, R[tag]["FOPID"]["x"], R[tag]["FOPID"]["f"], R[tag]["PID"]["x"], R[tag]["PID"]["f"])
    # second seed for reproducibility check
    xF, fF, hF = pp.pso(pp.cost, lb5, ub5, seed=7)
    R["seed7"] = dict(x=list(map(float, xF)), f=float(fF))
    json.dump(R, open(cache, "w"))
nums["pso"] = {k: (dict(FOPID=dict(x=R[k]["FOPID"]["x"], f=R[k]["FOPID"]["f"], it=len(R[k]["FOPID"]["h"])),
                        PID=dict(x=R[k]["PID"]["x"], f=R[k]["PID"]["f"], it=len(R[k]["PID"]["h"])))
                   if k != "seed7" else R[k]) for k in R}

# ---------------------------------------------------------------- PSO convergence
fig, ax = plt.subplots(figsize=(W, 1.9))
ax.plot(R["std"]["FOPID"]["h"], color=K, label="FOPID (5 parameters)")
ax.plot(R["std"]["PID"]["h"], color=C[1], ls="--", label="PID (3 parameters)")
ax.set_xlabel("iteration")
ax.set_ylabel("best cost $J$")
ax.legend(loc="upper right")
simlabel(ax)
save(fig, "fig06_pso")

# ---------------------------------------------------------------- cost map (Ki, Kd) at the FOPID optimum
xF = np.array([v.TUNE_FOPID[k] for k in ["Kp", "Ki", "Kd", "lam", "mu"]])
Ki = np.linspace(0.1, 6.0, 45)
Kd = np.linspace(0.0, 5.0, 41)
J = np.zeros((Kd.size, Ki.size))
for a, kd in enumerate(Kd):
    for b, ki in enumerate(Ki):
        J[a, b] = pp.cost([xF[0], ki, kd, xF[3], xF[4]])
fig, ax = plt.subplots(figsize=(W, 2.5))
Jr = J / pp.cost(list(xF))
im = ax.pcolormesh(Ki, Kd, np.log10(Jr), cmap="Greys", shading="auto", vmin=0, vmax=np.log10(5))
cs = ax.contour(Ki, Kd, Jr, levels=[1.02, 1.05, 1.2, 1.5, 2, 3], colors=[K], linewidths=0.6)
ax.clabel(cs, fmt=lambda x: f"{x:g}", fontsize=6)
ax.plot(xF[1], xF[2], "*", color=C[1], ms=8)
ax.axvline(3.0, color=C[0], ls="--", lw=0.8)
ax.text(3.08, 4.6, "search bound", color=C[0], fontsize=6.5)
cb = fig.colorbar(im, ax=ax)
cb.set_label(r"$\log_{10}(J/J^*)$")
ax.set_xlabel("$K_i$")
ax.set_ylabel("$K_d$")
ax.grid(False)
save(fig, "fig06_costmap")
nums["plateau_frac_5pct"] = float(np.mean(J / pp.cost(list(xF)) < 1.05))
nums["Jstar_map"] = float(pp.cost(list(xF)))
nums["Kd_slice"] = {str(kd): float(pp.cost([xF[0], xF[1], kd, xF[3], xF[4]])) for kd in [0.0, 0.5, 1.0, 1.75, 2.5, 4.0, 5.0]}

# ---------------------------------------------------------------- robustness band
scen = {"nominal": {}, "gain -30 %": dict(kP=0.084), "gain +30 %": dict(kP=0.156), "delay x2": dict(td=1.8),
        "ambient 28 °C": "amb28", "noise 0.05 °C": "noise"}
rob = {}
fig, ax = plt.subplots(1, 2, figsize=(W, 2.2), sharey=True)
for i, (name, x) in enumerate([("PID", v.TUNE_PID), ("FOPID", v.TUNE_FOPID)]):
    rob[name] = {}
    traces = []
    for s, p in scen.items():
        kw = dict(N=1801)
        if p == "amb28":
            kw.update(Tamb=28.0, Tref=18.2)
        elif p == "noise":
            kw.update(noise=0.05, seed=3)
        elif p:
            kw.update(p=p)
        o = v.simulate(v.FOPID(**x), **kw)
        m = v.metrics(o)
        rob[name][s] = dict(IAE=m["IAE"], undershoot=m["undershoot"], settle=m["settle"], cost=v.pso_cost(dict(o, Tref=o["Tref"], t=o["t"], Tc=o["Tc"][:3001])))
        traces.append((o["Tc"] - o["Tref"]) / (o["Tc"][0] - o["Tref"][0]))
    tr = np.array(traces)
    a = ax[i]
    a.fill_between(o["t"], tr.min(0), tr.max(0), color="0.8", label="envelope of 6 scenarios")
    a.plot(o["t"], tr[0], color=K, lw=0.9, label="nominal")
    a.axhline(0, color=G2, lw=0.6)
    a.set_xlim(0, 60)
    a.set_ylim(-0.05, 1.05)
    a.set_title(name)
    a.set_xlabel("time (s)")
ax[0].set_ylabel("normalised error")
ax[1].legend(loc="upper right")
for i, (name, x) in enumerate([("FOPID", v.TUNE_FOPID)]):
    pass
simlabel(ax[1])
save(fig, "fig06_robust")
nums["robust"] = rob

# widened-bound solution under delay x2
xw = R["wide"]["FOPID"]["x"]
for lab, x in [("std", R["std"]["FOPID"]["x"]), ("wide", xw)]:
    o = v.simulate(v.FOPID(*x), N=3001, p=dict(td=1.8))
    m = v.metrics(o)
    on = v.simulate(v.FOPID(*x), N=3001)
    mn = v.metrics(on)
    o2 = v.simulate(v.FOPID(*x), N=3001, noise=0.05, seed=3)
    nums[f"bounds_{lab}"] = dict(x=x, nominal_undershoot=mn["undershoot"], delay2_undershoot=m["undershoot"],
                                 delay2_settle=m["settle"], pwm_std_noise=float(np.std(o2["pwm"][1500:])))

# ---------------------------------------------------------------- fair comparison
def scen_runs(ctrl_factory, noise=0.03):
    out = {}
    # S1: tuning scenario (from ambient to 18.2 at 23.7), 300 s
    o = v.simulate(ctrl_factory(), N=3001, noise=noise, seed=11)
    out["S1"] = o
    # S2: unseen reference (from ambient to 14.0 at 23.7)
    out["S2"] = v.simulate(ctrl_factory(), N=3001, Tref=14.0, noise=noise, seed=12)
    # S3: disturbance: 2 W on the cold face for 60 s at steady state
    N = 3001
    t = np.arange(N) * v.TS
    Qd = np.where((t >= 150) & (t < 210), 2.0, 0.0)
    out["S3"] = v.simulate(ctrl_factory(), N=N, Qd=Qd, noise=noise, seed=13)
    return out


ctrls = {"Hysteresis": lambda: v.Hysteresis(h=0.3, on=255),
         "PID": lambda: v.FOPID(**v.TUNE_PID),
         "FOPID": lambda: v.FOPID(**v.TUNE_FOPID)}
runs = {n: scen_runs(f) for n, f in ctrls.items()}
tab = {}
for n, rr in runs.items():
    tab[n] = {}
    for s, o in rr.items():
        if s == "S3":
            m = v.metrics(o, t_from=150.0)
            seg = (o["t"] >= 150)
            tab[n][s] = dict(IAE=m["IAE"], peak=float(np.max(np.abs(o["e"][seg]))), duty=m["duty"])
        else:
            m = v.metrics(o)
            ss = o["t"] >= 200
            tab[n][s] = dict(IAE=m["IAE"], ITAE=m["ITAE"], undershoot=m["undershoot"], settle=m["settle"],
                             mean_abs_ss=float(np.mean(np.abs(o["e"][ss]))), duty=m["duty"],
                             pwm_std_ss=float(np.std(o["pwm"][ss])), J=v.pso_cost(o) if s == "S1" else None)
nums["fair"] = tab

fig, ax = plt.subplots(3, 1, figsize=(W, 4.0), sharex=False)
for n, ls, col in [("Hysteresis", ":", G1), ("PID", "--", C[1]), ("FOPID", "-", K)]:
    o = runs[n]["S1"]
    ax[0].plot(o["t"], o["Tc"], ls=ls, color=col, lw=0.9, label=n)
    o = runs[n]["S2"]
    ax[1].plot(o["t"], o["Tc"], ls=ls, color=col, lw=0.9)
    o = runs[n]["S3"]
    ax[2].plot(o["t"], o["Tc"], ls=ls, color=col, lw=0.9)
ax[0].set_title("S1: tuning scenario, $T_{\\mathrm{ref}}$ = 18.2 °C", fontsize=7)
ax[1].set_title("S2: unseen reference, $T_{\\mathrm{ref}}$ = 14.0 °C", fontsize=7)
ax[2].set_title("S3: 2 W load on the cold face from 150 to 210 s", fontsize=7)
ax[0].set_xlim(0, 120); ax[1].set_xlim(0, 120); ax[2].set_xlim(120, 300)
ax[0].set_ylim(15.5, 24); ax[1].set_ylim(11.5, 24); ax[2].set_ylim(15.5, 19.5)
for a in ax:
    a.set_ylabel("°C")
ax[2].set_xlabel("time (s)")
ax[0].legend(loc="upper right", ncol=3)
simlabel(ax[2])
fig.tight_layout()
save(fig, "fig06_compare")

# metrics bar chart (S1, S2, S3) for PID vs FOPID relative
fig, ax = plt.subplots(figsize=(W, 1.9))
labels = ["IAE S1", "ITAE S1", "IAE S2", "ITAE S2", "IAE S3", "peak S3", "J (S1)"]
pid = [tab["PID"]["S1"]["IAE"], tab["PID"]["S1"]["ITAE"], tab["PID"]["S2"]["IAE"], tab["PID"]["S2"]["ITAE"],
       tab["PID"]["S3"]["IAE"], tab["PID"]["S3"]["peak"], tab["PID"]["S1"]["J"]]
fo = [tab["FOPID"]["S1"]["IAE"], tab["FOPID"]["S1"]["ITAE"], tab["FOPID"]["S2"]["IAE"], tab["FOPID"]["S2"]["ITAE"],
      tab["FOPID"]["S3"]["IAE"], tab["FOPID"]["S3"]["peak"], tab["FOPID"]["S1"]["J"]]
rel = 100 * (np.array(fo) - np.array(pid)) / np.array(pid)
nums["rel_fopid_vs_pid"] = dict(zip(labels, rel.tolist()))
ax.bar(range(len(labels)), rel, color=[K if r < 0 else C[1] for r in rel], width=0.6)
ax.axhline(0, color=K, lw=0.6)
ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels, fontsize=6.5)
ax.set_ylabel("FOPID vs PID (%)")
for i, r in enumerate(rel):
    ax.text(i, r + (0.4 if r >= 0 else -0.4), f"{r:+.1f}", ha="center", va="bottom" if r >= 0 else "top", fontsize=6.5)
ax.set_ylim(min(rel.min() - 3, -3), max(rel.max() + 3, 3))
simlabel(ax)
save(fig, "fig06_metrics")
json.dump(nums, open(HERE + "/ch06_numbers.json", "w"), indent=1, default=float)
print(json.dumps({k: nums[k] for k in ["fair", "rel_fopid_vs_pid", "bounds_std", "bounds_wide", "plateau_frac_5pct"]}, indent=1, default=float))
print(json.dumps(nums["robust"], indent=1, default=float))
print(nums["pso"]["seed7"], nums["pso"]["std"]["FOPID"]["f"], nums["pso"]["std"]["PID"]["f"], nums["pso"]["wide"]["FOPID"])
