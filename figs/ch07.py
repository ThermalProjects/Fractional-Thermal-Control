"""Figures and numbers for Chapter 7 -- Bringing the Controller into the Microcontroller."""
import json
import numpy as np
from style import plt, save, W, K, G1, G2, C, simlabel
from style import HERE
import vplant as v

nums = {}
F = v.TUNE_FOPID

# ---------------------------------------------------------------- M sweep: cost and RAM
Ms = [2, 5, 10, 15, 20, 30, 50, 75, 100, 150, 200]
cost, iae, tailI = [], [], []
for M in Ms:
    o = v.simulate(v.FOPID(**F, M=M), N=3001)
    cost.append(v.pso_cost(o))
    iae.append(v.metrics(o)["IAE"])
    g = np.abs(v.gl_weights(1 - F["lam"], 20000))
    tailI.append(g[M + 1:].sum() / g[1:].sum())
ram = [12 * (M + 1) for M in Ms]
nums["Msweep"] = dict(M=Ms, cost=cost, IAE=iae, tailI=tailI, ram=ram)
fig, ax = plt.subplots(1, 2, figsize=(W, 2.0))
ax[0].plot(Ms, np.array(cost) / cost[Ms.index(30)], "o-", color=K, ms=3)
ax[0].axvline(30, color=C[1], lw=0.8, ls="--")
ax[0].set_xscale("log")
ax[0].set_xlabel("memory length $M$")
ax[0].set_ylabel(r"$J(M)/J(30)$")
ax[0].set_title("closed-loop cost", fontsize=7)
ax[1].plot(Ms, 100 * np.array(tailI), "s-", color=K, ms=3, label="weight beyond $M$ (%)")
ax[1].set_xscale("log")
ax[1].set_xlabel("memory length $M$")
ax[1].set_ylabel("weight dropped (%)")
ax2 = ax[1].twinx()
ax2.plot(Ms, ram, "^--", color=C[0], ms=3)
ax2.set_ylabel("RAM (bytes)", color=C[0])
ax2.grid(False)
ax[1].axvline(30, color=C[1], lw=0.8, ls="--")
ax[1].set_title("truncation and memory", fontsize=7)
fig.tight_layout()
simlabel(ax[0])
save(fig, "fig07_memory")


# ---------------------------------------------------------------- Ts sweep (controller period = n * 0.1 s)
def simulate_Ts(n, M, N=3001, noise=0.0, seed=0):
    rng = np.random.default_rng(seed)
    Ts = n * v.TS
    c = v.FOPID(**F, M=M, Ts=Ts)
    pl = v.Plant(23.7)
    fl = v.MedEMA(pl.Tc, alpha=min(0.2 * n, 1.0))
    out = dict(t=np.zeros(N), Tc=np.zeros(N), pwm=np.zeros(N), Tref=np.full(N, 18.2))
    pwm = 0
    for k in range(N):
        if k % n == 0:
            meas = pl.Tc + (rng.normal(0, noise) if noise else 0)
            Tf = fl(meas)
            pwm = c.step(18.2 - Tf, 18.2)
        pl.step(pwm, 23.7)
        out["t"][k] = k * v.TS
        out["Tc"][k] = pl.Tc
        out["pwm"][k] = pwm
    return out


Tsweep = {}
fig, ax = plt.subplots(figsize=(W, 2.1))
for n, ls in [(1, "-"), (2, "--"), (5, "-."), (10, ":")]:
    Ts = n * v.TS
    for M, keep in [(30, "fixed M"), (max(3, int(round(3.0 / Ts))), "horizon 3 s")]:
        o = simulate_Ts(n, M)
        Tsweep[f"Ts={Ts:.1f},M={M}"] = dict(cost=v.pso_cost(o), **{k: float(val) for k, val in v.metrics(o).items()})
    o = simulate_Ts(n, max(3, int(round(3.0 / Ts))))
    ax.plot(o["t"], o["Tc"], ls=ls, color=K, label=fr"$T_s$ = {Ts:.1f} s, $M$ = {max(3, int(round(3.0 / Ts)))}")
ax.axhline(18.2, color=C[1], lw=0.7, ls=":")
ax.set_xlim(0, 80)
ax.set_xlabel("time (s)")
ax.set_ylabel(r"$T_{\mathrm{sup}}$ (°C)")
ax.legend(loc="upper right")
simlabel(ax)
save(fig, "fig07_ts")
nums["Tsweep"] = Tsweep

# ---------------------------------------------------------------- filter speed in the loop
fil = {}
fig, ax = plt.subplots(figsize=(W, 2.0))
for a, ls in [(0.2, "-"), (0.05, "--"), (0.02, "-."), (0.005, ":")]:
    o = v.simulate(v.FOPID(**F), N=3001, alpha=a, noise=0.03, seed=2)
    m = v.metrics(o)
    fil[str(a)] = dict(IAE=m["IAE"], undershoot=m["undershoot"], settle=m["settle"])
    ax.plot(o["t"], o["Tc"], ls=ls, color=K, lw=0.8,
            label=fr"$\alpha$ = {a} ($\tau_f\approx$ {-0.1 / np.log(1 - a):.1f} s)")
ax.axhline(18.2, color=C[1], lw=0.7, ls=":")
ax.set_xlim(0, 200)
ax.set_xlabel("time (s)")
ax.set_ylabel(r"true $T_{\mathrm{sup}}$ (°C)")
ax.legend(loc="upper right")
simlabel(ax)
save(fig, "fig07_filter")
nums["filter"] = fil

# ---------------------------------------------------------------- ramp-aware anti-windup
aw = {}
for r in [8, 2, 1]:
    for mode in ["sat", "ramp"]:
        o = v.simulate(v.FOPID(**F, antiwindup=mode, ramp=r), N=3001)
        m = v.metrics(o)
        aw[f"ramp{r}_{mode}"] = dict(IAE=m["IAE"], undershoot=m["undershoot"], settle=m["settle"])
nums["aw"] = aw

# ---------------------------------------------------------------- float32 vs float64 replay
o = v.simulate(v.FOPID(**F), N=3001, noise=0.03, seed=9)
e = o["e"]
c64 = v.FOPID(**F, dtype=np.float64)
c32 = v.FOPID(**F, dtype=np.float32)
u64, u32, p64, p32 = [], [], [], []
for k in range(e.size):
    p64.append(c64.step(e[k], 18.2)); u64.append(c64.u)
    p32.append(c32.step(np.float32(e[k]), np.float32(18.2))); u32.append(c32.u)
du = np.array(u32) - np.array(u64)
dp = np.array(p32) - np.array(p64)
nums["float"] = dict(max_abs_du=float(np.max(np.abs(du))), rms_du=float(np.sqrt(np.mean(du ** 2))),
                     n_pwm_diff=int(np.sum(dp != 0)), max_dp=float(np.max(np.abs(dp))))
fig, ax = plt.subplots(2, 1, figsize=(W, 2.3), sharex=True)
ax[0].plot(o["t"], u64, color=K, lw=0.8, label="request, 64-bit (computer)")
ax[0].plot(o["t"], u32, color=C[1], lw=0.6, ls="--", label="request, 32-bit (microcontroller)")
ax[0].set_ylabel("request (counts)")
ax[0].set_ylim(-50, 300)
ax[0].legend(loc="upper right")
ax[1].plot(o["t"], du, color=K, lw=0.5)
ax[1].set_ylabel(r"$\Delta u$ (counts)")
ax[1].set_xlabel("time (s)")
ax[1].ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
simlabel(ax[1])
save(fig, "fig07_float")
json.dump(nums, open(HERE + "/ch07_numbers.json", "w"), indent=1, default=float)
print(json.dumps(nums, indent=1, default=float))
