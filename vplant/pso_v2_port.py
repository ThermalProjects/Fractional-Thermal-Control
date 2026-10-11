"""
pso_v2_port.py -- Python port of the MATLAB tuning set (sim_config_v2, simulate_fopid_v2,
cost_pso_v2, RUN_PSO_v2).  Same plant, same firmware law, same cost, same PSO settings.

Usage:
  python pso_v2_port.py check                 # reproduce the MATLAB reference costs
  python pso_v2_port.py run  [kp_ub] [ki_ub] [seed]
"""
import sys, json
import numpy as np
from numba import njit

# Virtual plant of the book (hot side recalibrated against the hot-face probe).
CFG = dict(dt=0.1, tFinal=300.0, Tamb=23.7, Tref=18.2, M=30, alpha=0.2, kP=0.12,
           R1=4.022803, R2=6.229903, Rconv=1.712813, Cc=10.305737, Cp=824.187344,
           Ch=67.943566, theta=0.302829, td=0.894863)
# Plant used by the original MATLAB tuning scripts (kept to reproduce their costs).
CFG_MATLAB = dict(CFG, R2=6.003481, Rconv=4.897706, Ch=56.890349)


def gl(a, M):
    w = np.empty(M + 1); w[0] = 1.0
    for k in range(1, M + 1):
        w[k] = w[k - 1] * ((k - 1) - a) / k
    return w


@njit(cache=True)
def _sim(Kp, Ki, Kd, wV, wD, TsV, TsD, dt, N, Tamb, Tref, alpha, kP,
         R1, R2, Rconv, Cc, Cp, Ch, theta, delay, noise):
    M = wV.shape[0] - 1
    beta = 1.0 - theta
    Tc = Tamb; Tm = Tamb; Th = Tamb
    eh = np.zeros(M + 1); Iacc = 0.0; pwm = 0; pwmT = 0
    m0 = Tc; m1 = Tc; m2 = Tc; Tf = Tc
    TL = np.zeros(N); PL = np.zeros(N)
    for k in range(N):
        m0 = m1; m1 = m2; m2 = Tc + noise[k]
        med = max(min(m0, m1), min(max(m0, m1), m2))
        Tf = Tf + alpha * (med - Tf)
        e = Tf - Tref
        for i in range(M, 0, -1):
            eh[i] = eh[i - 1]
        eh[0] = e
        sV = 0.0; sD = 0.0
        for i in range(M + 1):
            sV += wV[i] * eh[i]; sD += wD[i] * eh[i]
        v = sV / TsV
        uD = Kd * (sD / (1.0 + abs(sD))) / TsD
        uP = Kp * e
        uTry = uP + Ki * (Iacc + dt * v) + uD
        satH = ((uTry > 255.0) or (pwm < pwmT)) and (v > 0)
        satL = ((uTry < 0.0) or (pwm > pwmT)) and (v < 0)
        if (not satH) and (not satL):
            Iacc += dt * v
            if Ki > 1e-6:
                lim = 255.0 / Ki
                if Iacc > lim: Iacc = lim
                if Iacc < -lim: Iacc = -lim
        u = uP + Ki * Iacc + uD
        ui = int(u)                      # MATLAB fix(): truncation toward zero
        pwmT = min(max(ui, 0), 255)
        if pwm < pwmT:
            pwm = min(pwm + 8, pwmT)
        elif pwm > pwmT:
            pwm = max(pwm - 8, pwmT)
        PL[k] = pwm
        j = k - delay
        if j < 0: j = 0
        Pe = kP * PL[j]
        Pc = theta * Pe; Ph = beta * Pe
        dTc = (Tm - Tc) / (R1 * Cc) - Pc / Cc
        dTm = (Tc - Tm) / (R1 * Cp) + (Th - Tm) / (R2 * Cp)
        dTh = (Tm - Th) / (R2 * Ch) + Ph / Ch - (Th - Tamb) / (Rconv * Ch)
        Tc += dt * dTc; Tm += dt * dTm; Th += dt * dTh
        if Tc < -10: Tc = -10.0
        if Tc > 120: Tc = 120.0
        TL[k] = Tc
    return TL, PL


def simulate(x, cfg=CFG, noise_std=0.0, seed=0, **over):
    c = dict(cfg); c.update(over)
    Kp, Ki, Kd, lam, mu = x
    dt = c["dt"]; N = int(round(c["tFinal"] / dt)) + 1
    delay = max(1, int(round(c["td"] / dt)))
    return _sim(Kp, Ki, Kd, gl(1 - lam, c["M"]), gl(mu, c["M"]), dt ** (1 - lam), dt ** mu,
                dt, N, c["Tamb"], c["Tref"], c["alpha"], c["kP"], c["R1"], c["R2"], c["Rconv"],
                c["Cc"], c["Cp"], c["Ch"], c["theta"], delay,
                np.random.default_rng(seed).normal(0, noise_std, N) if noise_std > 0 else np.zeros(N))


def cost(x, cfg=CFG, **over):
    x = list(x)
    if len(x) == 3:
        x = x + [1.0, 1.0]
    c = dict(cfg); c.update(over)
    Tc, _ = simulate(x, c)
    N = len(Tc); dt = c["dt"]; t = np.arange(N) * dt; Tref = c["Tref"]
    if (not np.all(np.isfinite(Tc))) or np.any(np.abs(Tc) > 200):
        return 1e9
    e = Tc - Tref
    Jss = abs(np.mean(Tc[-21:]) - Tref)
    under = max(0.0, Tref - Tc.min())
    inband = np.abs(e) < 0.5
    out = np.where(~inband)[0]
    if len(out) == 0:
        settle = t[np.argmax(inband)]
    elif out[-1] < N - 1:
        settle = t[out[-1] + 1]
    else:
        settle = c["tFinal"]
    ITAE = np.sum(t * np.abs(e)) * dt
    return ITAE + 1000 * Jss ** 2 + 500 * under + settle


def pso(cfun, lb, ub, n=100, iters=300, wmax=0.9, wmin=0.4, c1=1.4962, c2=1.4962, seed=42):
    rng = np.random.default_rng(seed)
    lb = np.asarray(lb, float); ub = np.asarray(ub, float); D = len(lb)
    pos = lb + (ub - lb) * rng.random((n, D)); vel = np.zeros((n, D))
    fit = np.array([cfun(p) for p in pos])
    pb, pf = pos.copy(), fit.copy()
    gi = pf.argmin(); g, gf = pb[gi].copy(), pf[gi]
    vmax = 0.1 * (ub - lb); conv = []
    for it in range(1, iters + 1):
        w = wmax - (wmax - wmin) * it / iters
        for i in range(n):
            r1, r2 = rng.random(D), rng.random(D)
            vel[i] = w * vel[i] + c1 * r1 * (pb[i] - pos[i]) + c2 * r2 * (g - pos[i])
            vel[i] = np.clip(vel[i], -vmax, vmax)
            pos[i] = pos[i] + vel[i]
            lo = pos[i] < lb; hi = pos[i] > ub
            pos[i][lo] = lb[lo]; pos[i][hi] = ub[hi]
            vel[i][lo] *= -1; vel[i][hi] *= -1
            f = cfun(pos[i])
            if f < pf[i]:
                pf[i] = f; pb[i] = pos[i].copy()
        if pf.min() < gf:
            gi = pf.argmin(); gf = pf[gi]; g = pb[gi].copy()
        conv.append(gf)
        if it > 50 and abs(conv[-1] - conv[-51]) < 1e-6:
            break
    return g, gf, conv


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "check"
    if mode == "check":
        refs = [([27.138, 1.735, 2.791, 1.122, 1.186], 1681.5898),
                ([22.844, 2.479, 2.783, 0.960, 1.190], 2073.2342),
                ([44.739, 0.952, 1.187, 1.0, 1.0], 884.9879),
                ([49.99329999, 2.99575159, 1.75472605, 0.78085638, 1.24022096], 294.71797379),
                ([49.99709145, 1.59761882, 1.23985127], 313.66203562)]
        for x, r in refs:
            c = cost(x, CFG_MATLAB)
            print(f"{str(np.round(x, 3)):45s} python {c:10.4f}  MATLAB {r:10.4f}  diff {100*(c-r)/r:+.3f}%")
    else:
        kp_ub = float(sys.argv[2]) if len(sys.argv) > 2 else 100.0
        ki_ub = float(sys.argv[3]) if len(sys.argv) > 3 else 6.0
        seed = int(sys.argv[4]) if len(sys.argv) > 4 else 42
        lb5 = [5, 0.1, 0.0, 0.7, 0.7]; ub5 = [kp_ub, ki_ub, 5.0, 1.3, 1.3]
        xF, fF, hF = pso(cost, lb5, ub5, seed=seed)
        xP, fP, hP = pso(cost, lb5[:3], ub5[:3], seed=seed)
        res = dict(bounds=dict(lb=lb5, ub=ub5), seed=seed,
                   FOPID=dict(x=list(map(float, xF)), cost=float(fF), iters=len(hF)),
                   PID=dict(x=list(map(float, xP)), cost=float(fP), iters=len(hP)))
        print(json.dumps(res, indent=1))
        with open(f"pso_v2_kp{int(kp_ub)}_ki{int(ki_ub)}_s{seed}.json", "w") as f:
            json.dump(res, f, indent=1)
