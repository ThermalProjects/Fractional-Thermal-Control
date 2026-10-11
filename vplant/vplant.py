"""
vplant.py -- the book's virtual plant.

Three-node RC thermal network (cold face T_c, interior T_m, hot face T_h) with
transport delay and input gain k_P, plus the exact control law of the reference
firmware: median-of-3 + EMA filter, Grunwald-Letnikov FOPID (PID when
lambda = mu = 1), conditional anti-windup aware of the PWM ramp, symmetric ramp
of 8 counts per cycle, saturation 0..255 and integer truncation, fan rule.

Sign convention (book):  e = T_ref - T_sup.  The plant cools when PWM rises
(negative gain), so the controller output is   PWM = -C(e) + K_ff * T_ref.

Everything runs at T_s = 0.1 s.
"""
# Hot side (R2, Rconv, Ch) recalibrated against the hot-face probe of open-loop
# staircases, so that T_h respects the datasheet limit of the module (Chap. 9).
import json
import numpy as np

TS = 0.1
P = dict(R1=4.022803, R2=6.229903, Rconv=1.712813, Cc=10.305737, Cm=824.187344,
         Ch=67.943566, phi=0.302829, td=0.894863, kP=0.12, gam=0.0)
TUNE_FOPID = dict(Kp=49.99329999, Ki=2.99575159, Kd=1.75472605, lam=0.78085638, mu=1.24022096)
TUNE_PID = dict(Kp=49.99709145, Ki=1.59761882, Kd=1.23985127, lam=1.0, mu=1.0)


def load_params(path):
    """Read a parameter file (JSON) and return a dict usable as Plant(p=...)."""
    with open(path) as f:
        d = json.load(f)
    return {k: float(d[k]) for k in P if k in d}


def gl_weights(a, M=30):
    """Signed Grunwald-Letnikov weights g_0..g_M of order a."""
    g = np.empty(M + 1)
    g[0] = 1.0
    for j in range(1, M + 1):
        g[j] = g[j - 1] * ((j - 1) - a) / j
    return g


# ----------------------------------------------------------------- plant
class Plant:
    def __init__(self, Tamb=23.7, p=None, T0=None):
        self.p = dict(P)
        if p:
            self.p.update(p)
        T0 = Tamb if T0 is None else T0
        self.Tc = self.Tm = self.Th = T0
        if T0 != Tamb:          # start hot face / interior at ambient
            self.Tm = self.Th = Tamb
        self.delay = max(1, int(round(self.p["td"] / TS)))
        self.buf = []

    def step(self, pwm, Tamb, Qd=0.0, rconv_mult=1.0):
        """Advance one sampling period. pwm is the value applied now; the
        plant sees it after the transport delay. Qd: extra heat load on the
        cold face [W, model units]. rconv_mult: heat-sink resistance multiplier
        (fan failure > 1)."""
        p = self.p
        self.buf.append(float(pwm))
        u = self.buf.pop(0) if len(self.buf) > self.delay else self.buf[0]
        Pe = p["kP"] * u
        QJ = p["gam"] * u * u          # Joule heat released inside the module
        Pc, Ph = p["phi"] * Pe, (1 - p["phi"]) * Pe
        Rc = p["Rconv"] * rconv_mult
        dTc = (self.Tm - self.Tc) / (p["R1"] * p["Cc"]) - Pc / p["Cc"] + Qd / p["Cc"]
        dTm = (self.Tc - self.Tm) / (p["R1"] * p["Cm"]) + (self.Th - self.Tm) / (p["R2"] * p["Cm"]) + QJ / p["Cm"]
        dTh = (self.Tm - self.Th) / (p["R2"] * p["Ch"]) + Ph / p["Ch"] - (self.Th - Tamb) / (Rc * p["Ch"])
        self.Tc += TS * dTc
        self.Tm += TS * dTm
        self.Th += TS * dTh
        self.Tc = min(max(self.Tc, -10.0), 120.0)
        return self.Tc, self.Tm, self.Th


# ----------------------------------------------------------------- filter
class MedEMA:
    def __init__(self, x0, alpha=0.2):
        self.b = [x0, x0, x0]
        self.y = x0
        self.alpha = alpha

    def __call__(self, x):
        self.b = [self.b[1], self.b[2], x]
        a, b, c = self.b
        m = max(min(a, b), min(max(a, b), c))
        self.y += self.alpha * (m - self.y)
        return self.y


# ----------------------------------------------------------------- controller
class FOPID:
    """Firmware control law in the book convention e = T_ref - T_sup.
    dtype=np.float32 mimics the 32-bit float arithmetic of the MCU."""

    def __init__(self, Kp, Ki, Kd, lam=1.0, mu=1.0, Kff=0.0, M=30, ramp=8,
                 antiwindup="ramp", dtype=np.float64, Ts=TS):
        self.Kp, self.Ki, self.Kd, self.Kff = Kp, Ki, Kd, Kff
        self.M, self.ramp, self.aw, self.dt, self.Ts = M, ramp, antiwindup, dtype, Ts
        self.eh = np.zeros(M + 1, dtype=dtype)
        self.Iacc = dtype(0.0)
        self.pwm = 0
        self.pwmT = 0
        self.set_orders(lam, mu)

    def set_orders(self, lam, mu):
        self.lam, self.mu = lam, mu
        self.wV = gl_weights(1 - lam, self.M).astype(self.dt)
        self.wD = gl_weights(mu, self.M).astype(self.dt)
        self.TsV = self.dt(self.Ts ** (1 - lam))
        self.TsD = self.dt(self.Ts ** mu)

    def step(self, e, Tref=0.0):
        d = self.dt
        e = d(e)
        self.eh[1:] = self.eh[:-1]
        self.eh[0] = e
        v = d(np.dot(self.wV, self.eh)) / self.TsV
        S = d(np.dot(self.wD, self.eh))
        D = (S / (d(1) + abs(S))) / self.TsD
        # firmware works with -e: cooling demand grows when T_sup > T_ref
        uP = -self.Kp * e
        uD = -self.Kd * D
        vv = -v
        ff = self.Kff * Tref
        uTry = uP + self.Ki * (self.Iacc + d(self.Ts) * vv) + uD + ff
        if self.aw == "ramp":
            hi = ((uTry > 255) or (self.pwm < self.pwmT)) and (vv > 0)
            lo = ((uTry < 0) or (self.pwm > self.pwmT)) and (vv < 0)
        elif self.aw == "sat":
            hi = (uTry > 255) and (vv > 0)
            lo = (uTry < 0) and (vv < 0)
        else:
            hi = lo = False
        if not hi and not lo:
            self.Iacc = d(self.Iacc + d(self.Ts) * vv)
            if self.Ki > 1e-6 and self.aw != "none":
                lim = 255.0 / self.Ki
                self.Iacc = d(min(max(self.Iacc, -lim), lim))
        u = uP + self.Ki * self.Iacc + uD + ff
        self.u = float(u)
        self.pwmT = int(min(max(int(u), 0), 255))
        if self.ramp:
            if self.pwm < self.pwmT:
                self.pwm = min(self.pwm + self.ramp, self.pwmT)
            elif self.pwm > self.pwmT:
                self.pwm = max(self.pwm - self.ramp, self.pwmT)
        else:
            self.pwm = self.pwmT
        return self.pwm


class Hysteresis:
    def __init__(self, h=0.3, on=255, ramp=8):
        self.h, self.on, self.ramp = h, on, ramp
        self.state = 0
        self.pwm = 0

    def step(self, e, Tref=0.0):
        if e < -self.h:       # surface warmer than ref + h -> cool
            self.state = 1
        elif e > self.h:      # surface colder than ref - h -> stop
            self.state = 0
        tgt = self.on if self.state else 0
        if self.pwm < tgt:
            self.pwm = min(self.pwm + self.ramp, tgt)
        elif self.pwm > tgt:
            self.pwm = max(self.pwm - self.ramp, tgt)
        return self.pwm


class Supervisor:
    """Rule-based adaptive supervision (FOPID-AS): three regions by |e|,
    linear interpolation between the precision values (|e| <= e_low) and the
    far values (|e| >= e_high), hysteresis on leaving the precision region,
    aggressiveness factor A on K_p when e > 0 (surface colder than ref)."""

    def __init__(self, ctrl, near, far, e_low=0.3, e_high=3.0, hyst=0.1, A=1.2,
                 min_dorder=0.0):
        self.c, self.near, self.far = ctrl, near, far
        self.e_low, self.e_high, self.hyst, self.A = e_low, e_high, hyst, A
        self.inside = False
        self.f = 1.0

    def schedule(self, e):
        a = abs(e)
        if self.inside:
            if a > self.e_low + self.hyst:
                self.inside = False
        elif a < self.e_low:
            self.inside = True
        f = 0.0 if self.inside else min(max((a - self.e_low) / (self.e_high - self.e_low), 0.0), 1.0)
        self.f = f
        q = {k: self.near[k] + f * (self.far[k] - self.near[k]) for k in self.near}
        return q

    def step(self, e, Tref=0.0):
        q = self.schedule(e)
        c = self.c
        # A only outside the precision region: inside it, the sign of e is
        # mostly noise and a sign-dependent gain would chatter
        c.Kp = q["Kp"] * (self.A if (e > 0 and not self.inside) else 1.0)
        c.Kff = q.get("Kff", 0.0)
        if abs(q["lam"] - c.lam) > 1e-4 or abs(q["mu"] - c.mu) > 1e-4:
            c.set_orders(q["lam"], q["mu"])
        return c.step(e, Tref)


def magnus_dew(T, RH, a=17.67, b=243.5):
    g = np.log(np.maximum(RH, 1e-3) / 100.0) + a * T / (b + T)
    return b * g / (a - g)


def p_sat(T):
    """Saturation vapour pressure [hPa], Magnus form."""
    return 6.112 * np.exp(17.67 * T / (T + 243.5))


# ----------------------------------------------------------------- simulation
def simulate(ctrl, N=3001, Tref=18.2, Tamb=23.7, Qd=None, rconv=None, noise=0.0,
             seed=0, p=None, T0=None, filt=True, alpha=0.2, ctrl_on=0):
    """Closed loop. Tref, Tamb may be scalars or arrays of length N.
    Returns dict of arrays. Sample k: measure -> filter -> control -> apply -> plant."""
    rng = np.random.default_rng(seed)
    arr = lambda x: np.full(N, float(x)) if np.isscalar(x) else np.asarray(x, float)
    Tref, Tamb = arr(Tref), arr(Tamb)
    Qd = np.zeros(N) if Qd is None else arr(Qd)
    rconv = np.ones(N) if rconv is None else arr(rconv)
    pl = Plant(Tamb[0], p, T0)
    fl = MedEMA(pl.Tc, alpha)
    out = {k: np.zeros(N) for k in ["t", "Tc", "Tm", "Th", "Tsup", "pwm", "e", "Tref", "Tamb", "fan", "u"]}
    for k in range(N):
        meas = pl.Tc + (rng.normal(0, noise) if noise > 0 else 0.0)
        Ts_ = fl(meas) if filt else meas
        e = Tref[k] - Ts_
        pwm = ctrl.step(e, Tref[k]) if k >= ctrl_on else 0
        pl.step(pwm, Tamb[k], Qd[k], rconv[k])
        out["t"][k] = k * TS
        out["Tc"][k], out["Tm"][k], out["Th"][k] = pl.Tc, pl.Tm, pl.Th
        out["Tsup"][k], out["pwm"][k], out["e"][k] = Ts_, pwm, e
        out["Tref"][k], out["Tamb"][k] = Tref[k], Tamb[k]
        out["fan"][k] = 255 if pwm > 20 else 0
        out["u"][k] = getattr(ctrl, "u", np.nan) if hasattr(ctrl, "u") else np.nan
    return out


def open_loop(pwm, Tamb=23.7, p=None, noise=0.0, seed=0, rconv=None, Qd=None):
    rng = np.random.default_rng(seed)
    N = len(pwm)
    Tamb = np.full(N, float(Tamb)) if np.isscalar(Tamb) else np.asarray(Tamb, float)
    rconv = np.ones(N) if rconv is None else rconv
    Qd = np.zeros(N) if Qd is None else Qd
    pl = Plant(Tamb[0], p)
    o = {k: np.zeros(N) for k in ["t", "Tc", "Tm", "Th", "Tc_meas", "Th_meas"]}
    for k in range(N):
        pl.step(pwm[k], Tamb[k], Qd[k], rconv[k])
        o["t"][k] = k * TS
        o["Tc"][k], o["Tm"][k], o["Th"][k] = pl.Tc, pl.Tm, pl.Th
        o["Tc_meas"][k] = pl.Tc + (rng.normal(0, noise) if noise else 0)
        o["Th_meas"][k] = pl.Th + (rng.normal(0, noise) if noise else 0)
    return o


# ----------------------------------------------------------------- metrics
def metrics(o, band=0.5, t_from=0.0):
    t, Tc, Tref = o["t"], o["Tc"], o["Tref"]
    e = Tref - Tc
    m = t >= t_from
    dt = TS
    iae = np.sum(np.abs(e[m])) * dt
    itae = np.sum((t[m] - t_from) * np.abs(e[m])) * dt
    under = max(0.0, np.max(Tref[m] - Tc[m])) if np.any(m) else 0
    inb = np.abs(e) < band
    out = np.where(~inb & m)[0]
    settle = (t[out[-1] + 1] - t_from) if len(out) and out[-1] < len(t) - 1 else (0.0 if not len(out) else np.nan)
    last = t >= t[-1] - 2.0
    jss = abs(np.mean(Tc[-21:]) - Tref[-1])
    return dict(IAE=iae, ITAE=itae, undershoot=under, settle=settle, ess=jss,
                duty=np.mean(o["pwm"][m]) / 255 * 100,
                ripple=np.ptp(Tc[t >= t[-1] - 60]))


def pso_cost(o):
    """The book's single tuning cost (= MATLAB cost_pso_v2)."""
    t, Tc, Tref = o["t"], o["Tc"], o["Tref"][-1]
    if not np.all(np.isfinite(Tc)) or np.any(np.abs(Tc) > 200):
        return 1e9
    e = Tc - Tref
    Jss = abs(np.mean(Tc[-21:]) - Tref)
    under = max(0.0, Tref - Tc.min())
    inb = np.abs(e) < 0.5
    out = np.where(~inb)[0]
    if len(out) == 0:
        settle = t[np.argmax(inb)]
    elif out[-1] < len(t) - 1:
        settle = t[out[-1] + 1]
    else:
        settle = t[-1]
    ITAE = np.sum(t * np.abs(e)) * TS
    return ITAE + 1000 * Jss ** 2 + 500 * under + settle


# ----------------------------------------------------------------- book design values
# Supervisor regions designed around the single tuning (Chapter 8)
SUP_NEAR = dict(Kp=50.0, lam=1.00, mu=1.24)   # precision region, |e| <= e_low
SUP_FAR = dict(Kp=80.0, lam=0.78, mu=1.24)    # far region, |e| >= e_high
SUP_OPTS = dict(e_low=0.3, e_high=3.0, hyst=0.1, A=1.2)


# ----------------------------------------------------------------- observer (Chapter 9)
def model_matrices(p=None):
    """Discrete (forward Euler, T_s) three-node network: x = (Tc, Tm, Th),
    x+ = Ad x + bd u(k-d) + ed Tamb."""
    q = dict(P)
    if p:
        q.update(p)
    A = np.array([[-1 / (q["R1"] * q["Cc"]), 1 / (q["R1"] * q["Cc"]), 0],
                  [1 / (q["R1"] * q["Cm"]), -1 / (q["R1"] * q["Cm"]) - 1 / (q["R2"] * q["Cm"]), 1 / (q["R2"] * q["Cm"])],
                  [0, 1 / (q["R2"] * q["Ch"]), -1 / (q["R2"] * q["Ch"]) - 1 / (q["Rconv"] * q["Ch"])]])
    b = np.array([-q["phi"] * q["kP"] / q["Cc"], 0, (1 - q["phi"]) * q["kP"] / q["Ch"]])
    e = np.array([0, 0, 1 / (q["Rconv"] * q["Ch"])])
    return np.eye(3) + TS * A, TS * b, TS * e


def observer_gain(Cm, tau, p=None):
    """Luenberger gain: observer poles at continuous time constants tau, 1.25 tau, 1.5 tau."""
    from scipy.signal import place_poles
    Ad, _, _ = model_matrices(p)
    zp = np.exp(-TS / (tau * np.array([1.0, 1.25, 1.5])))
    return place_poles(Ad.T, np.atleast_2d(Cm).T, zp).gain_matrix.T


class Observer:
    """Predict-correct Luenberger observer of the three-node network."""

    def __init__(self, x0, L, p=None, Cm=((1.0, 0, 0), (0, 0, 1.0))):
        self.Ad, self.bd, self.ed = model_matrices(p)
        self.C = np.atleast_2d(np.array(Cm, float))
        self.x = np.array(x0, float)
        self.L = np.asarray(L, float)
        q = dict(P); q.update(p or {})
        self.d = max(1, int(round(q["td"] / TS)))
        self.ubuf = [0.0] * (self.d + 1)

    def step(self, y, u, Tamb):
        inn = np.atleast_1d(y) - self.C @ self.x
        xk = self.x + self.L @ inn
        self.ubuf.append(float(u))
        ud = self.ubuf.pop(0)
        self.x = self.Ad @ xk + self.bd * ud + self.ed * Tamb
        return xk, inn
