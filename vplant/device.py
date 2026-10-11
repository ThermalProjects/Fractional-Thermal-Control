"""
device.py -- the virtual device: the virtual plant behind the same serial
interface as the reference firmware (commands and CSV columns of Appendix B).

    from device import VirtualDevice
    dev = VirtualDevice()                      # reference parameters
    dev.command("SET_TREF,18.2"); dev.command("CONTROL_ON")
    for _ in range(3000):
        row = dev.step()                       # one 100 ms period -> CSV line

It runs the same law as the firmware (median-of-3 + EMA, FOPID with GL weights,
ramp-aware anti-windup, ramp, fan rule, supervisor, observer, dew-point
reference). Ambient conditions can be scripted with dev.Tamb / dev.RH, and
faults with dev.rconv_mult (fan failure) or dev.Qd (heat load on the plate).
"""
import math
import numpy as np
import vplant as v

COLUMNS = ("t_ms,Tamb,RH,Tdew,Tref,Tsup,Th,e,pwm,fan,Kp,Ki,Kd,lambda,mu,Kff,mode,"
           "Tsup_raw,Th_raw,state,loop_us,late_ms,req,sup,wfar,wrec,Tm_hat,inn_c,inn_h,"
           "alarm,Tdew_raw,bdew,ref,amb_flag,Qdew,n_below").split(",")


class VirtualDevice:
    def __init__(self, params=None, Tamb=23.7, RH=60.0, noise=0.03, seed=0):
        self.p = dict(v.P)
        if isinstance(params, str):
            self.p.update(v.load_params(params))
        elif params:
            self.p.update(params)
        self.rng = np.random.default_rng(seed)
        self.noise = noise
        self.Tamb, self.RH = Tamb, RH
        self.Qd, self.rconv_mult = 0.0, 1.0
        self.plant = v.Plant(Tamb, self.p)
        self.fs, self.fh = v.MedEMA(self.plant.Tc), v.MedEMA(self.plant.Th)
        self.k = 0
        # reference
        self.ref_mode, self.tset, self.bT, self.bdew = "FIXED", 18.0, 0.0, 0.0
        self.tdew = None; self.tdew_raw = float("nan"); self.tref = self.tset
        # controller
        self.mode, self.control, self.pwm_open = "FOPID", False, 0
        self.h, self.pwm_on, self.hyst = 0.3, 255, False
        self.gains = dict(Kp=50.0, Ki=3.00, Kd=1.75, lam=0.78, mu=1.24, Kff=0.0)
        self.ctrl = v.FOPID(**{k: self.gains[k] for k in ("Kp", "Ki", "Kd", "lam", "mu", "Kff")})
        self.sup_on = False
        self.sup = v.Supervisor(self.ctrl, near=dict(v.SUP_NEAR, Kff=0.0), far=dict(v.SUP_FAR, Kff=0.0), **v.SUP_OPTS)
        self.fan_mode = "AUTO"; self.fan = 0
        # observer and alarm
        self.obs = v.Observer([self.plant.Tc, self.plant.Tc, self.plant.Th], v.observer_gain([[1, 0, 0], [0, 0, 1]], 5.0, self.p), self.p)
        self.band = (0.04, 0.17); self.innbuf = []; self.out_count = 0; self.alarm = False
        self.qdew = 0.0; self.nbelow = 0

    # ------------------------------------------------------------ commands
    def command(self, s):
        s = s.strip()
        arg = lambda: float(s.split(",", 1)[1])
        g = self.gains
        if s == "CONTROL_ON":
            self.control = True; self._bumpless()
        elif s == "CONTROL_OFF":
            self.control = False
        elif s.startswith("SET_PWM,"):
            self.pwm_open = int(max(0, min(255, arg()))); self.mode = "OPEN"
        elif s in ("FAN_AUTO", "FAN_ON", "FAN_OFF"):
            self.fan_mode = s[4:]
        elif s.startswith("MODE,"):
            self.mode = s[5:]; self._bumpless()
        elif s.startswith("SET_H,"):
            self.h = arg()
        elif s.startswith("SET_PWMON,"):
            self.pwm_on = int(arg())
        elif s.split(",")[0] in ("SET_KP", "SET_KI", "SET_KD", "SET_LAMBDA", "SET_MU", "SET_KFF"):
            key = {"SET_KP": "Kp", "SET_KI": "Ki", "SET_KD": "Kd", "SET_LAMBDA": "lam", "SET_MU": "mu", "SET_KFF": "Kff"}[s.split(",")[0]]
            g[key] = arg()
            if key == "Kp":
                self.sup.near["Kp"] = g["Kp"]
        elif s == "SUP_ON":
            self.sup_on = True; self.sup.inside = True
        elif s == "SUP_OFF":
            self.sup_on = False
        elif s.startswith("SET_SUP,"):
            _, name, val = s.split(","); val = float(val)
            m = {"ELOW": "e_low", "EHIGH": "e_high", "HS": "hyst", "A": "A"}
            if name in m:
                setattr(self.sup, m[name], val)
            else:
                tgt = self.sup.near if name.endswith("N") else self.sup.far
                tgt[{"KP": "Kp", "LAM": "lam", "MU": "mu", "KFF": "Kff"}[name[:-1]]] = val
        elif s == "REF,DEW":
            self.ref_mode = "DEW"
        elif s == "REF,FIXED":
            self.ref_mode = "FIXED"
        elif s.startswith("SET_TREF,"):
            self.tset = arg(); self.ref_mode = "FIXED"; self.ctrl.Iacc = 0.0
        elif s.startswith("SET_BIAS,"):
            self.bT = arg()
        elif s.startswith("SET_BDEW,"):
            self.bdew = arg()
        elif s.startswith("SET_BAND,"):
            _, c, h = s.split(","); self.band = (float(c), float(h))
        elif s == "RESET_Q":
            self.qdew = 0.0; self.nbelow = 0
        elif s in ("CALIBRATE", "HEADER"):
            pass
        else:
            return "# unknown command: " + s
        return "# ok " + s

    def _bumpless(self):
        e = self.tref - self.fs.y
        self.ctrl.eh[:] = e; self.ctrl.Iacc = 0.0; self.hyst = False

    # ------------------------------------------------------------ one period
    def step(self):
        k, TS = self.k, v.TS
        pl = self.plant
        ts_raw = pl.Tc + self.rng.normal(0, self.noise); th_raw = pl.Th + self.rng.normal(0, self.noise)
        tsup, th = self.fs(ts_raw), self.fh(th_raw)
        if k % 20 == 0:                                   # ambient sensor every 2 s
            self.tdew_raw = float(v.magnus_dew(self.Tamb + self.rng.normal(0, 0.03), self.RH + self.rng.normal(0, 0.3)))
            self.tdew = self.tdew_raw if self.tdew is None else self.tdew + 0.1 * (self.tdew_raw - self.tdew)
        target = self.tdew - self.bdew if self.ref_mode == "DEW" else self.tset - self.bT
        stp = 0.05 * TS
        self.tref = target if abs(target - self.tref) <= stp else self.tref + math.copysign(stp, target - self.tref)
        e = self.tref - tsup
        c, g = self.ctrl, self.gains
        req = 0.0; wfar = 0.0
        if self.control:
            if self.mode == "OPEN":
                req = self.pwm_open; c.pwmT = self.pwm_open
            elif self.mode == "HYST":
                if e < -self.h: self.hyst = True
                elif e > self.h: self.hyst = False
                req = self.pwm_on if self.hyst else 0
            else:
                fop = self.mode == "FOPID"
                c.Ki = g["Ki"] if self.mode != "P" else 0.0
                c.Kd = g["Kd"] if self.mode in ("PID", "FOPID") else 0.0
                if self.sup_on and fop:
                    self.sup.schedule(e); wfar = self.sup.f
                    q = self.sup.schedule(e)
                    c.Kp = q["Kp"] * (self.sup.A if (e > 0 and not self.sup.inside) else 1.0)
                    c.Kff = q.get("Kff", 0.0); lam, mu = q["lam"], q["mu"]
                else:
                    c.Kp, c.Kff = g["Kp"], g["Kff"]
                    lam, mu = (g["lam"], g["mu"]) if fop else (1.0, 1.0)
                if abs(lam - c.lam) > 0.005 or abs(mu - c.mu) > 0.005:
                    c.set_orders(lam, mu)
                c.step(e, self.tref); req = c.u
        # output stage
        if not self.control:
            pwm = 0; c.pwm = 0; c.pwmT = 0
        else:
            tgt = int(max(0, min(255, int(req))))
            pwm = c.pwm
            pwm = min(pwm + 8, tgt) if pwm < tgt else max(pwm - 8, tgt) if pwm > tgt else pwm
            c.pwm, c.pwmT = pwm, tgt
        fan = 255 if (self.fan_mode == "ON" or (self.fan_mode == "AUTO" and pwm > 20)) else 0
        # observer and alarm
        x, inn = self.obs.step(np.array([tsup, th]), pwm, self.Tamb)
        self.innbuf.append(inn); self.innbuf = self.innbuf[-100:]
        mi = np.mean(self.innbuf, axis=0)
        out = abs(mi[0]) > self.band[0] or abs(mi[1]) > self.band[1]
        self.out_count = self.out_count + 1 if out else 0
        self.alarm = self.out_count >= 50
        # condensation counters
        self.qdew += max(0.0, self.tdew - tsup) * TS / 3600; self.nbelow += int(tsup < self.tdew)
        # plant advances one period (fan off -> worse heat sink)
        pl.step(pwm, self.Tamb, self.Qd, self.rconv_mult * (1.0 if fan else 3.5))
        self.k += 1
        mode = self.mode if self.control else "OFF"
        row = [k * 100, self.Tamb, self.RH, self.tdew, self.tref, tsup, th, e, pwm, int(fan > 0),
               c.Kp, g["Ki"], g["Kd"], c.lam, c.mu, c.Kff, mode,
               ts_raw, th_raw, int(self.hyst), 0, 0, req, int(self.sup_on), wfar, 0,
               x[1], inn[0], inn[1], int(self.alarm), self.tdew_raw, self.bdew, self.ref_mode,
               0, self.qdew, self.nbelow]
        return ",".join(f"{r:.4f}" if isinstance(r, float) else str(r) for r in row)

    @staticmethod
    def header():
        return ",".join(COLUMNS)


if __name__ == "__main__":
    # demo: 5 minutes of FOPID at 18.2 C, written to demo.csv
    dev = VirtualDevice()
    for c in ("SET_TREF,18.2", "CONTROL_ON"):
        print(dev.command(c))
    with open("demo.csv", "w") as f:
        f.write(dev.header() + "\n")
        for _ in range(3000):
            f.write(dev.step() + "\n")
    print("wrote demo.csv")
