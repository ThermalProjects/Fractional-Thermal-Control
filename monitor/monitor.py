"""
monitor.py -- minimal monitor and twin layer (B8, B33).

Live mode (device on a serial port):
    python monitor.py --port /dev/ttyACM0            # Windows: --port COM5
Simulation mode (the virtual device replaces the board):
    python monitor.py --sim --minutes 10 --cmd "SET_TREF,18.2" --cmd CONTROL_ON
With a twin (your fitted parameters):
    python monitor.py --port COM5 --twin ../vplant/params_reference.json

Every row is saved to run_YYYYMMDD_HHMMSS.csv. The plot refreshes once per
second: Tsup and Tref, the error with the running IAE, and (with --twin) the
residuals of the twin with their band. Type commands in the terminal and press
Enter to send them to the device.
"""
import argparse, datetime, os, sys, threading, queue, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "vplant"))
import vplant as v

ap = argparse.ArgumentParser()
ap.add_argument("--port"); ap.add_argument("--baud", type=int, default=115200)
ap.add_argument("--sim", action="store_true"); ap.add_argument("--minutes", type=float, default=10)
ap.add_argument("--cmd", action="append", default=[]); ap.add_argument("--twin")
ap.add_argument("--noplot", action="store_true")
a = ap.parse_args()

rows, lines = [], queue.Queue()
fname = datetime.datetime.now().strftime("run_%Y%m%d_%H%M%S.csv")
out = open(fname, "w")
header = None

# ---------------------------------------------------------------- source of rows
if a.sim:
    from device import VirtualDevice
    dev = VirtualDevice()
    for c in a.cmd:
        print(dev.command(c))
    def produce():
        lines.put(dev.header())
        for _ in range(int(a.minutes * 600)):
            lines.put(dev.step())
        lines.put(None)
    send = lambda s: print(dev.command(s))
else:
    import serial                                   # pip install pyserial
    ser = serial.Serial(a.port, a.baud, timeout=1)
    time.sleep(2)
    for c in a.cmd:
        ser.write((c + "\n").encode())
    def produce():
        while True:
            l = ser.readline().decode(errors="ignore").strip()
            if l:
                lines.put(l)
    send = lambda s: ser.write((s + "\n").encode())
threading.Thread(target=produce, daemon=True).start()
threading.Thread(target=lambda: [send(l.strip()) for l in sys.stdin], daemon=True).start()

# ---------------------------------------------------------------- twin layer
twin = None; band = None; res = []
if a.twin:
    tp = v.load_params(a.twin)

def twin_step(r):
    """Advance the twin one period with the applied PWM and measured Tamb."""
    global twin, band
    if twin is None:
        twin = v.Plant(r["Tamb"], tp); twin.Tc, twin.Th = r["Tsup"], r["Th"]
        twin.Tm = r["Tm_hat"] if np.isfinite(r["Tm_hat"]) else 0.5 * (r["Tsup"] + r["Th"])
    rc, rh = r["Tsup"] - twin.Tc, r["Th"] - twin.Th
    twin.step(r["pwm"], r["Tamb"])
    res.append((rc, rh))
    if band is None and len(res) == 12000:          # 20 minutes of healthy run
        R = np.array(res)
        band = (R.mean(0), 4 * R.std(0) + 0.3)
    flag = False
    if band is not None:
        recent = np.array(res[-600:])
        flag = bool(np.all(np.abs(recent - band[0]) > band[1], axis=0).any())   # 60 s outside
    return rc, rh, flag

# ---------------------------------------------------------------- loop
if not a.noplot:
    import matplotlib.pyplot as plt
    plt.ion(); fig, ax = plt.subplots(3 if a.twin else 2, 1, sharex=True, figsize=(8, 7))
last = time.time(); done = False
while not done:
    try:
        l = lines.get(timeout=0.2)
    except queue.Empty:
        l = ""
    if l is None:
        done = True
    elif l.startswith("#"):
        print(l)
    elif l and header is None and l.startswith("t_ms"):
        header = l.split(","); out.write(l + (",Tc_twin,Th_twin,rho_c,rho_h,twin_flag" if a.twin else "") + "\n")
    elif l and header:
        f = l.split(",")
        if len(f) != len(header):
            continue                                 # ignore malformed lines
        r = {h: (float(x) if x.replace(".", "", 1).replace("-", "", 1).replace("e", "", 1).isdigit() or x == "nan" else x)
             for h, x in zip(header, f)}
        extra = ""
        if a.twin and isinstance(r.get("Tamb"), float) and np.isfinite(r["Tamb"]):
            rc, rh, flag = twin_step(r)
            r["rho_c"], r["rho_h"], r["twin_flag"] = rc, rh, flag
            extra = f",{twin.Tc:.2f},{twin.Th:.2f},{rc:.3f},{rh:.3f},{int(flag)}"
        out.write(l + extra + "\n"); rows.append(r)
    if not a.noplot and rows and (time.time() - last > 1.0 or done):
        last = time.time()
        t = np.array([r["t_ms"] for r in rows]) / 1000
        e = np.array([r["e"] for r in rows], float)
        for x in ax: x.cla()
        ax[0].plot(t, [r["Tsup"] for r in rows], label="Tsup"); ax[0].plot(t, [r["Tref"] for r in rows], label="Tref")
        ax[0].set_ylabel("°C"); ax[0].legend(loc="upper right")
        ax[1].plot(t, e); ax[1].set_ylabel("e (°C)")
        ax[1].set_title(f"IAE = {np.nansum(np.abs(e)) * 0.1:.1f} °C·s")
        if a.twin and res:
            R = np.array(res); ax[2].plot(t[-len(R):], R[:, 0], label="rho_c"); ax[2].plot(t[-len(R):], R[:, 1], label="rho_h")
            if band is not None:
                for i in (0, 1):
                    ax[2].axhline(band[0][i] + band[1][i], ls=":"); ax[2].axhline(band[0][i] - band[1][i], ls=":")
            ax[2].set_ylabel("residual (°C)"); ax[2].legend(loc="upper left")
        ax[-1].set_xlabel("time (s)")
        plt.pause(0.01)
out.close()
print("saved", fname, len(rows), "rows")
if not a.noplot:
    plt.ioff(); plt.show()
