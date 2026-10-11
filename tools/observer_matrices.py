"""
observer_matrices.py -- matrices of the three-node observer for the firmware (B29).

    python observer_matrices.py ../vplant/params_reference.json          # tau_o = 5 s
    python observer_matrices.py my_params.json --tau 10

Prints Ad, bd, ed, the Luenberger gain L (two probes) and the delay in periods
as C arrays, ready to paste into firmware/FTC_reference/FTC_reference.ino.
"""
import argparse, os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "vplant"))
import vplant as v

ap = argparse.ArgumentParser()
ap.add_argument("params"); ap.add_argument("--tau", type=float, default=5.0)
a = ap.parse_args()
p = v.load_params(a.params)
Ad, bd, ed = v.model_matrices(p)
C = np.array([[1, 0, 0], [0, 0, 1.0]])
O = np.vstack([C, C @ Ad, C @ Ad @ Ad])
print("// rank of the observability matrix:", np.linalg.matrix_rank(O))
L = v.observer_gain(C, a.tau, p)
def f(x):
    t = "%.9g" % x
    return (t if ("." in t or "e" in t) else t + ".0") + "f"
print("float Ad[3][3] = {" + ", ".join("{" + ", ".join(f(x) for x in r) + "}" for r in Ad) + "};")
print("float bd[3] = {" + ", ".join(f(x) for x in bd) + "};")
print("float ed[3] = {" + ", ".join(f(x) for x in ed) + "};")
print("float Lg[3][2] = {" + ", ".join("{" + ", ".join(f(x) for x in r) + "}" for r in L) + "};")
print("#define DELAY_STEPS %d" % max(1, int(round(p["td"] / v.TS))))
