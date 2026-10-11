# B29 — An observer beside the controller

I have a three-node thermal model (cold face, interior, hot face) with R1, R2, Rconv, Cc, Cm, Ch, φ, kP, td fitted to my device, sampled at Ts = 0.1 s. My board measures both faces every period and logs the applied PWM and the ambient temperature.

Can you write:
1. a Python script that builds the matrices, discretises them with forward Euler (Ad = I + Ts·A), checks the rank of the observability matrix with both faces measured, and computes a 3×2 Luenberger gain by pole placement with time constants τo, 1.25τo, 1.5τo (τo = 5 s), printing everything as C arrays (the repo has `tools/observer_matrices.py` if you want to compare);
2. the firmware part: each period, innovation r = y − C·x̂⁻, correct x̂ = x̂⁻ + L·r, then predict x̂⁻ = Ad·x̂ + bd·u[k−d] + ed·Tamb, with u the PWM *applied* d periods ago?

Start x̂ at the first measured faces and their mean for the interior. No matrix library, and the observer must not change the control law. Log T̂m and both innovations.

Explain why the observer uses the applied PWM rather than the requested one, and how the cold-face innovation can correct the interior temperature.
