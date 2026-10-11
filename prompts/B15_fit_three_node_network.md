# B15 — Fit the three-node network

I logged an open-loop staircase on my Peltier device (0.1 s sampling; columns `t_ms, Tamb, Tsup, Th, pwm`). My model is a three-node thermal network (cold face Tc, interior Tm, hot face Th):

- Cc·dTc/dt = (Tm − Tc)/R1 − φ·kP·u(t − td)
- Cm·dTm/dt = (Tc − Tm)/R1 + (Th − Tm)/R2
- Ch·dTh/dt = (Tm − Th)/R2 + (1 − φ)·kP·u(t − td) − (Th − Tamb)/Rconv

Can you write a Python script that:
1. simulates it with forward Euler at 0.1 s, using the logged PWM delayed by td and the logged Tamb;
2. fits R1, R2, Rconv, Cc, Cm, Ch, φ and td by least squares on *both* faces, with kP = 0.12 fixed;
3. prints the parameters, the RMSE of each face and whether any parameter ended up on a bound;
4. plots data, fit and residuals.

Bounds: all R and C positive, 0 < φ < 1, 0.1 ≤ td ≤ 3 s. I'll give you a starting guess. It has to be fast (numba or vectorised), because the fit calls it thousands of times.

Explain why I shouldn't fit kP too, and how to recognise a parameter my data can't pin down.
