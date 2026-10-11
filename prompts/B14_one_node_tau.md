# B14 — Estimate τ with a one-node model

I have a CSV log of an open-loop staircase on my Peltier-cooled plate (0.1 s sampling; `pwm` is the applied drive, `Tsup` the cold-face temperature).

Can you write a Python script that cuts out the first step (PWM 0 → 25), fits T(t) = T0 + K·Δu·(1 − e^(−(t−L)/τ)) for t > L by least squares, prints K, τ, L and the RMSE, and plots data, fit and the 63 % point?

Use `scipy.optimize.least_squares` with bounds (τ > 0, 0 ≤ L ≤ 5 s), and find the step times from the `pwm` column, not from hard-coded indices.

Explain why the 63 % point gives τ, and what the RMSE would be if the model were perfect.
