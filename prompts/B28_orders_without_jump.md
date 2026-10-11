# B28 — Change λ and μ online without a jump

My FOPID (B20) computes Grünwald–Letnikov weights from λ and μ, and my supervisor (B27) now changes them every period.

Can you:
1. remember the orders used for the last recomputation and only recompute the two weight arrays and the two powers of Ts if either order moved by more than 0.005;
2. never clear the error buffer or the integral when the weights change;
3. log a flag (`wrec`) in the period when the weights were recomputed, and measure how long one recomputation takes with my `loop_us` column?

Everything together must stay under 50 % of the 100 ms period.

Explain why resetting the integral at a change of order makes the drive jump, and why keeping the error history doesn't.
