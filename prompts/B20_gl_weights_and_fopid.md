# B20 — Grünwald–Letnikov weights and the FOPID

My board runs a PID with ramp-aware anti-windup every 100 ms (B19). I want to turn it into a fractional-order PID with Grünwald–Letnikov weights and memory M = 30.

Can you write:
1. a function that fills an array with g0 = 1, gj = g(j−1)·(j − 1 − a)/j for j = 1..M;
2. a circular buffer with the last 31 errors;
3. each period: v = Ts^−(1−λ)·Σ gj^(1−λ)·e[k−j], S = Σ gj^(μ)·e[k−j], D = Ts^−μ·S/(1+|S|); integrate I += Ts·v under my existing anti-windup rule (using v instead of e); request −(Kp·e + Ki·I + Kd·D).

Precompute the weights and the powers of Ts only when λ or μ change (`SET_LAMBDA`, `SET_MU`). Floats only, no dynamic memory. Log λ and μ.

Check for me: with λ = μ = 1 it must give exactly the PID of B19.

Explain why the integral is an accumulator of a fractional derivative of order 1−λ, and what would happen to the steady error if I just truncated a fractional integral of order λ.
