# B32 — Condensation indicator in the log

My board logs, every 100 ms, the filtered dew point Tdew and the surface temperature Tsup, and runs in dew-point mode with margin bdew.

Can you:
1. in the firmware, accumulate Q += max(0, Tdew − Tsup)·Ts/3600 and count the periods with Tsup < Tdew; log both and reset them with `RESET_Q` (don't accumulate while `amb_flag` is up);
2. in my monitor, show Q in °C·h and the % of time below the dew point;
3. a Python script that recomputes both from a finished log and also estimates the mass with ṁ = hm·A·Δp/(Rv·Tf), with hm and A as parameters (start with hm = 0.006 m/s, A = 0.04 × 0.04 m²), and plots the cumulative mass?

Explain why Q is only roughly proportional to the condensed mass, and which assumptions of that formula I should trust least on my bench.
