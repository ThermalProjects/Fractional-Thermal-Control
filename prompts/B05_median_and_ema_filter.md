# B5 — Median of three + EMA

I read two temperatures every 100 ms. The readings tremble a bit and once in a while there's a single crazy spike (a loose contact, the PWM switching…).

Can you add two filters in series to each channel? First a median of the last three samples (kills isolated spikes), then an exponential moving average y += 0.2·(m − y). Initialise both with the first valid reading so the filter doesn't start from zero. Keep the raw value too, because I log both.

Explain why the median goes first and the EMA second, and roughly how much delay the EMA adds with α = 0.2 at 100 ms (τ ≈ Ts/α).
