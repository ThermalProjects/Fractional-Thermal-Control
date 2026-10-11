# B13 — Hysteresis controller

My firmware computes e = Tref − Tsup every 100 ms and has an output stage with a ramp (B9). Careful: the plant cools when PWM goes up.

Can you add a hysteresis controller? Cooling turns on when e < −h and off when e > +h; it requests PWM_ON while on and 0 while off. Make h and PWM_ON settable over serial (`SET_H`, `SET_PWMON`) and select it with `MODE,HYST`. The controller only produces a request — the ramp and the fan rule stay in the output stage. Log the on/off state in its own column.

Explain what sets the period and the amplitude of the oscillation I'll see, and why it doesn't stay inside the ±h band.
