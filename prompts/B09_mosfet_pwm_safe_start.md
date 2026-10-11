# B9 — MOSFET + PWM with a safe start

My 100 ms loop is running. A Peltier cell (12 V, up to about 2 A) is switched on the low side by a MOSFET (a ready-made module or a logic-level N-MOSFET with 220 Ω on the gate and a 10 kΩ pull-down), driven from pin D6.

Can you write an output stage `applyPWM(int requested)` that clamps the request to 0–255, moves the applied value towards it by at most 8 counts per call, and writes it with `analogWrite`? Plus the serial commands `SET_PWM,<n>` (open-loop tests), `CONTROL_ON` and `CONTROL_OFF`.

Safety first: in `setup()`, D6 and D5 must be set as outputs and written LOW before anything else. With control off, the applied PWM goes to 0 immediately (no ramp down). Log the applied value every period in the `pwm` column.

Explain why the pull-down matters during the first milliseconds after power-up, and how long the ramp takes from 0 to 255.
