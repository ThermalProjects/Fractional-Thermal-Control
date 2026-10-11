# B18 — PI with anti-windup

My firmware runs every 100 ms with e = Tref − Tsup and an output stage that clamps the request to 0–255 and ramps the applied PWM by at most 8 counts per period. The plant cools when PWM rises, so the request is −(Kp·e + Ki·I).

Can you implement the PI with conditional anti-windup? Compute the tentative request with the updated integral, then don't update the integral if
(a) the tentative request is above 255 *or* the applied PWM is still below the previous request, and the new increment would push the request higher; or
(b) the tentative request is below 0 *or* the applied PWM is still above the previous request, and the increment would push it lower.
Also clamp |Ki·I| ≤ 255. Reset I to zero when the controller is switched on or the reference is changed.

Keep the controller and the output stage in separate functions; the controller reads the applied and previously requested PWM from the output stage, it doesn't keep its own copies.

Explain why the ramp alone can cause windup even if the request never goes above 255.
