# B25 — Ramp and ramp-aware anti-windup, checked

My firmware's output stage ramps the PWM; the FOPID's anti-windup reads the applied and requested PWM of the previous period. I want to prove the rule works.

Can you write a Python script for the virtual plant that runs scenario S1 (23.7 → 18.2 °C, 300 s) with ramp steps of 8, 2 and 1 counts per period, each with the ramp-aware rule on and off, and prints a small table of undershoot and settling time, plus one plot with the six responses?

Explain why the undershoot grows so much with a slow ramp when the rule is off.
