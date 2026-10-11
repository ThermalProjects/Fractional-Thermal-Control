# B19 — The complete PID

I have a PI with ramp-aware anti-windup (B18) in a 100 ms loop; the plant cools when PWM rises.

Can you add the regularised derivative D = (1/Ts)·S/(1+|S|) with S = e[k] − e[k−1], so the request is −(Kp·e + Ki·I + Kd·D)? Add `SET_KP`, `SET_KI`, `SET_KD` and `MODE,PID`, and log the three gains.

No derivative kick: when the controller starts or the reference changes, set e[k−1] = e[k].

Explain the largest kick the derivative can ever give with Kd = 1.24, and why that limit doesn't depend on how noisy the signal is.
