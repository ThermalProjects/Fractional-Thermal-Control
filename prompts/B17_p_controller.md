# B17 — P controller

My firmware has an output stage (clamp 0–255, ramp 8 counts per period) and computes e = Tref − Tsup every 100 ms. The plant cools when PWM rises.

Can you add a proportional controller as a selectable mode (`MODE,P`) that requests PWM* = −Kp·e, with Kp settable by `SET_KP,<v>`? The request goes through the existing output stage.

Then tell me what steady-state error I should expect for a given steady PWM, and how it changes if I double Kp. I'll check it on the virtual plant first and then on the device, 5 °C below ambient.
