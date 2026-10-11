# B27 — Region scheduler with hysteresis

My board runs a FOPID every 100 ms with e = Tref − Tsup from the filtered temperature. Kp, Ki, Kd, λ, μ and Kff are global variables I can change over serial.

Can you write a supervisor function, called once per period before the controller?
1. A boolean "inside": enter when |e| < e_low, leave when |e| > e_low + hs.
2. w = 0 when inside; otherwise w = clamp((|e| − e_low)/(e_high − e_low), 0, 1).
3. For Kp, λ, μ, Kff: value = near + w·(far − near).
4. Multiply Kp by A when e > 0 and not inside.
5. Commands `SUP_ON` / `SUP_OFF` and `SET_SUP,<name>,<value>` for the eleven values (ELOW, EHIGH, HS, A, KPN, LAMN, MUN, KFFN, KPF, LAMF, MUF, KFFF). When off, the FOPID keeps the single tuning.

Defaults: e_low = 0.3, e_high = 3.0, hs = 0.1, A = 1.2; near = (Kp 50, λ 1.00, μ 1.24, Kff 0); far = (Kp 80, λ 0.78, μ 1.24, Kff 0). No dynamic memory, no blocking, and don't touch the integral or the error history. Log w, Kp, λ and μ.

Explain what would happen to the drive if I just switched between the two sets at a single threshold.
