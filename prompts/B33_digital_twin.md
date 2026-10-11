# B33 — A digital twin synchronised with the device

My board sends one CSV row every 100 ms with, among others, the time, Tamb, the applied PWM, Tsup, Th and T̂m. I have the book's Python virtual plant: `step(pwm, Tamb)` advances one period and the parameters come from a JSON file.

Can you write a program that:
1. reads the serial port line by line, saves every row to a dated CSV, and for each row advances the virtual plant once with that row's PWM and Tamb;
2. computes the residuals ρc = Tsup − T̃c and ρh = Th − T̃h;
3. after 20 minutes freezes a band of mean ± (4σ + 0.3 °C) for each residual and raises a flag when a residual stays outside for 60 s;
4. shows a live plot: measured vs twin temperatures, residuals with bands, and the flags;
5. has a simulation mode where the virtual plant also replaces the device?

Start the twin from the first measured faces and the first T̂m; if a row is missing, still advance one step using the last PWM. (The repo's `monitor/monitor.py` does this — compare with yours.)

Explain why the twin gets the applied PWM and not the reference, and why it isn't corrected with the measurements like the observer.
