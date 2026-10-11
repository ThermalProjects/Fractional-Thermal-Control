# B8 — CSV log and live plot

My board prints one CSV row every 100 ms at 115200 baud with these columns:
`t_ms, Tamb, RH, Tdew, Tref, Tsup, Th, e, pwm, fan, Kp, Ki, Kd, lambda, mu, Kff, mode`

Can you write a Python program that opens the serial port, saves every complete row to `run_YYYYMMDD_HHMMSS.csv` (with the header), and shows a live plot with two panels: Tsup and Tref vs time, and the error e vs time with the running IAE = Σ|e|·Ts in the title?

Ignore incomplete or garbled lines without crashing, don't lose rows while it's drawing, and refresh the plot once per second.

Explain how the program keeps reading the port while the plot is redrawn, and how I can check afterwards that no row was lost.
