# B12 — Open-loop staircase, logged

My device takes `SET_PWM,<n>`, `FAN_ON` and `FAN_AUTO` over serial and logs one CSV row every 100 ms.

Can you extend my Python monitor with an experiment mode? It reads a small text file with lines `time_s,command`, sends each command at its time, and saves the log with the schedule file name in the header. And write me the staircase schedule: fan on; PWM 0 for 20 s; then 25, 50, 75, 100, 75, 50, 25 and 0, each for 200 s; then fan back to auto.

The schedule runs on the computer's clock, but every row keeps the board's `t_ms`. If the serial link drops, stop the experiment and send `SET_PWM,0`.

Explain how I can check in the log that each command arrived on time, using only the `pwm` column.
