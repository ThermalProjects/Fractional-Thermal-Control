# B24 — Measure the time of each loop

My board runs a 100 ms loop scheduled with `tLast += 100`. Every period it reads two NTCs (10 conversions each), runs a FOPID, writes PWM and prints a CSV row; every 2 s it reads an I2C temperature/humidity sensor that needs 15–60 ms to convert.

Can you:
1. add two log columns: `loop_us` (how long the work of the period took) and `late_ms` (how late the period started);
2. turn the sensor read into a two-state machine: start the conversion in one period, read the six bytes in the first period at least 60 ms later?

The measurement mustn't change the timing it measures by more than a few tens of µs, and no `delay()` anywhere in `loop()`.

Explain which part of the period takes longest, and how I can see in the log whether the I2C read was blocking.
