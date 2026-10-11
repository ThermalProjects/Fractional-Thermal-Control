# B4 — A timed loop at Ts = 0.1 s

My sketch reads two NTC channels and converts them to °C with a function I already have. Now I need a fixed sampling period of exactly 100 ms, because a controller will run on it later.

Can you restructure my `loop()` so measuring and printing happen once every 100 ms using `millis()`, with no `delay()` at all? Each sample should be the average of 10 ADC conversions per channel. Print one line per period: `t_ms,Tsup,Th`. Serial at 115200 baud.

Important: the period must not drift if one iteration runs late — use `tLast += 100`, not `tLast = millis()`.

Explain the difference between those two lines and what happens to the sampling grid in each case when one cycle is late.
