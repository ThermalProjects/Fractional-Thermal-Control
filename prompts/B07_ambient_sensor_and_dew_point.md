# B7 — Ambient sensor and dew point

My 100 ms loop is already running. I have an I2C temperature/humidity sensor at address 0x44. It starts one measurement when I send the bytes 0x24 0x00, and 15–60 ms later it gives six bytes: temperature (2 bytes), checksum, humidity (2 bytes), checksum. T = −45 + 175·sT/65535 and RH = 100·sRH/65535.

Can you write:
1. a non-blocking state machine that starts a measurement every 2 s and reads the result in a later loop period (never wait inside a period!);
2. a function that gives the dew point with the Magnus formula, a = 17.67, b = 243.5;
3. the reference: Tref = Tset − bT in fixed mode, Tdew − bdew in dew-point mode.

If the sensor doesn't answer, keep the last good values and raise a flag.

Explain why the formula uses ln(RH/100), and what happens to Tdew if RH comes back as 0.
