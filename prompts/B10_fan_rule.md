# B10 — The fan rule

My output stage (B9) sets the Peltier PWM on D6. The fan is on D5 through a second low-side MOSFET with a flyback diode across the fan.

Can you add the fan rule after the ramp: fan = 255 if the applied PWM > 20, otherwise 0? Plus the commands `FAN_AUTO`, `FAN_ON` and `FAN_OFF`, and a `fan` column (0/1) in the log.

Explain why the fan is either fully on or off instead of proportional to the drive, and what the flyback diode protects.
