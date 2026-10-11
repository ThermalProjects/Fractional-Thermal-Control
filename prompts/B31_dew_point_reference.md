# B31 — Dew-point reference with bdew over serial

My board reads an ambient temperature/humidity sensor every 2 s without blocking and computes the dew point with Magnus. The controller runs every 100 ms on e = Tref − Tsup.

Can you:
1. filter the dew point with an EMA whose time constant is 20 s at an update every 2 s, starting from the first valid reading;
2. dew-point mode: Tref = Tdew_filtered − bdew; fixed mode: Tref = Tset − bT;
3. add `REF,DEW`, `REF,FIXED`, `SET_BDEW,<b>`, `SET_TREF,<T>`, each echoed back;
4. if the sensor doesn't answer for 10 s, hold the last reference and raise `amb_flag`?

When I switch mode, the reference must not jump: start from the current one and move at no more than 0.05 °C/s. Log the raw and filtered dew point, bdew and the mode.

Explain why I filter the dew point but not the surface temperature any harder than in B5.
