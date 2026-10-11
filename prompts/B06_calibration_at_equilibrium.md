# B6 — Calibration of offsets, only at equilibrium

My device reads two NTCs (Tsup on A0, Th on A1) and an ambient temperature Tamb from an I2C sensor, every 100 ms. I want a one-point offset calibration of each thermistor against Tamb.

Can you write a calibration routine that:
1. checks equilibrium: over the last 30 s both thermistors changed by less than 0.05 °C and both are within 1 °C of Tamb;
2. only then averages 150 samples per channel and stores the offsets δ = Tamb − mean(T);
3. otherwise keeps waiting and prints why.

It runs at start-up and when I send `CALIBRATE`. The Peltier and the fan must be off while it runs, and it must not block the 100 ms loop.

Explain why calibrating with a warm heat sink gives a wrong offset, and how the equilibrium test catches it.
