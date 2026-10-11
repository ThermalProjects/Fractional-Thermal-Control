# B1 — Wire the two NTC dividers

Hi! I'm wiring two NTC thermistors to a microcontroller board (a Nano 33 IoT in my case: 3.3 V logic, 12-bit ADC). Each one goes in a voltage divider: NTC from 3.3 V to the analog pin, a 10 kΩ resistor from the pin to GND. Cold face on A0, hot face on A1. The probes are 10 kΩ at 25 °C, B = 3435 K.

Can you give me a short checklist to verify the wiring with a multimeter *before* I power anything? Say what resistance I should read between which points, with the board unplugged, and what changes when I hold the probe in my fingers.

Also tell me what would go wrong if I swapped the NTC and the fixed resistor, and what the pin voltage should be at 25 °C.
