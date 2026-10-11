# B2 — Read the raw ADC

My board is a Nano 33 IoT (3.3 V, 12-bit ADC). I've got two NTC dividers: NTC between 3.3 V and the pin, 10 kΩ from the pin to GND. Cold face on A0, hot face on A1.

Can you write me a tiny sketch that sets the ADC to 12 bits and prints the raw counts of A0 and A1, comma-separated, five times per second? No libraries and no conversion to temperature yet — I just want to see the numbers move.

Please explain which function sets the ADC resolution on this board (and what to do on an Uno, which is 10-bit), and what count I should expect when the thermistor is exactly at 25 °C.
