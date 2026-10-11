# B3 — From counts to temperature

I've got two NTC thermistors on my board (Nano 33 IoT: 3.3 V, 12-bit ADC). Each one sits on the top side of a voltage divider with a 10 kΩ fixed resistor to GND, read on A0 and A1. The NTC is 10 kΩ at 25 °C with B = 3435 K.

Can you write me a function `float countsToCelsius(int N)`? First get the NTC resistance from the divider, R = Rf·(4095/N − 1), then use the beta equation 1/T = 1/T0 + ln(R/R0)/B with T0 = 298.15 K, and return °C. Keep it simple: plain floats, no libraries.

Watch out for N = 0 and N = 4095 so it doesn't divide by zero — return something obviously wrong, like −99, so I notice.

And explain why the 3.3 V doesn't show up anywhere in the formula, and which line I'd change if the NTC were on the bottom side of the divider.
