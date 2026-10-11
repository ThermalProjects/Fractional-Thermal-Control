# B26 — Arduino–Python consistency check (Δu)

I have a CSV log of my board running the FOPID (columns include `e`, `req`, `pwm`, the five parameters and M), and a Python function with the same law, including the anti-windup rule that reads the applied and requested PWM of the previous period.

Can you write a script that feeds the logged error of each row into the Python law, in order, using the logged applied PWM for the anti-windup rule, and computes Δu = u_Python − u_board for every row? Report max |Δu|, plot Δu vs time, and list the first ten rows where |Δu| > 1.

Don't refit or tweak anything, and start the Python controller from the same initial state as the board.

Explain what kind of bug gives a constant Δu, one that grows with time, and one that only shows up during transients.
