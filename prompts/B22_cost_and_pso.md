# B22 — Cost function + PSO on the virtual plant

I have a Python virtual plant that simulates my device with the exact firmware law of my FOPID (λ = μ = 1 gives the PID) and returns the cold-face temperature every 0.1 s.

Can you write:
1. `cost(x)`: run the scenario (Tamb = 23.7 °C, Tref = 18.2 °C, 300 s, starting at ambient) and return J = ITAE + 1000·Jss² + 500·U + ts, where ITAE = Σ t·|e|·Ts, Jss = |mean of the last 2 s − Tref|, U = how far it went below Tref, ts = time after which |e| stays below 0.5 °C;
2. a PSO: 100 particles, up to 300 iterations, inertia from 0.9 to 0.4, c1 = c2 = 1.4962, velocity limit 10 % of the range, bounce off the bounds, stop after 50 iterations without improvement;
3. run it for x = (Kp, Ki, Kd, λ, μ) and for x = (Kp, Ki, Kd) with my bounds, and save the best point, its cost and the convergence history.

Use a seeded random generator, and make the cost return a big finite number (not an exception) if a run blows up.

Explain what it means if the best point sits on a bound, and how to tell a real optimum from a plateau.
