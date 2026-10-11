# B23 — Load the gains and compare

I ran scenario S1 on my device three times with the PID and three times with the FOPID (CSV logs with my column dictionary), and I have a prediction of each from my virtual plant.

Can you write a script that computes for every run: IAE, ITAE, undershoot, settling time to ±0.5 °C, mean |e| after 200 s and mean duty cycle; prints a table with mean ± standard deviation over the three runs of each controller next to the prediction; and plots the six runs and the two predictions?

Align the runs on the moment the controller was switched on, and use the logged Tref of each row.

Explain how to decide whether the difference between the two controllers is bigger than the run-to-run spread.
