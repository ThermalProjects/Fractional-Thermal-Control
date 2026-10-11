# B16 — Validate on another run and export your virtual plant

I fitted my three-node network to run 1 of my staircase (parameters below). The book's virtual plant reads its parameters from a JSON file (`vplant/params_reference.json` is the template) and runs the same serial commands as my firmware.

Can you write a script that:
1. simulates run 2 with my parameters, *without refitting*, and prints the RMSE of each face;
2. writes my parameters into a new JSON file for the virtual plant;
3. runs the staircase schedule on the virtual plant and plots simulated and measured Tsup and Th on the same axes?

Use each run's own logged Tamb, and start each node at the first measured value (interior = average of the two faces).

Explain what a validation RMSE much bigger than the fitting RMSE would tell me.
