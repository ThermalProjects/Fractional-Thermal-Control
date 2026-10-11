# B30 — Fault warning from the innovation

My observer (B29) logs two innovations every 100 ms. I have a healthy log of 15+ minutes with a few reference changes.

Can you:
1. write a Python script that computes the 10 s moving average of each innovation in the healthy log and prints five-sigma bands (|mean| + 5·std) for each;
2. add to the firmware a flag `alarm` that rises when either 10 s average stays outside its band for 5 s, with the bands settable by `SET_BAND,<cold>,<hot>`?

Then I'll test it on the virtual plant with a simulated fan failure and on the device by unplugging the fan for two minutes.

Explain why I watch the innovation and not just a threshold on Th, and why a very fast observer would hide the fault.
