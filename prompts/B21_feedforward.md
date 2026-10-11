# B21 — Feedforward over serial

My FOPID (B20) requests −C(e). I want to add a feedforward term Kff·Tref, settable with `SET_KFF,<v>` and logged, default 0.

Can you add it to the request *and* to the tentative request used by the anti-windup rule? Keep everything else identical.

Explain why the feedforward must also enter the anti-windup test, and what a too-large Kff would look like in the log.
