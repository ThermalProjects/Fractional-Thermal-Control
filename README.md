# Fractional Thermal Control in Embedded Systems — book repository

Companion repository of the book *Fractional Thermal Control in Embedded Systems:
Build, Measure and Regulate a Low-Cost Thermoelectric Experiment* (Springer Nature).

Everything you need to follow the build steps B1–B34 is here: the reference firmware,
the virtual plant, a monitor with a digital twin, the prompts of every build step and
the scripts that produce every simulated figure of the book.

```
requirements.txt          Python packages (pip install -r requirements.txt)
firmware/FTC_reference/   reference firmware (one sketch, all controllers)
vplant/                   virtual plant, virtual device, reference parameters
monitor/                  live monitor + twin layer (serial or simulation)
tools/                    observer matrices for your fitted model (B29)
prompts/                  one prompt per build step (B01 … B34)
figs/                     one script per chapter for the simulated figures
```

## Installation

You need **Python 3.10 or newer** (from python.org; on Windows tick *Add python.exe to PATH*
during installation) and the packages listed in `requirements.txt`:

| Package | What for |
|---|---|
| numpy, scipy, matplotlib | virtual plant, analysis and plots (required) |
| pandas | reading the CSV logs (required) |
| pyserial | talking to the real board (only with hardware) |
| numba | speeds up the tuning (PSO) and the Chapter 4 fits (optional, but much faster) |

### Windows (PowerShell)

Open PowerShell in the unzipped `FTC_repo` folder and run:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

If `Activate.ps1` is refused with a permissions error, run this once and activate again:

```powershell
Set-ExecutionPolicy -Scope CurrentUser RemoteSigned
```

Each time you open a new PowerShell window, activate the environment again with
`.\.venv\Scripts\Activate.ps1` before running the scripts.

### macOS / Linux

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Quick start without hardware

```
cd monitor
python monitor.py --sim --minutes 5 --cmd "SET_TREF,18.2" --cmd CONTROL_ON
```

A live plot opens and a `run_YYYYMMDD_HHMMSS.csv` is saved. While it runs you can type
commands in the terminal and press Enter (`MODE,PID`, `SUP_ON`, `REF,DEW`,
`SET_BDEW,1.5`, …).

## With your device

1. Open `firmware/FTC_reference/FTC_reference.ino`, check the BOARD and PINS sections
   for your board, and upload it.
2. Wait for `# calibrated …` (the device must be at room temperature; control is refused
   until the calibration is done).
3. With the environment activated: `python monitor/monitor.py --port COM5`, then type
   `SET_TREF,18.2` and `CONTROL_ON`. Your port name is shown in the Arduino IDE
   (Tools > Port) or in the Windows Device Manager; on macOS/Linux it looks like
   `/dev/ttyACM0` or `/dev/cu.usbmodem…`.
4. With your own fitted model (B16): copy `vplant/params_reference.json`, put your values
   in it, and add `--twin my_params.json` to run the digital twin beside the device.

## Boards

The numbers in the book are computed for an Arduino Nano 33 IoT (3.3 V, 12-bit ADC),
but the firmware runs on any board with two analog inputs, I²C and two PWM outputs.

- 3.3 V boards, 12-bit ADC (Nano 33 IoT, ESP32, Raspberry Pi Pico): `ADC_BITS 12`.
- 5 V AVR boards (Uno, Nano, Leonardo): `ADC_BITS 10`, wire the dividers to 5 V, and remove
  the `analogReadResolution` line. An Uno has only 2 KB of RAM and cannot hold the whole
  firmware (calibration history, GL buffers, observer); it is fine for Part I, and for the
  rest use a Nano 33 IoT, ESP32 or Pico.
- A Raspberry Pi has no ADC: use an external I²C or SPI converter.

## Reference parts

| Part | Reference |
|---|---|
| Board | Arduino Nano 33 IoT |
| NTC probes (2) | Littelfuse USP16673, 10 kΩ ±1 %, B25/85 = 3435 K |
| Ambient sensor | DFRobot SEN0546 (I²C, 0x44), ±0.5 °C, ±3 % RH |
| Peltier module | European Thermodynamics ETH-241-10-13-E-H1 (40 × 40 mm) |
| Switches (2) | XY-MOS low-side MOSFET modules (or discrete logic-level N-MOSFETs) |

## Serial commands and log columns

See Appendix B of the book. In short: `CALIBRATE`, `CONTROL_ON/OFF`, `SET_PWM,n`,
`FAN_AUTO/ON/OFF`, `MODE,HYST|P|PI|PID|FOPID`, `SET_H`, `SET_PWMON`, `SET_KP/KI/KD`,
`SET_LAMBDA`, `SET_MU`, `SET_KFF`, `SUP_ON/OFF`, `SET_SUP,<name>,<v>`, `REF,DEW|FIXED`,
`SET_TREF`, `SET_BIAS`, `SET_BDEW`, `SET_BAND,<cold>,<hot>`, `RESET_Q`, `HEADER`.

The log has 36 columns; columns are only ever appended, never renamed or reordered,
so every script written for an early log also reads a later one.

## Regenerating the figures

```
cd figs
python ch08.py          # figures go to figs/out/ (or set FIGOUT=...)
```

Each script also writes `chXX_numbers.json` with the numbers quoted in the text.
Change the parameters in `vplant/params_reference.json` (or load yours) to see how
your device differs from the reference one.

## Safety

The power path carries up to about 2 A at 12 V. Never run the Peltier cell without its
heat sink, keep the 12 V wiring away from the board's pins, and keep the electronics
above and away from the cold plate: it makes water. See Appendix A of the book.

## License and citation

The code (firmware, virtual plant, monitor, tools and figure scripts) is released under
the MIT License (`LICENSE`). The prompts in `prompts/` are released under Creative Commons
Attribution 4.0 (`prompts/LICENSE`). If you use this repository in your teaching or
research, please cite it: GitHub shows a *Cite this repository* button built from
`CITATION.cff`.

