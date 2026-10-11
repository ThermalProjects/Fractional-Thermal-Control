/*
  FTC_reference.ino -- reference firmware of the book
  "Fractional Thermal Control in Embedded Systems".

  One sketch for the whole book. Every controller (hysteresis, P, PI, PID,
  FOPID), the supervisor, the observer and the dew-point reference are here,
  selected and tuned over the serial port. Commands and CSV columns are listed
  in Appendix B of the book and in README.md.

  Reference board: Arduino Nano 33 IoT (3.3 V, 12-bit ADC). For other boards
  change only the BOARD section below.

  Sign convention: e = Tref - Tsup. The plant cools when PWM rises, so the
  controller requests PWM* = -C(e) + Kff * Tref.
*/

#include <Wire.h>
#include <math.h>

// ============================================================ BOARD
// 3.3 V boards with 12-bit ADC (Nano 33 IoT, ESP32, Pico): ADC_BITS 12
// 5 V AVR boards (Uno, Nano, Leonardo): ADC_BITS 10 and wire the dividers to 5 V
#define ADC_BITS 12
const float ADC_MAX = (float)((1UL << ADC_BITS) - 1);

// ============================================================ PINS
#define PIN_TSUP   A0      // cold-face NTC divider
#define PIN_TH     A1      // hot-face NTC divider
#define PIN_CELL   6       // Peltier switch (PWM)
#define PIN_FAN    5       // fan switch
#define AMB_ADDR   0x44    // I2C temperature-humidity sensor

// ============================================================ CONSTANTS
const unsigned long TS_MS = 100;          // sampling period
const float TS = 0.1f;
const float RF = 10000.0f, R0 = 10000.0f, BETA = 3435.0f, T0K = 298.15f;
const int   N_AVG = 10;                   // conversions averaged per sample
const float ALPHA = 0.2f;                 // EMA after the median of three
const int   RAMP = 8;                     // PWM counts per period
const int   FAN_THRESHOLD = 20;           // fan on above this drive
const float MAG_A = 17.67f, MAG_B = 243.5f;
const float TAU_REF = 20.0f;              // dew-point filter time constant (s)
const float REF_SLEW = 0.05f;             // C/s when the reference mode changes
const float T_MIN_OK = -30.0f, T_MAX_OK = 90.0f;   // probe rating

#define M 30                              // GL memory length

// ============================================================ STATE: measurement
float tsupRaw, thRaw, tsup, th;           // raw and filtered temperatures
float medS[3], medH[3];
float offS = 0, offH = 0;                 // calibration offsets
float tamb = NAN, rh = NAN, tdewRaw = NAN, tdew = NAN;
bool  ambFlag = false;
unsigned long ambLastOk = 0;

// ============================================================ STATE: reference
enum RefMode { REF_FIXED, REF_DEW };
RefMode refMode = REF_FIXED;
float tset = 18.0f;                       // fixed set point
float bT = 0.0f;                          // THERMAL_BIAS (fixed mode)
float bdew = 0.0f;                        // margin below the dew point
float tref = 18.0f;                       // reference in force (slew-limited)

// ============================================================ STATE: controller
enum Mode { CTRL_OPEN, CTRL_HYST, CTRL_P, CTRL_PI, CTRL_PID, CTRL_FOPID };
Mode mode = CTRL_FOPID;
bool control = false;
int  pwmOpen = 0;
float h = 0.3f; int pwmOn = 255; bool hystState = false;
float Kp = 50.0f, Ki = 3.00f, Kd = 1.75f, lam = 0.78f, mu = 1.24f, Kff = 0.0f;
float eh[M + 1];                          // error history, eh[0] = newest
float wV[M + 1], wD[M + 1], tsV, tsD;     // GL weights and powers of Ts
float lamW = -1, muW = -1;                // orders used for the weights
bool  wrec = false;
float I = 0;                              // integral accumulator
float req = 0;                            // requested drive (full precision)
int   pwmTarget = 0, pwm = 0, fan = 0;
enum FanMode { FAN_AUTO_M, FAN_ON_M, FAN_OFF_M };
FanMode fanMode = FAN_AUTO_M;

// ============================================================ STATE: supervisor
bool  supOn = false, inside = true;
float eLow = 0.3f, eHigh = 3.0f, hs = 0.1f, aggr = 1.2f;
float KpN = 50, lamN = 1.00f, muN = 1.24f, KffN = 0;      // near set
float KpF = 80, lamF = 0.78f, muF = 1.24f, KffF = 0;      // far set
float wfar = 0;

// ============================================================ STATE: observer
// Three-node network discretised with Ts (forward Euler) and the Luenberger
// gain for tau_o = 5 s. Defaults: reference virtual plant. Replace with the
// output of repo/tools/observer_matrices.py for your fitted model (B29).
float Ad[3][3] = {{0.997587917f, 0.00241208262f, 0.0f}, {3.01609692e-05f, 0.999950363f, 1.94756853e-05f}, {0.0f, 0.000236249204f, 0.998904457f}};
float bd[3] = {-0.000352614083f, 0.0f, 0.000123132366f};
float ed[3] = {0.0f, 0.0f, 0.000859294033f};
float Lg[3][2] = {{0.0304203876f, 0.00167441291f}, {0.107045237f, 0.0105945379f}, {0.00166172428f, 0.0149411947f}};
#define DELAY_STEPS 9
float xp[3];                              // predicted state
float innC = 0, innH = 0, tmHat = NAN;
int   ubuf[DELAY_STEPS + 1]; int ubufIdx = 0;
bool  obsInit = false;
// innovation alarm: 10 s mean outside +-band for 5 s
float bandC = 0.04f, bandH = 0.17f;
float innCbuf[100], innHbuf[100]; int innIdx = 0; float innCsum = 0, innHsum = 0;
int   outCount = 0; bool alarm = false;

// ============================================================ STATE: condensation
double qdew = 0; unsigned long nBelow = 0;

// ============================================================ STATE: timing
unsigned long tStart, tLast, loopUs = 0; long lateMs = 0;

// ============================================================ calibration
bool calibrating = false;
float histS[300], histH[300]; int histN = 0, histIdx = 0;   // last 30 s
double sumS = 0, sumH = 0; int calN = 0;

// =================================================================
// MEASUREMENT
// =================================================================
float countsToC(float n) {
  if (n < 1.0f || n > ADC_MAX - 1.0f) return -99.0f;   // open or shorted probe
  float r = RF * (ADC_MAX / n - 1.0f);                  // NTC on the high side
  return 1.0f / (1.0f / T0K + logf(r / R0) / BETA) - 273.15f;
}

float readAvg(int pin) {
  long s = 0;
  for (int i = 0; i < N_AVG; i++) { s += analogRead(pin); delayMicroseconds(50); }
  return (float)s / N_AVG;
}

float median3(float *b, float x) {
  b[0] = b[1]; b[1] = b[2]; b[2] = x;
  float a = b[0], c = b[1], d = b[2];
  return fmaxf(fminf(a, c), fminf(fmaxf(a, c), d));
}

void measure() {
  tsupRaw = countsToC(readAvg(PIN_TSUP)) + offS;
  thRaw   = countsToC(readAvg(PIN_TH)) + offH;
  tsup += ALPHA * (median3(medS, tsupRaw) - tsup);
  th   += ALPHA * (median3(medH, thRaw) - th);
}

// Ambient sensor: start a conversion, collect it at least 60 ms later.
enum AmbState { AMB_IDLE, AMB_WAIT };
AmbState ambState = AMB_IDLE; unsigned long ambT = 0;

float magnusDew(float T, float RH) {
  float g = logf(RH / 100.0f) + MAG_A * T / (MAG_B + T);
  return MAG_B * g / (MAG_A - g);
}

void serviceAmbient(unsigned long now) {
  if (ambState == AMB_IDLE && now - ambT >= 2000) {
    Wire.beginTransmission(AMB_ADDR); Wire.write(0x24); Wire.write(0x00); Wire.endTransmission();
    ambState = AMB_WAIT; ambT = now;
  } else if (ambState == AMB_WAIT && now - ambT >= 60) {
    ambState = AMB_IDLE;
    if (Wire.requestFrom(AMB_ADDR, 6) == 6) {
      uint16_t st = ((uint16_t)Wire.read() << 8) | Wire.read(); Wire.read();
      uint16_t sh = ((uint16_t)Wire.read() << 8) | Wire.read(); Wire.read();
      tamb = -45.0f + 175.0f * st / 65535.0f;
      rh = 100.0f * sh / 65535.0f;
      tdewRaw = magnusDew(tamb, rh);
      const float a = 2.0f / TAU_REF;                     // update every 2 s
      tdew = isnan(tdew) ? tdewRaw : tdew + a * (tdewRaw - tdew);
      ambLastOk = now;
    }
  }
  ambFlag = (now - ambLastOk > 10000);
}

// Calibration: offsets against Tamb, only at equilibrium (B6).
void serviceCalibration() {
  histS[histIdx] = tsup - offS; histH[histIdx] = th - offH;
  histIdx = (histIdx + 1) % 300; if (histN < 300) histN++;
  if (!calibrating || isnan(tamb) || histN < 300) return;
  float mnS = 1e9, mxS = -1e9, mnH = 1e9, mxH = -1e9;
  for (int i = 0; i < 300; i++) {
    mnS = fminf(mnS, histS[i]); mxS = fmaxf(mxS, histS[i]);
    mnH = fminf(mnH, histH[i]); mxH = fmaxf(mxH, histH[i]);
  }
  bool eq = (mxS - mnS < 0.05f) && (mxH - mnH < 0.05f) &&
            fabsf(histS[(histIdx + 299) % 300] - tamb) < 1.0f &&
            fabsf(histH[(histIdx + 299) % 300] - tamb) < 1.0f;
  if (!eq) return;                                        // keep waiting
  sumS += tsup - offS; sumH += th - offH; calN++;
  if (calN >= 150) {
    offS = tamb - (float)(sumS / calN); offH = tamb - (float)(sumH / calN);
    calibrating = false; sumS = sumH = 0; calN = 0;
    Serial.print(F("# calibrated offS=")); Serial.print(offS, 3);
    Serial.print(F(" offH=")); Serial.println(offH, 3);
  }
}

// =================================================================
// REFERENCE
// =================================================================
void updateReference() {
  float target = tset - bT;
  if (refMode == REF_DEW && !isnan(tdew) && !ambFlag) target = tdew - bdew;
  if (refMode == REF_DEW && (isnan(tdew) || ambFlag)) target = tref;   // hold
  float step = REF_SLEW * TS;
  if (target > tref + step) tref += step;
  else if (target < tref - step) tref -= step;
  else tref = target;
}

// =================================================================
// CONTROLLER
// =================================================================
void glWeights(float *g, float a) {
  g[0] = 1.0f;
  for (int j = 1; j <= M; j++) g[j] = g[j - 1] * ((j - 1) - a) / j;
}

void updateWeights(float l, float m) {
  wrec = false;
  if (fabsf(l - lamW) > 0.005f || fabsf(m - muW) > 0.005f) {
    glWeights(wV, 1.0f - l); glWeights(wD, m);
    tsV = powf(TS, 1.0f - l); tsD = powf(TS, m);
    lamW = l; muW = m; wrec = true;              // history and I untouched
  }
}

void bumplessStart(float e) {
  for (int j = 0; j <= M; j++) eh[j] = e;        // no derivative kick
  I = 0; hystState = false;
}

void supervise(float e) {
  if (!supOn) { wfar = 0; return; }
  float a = fabsf(e);
  if (inside) { if (a > eLow + hs) inside = false; }
  else if (a < eLow) inside = true;
  wfar = inside ? 0.0f : constrain((a - eLow) / (eHigh - eLow), 0.0f, 1.0f);
  Kp  = KpN + wfar * (KpF - KpN);
  lam = lamN + wfar * (lamF - lamN);
  mu  = muN + wfar * (muF - muN);
  Kff = KffN + wfar * (KffF - KffN);
  if (e > 0 && !inside) Kp *= aggr;
}

float controller(float e) {
  if (mode == CTRL_OPEN) return pwmOpen;
  if (mode == CTRL_HYST) {
    if (e < -h) hystState = true; else if (e > h) hystState = false;
    return hystState ? pwmOn : 0;
  }
  // P, PI, PID and FOPID share one law; P/PI/PID use integer orders.
  float l = (mode == CTRL_FOPID) ? lam : 1.0f;
  float m = (mode == CTRL_FOPID) ? mu : 1.0f;
  updateWeights(l, m);
  for (int j = M; j > 0; j--) eh[j] = eh[j - 1];
  eh[0] = e;
  float sV = 0, sD = 0;
  for (int j = 0; j <= M; j++) { sV += wV[j] * eh[j]; sD += wD[j] * eh[j]; }
  float v = -sV / tsV;                           // sign: plant gain is negative
  float D = -(sD / (1.0f + fabsf(sD))) / tsD;
  bool useI = (mode != CTRL_P), useD = (mode == CTRL_PID || mode == CTRL_FOPID);
  float ki = useI ? Ki : 0, kd = useD ? Kd : 0;
  float uP = -Kp * e, ff = Kff * tref;
  float uTry = uP + ki * (I + TS * v) + kd * D + ff;
  // ramp-aware conditional integration
  bool hi = ((uTry > 255) || (pwm < pwmTarget)) && (v > 0);
  bool lo = ((uTry < 0) || (pwm > pwmTarget)) && (v < 0);
  if (useI && !hi && !lo) {
    I += TS * v;
    float lim = 255.0f / fmaxf(ki, 1e-6f);
    I = constrain(I, -lim, lim);
  }
  return uP + ki * I + kd * D + ff;
}

void applyOutput(float r) {
  if (!control) { pwmTarget = 0; pwm = 0; }      // off: no ramp down
  else {
    pwmTarget = constrain((int)r, 0, 255);
    if (pwm < pwmTarget) pwm = min(pwm + RAMP, pwmTarget);
    else if (pwm > pwmTarget) pwm = max(pwm - RAMP, pwmTarget);
  }
  if (fanMode == FAN_ON_M) fan = 255;
  else if (fanMode == FAN_OFF_M) fan = 0;
  else fan = (pwm > FAN_THRESHOLD) ? 255 : 0;
  analogWrite(PIN_CELL, pwm);
  analogWrite(PIN_FAN, fan);
}

// =================================================================
// OBSERVER (B29) and ALARM (B30)
// =================================================================
void observe() {
  if (!obsInit) {
    xp[0] = tsup; xp[2] = th; xp[1] = 0.5f * (tsup + th);
    for (int i = 0; i <= DELAY_STEPS; i++) ubuf[i] = 0;
    obsInit = true;
  }
  innC = tsup - xp[0]; innH = th - xp[2];
  float x[3];
  for (int i = 0; i < 3; i++) x[i] = xp[i] + Lg[i][0] * innC + Lg[i][1] * innH;
  tmHat = x[1];
  int ud = ubuf[ubufIdx];                  // PWM applied DELAY_STEPS periods ago
  ubuf[ubufIdx] = pwm; ubufIdx = (ubufIdx + 1) % (DELAY_STEPS + 1);
  float ta = isnan(tamb) ? x[2] : tamb;
  for (int i = 0; i < 3; i++)
    xp[i] = Ad[i][0] * x[0] + Ad[i][1] * x[1] + Ad[i][2] * x[2] + bd[i] * ud + ed[i] * ta;
  // 10 s means of the innovations
  innCsum += innC - innCbuf[innIdx]; innHsum += innH - innHbuf[innIdx];
  innCbuf[innIdx] = innC; innHbuf[innIdx] = innH; innIdx = (innIdx + 1) % 100;
  bool out = fabsf(innCsum / 100) > bandC || fabsf(innHsum / 100) > bandH;
  outCount = out ? outCount + 1 : 0;
  alarm = outCount >= 50;                  // 5 s outside
}

// =================================================================
// SAFETY
// =================================================================
void safetyCheck() {
  bool bad = tsupRaw < T_MIN_OK || tsupRaw > T_MAX_OK || thRaw < T_MIN_OK || thRaw > T_MAX_OK;
  if (bad && control) {
    control = false; applyOutput(0);
    Serial.println(F("# probe out of range: control disabled"));
  }
}

// =================================================================
// LOG -- column dictionary (Appendix B); columns are only ever appended
// =================================================================
const char *modeName() {
  switch (mode) { case CTRL_OPEN: return "OPEN"; case CTRL_HYST: return "HYST"; case CTRL_P: return "P";
    case CTRL_PI: return "PI"; case CTRL_PID: return "PID"; default: return "FOPID"; }
}

void logRow(unsigned long now, float e) {
  Serial.print(now - tStart); Serial.print(',');
  Serial.print(tamb, 2); Serial.print(','); Serial.print(rh, 2); Serial.print(',');
  Serial.print(tdew, 2); Serial.print(','); Serial.print(tref, 2); Serial.print(',');
  Serial.print(tsup, 2); Serial.print(','); Serial.print(th, 2); Serial.print(',');
  Serial.print(e, 3); Serial.print(','); Serial.print(pwm); Serial.print(',');
  Serial.print(fan > 0 ? 1 : 0); Serial.print(',');
  Serial.print(Kp, 4); Serial.print(','); Serial.print(Ki, 4); Serial.print(',');
  Serial.print(Kd, 4); Serial.print(','); Serial.print(lam, 4); Serial.print(',');
  Serial.print(mu, 4); Serial.print(','); Serial.print(Kff, 4); Serial.print(',');
  Serial.print(control ? modeName() : "OFF"); Serial.print(',');
  // appended columns
  Serial.print(tsupRaw, 2); Serial.print(','); Serial.print(thRaw, 2); Serial.print(',');
  Serial.print(hystState ? 1 : 0); Serial.print(',');
  Serial.print(loopUs); Serial.print(','); Serial.print(lateMs); Serial.print(',');
  Serial.print(req, 4); Serial.print(',');
  Serial.print(supOn ? 1 : 0); Serial.print(','); Serial.print(wfar, 3); Serial.print(',');
  Serial.print(wrec ? 1 : 0); Serial.print(',');
  Serial.print(tmHat, 2); Serial.print(','); Serial.print(innC, 3); Serial.print(',');
  Serial.print(innH, 3); Serial.print(','); Serial.print(alarm ? 1 : 0); Serial.print(',');
  Serial.print(tdewRaw, 2); Serial.print(','); Serial.print(bdew, 2); Serial.print(',');
  Serial.print(refMode == REF_DEW ? "DEW" : "FIXED"); Serial.print(',');
  Serial.print(ambFlag ? 1 : 0); Serial.print(',');
  Serial.print(qdew, 4); Serial.print(','); Serial.println(nBelow);
}

void printHeader() {
  Serial.println(F("t_ms,Tamb,RH,Tdew,Tref,Tsup,Th,e,pwm,fan,Kp,Ki,Kd,lambda,mu,Kff,mode,"
                   "Tsup_raw,Th_raw,state,loop_us,late_ms,req,sup,wfar,wrec,Tm_hat,inn_c,inn_h,"
                   "alarm,Tdew_raw,bdew,ref,amb_flag,Qdew,n_below"));
}

// =================================================================
// SERIAL COMMANDS (Appendix B)
// =================================================================
float argf(const String &s) { int c = s.indexOf(','); return s.substring(c + 1).toFloat(); }

void setSup(const String &s) {             // SET_SUP,<name>,<value>
  int c1 = s.indexOf(','), c2 = s.indexOf(',', c1 + 1);
  String n = s.substring(c1 + 1, c2); float v = s.substring(c2 + 1).toFloat();
  if (n == "ELOW") eLow = v; else if (n == "EHIGH") eHigh = v; else if (n == "HS") hs = v;
  else if (n == "A") aggr = v;
  else if (n == "KPN") KpN = v; else if (n == "LAMN") lamN = v; else if (n == "MUN") muN = v;
  else if (n == "KFFN") KffN = v;
  else if (n == "KPF") KpF = v; else if (n == "LAMF") lamF = v; else if (n == "MUF") muF = v;
  else if (n == "KFFF") KffF = v;
  else { Serial.println(F("# unknown supervisor value")); return; }
  Serial.print(F("# SUP ")); Serial.print(n); Serial.print('='); Serial.println(v, 4);
}

void commands() {
  if (!Serial.available()) return;
  String s = Serial.readStringUntil('\n'); s.trim();
  if (s.length() == 0) return;
  float e = tref - tsup;
  if (s == "CALIBRATE") { calibrating = true; control = false; }
  else if (s == "CONTROL_ON") {
    if (calibrating) { Serial.println(F("# refused: calibration pending (device not yet at equilibrium)")); return; }
    control = true; bumplessStart(e);
  }
  else if (s == "CONTROL_OFF") { control = false; applyOutput(0); }
  else if (s.startsWith("SET_PWM,")) { pwmOpen = constrain((int)argf(s), 0, 255); mode = CTRL_OPEN; }
  else if (s == "FAN_AUTO") fanMode = FAN_AUTO_M;
  else if (s == "FAN_ON") fanMode = FAN_ON_M;
  else if (s == "FAN_OFF") fanMode = FAN_OFF_M;
  else if (s.startsWith("MODE,")) {
    String m = s.substring(5);
    if (m == "HYST") mode = CTRL_HYST; else if (m == "P") mode = CTRL_P; else if (m == "PI") mode = CTRL_PI;
    else if (m == "PID") mode = CTRL_PID; else if (m == "FOPID") mode = CTRL_FOPID;
    else { Serial.println(F("# unknown mode")); return; }
    bumplessStart(e);
  }
  else if (s.startsWith("SET_H,")) h = argf(s);
  else if (s.startsWith("SET_PWMON,")) pwmOn = constrain((int)argf(s), 0, 255);
  else if (s.startsWith("SET_KP,")) Kp = KpN = argf(s);
  else if (s.startsWith("SET_KI,")) Ki = argf(s);
  else if (s.startsWith("SET_KD,")) Kd = argf(s);
  else if (s.startsWith("SET_LAMBDA,")) lam = argf(s);
  else if (s.startsWith("SET_MU,")) mu = argf(s);
  else if (s.startsWith("SET_KFF,")) Kff = argf(s);
  else if (s == "SUP_ON") { supOn = true; inside = true; }
  else if (s == "SUP_OFF") supOn = false;
  else if (s.startsWith("SET_SUP,")) setSup(s);
  else if (s == "REF,DEW") refMode = REF_DEW;
  else if (s == "REF,FIXED") refMode = REF_FIXED;
  else if (s.startsWith("SET_TREF,")) { tset = argf(s); refMode = REF_FIXED; I = 0; }
  else if (s.startsWith("SET_BIAS,")) bT = argf(s);
  else if (s.startsWith("SET_BDEW,")) bdew = argf(s);
  else if (s.startsWith("SET_BAND,")) {   // SET_BAND,<cold>,<hot>
    int c1 = s.indexOf(','), c2 = s.indexOf(',', c1 + 1);
    bandC = s.substring(c1 + 1, c2).toFloat(); bandH = s.substring(c2 + 1).toFloat();
  }
  else if (s == "RESET_Q") { qdew = 0; nBelow = 0; }
  else if (s == "HEADER") printHeader();
  else { Serial.print(F("# unknown command: ")); Serial.println(s); return; }
  Serial.print(F("# ok ")); Serial.println(s);
}

// =================================================================
void setup() {
  pinMode(PIN_CELL, OUTPUT); digitalWrite(PIN_CELL, LOW);   // outputs off first
  pinMode(PIN_FAN, OUTPUT);  digitalWrite(PIN_FAN, LOW);
  Serial.begin(115200);
  analogReadResolution(ADC_BITS);         // remove this line on AVR boards
  Wire.begin();
  float s0 = countsToC(readAvg(PIN_TSUP)), h0 = countsToC(readAvg(PIN_TH));
  tsup = tsupRaw = s0; th = thRaw = h0;
  for (int i = 0; i < 3; i++) { medS[i] = s0; medH[i] = h0; }
  tref = tset;
  updateWeights(lam, mu);
  bumplessStart(0);
  calibrating = true;                     // calibrates as soon as it is at equilibrium
  printHeader();
  tStart = tLast = millis();
}

void loop() {
  unsigned long now = millis();
  serviceAmbient(now);                    // never blocks
  commands();
  if (now - tLast < TS_MS) return;
  lateMs = (long)(now - tLast) - (long)TS_MS;
  tLast += TS_MS;                         // the period does not drift
  unsigned long t0 = micros();

  measure();
  serviceCalibration();
  updateReference();
  float e = tref - tsup;
  supervise(e);
  req = control ? controller(e) : 0;
  applyOutput(req);
  safetyCheck();
  observe();
  if (!ambFlag && !isnan(tdew)) {
    qdew += fmaxf(0.0f, tdew - tsup) * TS / 3600.0f;
    if (tsup < tdew) nBelow++;
  }
  loopUs = micros() - t0;
  logRow(now, e);
}
