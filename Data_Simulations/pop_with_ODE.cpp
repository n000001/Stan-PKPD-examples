$PARAM CL = 10, Q = 15, VC = 35, VP = 105, KA = 2.0, WT = 70

$SET delta=0.1 // simulation grid

$CMT GUT CENT PERI 

$MAIN
// Individual PK parameters
double CLi = exp(log(CL) + 0.75*log(WT/70) + ETA(1));
double Qi = exp(log(Q) + 0.75*log(WT/70) + ETA(2));
double VCi = exp(log(VC) + log(WT/70) + ETA(3));
double VPi = exp(log(VP) + log(WT/70) + ETA(4));

// Reparametrization
double k10 = CLi / VCi;
double k12 = Qi / VCi;
double k21 = Qi / VPi;

$ODE 
dxdt_GUT = -KA * GUT;
dxdt_CENT = KA * GUT - (k10 + k12) * CENT + k21 * PERI;
dxdt_PERI = k12 * CENT - k21 * PERI;

$OMEGA name="IIV"
0.0625 0.0625 0.0625 0.0625 

$SIGMA 0.01 0.01

$TABLE
double CP = CENT/VCi;
double DV = CENT/VCi * exp(EPS(1));
double WEIGHT = WT;

$CAPTURE CP DV WEIGHT

