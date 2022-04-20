$PROB
# Example population PK model

$PARAM TVKA = 2.0, TVQ = 15, TVCL = 10, TVV2 = 35, TVV3 = 105

$PARAM @covariates
WT = 70

$PKMODEL cmt="GUT CENT PERI", depot=TRUE

$SET end=240, delta=0.5

$MAIN
double CL = exp(log(TVCL) + 0.75*log(WT/70) + ECL);
double Q = exp(log(TVQ) + 0.75*log(WT/70) + ECQ);
double V2  = exp(log(TVV2)  +      log(WT/70) + EV2 );
double V3  = exp(log(TVV3)  +      log(WT/70) + EV3 );
double KA = exp(log(TVKA));

$OMEGA @labels ECL ECQ EV2 EV3 
0.3 0.1 0.1 0.4 

$SIGMA 0

$TABLE
capture IPRED = CENT/V3;
capture DV = IPRED*exp(EPS(1));

$CAPTURE CL Q V2 V3 KA ECL

  
