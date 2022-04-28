## Template to simulate PKPD data
## Fribgerg-Karlsson population model
## Simulate ME-2 plasma concentrations and ANC values
## using mrgsolve.

library(mrgsolve)
library(rstan)
library(rjson)
library(dplyr)

modelName <- "EffectData_5sub"

nSub <- 5; # number of subjects
nIIV <- 5; # number of parameters with inter-individual variations

code <- '
$PARAM CL = 10, Q = 15, VC = 35, VP = 105, KA = 2.0, ke0 = 1.5, EC50 = 120,
WT = 70

$SET delta=0.1 // simulation grid

$CMT GUT CENT PERI EFFECT

$MAIN
// Individual parameters
double CLi = exp(log(CL) + 0.75*log(WT/70) + ETA(1));
double Qi = exp(log(Q) + 0.75*log(WT/70) + ETA(2));
double VCi = exp(log(VC) + log(WT/70) + ETA(3));
double VPi = exp(log(VP) + log(WT/70) + ETA(4));
double ke0i = exp(log(ke0) + ETA(5));

// Reparametrization
double k10 = CLi / VCi;
double k12 = Qi / VCi;
double k21 = Qi / VPi;

$ODE 
dxdt_GUT = -KA * GUT;
dxdt_CENT = KA * GUT - (k10 + k12) * CENT + k21 * PERI;
dxdt_PERI = k12 * CENT - k21 * PERI;
dxdt_EFFECT = ke0 * CENT - ke0 * EFFECT;

$OMEGA name="IIV"
0.0625 0.0625 0.0625 0.0625 0.125

$SIGMA 0.01 0.005

$TABLE
double CP = CENT/VCi;
double DV1 = CENT/VCi * exp(EPS(1));
double DV2 = 100 * ((EFFECT/VCi) / (EC50 + (EFFECT/VCi))) * exp(EPS(2));
double WEIGHT = WT;

$CAPTURE CP DV1 DV2 WEIGHT
'

mod <- mread("acum", tempdir(), code)
e1 <- expand.ev(amt = 80 * 1000, addl = 14, ii = 12, WT = rnorm(nSub, 70, 10))
out <- mod %>% data_set(e1) %>% carry.out(dose) %>% Req(CP,DV1,DV2) %>% mrgsim(end=500)


## Observation and dosing times
doseTimes <- seq(0.0, 168.0, by = 12)
xpk <- c(0, 0.083, 0.167, 0.25, 0.5, 0.75, 1, 1.5, 2,3,4,6,8)
xpk <- c(xpk, xpk + 12, seq(24, 156, by = 12), 
         c(xpk, 12, 18, 24, 30, 36) + 168)
time <- xpk

mod %>% data_set(e1) %>%
  carry.out(cmt, ii, addl, rate, amt, evid, ss) %>%
  mrgsim(Req = "CP, DV1, DV2, WEIGHT", end = -1, add = time, rescort = 3) %>%
  plot(DV1 ~ time)

mod %>% data_set(e1) %>%
  carry.out(cmt, ii, addl, rate, amt, evid, ss) %>%
  mrgsim(Req = "CP, DV1, DV2, WEIGHT", end = -1, add = time, rescort = 3) %>%
  plot(DV2 ~ time)

## Assemble data set for Stan
xdata <- mod %>% data_set(e1) %>%
  carry.out(cmt, ii, addl, rate, amt, evid, ss) %>%
  mrgsim(Req = "DV1, DV2, WEIGHT", end = -1, add = time, rescort = 3) %>%
  as.data.frame

xdata <- xdata %>%
  mutate(DV1 = ifelse(time %in% xpk & time != 0 & evid == 0, DV1, NA),
         DV2 = ifelse(time %in% xpk & time != 0 & evid == 0, DV2, NA))


xdata$cmt[xdata$cmt == 0] <- 2 #switch from mrgsolve to NONMEM convention

nt <- length(time)
N <- nrow(xdata) #total number of observations

## Added by Nina 
timePK <- xdata$time[iObsPK]
startPK <- rep(0, nSub) ;endPK <- rep(0, nSub)
j = 2 ; startPK[1] = 1 ; endPK[nSub] = length(timePK);
for (i in 2:length(timePK)) {
  if (timePK[i] < timePK[i-1]) { 
    startPK[j] <- i
    if (j >= 2) endPK[j-1] <- i-1
    j = j + 1
  }
}

dose <- xdata$amt 

## Subject specific data
start <- (1:N)[!duplicated(xdata$ID)] 
end <- c(start[-1] -1, N)
xsub <- subset(xdata, !duplicated(ID))
weight <- xsub$WEIGHT

## Look at simulated data using plots
p1 <- ggplot(xdata %>% filter(!is.na(DV1)), aes(x = time, y = DV1))
p1 + geom_point() + geom_line() +
  labs(x = "time (h)", y = "ME-2 plasma concentration (ng/mL)") +
  theme(text = element_text(size = 12), axis.text = element_text(size = 12),
        legend.position = "none", strip.text = element_text(size = 8)) +
  facet_wrap(~ID)

p1 <- ggplot(xdata %>% filter(!is.na(DV2)), aes(x = time, y = DV2))
p1 + geom_point() + geom_line() +
  labs(x = "time (h)",
       y = "effect") +
  theme(text = element_text(size = 12), axis.text = element_text(size = 12),
        legend.position = "none", strip.text = element_text(size = 8)) +
  facet_wrap(~ID)

## Indices of records containing observed concentrations
iObsPK <- with(xdata, (1:nrow(xdata))[!is.na(DV1) & evid == 0])
nObsPK <- length(iObsPK)



## Parameters for priors (only used for informed initial estimates)
CLPrior = 10
QPrior = 15
V1Prior = 35
V2Prior = 105
kaPrior = 2
ke0Prior = 2
EC50Prior = 120 
CLPriorCV = 0.10
QPriorCV = 0.18
V1PriorCV = 0.14
V2PriorCV = 0.17
kaPriorCV = 0.16
ke0PriorCV = 0.15
EC50PriorCV = 0.18


## create data set
data <- with(xdata,
             list(
               nt = nt,
               N= N,
               nObsPK = nObsPK,
               iObsPK = iObsPK,
               amt = amt,
               cmt = cmt,
               evid = evid,
               time = time,
               timePK = timePK,
               startPK = startPK,
               endPK = endPK,
               ii = ii,
               addl = addl,
               ss = ss,
               rate = rate,
               dose = dose,
               cObs = DV1[iObsPK],
               effObs = DV2[iObsPK],
  
               nSubjects = nSub,
               nIIV = nIIV,
               start = start,
               end = end,
               weight = weight,
               
               # Priors for effect
               ke0Prior = ke0Prior,
               ke0PriorCV = ke0PriorCV,
               EC50Prior = EC50Prior,
               EC50PriorCV = EC50PriorCV
             ))

## create initial estimates
init <- function(){
  list(CL_pop = abs(rnorm(1, 0, 20)),
       Q_pop = abs(rnorm(1, 0, 20)),
       VC_pop = abs(rnorm(1, 0, 100)),
       VP_pop = abs(rnorm(1, 0, 200)),
       ka_pop = abs(rnorm(1, 0, 5)),
       sigma_c = 0.2,
       sigma_eff = 0.2,
       ke0_pop = exp(rnorm(1, log(ke0Prior), ke0PriorCV)),
       EC50_pop = exp(rnorm(1, log(EC50Prior), EC50PriorCV)),
       L = diag(nIIV),
       omega = exp(rnorm(nIIV, log(0.05), 0.5)),
       etaStd = matrix(rep(0, nIIV * nSub), nrow = nIIV))
}

with(data, stan_rdump(ls(data), file = paste0(modelName,".data.R")))
jsonData <- toJSON(data)
write(jsonData, paste0(modelName, ".data", ".json"))

inits <- init()
with(inits, stan_rdump(ls(inits), file = paste0(modelName,".init.R")))
jsonInits <- toJSON(inits)
write(jsonInits, paste0(modelName, ".init", ".json"))
