set.seed(0)
rm(list = ls())
gc()

modelName <- "sim_"

library(rstan)
library(mrgsolve)
library(rjson)
library(dplyr)



code <- '
$PARAM CL = 7, V1 = 65, KA = 2.5

$SET delta = 1  // simulation grid

$CMT GUT CENT 

$GLOBAL
#define CP (CENT/V1)

$MAIN
double CLi = exp(log(CL) + ETA(1));
double V1i = exp(log(V1) + ETA(2));

// Reparametrization
double k10 = CLi / V1i;

$ODE 
dxdt_GUT = -KA * GUT;
dxdt_CENT = KA * GUT - k10 * CENT;

$OMEGA  name="IIV"
0.3 0.3 0.3

$SIGMA 0.1

$TABLE
capture DV = CP * exp(EPS(1));

$CAPTURE DV CP
'
mod <- mread("popex", tempdir(), code)

pk_data <- read.table('pk_data.csv', header=T, sep=',')
pk_data['CMT'] = 1
pk_data <- subset(pk_data, select=-c(DV, DATE, X))

out <- mod %>%
  data_set(pk_data) %>%
  mrgsim(time=pk_data$TIME) 
dev.new(width = 550, height = 330, unit = "px")
plot(out, DV~time | factor(ID))

data <- with(out, list(ID = out$ID,
                       TIME = out$TIME,
                       DV = out$DV))
#Save
data_pop_dir = "pop_data/"
with(data, stan_rdump(ls(data), file = paste0(data_pop_dir, modelName,"pop.data.R")))
jsonData <- toJSON(data)
write(jsonData, paste0(data_pop_dir, modelName, "pop", ".json"))

