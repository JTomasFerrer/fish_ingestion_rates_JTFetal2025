
2 de 10
Esborrany article Ingestion Rate
Safata d'entrada

Joaquim Tomàs Ferrer
Fitxers adjunts9 de juny 2025 10:46
Bon dia, A la fi tenc un esborrany de s'article d'ingesta lo suficientment madur com perquè pugui passar una ronda de revisions entre tots es coautors. Aquí vos
14

Joaquim Tomàs Ferrer
23 de juny 2025 13:11 (fa 5 dies)
Ok, afegesc en Pablo a sa llista de suggested reviewers. Ho pujam des des teu compte (ho puc fer jo si me dones ses teves credencials) o me faig un compte jo? M

MIGUEL PALMER VIDAL
Fitxers adjunts
23 de juny 2025 13:40 (fa 5 dies)
per a mi

Hola,

Adjunt petites modificacions al supp material (ok amb el nous 
exemples!) i al readme del github (alguns typos).

El codi de R (jtf_feeding_analysis.R) no funcionava perquè hi havia 
NAs. Ho he modificat, però has de fer un run i comprovar que te surt 
el mateix que el paper (descomptades petites desviacions degudes a que 
el procés d'estima de paràmetres és heurístic).

Respecte a la submissió, revisa el link al seminari de PROA que t'he 
enviat abans. Com ara ets csic, crec que tu pots ser corresponding 
author. Si és així tira endavant com a corresponding author. Però molt 
alerta en seguir les instruccions específiques per a plosone (al 
seminari PROA i a la plana de submissio de plosone) per tal de no 
haver de pagar.

m

 3 fitxers adjunts
  •  Escanejat per Gmail
#---------------------------------------
# Size-independent, between-individual variability in feed ingestion rate in 
# European seabass (Dicentrarchus labrax)
#
# Tomas-Ferrer, J., Moro-Martinez,I., Massuti-Pascual, E.,
# Grau, A., Palmer, M.
#---------------------------------------

#-----
# 1: loading libraries
#-----
remove(list=ls())
library(cmdstanr)

#-----
# 2: Loading data
#-----
load("input_jtf.RData")
ls()

# NF: number of fish
# NR: number of replicates per fish
# NC: number of cages
# CONSUMED: [NF+1,NR,NC] total number of pellets consumed by fish and by trial (cage*replicate)
#           last NF accounts for non-consumed pellets
# NP: maximum number of pellets delivered in any trial
# consumed [NF+1,NP,NR,NC] accumulated number of pellets consumed at each time step, by each fish, cage and replicate
# ID.list: [NF] ID for fish (code color)
# ID.cage: [NC] ID for cage (1 to 6)
# ID.replicate: [NR] Id for replicate (1 to 8)
# L2: [NF,NC] quadratic structural length (cm2)
# weight: [NF,NC] (gr)
# temperature: [NR,NC] (Celsius) 
# date: [NR,NC]
# food: [NR,NC] Ratio pellets monitored/pellets to be delivered
# TRAFFIC: [NR,NC] Number of boats passing through the port channel (next to the cages) during the trial
# diet [NC] diet level at each cage (1,2,3) (low, medium, large)
# included [NR,NC] valid replicates (excluding 2 non-valid trials) 
# last_row [NR,NC] number of pellets monitored at each trial

# Changing NAs by -999
NAs=which(is.na(CONSUMED),arr.ind = T)
CONSUMED[NAs]=-999

#-----       
# 3: STAN model
#-----
sink("model.stan")
cat(" // first line
functions {
}

data {

  int NF; // number of fish in a cage
  int NR; // number of rplicated trials per cage
  int NC; // number of cages
  array [(NF+1),NR,NC] int consumed; //cummulated number of pellets consumed by each fish
  array [NF,NC] real L2;  // quadratic structural length (or fish weight)
  array [NR,NC] real temperature;  // temperature
  array [NR,NC] real traffic;  // temperature
  array [NC] int diet;    //long term diet, as category
  array [NR,NC] real food;  // food actually delivered in a trial
  array [NR,NC] int included;  // Indicator for valid (1) or outlier (0) samples
  
  array [2] real prior_beta;
  array [2] real prior_beta_F_sd;
  array [2] real prior_beta_R_sd;
  array [2] real prior_beta_L2;
  array [2] real prior_beta_T;
  array [2] real prior_beta_T2;
  array [2] real prior_beta_TRAFFIC;
  array [2] real prior_beta_diet2;
  array [2] real prior_beta_diet3;
  array [2] real prior_beta_food;
}

transformed data {
}

parameters {
  real beta;
  real beta_diet2;
  real beta_diet3;
  
  array [NF,NC] real beta_F;
  real  <lower=0.0> beta_F_sd;
  array [NR,NC] real beta_R;
  real  <lower=0.0> beta_R_sd;
  real beta_L2;
  //real beta_T;
  //real beta_T2;
  real beta_TRAFFIC;
  real beta_food;
}

transformed parameters {
  array [3] real beta_delta;
  beta_delta[1]=0.0;
  beta_delta[2]=beta_diet2;
  beta_delta[3]=beta_diet3;
}   

model {

  beta ~ normal(prior_beta[1], prior_beta[2]);
  beta_diet2 ~ normal(prior_beta_diet2[1], prior_beta_diet2[2]);
  beta_diet3 ~ normal(prior_beta_diet3[1], prior_beta_diet3[2]);
  beta_F_sd ~ normal(prior_beta_F_sd[1], prior_beta_F_sd[2]);
  beta_R_sd ~ normal(prior_beta_R_sd[1], prior_beta_R_sd[2]);
  beta_L2 ~ normal(prior_beta_L2[1], prior_beta_L2[2]);
  //beta_T ~ normal(prior_beta_T[1], prior_beta_T[2]);
  //beta_T2 ~ normal(prior_beta_T2[1], prior_beta_T2[2]);
  beta_TRAFFIC ~ normal(prior_beta_TRAFFIC[1], prior_beta_TRAFFIC[2]);
  beta_food ~ normal(prior_beta_food[1], prior_beta_food[2]);
  
  for (cage in 1:NC){ // fish random effects 
    for (i in 1:NF){
      beta_F[i,cage] ~ normal(0.0, beta_F_sd);
    }
  }
  for (cage in 1:NC){ // replicate radom effects
    for (j in 1:NR){
      beta_R[j,cage] ~ normal(0.0, beta_R_sd);
    }
  }
  
  {//local
    vector [(NF+1)] prob;
    vector [(NF+1)] temp;
    
    for (cage in 1:NC){
      for (j in 1:NR){
      if (included[j, cage] == 1) {  // Only process valid samples
        for (i in 1:NF){
          // linear combination
          temp[i]=beta+beta_delta[diet[cage]]+
                  beta_F[i,cage]+
                  beta_R[j,cage]+
                  beta_L2*L2[i,cage]+
                  //beta_T*temperature[j,cage]+
                  //beta_T2*((temperature[j,cage])^2)+
                  beta_TRAFFIC*traffic[j,cage]+
                  beta_food*food[j,cage]
                  ;
        }
        temp[NF+1]=-sum(temp[1:NF]); // setting non consumed score
        prob=softmax(temp);
        consumed[1:(NF+1),j,cage] ~ multinomial(prob[1:(NF+1)]);
      }
      }
    }
  }//end local
}

generated quantities {
}

" # end of model code
,fill = TRUE)
sink()
#-----

#-----
# 4: Priors and initial values
#-----
# priors
prior_beta = c(0.0, 100.0);
prior_beta_F_sd = c(0.5,5.0);
prior_beta_R_sd = c(0.2,5.0);
prior_beta_L2 = c(0.0, 5.0);
prior_beta_T = c(0.0, 5.0);
prior_beta_T2 = c(0.0,5.0);
prior_beta_N = c(0.0,5.0);
prior_beta_TRAFFIC = c(0.0,5.0);
prior_beta_diet2 = c(0.0,5.0);
prior_beta_diet3 = c(0.0,5.0);
prior_beta_food = c(0.0,5.0);

n.chains = 2
initializer = function() list(
  "beta"=prior_beta[1],
  "beta_F_sd"=prior_beta_F_sd[1],
  "beta_F"=array(0,dim=c(NF,NC)),
  "beta_R_sd"=prior_beta_R_sd[1],
  "beta_R"=array(0,dim=c(NR,NC)),
  "beta_L2"=prior_beta_L2[1],
  #"beta_T"=prior_beta_T[1],
  #"beta_T2"=prior_beta_T2[1],
  "beta_TRAFFIC"=prior_beta_TRAFFIC[1],
  "beta_diet2"=prior_beta_diet2[1],
  "beta_diet3"=prior_beta_diet3[1],
  "beta_food"=prior_beta_food[1]
)
inits = list()
for (chain in 1:n.chains) inits[[chain]] = initializer()

#-----
# 5: compiling and running the model
#-----
mod = cmdstan_model("model.stan")

fit = mod$sample(
  data =list (
    NF=NF,
    NR=NR,
    NC=NC,
    consumed=CONSUMED,
    L2=(L2-mean(L2))/sd(L2),
    #L2=(weight-mean(weight))/sd(weight),
    temperature=(temperature-mean(temperature))/sd(temperature),
    traffic=(TRAFFIC-mean(TRAFFIC))/sd(TRAFFIC),
    diet=diet,
    #food=(last_row-mean(last_row))/sd(last_row),
    food=(food-mean(food))/sd(food),
    included=included,
    prior_beta = prior_beta,
    prior_beta_F_sd = prior_beta_F_sd,
    prior_beta_R_sd = prior_beta_R_sd,
    prior_beta_L2 = prior_beta_L2,
    prior_beta_T = prior_beta_T,
    prior_beta_T2 = prior_beta_T2,
    prior_beta_TRAFFIC = prior_beta_TRAFFIC,
    prior_beta_diet2 = prior_beta_diet2,
    prior_beta_diet3 = prior_beta_diet3,
    prior_beta_food = prior_beta_food
  ),
  chains = n.chains,
  parallel_chains = n.chains,
  iter_warmup = 2000,
  iter_sampling = 2000,
  init = inits,
  max_treedepth = 12,
  adapt_delta = 0.8
)

#-----
# 6: saving output
#-----
fit$save_object(file = paste("analysis_jtf",".RDS",sep=""))
sink(file = paste("analysis_jtf",".txt",sep=""))
fit$cmdstan_diagnose()
sink()
file.rename("model.stan",paste("analysis_jtf","model.R",sep=""))
#-----

#-----
# 7: Exploring results:
#-----
# remove(list=ls())
fit = readRDS("analysis_jtf.RDS")  # Results from Analysis.R script
fit$time()
draws=fit$draws(format="matrix")
