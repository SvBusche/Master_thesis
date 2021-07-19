library(EcoAgents)
library(snow)
library(parallel)

# Create a simulation object
sim <- init_simulation(cbind(c(-100, -100, 100, 100), c(-100, 100, 100, -100)),
                       gridFieldSize = 1,
                       gridFieldLayers = 5)

# add organisms (15 cells of each as 'starter culture')
models <- list()
models[['V.parv']] <- readRDS("/Users/Svenja/Desktop/gapseq/V.parv/V.parv_1.RDS")
models[['F.peri']] <- readRDS("/Users/Svenja/Desktop/gapseq/F.peri/F.peri1.RDS")
models[['P.mela']] <- readRDS("/Users/Svenja/Desktop/gapseq/P.mela/P.mela1.RDS")
models[['R.muci']] <- readRDS("/Users/Svenja/Desktop/gapseq/Rothia/R.mucila_1.RDS")
models[['H.para']] <- readRDS("/Users/Svenja/Desktop/gapseq/H.parainfluenzae/H.para1.RDS")
models[['S.vest']] <- readRDS("/Users/Svenja/Desktop/gapseq/S.vestibularis/S.vestibularis_1.RDS")


#composition of the microbiome
#colonization experiment
sim <- add_organism(sim, model = models[["S.vest"]], name = "S.vest", ncells = 1, distribution.radius = 10, open.bounds = 0.5)
sim <- add_organism(sim, model = models[["P.mela"]], name = "P.mela", ncells = 1, distribution.radius = 10, open.bounds = 0.5)
sim <- add_organism(sim, model = models[["V.parv"]], name = "V.parv", ncells = 1, distribution.radius = 10, open.bounds = 0.5)
sim <- add_organism(sim, model = models[["F.peri"]], name = "F.peri", ncells = 1, distribution.radius = 10, open.bounds = 0.5)
sim <- add_organism(sim, model = models[["R.muci"]], name = "R.muci", ncells = 1, distribution.radius = 10, open.bounds = 0.5)
sim <- add_organism(sim, model = models[["H.para"]], name = "H.para", ncells = 1, distribution.radius = 10, open.bounds = 0.5)

#influence of the diet with a composition found in a study
sim <- add_organism(sim, model = models[["S.vest"]], name = "S.vest", ncells = 46, distribution.radius = 10, open.bounds = 0.5)
sim <- add_organism(sim, model = models[["P.mela"]], name = "P.mela", ncells = 21, distribution.radius = 10, open.bounds = 0.5)
sim <- add_organism(sim, model = models[["V.parv"]], name = "V.parv", ncells = 11, distribution.radius = 10, open.bounds = 0.5)
sim <- add_organism(sim, model = models[["F.peri"]], name = "F.peri", ncells = 9, distribution.radius = 10, open.bounds = 0.5)
sim <- add_organism(sim, model = models[["R.muci"]], name = "R.muci", ncells = 4, distribution.radius = 10, open.bounds = 0.5)
sim <- add_organism(sim, model = models[["H.para"]], name = "H.para", ncells = 9, distribution.radius = 10, open.bounds = 0.5)

# Plot distribution
p <- plot_cells(sim, xlim = c(-100, 100), ylim= c(-100, 100), iter=144)
p


#step to make the simulation run on Mac
cl <- parallel::makeCluster(2, setup_strategy = "sequential")
if (Sys.getenv("RSTUDIO") == "1" && !nzchar(Sys.getenv("RSTUDIO_TERM")) && 
    Sys.info()["sysname"] == "Darwin" && getRversion() >= "4.0.0") {
  parallel:::setDefaultClusterOptions(setup_strategy = "sequential")
}

#starting with the mealplan 
#BREAKFAST
#mediterranean diet
dt_med.diet_2.5l <- fread("med.diet_2500ml_3.csv")
dt_med.diet_2.5l
sim <- add_compounds(sim,
                     compounds = dt_med.diet_2.5l$cpd.id,
                     concentrations = dt_med.diet_2.5l$mM,
                     compound.names = dt_med.diet_2.5l$cpd.name,
                     is.constant = dt_med.diet_2.5l$is.constant)
#or lowcarb diet
dt_lowcarb_2.5l <- fread("lowcarbDiet_2500ml_3.csv")
dt_lowcarb_2.5l
sim <- add_compounds(sim,
                     compounds = dt_lowcarb_2.5l$cpd.id,
                     concentrations = dt_lowcarb_2.5l$mM,
                     compound.names = dt_lowcarb_2.5l$cpd.name,
                     is.constant = dt_lowcarb_2.5l$is.constant)

update_sim <- function(sim) {
  if(sim@n_rounds == 27) { # "%%" ist der Modulo operator
    sim <- dilute_compounds(sim, dilution.factor = 0.5)
  }
  return(sim)
}
sim <- run_simulation(sim, niter = 36, lim_time = 300,
                      n.cores = 4, verbose = 1, on.iteration = update_sim)

#LUNCH
#mediterranean diet
sim <- add_compounds(sim,
                     compounds = dt_med.diet_2.5l$cpd.id,
                     concentrations = dt_med.diet_2.5l$mM,
                     compound.names = dt_med.diet_2.5l$cpd.name,
                     is.constant = dt_med.diet_2.5l$is.constant)
#or lowcarb diet
sim <- add_compounds(sim,
                     compounds = dt_lowcarb_2.5l$cpd.id,
                     concentrations = dt_lowcarb_2.5l$mM,
                     compound.names = dt_lowcarb_2.5l$cpd.name,
                     is.constant = dt_lowcarb_2.5l$is.constant)

update_sim <- function(sim) {
  if(sim@n_rounds == 63) { # "%%" ist der Modulo operator
    sim <- dilute_compounds(sim, dilution.factor = 0.5)
  }
  return(sim)
}
sim <- run_simulation(sim, niter = 36, lim_time = 300,
                      n.cores = 4, verbose = 1, on.iteration = update_sim)

#DINNER
#mediterranean diet
sim <- add_compounds(sim,
                     compounds = dt_med.diet_2.5l$cpd.id,
                     concentrations = dt_med.diet_2.5l$mM,
                     compound.names = dt_med.diet_2.5l$cpd.name,
                     is.constant = dt_med.diet_2.5l$is.constant)

#or lowcarb diet
sim <- add_compounds(sim,
                     compounds = dt_lowcarb_2.5l$cpd.id,
                     concentrations = dt_lowcarb_2.5l$mM,
                     compound.names = dt_lowcarb_2.5l$cpd.name,
                     is.constant = dt_lowcarb_2.5l$is.constant)
update_sim <- function(sim) {
  if(sim@n_rounds == 90) { # "%%" ist der Modulo operator
    sim <- dilute_compounds(sim, dilution.factor = 0.5)
  }
  return(sim)
}
sim <- run_simulation(sim, niter = 72, lim_time = 500,
                      n.cores = 4, verbose = 1, on.iteration = update_sim)

#evaluation growth
show(sim)
plot_growth(sim, tlim = NULL)
p <- plot_cells(sim, xlim = c(-100, 100), ylim= c(-100, 100))
p

#change of the organic acids/SCFA concentrations over time
plot_compounds(sim, compounds = c("cpd00029_e0","cpd00159_e0", "cpd00047_e0",
                                  "cpd00141_e0","cpd00036_e0","cpd00106_e0"))

plot_compounds(sim, compounds = c("cpd00263_e0"), tlim = c(0,6))
plot_compounds(sim, compounds = c("cpd00305_e0"))

#global metabolite concentrations over time
history28 <- sim@history[[28]]$global_compounds
View(history28)
fwrite(history, file = "sim_concentrations_med.diet_28.csv")
history29 <- sim@history[[29]]$global_compounds
fwrite(history, file = "sim_concentrations_med.diet_29.csv")
View(history29)
history36 <- sim@history[[36]]$global_compounds
View(history36)
fwrite(history, file = "sim_concentrations_med.diet_36.csv")
history37 <- sim@history[[37]]$global_compounds
View(history37)
fwrite(history, file = "sim_concentrations_med.diet_37.csv")
history64 <- sim@history[[64]]$global_compounds
View(history64)
fwrite(history, file = "sim_concentrations_med.diet_64.csv")
history65 <- sim@history[[65]]$global_compounds
View(history65)
fwrite(history, file = "sim_concentrations_med.diet_65.csv")
history72 <- sim@history[[72]]$global_compounds
View(history72)
fwrite(history, file = "sim_concentrations_med.diet_72.csv")
history73 <- sim@history[[73]]$global_compounds
View(history73)
fwrite(history, file = "sim_concentrations_med.diet_73.csv")
history91 <- sim@history[[91]]$global_compounds
View(history91)
fwrite(history, file = "sim_concentrations_med.diet_91.csv")
history92 <- sim@history[[92]]$global_compounds
View(history92)
fwrite(history, file = "sim_concentrations_med.diet_92.csv")
history144 <- sim@history[[144]]$global_compounds
View(history144)
fwrite(history, file = "sim_concentrations_med.diet_144.csv")
# plot the environment with the different B-vitamins
plot_environment(sim, compounds = c("cpd00220_e0","cpd00393_e0", "cpd00218_e0","cpd00644_e0"),
                 xlim = c(-50,50), ylim= c(-50,50))

plot_environment(sim, compounds = c("cpd00263_e0"),
                 xlim = c(-100,100), ylim= c(-100,100))
plot_environment(sim, compounds = c("cpd00644_e0"),
                 xlim = c(-100,100), ylim= c(-100,100))
plot_environment(sim, compounds = c("cpd00218_e0"),
                 xlim = c(-100,100), ylim= c(-100,100))
plot_environment(sim, compounds = c("cpd00393_e0"),
                 xlim = c(-100,100), ylim= c(-100,100))
plot_environment(sim, compounds = c("cpd00305_e0"),
                 xlim = c(-100,100), ylim= c(-100,100))
plot_environment(sim, compounds = c("cpd00220_e0"),
                 xlim = c(-100,100), ylim= c(-100,100))


plot_environment(sim, compounds = c("cpd00218_e0"),
                 xlim = c(-50,50), ylim= c(-50,50))


# plot lactate, acetate, formate, Propionate
plot_environment(sim, compounds = c("cpd00029_e0","cpd00159_e0", "cpd00047_e0",
                                    "cpd00141_e0"),
                 xlim = c(-50,50), ylim= c(-50,50))
plot_environment(sim, compounds = c("cpd00159_e0"
                                    ),
                 xlim = c(-50,50), ylim= c(-50,50))


#get the exchanges at different times depending on the diet chosen for the simulation
#preset for mediterranean diet in a file
sum <- summary_exchanges(sim, iter= 37)
fwrite(sum, file = "sim_Meddiet2.5l_after_lunch_sim4.csv")
sum <- summary_exchanges(sim, iter= 73)
fwrite(sum, file = "sim_Meddiet2.5l_after_dinner_sim4.csv")
sum <- summary_exchanges(sim, iter= 144)
fwrite(sum, file = "sim_Meddiet2.5l_after_24h_sim4.csv")

#preset for lowcarb diet in a file
sum <- summary_exchanges(sim, iter= 37)
fwrite(sum, file = "sim_lowcarb2.5l_after_lunch_sim2.csv")
sum <- summary_exchanges(sim, iter= 73)
fwrite(sum, file = "sim_lowcarb2.5l_after_dinner_sim2.csv")
sum <- summary_exchanges(sim, iter= 144)
fwrite(sum, file = "sim_lowcarb2.5l_after_24_sim2.csv")

# get uptake/production summary
dt <- summary_exchanges(sim, iter = NULL)
dt[grepl("Lactate", compound.name)]
dt[grepl("Butyrate", compound.name)]
dt[grepl("Acetate", compound.name)]
dt[grepl("Propionate", compound.name)]
dt[grepl("Formate", compound.name)]
dt[grepl("Succinate", compound.name)]
dt[grepl("Fumarate", compound.name)]
dt[grepl("citrate", compound.name)]
dt[grepl("malate", compound.name)]
dt[grepl("cpd00137", compound)]
dt[grepl("cpd00130", compound)]
dt[grepl("ethanol", compound.name)]

#b-vitamins
dt[grepl("thiamine", compound.name)]
dt[grepl("riboflavin", compound.name)]
dt[grepl("niacin", compound.name)]
dt[grepl("pantothenic acid", compound.name)]
dt[grepl("pyridoxine", compound.name)]
dt[grepl("adenosylcobalamine", compound.name)]
dt[grepl("biotin", compound.name)]
dt[grepl("folic acid", compound.name)]
dt[grepl("cpd11606_e0", compound)]
dt[grepl("cpd03226_e0", compound)]
dt[grepl("cpd01628_e0", compound)]
dt[grepl("retinol", compound.name)]
dt[grepl("nicotinamide", compound.name)]
dt[grepl("cpd00300_e0", compound)]

#amino acids
dt[grepl("alanine", compound.name)]
dt[grepl("arginine", compound.name)]
dt[grepl("asparagine", compound.name)]
dt[grepl("aspartate", compound.name)]
dt[grepl("cysteine", compound.name)]
dt[grepl("glutamate", compound.name)]
dt[grepl("glutamine", compound.name)]
dt[grepl("cpd00053", compound)]
dt[grepl("cpd00033", compound)]
dt[grepl("glycin", compound.name)]
dt[grepl("histidine", compound.name)]
dt[grepl("isoleucine", compound.name)]
dt[grepl("leucine", compound.name)]
dt[grepl("lysine", compound.name)]
dt[grepl("methionine", compound.name)]
dt[grepl("phenylalanine", compound.name)]
dt[grepl("proline", compound.name)]
dt[grepl("serine", compound.name)]
dt[grepl("threonine", compound.name)]
dt[grepl("tryptophan", compound.name)]
dt[grepl("tyrosine", compound.name)]
dt[grepl("valine", compound.name)]

#sugars
dt[grepl("cpd00105_e0", compound)]
dt[grepl("cpd11746_e0", compound)]
dt[grepl("fructose", compound.name)]
dt[grepl("glucose", compound.name)]
dt[grepl("maltose", compound.name)]
dt[grepl("sucrose", compound.name)]
dt[grepl("galactose", compound.name)]
dt[grepl("lactose", compound.name)]
dt[grepl("cellulose", compound.name)]
dt[grepl("starch", compound.name)]
dt[grepl("sorbitol", compound.name)]
dt[grepl("D-ribose", compound.name)]


#fatty acids
dt[grepl("cpd01107_e0", compound)]
dt[grepl("cpd01741_e0	", compound)]
dt[grepl("cpd03847_e0", compound)]
dt[grepl("cpd05237_e0", compound)]
dt[grepl("cpd16351_e0", compound)]
dt[grepl("cpd00214_e0", compound)]
dt[grepl("cpd24916_e0", compound)]
dt[grepl("cpd01080_e0", compound)]
dt[grepl("cpd00536_e0", compound)]
dt[grepl("cpd03850_e0", compound)]
dt[grepl("cpd15016_e0", compound)]
dt[grepl("cpd01122_e0", compound)]
dt[grepl("cpd03848_e0", compound)]
dt[grepl("ccpd16341_e0", compound)]
dt[grepl("cpd16340_e0", compound)]
dt[grepl("cpd00188_e0", compound)]
dt[grepl("cpd05196_e0", compound)]
dt[grepl("cpd05231_e0", compound)]
dt[grepl("cpd16342_e0", compound)]
dt[grepl("cpd16301_e0", compound)]
dt[grepl("cpd03852_e0", compound)]
dt[grepl("cpd03846_e0", compound)]
dt[grepl("cpd00160_e0", compound)]
dt[grepl("cpd05274_e0", compound)]
dt[grepl("cpd05235_e0", compound)]
dt[grepl("cholesterol", compound.name)]
dt[grepl("Palmitate", compound.name)]
dt[grepl("Dodeconate", compound.name)]

dt[grepl("nicotinamide", compound.name)]
dt[grepl("GABA", compound.name)]

