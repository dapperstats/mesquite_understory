library(lme4)
library(MCMCglmm)
library(viridis)
library(coda)
source("functions.R")

dens <- read.csv("exclosure_density.csv", stringsAsFactors = FALSE)
mesq <- read.csv("mesquite_cover.csv", stringsAsFactors = FALSE)

mdata <- prep_data(dens, mesq)
prior <- set_prior()

# note: running the model will take a few hours -> day
run_model(mdata, prior, reps = 3)

modlist <- load_and_list_runs()

fig1(mdata)
fig2(modlist, mdata)
fig3(modlist, mdata)

table1(mdata)
table2(modlist)
table3(modlist)
table4(modlist)
table5(modlist, mdata)

