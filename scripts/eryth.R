#  eryth.R
#  
#  Copyright 2015 Christian Diener <ch.diener@gmail.com>
#  
#  MIT license. See LICENSE for more information.

# Script to analyze the k-cone basis for the Human Red Blood Cell Model

devtools::load_all("..")
library(ggplot2)

data(eryth)
species = get_species(eryth)
ss_concs = read.csv("../stuff/eryth_ic.csv")
ss_ids = ss_concs$id
ss_concs = ss_concs$initial
names(ss_concs) = ss_ids
no_c = !(species %in% ss_ids)
constants = species[no_c]
ss_concs = ss_concs[species[!no_c]]

S = get_stochiometry(eryth, const=constants) 

linB = get_kcone_basis(S, get_ma_terms(S, ss_concs))
polB = get_polytope_basis(S, get_ma_terms(S, ss_concs))

# plot basis vectors
plin = plot_basis(linB)
ppol = plot_basis(polB)

stab_lin = stability_analysis(linB, S, ss_concs)
stab_pol = stability_analysis(polB, S, ss_concs)
write("\nBasis stability:", file="")
print( table(stab_pol$what) )

f_pol = get_fluxes(polB, S, eryth, ss_concs)
write("\nLargest mean fluxes in steady state:", "")
print( head(f_pol,10) )

# Calculate basis dynamics for first vector
tc = timecourse( ss_concs+0.5, seq(0,10,0.1), 2e4*polB[,1],S )
pl = plot(tc) + xlab("Time")

# Save plots
ggsave("linB.pdf", plin, width=8, height=6)
ggsave("polB.pdf", ppol, width=8, height=6)
ggsave("tc1.pdf", pl, width=12, height=8)
