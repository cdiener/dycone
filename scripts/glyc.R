#  glyc.R
#  
#  Copyright 2015 Christian Diener <ch.diener@gmail.com>
#  
#  MIT license. See LICENSE for more information.

# This is an example studying a small sample of glycolysis. We extract 
# all the required information from a Biomodels SBML file and execute
# a differential k-cone analysis for a version of the system where we
# disturbed the reaction responsible for ethanol production.

library(sybilSBML)
load_all("..")

# Extracting the stochiometric matrix and steady state concentrations

sbml = readSBMLmod("../stuff/BIOMD0000000042.xml")
rownames(sbml@S) = sbml@met_id
reacts = as.reactions(sbml@S, reversible=sbml@react_rev, r_names=sbml@react_id)
S = get_stochiometry(reacts) # to get the irreversible stochiometric matrix
species_data = get_sbml_species("../stuff/BIOMD0000000042.xml") 
concs = species_data$initialAmount
names(concs) = species_data$id

# Calculating the k-cone and its eigenpathway

kc = get_polytope_basis(S, get_ma_terms(S, concs))
ep = eigenpathways(kc)[,1]

# As parameters we choose a faster version of the eigenpathway
k = 1e3*ep

# Now we generate several new steady states coming from a disturbed ethanol 
# production 

scales = 10^c(-2,-4)
d_concs = concs
d_eps = ep 
d_basis = list(kc)
for( i in 1:length(scales) ) {
	d_k = k
	d_k[30] = d_k[30]*scales[i]
	tc = timecourse(concs, c(0,1e5), d_k, S)
	d_concs = rbind(d_concs, tc[2,-1])
	basis = get_polytope_basis(S, get_ma_terms(S, tc[2,-1]))
	d_basis[[i+1]] = basis
	d_eps = rbind(d_eps, eigenpathways(basis)[,1])
}

# Lets plot a projection of the respective k-cones
plot_red(d_basis)
readline("Press [Enter] to continue...")

# And the differnce between the last and initial basis
plot( d(d_basis[[3]], kc) )
readline("Press [Enter] to continue...")

# A similar picture can be obtained using only the eigenpathways
barplot( d_eps, col=TRANSCOL(3), names=1:ncol(d_eps), xlab="reactions", beside=T )

# And let's see the hypothesis
h = hyp(kc, d_basis[[3]], reacts)
print(h)
