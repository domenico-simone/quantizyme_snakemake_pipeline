# This file is sourced by an R script, i.e. all variables must be set according to R syntax rules.

# If a filename is enterd here: load the clustalw alignment for the whole fasta from this file, instead of calculating it.
# (e.g. when analysis is to be resumed) 
# if a filename is entered here, it must reside in the folder where the script is started
# if you do not want to load the clustalw alignment, outcomment this command or enter an empty string ("") 
# phylip.all = "GT48_align1.phy"	# for GT48_enzyme.fasta  1233 contigs
# phylip.all = "GT48small_align1.phy"			# for the small test.fasta with 24 sequences
phylip.all = ""


# number of boots in clustalw 
nr.boot = 100  		# clustalw TODO (not in source code parametrized) 


# Each subtree will be subdivided in a number of groups.
# A sequence is randomly picked from each group.
# The random component can be critical.
# It is therefore possible to repeat the random picking, in order to see if repeated picking delievers different results
# NOER This variable has to be defined. (default is 1) 
# nr.trials.random.picking = 1
nr.trials.random.picking = 10   


# Thinning: Approximate size of a subgroup in a subtree, in percent. Default = 4
subgroup.percent = 30 


# Number of CPUs used by clustalO
clustalo.nr.cpu = 20

# Number of CPU's used when creating the profile HMM (hmmbuild command):
hmmbuild.nr.cpu = 20

# Number of CPU's used when searching with the profile HMM (nhmmer command):
nhmmer.nr.cpu = 20


# E-value for the nhmmer search 
nhmmer.evalue = 1   

