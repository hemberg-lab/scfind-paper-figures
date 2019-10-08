using StatsBase
include("parsePubmedAbstractFunctions.jl")

#######################################################
# This first bit should take a few minutes to run
phrases, pubmedMap = readPhrases();
println("Read " * string(length(pubmedMap)) * " phrases from PubMed.")


#######################################################
# Parse the data mapping PMID -> gene to allow for faster access. It should take ~5 mins to create this map
pmid2gene, uniquePMIDs, allNCBIIDs = buildPMID2geneMap();
println("Built index for PMID -> gene with " * string(length(keys(pmid2gene))) * " keys.")

#######################################################
# Build a dictionary for NCBIIDs to gene names, for some reason, this one-liner only works for mm and not for hs. This should be done in less than one minute


# genename2NCBIID_hs = readdlm("/nfs/team218/MH/PubMedPhrases/genename2NCBIID_hs.tsv", '\t', header=true)[1]
NCBIID_hs2genename_mm = buildNCBIID2genenameMmMap();
NCBIID2genename_hs = buildNCBIID2genenameHsMap();
