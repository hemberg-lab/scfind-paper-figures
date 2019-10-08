using StatsBase
include("parsePubmedAbstractFunctions.jl")

#######################################################
# Parse the data mapping PMID -> gene to allow for faster access. It should take ~5 mins to create this map
@time pmid2gene, uniquePMIDsGenes, allNCBIIDs = buildPMID2geneMap();
println("Built index for PMID -> NCBIID with " * string(length(keys(pmid2gene))) * " keys.")

# Create a map to take us from NCBI ID to Mm gene name
@time NCBIID_hs2genename_hs = buildNCBIID2genenameHsMap();
println("Built index for NCBIID -> genename with " * string(length(keys(NCBIID_hs2genename_hs))) * " keys.")

# Build an index mapping variants to PMID
pmid2variant, uniquePMIDsVariants, allVariants = buildPMID2mutationMap("../data/pmid2mutation.tsv")

# Now we can build a map from variant to Mm gene name
variant2genename_hs = Dict{String, Array{String, 1}}(i => String[] for i in allVariants);
for i in 1:length(uniquePMIDsVariants)
    if haskey(pmid2gene, uniquePMIDsVariants[i])
        ncbiids = pmid2gene[uniquePMIDsVariants[i]];
        vars = pmid2variant[uniquePMIDsVariants[i]];
        for j in 1:length(vars)
            for k in 1:length(ncbiids)
                if haskey(NCBIID_hs2genename_hs, ncbiids[k])
                    push!(variant2genename_hs[vars[j]], NCBIID_hs2genename_hs[ncbiids[k]]);
                end
            end
        end
    end
end
allVariantsFound = collect(keys(variant2genename_hs));
variant2genename_hs_count = Dict{String, Array{Int64, 1}}(i => Int64[] for i in allVariantsFound);
for i in 1:length(allVariantsFound)
    if length(variant2genename_hs[allVariantsFound[i]])==0
        delete!(variant2genename_hs, allVariantsFound[i]);
        delete!(variant2genename_hs_count, allVariantsFound[i]);
    else
        genenamesUnique_hs, counts = rle(sort(variant2genename_hs[allVariantsFound[i]]));
        variant2genename_hs[allVariantsFound[i]] = genenamesUnique_hs;
        variant2genename_hs_count[allVariantsFound[i]] = counts;
    end
end
saveMapFile("~/variant2genename_hs.tsv", variant2genename_hs, variant2genename_hs_count);

