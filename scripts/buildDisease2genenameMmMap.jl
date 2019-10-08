using StatsBase, DelimitedFiles
include("parsePubmedAbstractFunctions.jl")

#######################################################
# Parse the data mapping PMID -> gene to allow for faster access. It should take ~5 mins to create this map
@time pmid2gene, uniquePMIDsGenes, allNCBIIDs = buildPMID2geneMap();
println("Built index for PMID -> NCBIID with " * string(length(keys(pmid2gene))) * " keys.")

# Create a map to take us from NCBI ID to Mm gene name
@time NCBIID_hs2genename_mm = buildNCBIID2genenameMmMap();
println("Built index for NCBIID -> genename with " * string(length(keys(NCBIID_hs2genename_mm))) * " keys.")

@time pmid2meshID, meshID2pmid = buildMeshID2PMIDMap();

#Build an index mapping disease names to PMID
@time disease2meshID, meshID2disease = buildMeshID2DiseaseMap();
#save these files
writedlm()

@time disease2meshID, meshID2disease, pmid2meshID, meshID2pmid = buildMeshID2DiseaseMap();


pmid2meshID, meshID2pmid, uniquePMIDsDisease, uniqueDiseaseNames, allMESHIDs, disease2meshID, meshID2disease = buildPMID2diseaseMap();

#First build a map from MESHIDs to gene names since so many of the disease names map to the same terms
uniqueMESHIDs = keys(meshID2pmid);
meshID2genename_mm = Dict{String, Array{String, 1}}(i => String[] for i in uniqueMESHIDs);
for i in 1:length(uniqueMESHIDs)
    if haskey(meshID2pmid, uniqueMESHIDs[i])
        pmids = meshID2pmid[uniqueMESHIDs[i]]
        ncbiids = Int64[];
        for j in 1:length(pmids)
            if haskey(pmid2gene, pmids[j])
                append!(ncbiids, pmid2gene[pmids[j]])
            end
        end
        for j in 1:length(ncbiids)
            if haskey(NCBIID_hs2genename_mm, ncbiids[j])
                push!(meshID2genename_mm[uniqueMESHIDs[i]], NCBIID_hs2genename_mm[ncbiids[j]]);
            end
        end
    end
end
meshID2genename_mm_count = Dict{String, Array{Int64, 1}}(i => Int64[] for i in uniqueMESHIDs);
for i in 1:length(uniqueMESHIDs)
    if length(meshID2genename_mm[uniqueMESHIDs[i]])==0
        delete!(meshID2genename_mm, uniqueMESHIDs[i]);
        delete!(meshID2genename_mm_count, uniqueMESHIDs[i]);
    else
        genenamesUnique_mm, counts = rle(sort(meshID2genename_mm[uniqueMESHIDs[i]]));
        meshID2genename_mm[uniqueMESHIDs[i]] = genenamesUnique_mm;
        meshID2genename_mm_count[uniqueMESHIDs[i]] = counts;
    end
end
saveMapFile("~/meshID2genename_mm.tsv", meshID2genename_mm, meshID2genename_mm_count);
    
#Filter to retain only high confidence interactions
#Now we can build a map from disease to Mm gene name
disease2genename_mm = Dict{String, Array{String, 1}}(i => String[] for i in uniqueDiseaseNames);
disease2genename_mm_count = Dict{String, Array{Int64, 1}}(i => Int64[] for i in uniqueDiseaseNames);
for i in 1:length(uniqueDiseaseNames)
    meshID = disease2meshID[uniqueDiseaseNames[i]];
    if haskey(meshID2genename_mm, meshID)
        disease2genename_mm[uniqueDiseaseNames[i]] = meshID2genename_mm[meshID]
        disease2genename_mm_count[uniqueDiseaseNames[i]] = meshID2genename_mm_count[meshID]
    else
        delete!(disease2genename_mm, uniqueDiseaseNames[i]);
        delete!(disease2genename_mm_count, uniqueDiseaseNames[i]);        
    end
end
saveMapFile("~/disease2genename_mm.tsv", disease2genename_mm, disease2genename_mm_count);
disease2genename_mm_filter
#Invert the mapping to get diseases associated with each gene list