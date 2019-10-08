
function readPhrases(fileName::String="../data/all_dictionary.pmid")
    @time pubmedMap = readdlm(fileName, '\t');
    println("Read " * string(length(pubmedMap)) * " phrases from PubMed.")
    phrases = String[]; # parse the keys
    for i in 1:length(pubmedMap)
        push!(phrases, split(pubmedMap[i], '|')[1]);
    end
    return phrases, pubmedMap;
end

function buildPMID2geneMap(fileName::String="../data/pmid2gene_unique.tsv")
    @time tmp = readdlm(fileName, '\t', quotes=false);
    println("Read " * string(size(tmp, 1)) * " PMID -> gene relations.")
    allPMIDs = tmp[2:end,1];
    order = sortperm(allPMIDs);
    uniquePMIDs = unique(allPMIDs[order]);
    allNCBIIDs = tmp[2:end,2][order]
    pmid2gene = Dict{Int64, Array{Int64, 1}}(i => zeros(0) for i in uniquePMIDs)
    ncbiInd = 1;
    for i in 1:length(uniquePMIDs)
        push!(pmid2gene[uniquePMIDs[i]], allNCBIIDs[ncbiInd])
        while (ncbiInd<length(allNCBIIDs)) && (allPMIDs[order[ncbiInd]]==allPMIDs[order[ncbiInd+1]])
            ncbiInd = ncbiInd + 1;
            push!(pmid2gene[uniquePMIDs[i]], allNCBIIDs[ncbiInd])
        end
        ncbiInd = ncbiInd + 1;
        if ncbiInd>=length(allNCBIIDs)
            break
        end
    end
    return pmid2gene, uniquePMIDs, allNCBIIDs;
end

function buildMeshID2DiseaseMap(fileName::String="../data/pmid2disease.tsv", stopwordsFileName::String="../data/stopwords.txt")
    stopwords = vec(convert(Array{String, 2}, readdlm(stopwordsFileName)));
    
    # This file is large and reading can take >20 mins
    @time tmp = readdlm(fileName, '\t', quotes=false);
    println("Read " * string(size(tmp, 1)) * " PMID -> disease relations.")
    
    # Build a map from disease name to MESHID
    allMESHIDs = tmp[2:end,2];
    uniqueMESHIDs = unique(allMESHIDs);
    allDiseaseNames = tmp[2:end,3];
    
    tmp2 = convert(Array{String, 1}, map(string, allDiseaseNames));
    diseaseNames = String[];
    meshIDsKeepInds = Int64[];
    for i in 1:length(tmp2)
        tokens = convert(Array{String, 1}, split(string(tmp2[i]), "|"));
        for j in 1:length(tokens)
            tokens[j] = replace(tokens[j], r",|;" => s"")
            
            # Remove stop words and sort the tokens in alphabetical order to build the index
            tokens2 = map(lowercasefirst, split(tokens[j], ' '));
            tokens3 = sort(tokens2[findall(indexin(tokens2, stopwords).==nothing)]);
            if length(tokens3)>0
                t = try
                    join(tokens3, ' ');
                catch
                    println(repr(tmp[i]))
                    "";
                end
                if ((length(repr(t))>0) .& (repr(t)!="nothing"))
                    push!(diseaseNames, lowercase(replace(repr(t), "\""=>"")));
                    push!(meshIDsKeepInds, i);
                end
            end
        end
    end
    println("Found " * repr(length(diseaseNames)) * " different names and " * repr(length(uniqueMESHIDs)) * " MESHIDs")

    meshID2disease = Dict{String, Array{String, 1}}(i => String[] for i in uniqueMESHIDs);
    meshID2diseaseCounts = Dict{String, Array{Int64, 1}}(i => Int64[] for i in uniqueMESHIDs);
    for i in 1:length(uniqueMESHIDs)
        diseases = sort(diseaseNames[findall(allMESHIDs[meshIDsKeepInds].==uniqueMESHIDs[i])]);
        diseaseNamesUnique, counts = rle(diseases);
            
        # sort in ascending order
        order = sortperm(counts);
        meshID2disease[uniqueMESHIDs[i]] = diseaseNamesUnique[order];
        meshID2diseaseCounts[uniqueMESHIDs[i]] = counts[order];
    end

        
    return meshID2disease, meshID2diseaseCounts, diseaseNames;
end

function buildMeshID2PMIDMap(fileName::String="../data/pmid2disease.tsv")
    @time tmp = readdlm(fileName, '\t', quotes=false);
    println("Read " * string(size(tmp, 1)) * " PMID -> disease relations.")
    allPMIDs = tmp[2:end,1];
    order = sortperm(allPMIDs);
    uniquePMIDs = unique(allPMIDs[order]);
    allMESHIDs = tmp[2:end,2][order];
    uniqueMESHIDs = unique(allMESHIDs);
        
    # Now build maps from MESHID to PMID and vice versa
    pmid2meshID = Dict{Int64, Array{String, 1}}(i => zeros(0) for i in uniquePMIDs);
    meshID2pmid = Dict{String, Array{Int64, 1}}(i => zeros(0) for i in uniqueMESHIDs);
    index = 1;
    for i in 1:length(uniquePMIDs)
            
        # find the MESHIDs that are associated with this PMID
        push!(pmid2meshID[uniquePMIDs[i]], allMESHIDs[index])
        push!(meshID2pmid[allMESHIDs[index]], uniquePMIDs[i])
        while (index<length(allMESHIDs)) && (allPMIDs[order[index]]==allPMIDs[order[index+1]])
            index = index + 1;
            push!(pmid2meshID[uniquePMIDs[i]], allMESHIDs[index])
            push!(meshID2pmid[allMESHIDs[index]], uniquePMIDs[i])
        end
        index = index + 1;
        if index>=length(allMESHIDs)
            break
        end
    end
    return pmid2meshID, meshID2pmid;
end

function buildPMID2diseaseMap()
    @time tmp = readdlm("../data/pmid2disease.tsv", '\t', quotes=false);
    println("Read " * string(size(tmp, 1)) * " PMID -> disease relations.")
    allPMIDs = tmp[2:end,1];
    order = sortperm(allPMIDs);
    uniquePMIDs = unique(allPMIDs[order]);
    allMESHIDs = tmp[2:end,2][order]
    allDiseaseNames = tmp[2:end,3][order]
        
    # Also build a map for disease name to MeSHID
    tmp = convert(Array{String, 1}, map(string, unique(allDiseaseNames)));
    uniqueDiseaseNames = String[];
    for i in 1:length(tmp)
        tokens = convert(Array{String, 1}, split(string(tmp[i]), "|"));
        for j in 1:length(tokens)
            push!(uniqueDiseaseNames, tokens[j]);
        end
    end
    uniqueDiseaseNames = unique(uniqueDiseaseNames);
    uniqueMESHIDs = unique(allMESHIDs);
    pmid2meshID = Dict{Int64, Array{String, 1}}(i => zeros(0) for i in uniquePMIDs)
    meshID2pmid = Dict{String, Array{Int64, 1}}(i => zeros(0) for i in uniqueMESHIDs)
    index = 1;
    for i in 1:length(uniquePMIDs)

        push!(pmid2meshID[uniquePMIDs[i]], allMESHIDs[index])
        push!(meshID2pmid[allMESHIDs[index]], uniquePMIDs[i])
        while (index<length(allMESHIDs)) && (allPMIDs[order[index]]==allPMIDs[order[index+1]])
            index = index + 1;
            push!(pmid2meshID[uniquePMIDs[i]], allMESHIDs[index])
            push!(meshID2pmid[allMESHIDs[index]], uniquePMIDs[i])
        end
        index = index + 1;
        if index>=length(allMESHIDs)
            break
        end
    end

    return pmid2meshID, meshID2pmid, uniquePMIDs, uniqueDiseaseNames, allMESHIDs, disease2meshID, meshID2disease;
end

function buildPMID2mutationMap()
    @time tmp = readdlm("../data/pmid2mutation_edit.tsv", '\t', quotes=false);
    println("Read " * string(size(tmp, 1)) * " PMID -> disease relations.")
    allPMIDs = tmp[2:end,1];
    order = sortperm(allPMIDs);
    uniquePMIDs = unique(allPMIDs[order]);
    allVars = tmp[2:end,2][order]
    pmid2var = Dict{Int64, Array{String, 1}}(i => zeros(0) for i in uniquePMIDs)
    ind = 1;
    for i in 1:length(uniquePMIDs)
        push!(pmid2var[uniquePMIDs[i]], allVars[ind])
        while (ind<length(allVars)) && (allPMIDs[order[ind]]==allPMIDs[order[ind+1]])
            ind = ind + 1;
            push!(pmid2var[uniquePMIDs[i]], allVars[ind])
        end
        ind = ind + 1;
        if ind>=length(allVars)
            break
        end
    end
    return pmid2var, uniquePMIDs, allVars;
end


function buildNCBIID2genenameMmMap()
    genename2NCBIID_mm_array = readdlm("../data/genename2NCBIID_mm.tsv", '\t', header=true)[1];
    NCBIID2genename_mm = Dict(zip(genename2NCBIID_mm_array[:,2], genename2NCBIID_mm_array[:,1]));

    ensembl2NCBIID_hs = readdlm("../data/ensemblID2NCBIID_hs.tsv", '\t', header=true)[1];
    ensembl2NCBIID_hs = ensembl2NCBIID_hs[findall(ensembl2NCBIID_hs[:,2].!=""),:];

    ensembl2genename_mm = readdlm("../data/ensemblID2genename_mm.tsv", '\t', header=true)[1]
    ensembl2genename_mm = ensembl2genename_mm[findall(ensembl2genename_mm[:,2].!=""),:];
    ensembl2genename_mm = Dict(zip(ensembl2genename_mm[:,1], ensembl2genename_mm[:,2]));
    ensembl_hs2mm = readdlm("../data/ensembl_hs2mm.tsv", '\t', header=true)[1]
    ensembl_hs2mm = ensembl_hs2mm[findall(ensembl_hs2mm[:,2].!=""),:];
    ensembl_hs2mm = Dict(zip(ensembl_hs2mm[:,1], ensembl_hs2mm[:,2]))
        
    # now build a dictionary to map from human NCBI to mouse gene name
    NCBIID_hs2genename_mm = Dict{Int64, String}(i => "" for i in ensembl2NCBIID_hs[:,2]);
    for i in 1:size(ensembl2NCBIID_hs, 1)
        foundKey = false;
        if haskey(ensembl_hs2mm, ensembl2NCBIID_hs[i,1]) # find the mouse ensembl ID
            ensembl_mm = ensembl_hs2mm[ensembl2NCBIID_hs[i,1]]
            if haskey(ensembl2genename_mm, ensembl_mm)
                NCBIID_hs2genename_mm[ensembl2NCBIID_hs[i,2]] = ensembl2genename_mm[ensembl_mm]
                foundKey = true;
            end
        end
        if ~foundKey
            delete!(NCBIID_hs2genename_mm, ensembl2NCBIID_hs[i,2])
        end
    end
    return NCBIID_hs2genename_mm;
end

function buildNCBIID2genenameHsMap()

    ensembl2NCBIID_hs = readdlm("../data/ensemblID2NCBIID_hs.tsv", '\t', header=true)[1];
    ensembl2NCBIID_hs = ensembl2NCBIID_hs[findall(ensembl2NCBIID_hs[:,2].!=""),:];
    NCBIID2ensembl_hs = Dict(zip(ensembl2NCBIID_hs[:,2], ensembl2NCBIID_hs[:,1]));
    ncbiids = collect(keys(NCBIID2ensembl_hs));

    ensembl2genename_hs = readdlm("../data/ensemblID2genename_hs.tsv", '\t', header=true)[1]
    ensembl2genename_hs = ensembl2genename_hs[findall(ensembl2genename_hs[:,2].!=""),:];
    ensembl2genename_hs = Dict(zip(ensembl2genename_hs[:,1], ensembl2genename_hs[:,2]));

    # now build a dictionary to map from NCBI IDs to human gene name
    NCBIID2genename_hs = Dict{Int64, String}(i => "" for i in ncbiids);
    for i in 1:length(ncbiids)
        foundKey = false;
        if haskey(NCBIID2ensembl_hs, ncbiids[i])
            if haskey(ensembl2genename_hs, NCBIID2ensembl_hs[ncbiids[i]])
                NCBIID2genename_hs[ncbiids[i]] = ensembl2genename_hs[NCBIID2ensembl_hs[ncbiids[i]]];
                foundKey = true;
            end

        end
        if ~foundKey
            delete!(NCBIID2genename_hs, ncbiids[i])
        end
    end
    return NCBIID2genename_hs;
end

function buildPhrase2GeneMap(saveFileStart)
    phrase2genes = Dict{String, Array{Int64, 1}}(i => zeros(0) for i in phrases);
    phrase2counts = Dict{String, Array{Int64, 1}}(i => zeros(0) for i in phrases);
    phrase2genenames_mm = Dict{String, Array{String, 1}}(i => zeros(0) for i in phrases);
    phrase2genenamescounts_mm = Dict{String, Array{Int64, 1}}(i => zeros(0) for i in phrases);
    fileHandle = open("~/phrase2NCBIID.tsv", "w")
    fileHandle_mm = open("~/phrase2genename_mm.tsv", "w")
    for i in 1:length(pubmedMap)
        tokens = split(pubmedMap[i], '|');
        pmids = map(parse, tokens[2:end]);
        geneIDs = Int64[];
        for j in 1:length(pmids)
            if haskey(pmid2gene, pmids[j])
                append!(geneIDs, pmid2gene[pmids[j]])
            end
        end
        if length(geneIDs)>0
            geneIDsUnique, counts = rle(sort(geneIDs));
            phrase2genes[tokens[1]] = geneIDsUnique;
            phrase2counts[tokens[1]] = counts;
            #save information to file
            write(fileHandle, tokens[1] * "\t")
            for j in 1:(length(geneIDsUnique)-1)
                write(fileHandle, string(geneIDsUnique[j]) * ",");
            end
            write(fileHandle, string(geneIDs[end]) * "\t");
            for j in 1:(length(geneIDsUnique)-1)
                write(fileHandle, string(counts[j]) * ",");
            end
            write(fileHandle, string(counts[end]) * "\n");
            #map the NCBI gene IDs to mm gene names
            genenames_mm = String[]
            for j in 1:length(geneIDs)
                if haskey(NCBIID_hs2genename_mm, geneIDs[j])
                    push!(genenames_mm, NCBIID_hs2genename_mm[geneIDs[j]])
                elseif haskey(NCBIID2genename_mm, geneIDs[j])
                    push!(genenames_mm, NCBIID2genename_mm[geneIDs[j]])
                end
            end
            if length(genenames_mm)>0
                genenamesUnique_mm, counts = rle(sort(genenames_mm));
                phrase2genenames_mm[tokens[1]] = genenamesUnique_mm;
                phrase2genenamescounts_mm[tokens[1]] = counts;
                #save information to file
                write(fileHandle_mm, tokens[1] * "\t")
                for j in 1:(length(genenamesUnique_mm)-1)
                    write(fileHandle_mm, string(genenamesUnique_mm[j]) * ",");
                end
                write(fileHandle_mm, string(genenamesUnique_mm[end]) * "\t");
                for j in 1:(length(genenamesUnique_mm)-1)
                    write(fileHandle_mm, string(counts[j]) * ",");
                end
                write(fileHandle_mm, string(counts[end]) * "\n");
            else
                delete!(phrase2genenames_mm, tokens[1]);
                delete!(phrase2genenamescounts_mm, tokens[1]);            
            end
        else
            delete!(phrase2genes, tokens[1]);
            delete!(phrase2counts, tokens[1]);
            delete!(phrase2genenames_mm, tokens[1]);
            delete!(phrase2genenamescounts_mm, tokens[1]);
        end
    end
    close(fileHandle)
    close(fileHandle_mm)
    println("Built index for phrase -> gene with " * string(length(keys(phrase2genes))) * " keys.")
end

function buildPhrase2GeneMapHs(ncbiid2genename_hs)
    phrase2ncbiid, phrase2ncbiid_count = readMapFile("../data/phrase2NCBIID.tsv");
    phrases = collect(keys(phrase2ncbiid));
    fileHandle_hs = open("~/phrase2genename_hs.tsv", "w");
    for i in 1:length(phrases)
        ncbiids = phrase2ncbiid[phrases[i]];
        genenames_hs = String[];
        cs = phrase2ncbiid_count[phrases[i]];
        counts_hs = Int64[];
        for j in 1:length(ncbiids)
            if haskey(ncbiid2genename_hs, parse(Int64, ncbiids[j]))
                push!(genenames_hs, ncbiid2genename_hs[parse(Int64, ncbiids[j])]);
                push!(counts_hs, cs[j]);
            end
        end
        if length(genenames_hs)>0
            genenamesUnique_hs, counts2 = rle(sort(genenames_hs));
            write(fileHandle_hs, phrases[i] * "\t")
            for j in 1:(length(genenamesUnique_hs)-1)
                write(fileHandle_hs, string(genenamesUnique_hs[j]) * ",");
            end
            write(fileHandle_hs, string(genenamesUnique_hs[end]) * "\t");
            for j in 1:(length(genenamesUnique_hs)-1)
                inds = findall(genenames_hs.==genenamesUnique_hs[j]);
                write(fileHandle_hs, string(sum(counts_hs[inds])) * ",");
            end
            inds = findall(genenames_hs.==genenamesUnique_hs[end]);
            write(fileHandle_hs, string(sum(counts_hs[inds])) * "\n");
        end
    end
    close(fileHandle_hs)
end


function filterMap(wordmap::Dict{String, Array{String, 1}}, countmap::Dict{String, Array{Int64, 1}}; minCount::Int64=2, maxPhrase::Int64=49, minPhrase::Int64=2, minFraction::Float64=0.01)
    # Use this method to filter out mappings with few hits
    ks = collect(keys(wordmap));
    wordmap_filter = Dict{String, Array{String, 1}}(i => String[] for i in ks);
    countmap_filter = Dict{String, Array{Int64, 1}}(i => Int64[] for i in ks);
    for i in 1:length(ks)
        indsKeep = findall(countmap[ks[i]].>minCount);
        indsKeep = indsKeep[findall(countmap[ks[i]][indsKeep].>minFraction*sum(countmap[ks[i]][indsKeep]))];
        if length(indsKeep)>minPhrase # keep the maxGenes best matches
            nKeep = maxPhrase;
            if length(indsKeep)<=maxPhrase
                nKeep = length(indsKeep) - 1;
            end
            inds = indsKeep[sortperm(countmap[ks[i]][indsKeep])[end-nKeep:end]];
            wordmap_filter[ks[i]] = wordmap[ks[i]][inds];
            countmap_filter[ks[i]] = countmap[ks[i]][inds];
        else
            delete!(wordmap_filter, ks[i]);
            delete!(countmap_filter, ks[i]);
        end
    end
    return wordmap_filter, countmap_filter;
end

function saveMapFile(fileName::String, wordmap::Dict{Int, Array{Int, 1}})
    fileHandle = open(fileName, "w");
    ks = collect(keys(wordmap));
    for i in 1:length(ks)
        write(fileHandle, repr(ks[i]) * "\t")
        for j in 1:(length(wordmap[ks[i]])-1)
            write(fileHandle, repr(wordmap[ks[i]][j]) * ",");
        end
        write(fileHandle, repr(wordmap[ks[i]][end]) * "\n");
    end
    close(fileHandle);
end

function saveMapFile(fileName::String, wordmap::Dict{String, Array{String, 1}})
    fileHandle = open(fileName, "w");
    ks = collect(keys(wordmap));
    for i in 1:length(ks)
        write(fileHandle, ks[i] * "\t")
        for j in 1:(length(wordmap[ks[i]])-1)
            write(fileHandle, wordmap[ks[i]][j] * ",");
        end
        write(fileHandle, wordmap[ks[i]][end] * "\n");
    end
    close(fileHandle);
end

function saveMapFile(fileName::String, wordmap::Dict{String, String})
    fileHandle = open(fileName, "w");
    ks = collect(keys(wordmap));
    for i in 1:length(ks)
        write(fileHandle, ks[i] * "\t" * wordmap[ks[i]] * "\n")
    end
    close(fileHandle);
end

function readMap(fileName::String="W2Vdata/pmid2NCBIID.tsv")
    @time tmp = readdlm(fileName, '\t', header=false);
    pmids = tmp[:,1];    
    pmid2gene = Dict{Int, Array{Int, 1}}(i => Int[] for i in pmids);
    for j in 1:length(pmids)
        if occursin(",", string(tmp[j,2]))
            append!(pmid2gene[phrases[j]], convert(Array{String, 1}, split(tmp[j,2], ',')));
            tmp2 = convert(Array{String, 1}, split(string(tmp[j,3]), ','));
            for i in 1:length(tmp2)
                push!(pmid2gene[pmids[j]], parse(Int64, tmp2[i]));
            end
        else
            push!(pmid2gene[pmids[j]], parse(Int64, string(tmp[j,2])));
            push!(pmid2gene[pmids[j]], tmp[j,3]);
        end
    end
    return pmid2gene;
end
                        
                        
function saveMapFile(fileName::String, wordmap::Dict{String, Array{String, 1}}, countmap::Dict{String, Array{Int64, 1}}; sortkeys::Bool=false, sep::String=",")
    fileHandle = open(fileName, "w");
    ks = collect(keys(wordmap));
    for i in 1:length(ks)
        if sortkeys
            tokens = sort(split(ks[i], [' ', ':']));
            for j in 1:(length(tokens)-1)
                write(fileHandle, tokens[j] * " ")
            end
            write(fileHandle, tokens[end] * "\t")
        else
            write(fileHandle, ks[i] * "\t")
        end
        if length(wordmap[ks[i]])>0
            for j in 1:(length(wordmap[ks[i]])-1)
                write(fileHandle, wordmap[ks[i]][j] * sep);
            end
            write(fileHandle, wordmap[ks[i]][end] * "\t");
            for j in 1:(length(countmap[ks[i]])-1)
                write(fileHandle, string(countmap[ks[i]][j]) * sep);
            end
            write(fileHandle, string(countmap[ks[i]][end]) * "\n");
        end
    end
    close(fileHandle);
end
    
    
function readMapFile(fileName::String="../data/phrase2genename_mm.tsv")
    # Use this function to parse the files mapping phrases or gene names across. This takes ~30 minutes to run, suggesting that we also want to save in a binary format to speed up loading of the index
    @time tmp = readdlm(fileName, '\t', header=false);
    phrases = map(repr, tmp[:,1]); 
    for j in 1:length(phrases)
        phrases[j] = replace(phrases[j], r"\"" => "");                            
    end
    phrase2genename_mm = Dict{String, Array{String, 1}}(i => String[] for i in phrases);
    phrase2genenamescounts_mm = Dict{String, Array{Int64, 1}}(i => Int64[] for i in phrases);
    for j in 1:length(phrases)
        if occursin(",", string(tmp[j,2]))
            append!(phrase2genename_mm[phrases[j]], convert(Array{String, 1}, split(tmp[j,2], ',')));
            tmp2 = convert(Array{String, 1}, split(string(tmp[j,3]), ','));
            for i in 1:length(tmp2)
                push!(phrase2genenamescounts_mm[phrases[j]], parse(Int64, tmp2[i]));
            end
        else
            push!(phrase2genename_mm[phrases[j]], string(tmp[j,2]));
            push!(phrase2genenamescounts_mm[phrases[j]], tmp[j,3]);
        end
    end
    return phrase2genename_mm, phrase2genenamescounts_mm;
end

function mergeMapFile(fileStart::String="../datagenename_hs2words_", genenamesfile::String="../data/ensemblID2genename_hs.tsv")
    # Merge several different maps
    dir = Base.Filesystem.dirname(fileStart)
    namestart = Base.Filesystem.basename(fileStart)
    files = Base.Filesystem.readdir(dir)
    targets = nothing;
    targetcounts = nothing; 
    genenames = readdlm(genenamesfile, header=true)[1][:,2];
    genenames = unique(genenames[findall(!isempty, genenames)]);
    for f in files
        if startswith(f, namestart)
            println("Opening " * dir * "/" * f)
            @time tmp = readdlm(dir * "/" * f, '\t', header=false);
            ks = tmp[:,1];
            if targets==nothing
                targets = Dict{String, Array{String, 1}}(i => String[] for i in genenames);
                targetcounts = Dict{String, Array{Int64, 1}}(i => Int64[] for i in genenames);
            end
            for j in 1:length(ks)
                if haskey(targets, ks[j])
                    if occursin(",", string(tmp[j,2]))
                        ws = convert(Array{String, 1}, split(tmp[j,2], ','))
                        tmp2 = convert(Array{String, 1}, split(string(tmp[j,3]), ','));
                        for i in 1:length(ws)
                            ind = findall(ws[i].==targets[ks[j]])
                            if length(ind)==1
                                targetcounts[ks[j]][ind[1]] += parse(Int64, tmp2[i])
                            elseif isempty(ind)
                                push!(targets[ks[j]], ws[i]);
                                push!(targetcounts[ks[j]], parse(Int64, tmp2[i]));
                            end
                        end
                    else
                        ind = findall(string(tmp[j,2]).==targets[ks[j]])
                        if length(ind)==1
                            targetcounts[ks[j]][ind[1]] += tmp[j,3]
                        elseif isempty(ind)
                            push!(targets[ks[j]], string(tmp[j,2]));
                            push!(targetcounts[ks[j]], tmp[j,3]);
                        end
                    end
                else
                    println("Could not find gene " * ks[j])
                end
            end
        end
    end
    for i in genenames
        if length(targets[i])==0
            delete!(targets, i)
            delete!(targetcounts, i)
        end
    end
    return targets, targetcounts;
end

    
function invertMapFile(keyList::Array{String, 1}; fileName::String="../data/phrase2genename_mm.tsv")
    # Use this function to parse and invert the files mapping phrases or gene names across
    @time tmp = readdlm(fileName, '\t', header=false);

    geneName2phrase = Dict{String, Array{String, 1}}(i => String[] for i in keyList);
    geneName2phraseCount = Dict{String, Array{Int64, 1}}(i => Int64[] for i in keyList);
    phrases = tmp[:,1];
    for j in 1:length(phrases)
        if occursin(",", repr(tmp[j,2]))
            geneNames = convert(Array{String, 1}, split(tmp[j,2], ','));
            tmp2 = convert(Array{String, 1}, split(string(tmp[j,3]), ','));

            for i in 1:length(geneNames)
                if haskey(geneName2phrase, geneNames[i])
                    push!(geneName2phrase[geneNames[i]], phrases[j])
                    push!(geneName2phraseCount[geneNames[i]], parse(Int64, tmp2[i]))
                end
            end
        else
            if haskey(geneName2phrase, tmp[j,2])
                push!(geneName2phrase[tmp[j,2]], phrases[j])
                push!(geneName2phraseCount[tmp[j,2]], tmp[j,3])
            end
        end
    end
    for i in keyList
        if haskey(geneName2phrase, i)
            if length(geneName2phrase[i])==0
                delete!(geneName2phrase, i);
                delete!(geneName2phraseCount, i);
            end
        end
    end
    return geneName2phrase, geneName2phraseCount;
 end
    

function parseImplicitVectors(l::Int64)
    # Get the list of words that were indexed and re-order them to allow for fast access
    fullModel = readdlm("/lustre/scratch117/cellgen/team218/MH/ImplicitVectors/full_model.csv", '\t', skipstart=3);
    inds = convert(Array{Int64, 1}, fullModel[:,1]);
    fullModelID2Word = Dict(zip(fullModel[:,1], fullModel[:,2])); fullModel[sortperm(inds),2];
                                                        
    # read the files that allow us to map PMIDs to NCBIIDs and then NCBIIDs to genenames for Mm or Hs
    NCBIID2genename_hs = buildNCBIID2genenameHsMap();
    genenames_hs = collect(values(NCBIID2genename_hs));
    NCBIID2genename_mm = buildNCBIID2genenameMmMap();
    genenames_mm = collect(values(NCBIID2genename_mm));
    @time pmid2gene, uniquePMIDs, allNCBIIDs = buildPMID2geneMap();

    pmid2wordsFileStart="/lustre/scratch117/cellgen/team218/MH/ImplicitVectors/implicit_vectors_";

        genename_mm2words = Dict{String, Array{String, 1}}(i => String[] for i in genenames_mm);
        genename_hs2words = Dict{String, Array{String, 1}}(i => String[] for i in genenames_hs);
        ncbiid2words = Dict{Int64, Array{String, 1}}(i => String[] for i in allNCBIIDs);

        fileName = pmid2wordsFileStart * string(l) * ".csv"
        println("Opening new file, " * fileName);
        @time run(`gunzip $fileName.gz`)
        fileHandle = open("../data/implicit_vectors_" * string(l) * ".csv");
        line = readline(fileHandle);
        ncbiids = Int64[];
        wordInds = String[];
        while length(line)>0
            try
            tokens = split(line, '\t');
            pmid = parse(Int64, tokens[1]);
            if haskey(pmid2gene, pmid)
                ncbiids = pmid2gene[pmid];
                wordInds = split(tokens[2], " ");
                for j in 1:length(wordInds)
                    wInd = parse(Int64, wordInds[j]);
                    if haskey(fullModelID2Word, wInd)
                        for k in 1:length(ncbiids)
                            if haskey(NCBIID2genename_mm, ncbiids[k])
                                push!(genename_mm2words[NCBIID2genename_mm[ncbiids[k]]], string(fullModelID2Word[wInd]));
                            end
                            if haskey(NCBIID2genename_hs, ncbiids[k])
                                push!(genename_hs2words[NCBIID2genename_hs[ncbiids[k]]], string(fullModelID2Word[wInd]));
                            end
                            push!(ncbiid2words[ncbiids[k]], string(fullModelID2Word[wInd]));
                        end
                    end
                end
            end
            catch
                println("Error at: " * line)
                break;
            end
            line = readline(fileHandle);
        end
        close(fileHandle);
        @time run(`gzip $fileName`)
                                                            
        # Now we can compress by counting how many times each gene was associated with each word
        genename_hs2words_counts = Dict{String, Array{Int64, 1}}(i => Int64[] for i in genenames_hs);
        ncbiid2words_counts = Dict{Int64, Array{Int64, 1}}(i => Int64[] for i in allNCBIIDs);
        for i in 1:length(genenames_hs)
            if haskey(genename_hs2words, genenames_hs[i])
                if length(genename_hs2words[genenames_hs[i]])>0
                    words_hs, counts_hs = rle(sort(genename_hs2words[genenames_hs[i]]));
                    genename_hs2words[genenames_hs[i]] = words_hs;
                    genename_hs2words_counts[genenames_hs[i]] = counts_hs;
                else
                    delete!(genename_hs2words, genenames_hs[i]);
                    delete!(genename_hs2words_counts, genenames_hs[i]);
                end
            end
        end
        saveMapFile("~/genename_hs2words_" * string(l) * ".tsv", genename_hs2words, genename_hs2words_counts);
        genename_mm2words_counts = Dict{String, Array{Int64, 1}}(i => Int64[] for i in genenames_mm);
        for i in 1:length(genenames_mm)
            if haskey(genename_mm2words, genenames_mm[i])
                if length(genename_mm2words[genenames_mm[i]])>0
                    words_mm, counts_mm = rle(sort(genename_mm2words[genenames_mm[i]]));
                    genename_mm2words[genenames_mm[i]] = words_mm;
                    genename_mm2words_counts[genenames_mm[i]] = counts_mm;
                else
                    delete!(genename_mm2words, genenames_mm[i]);
                    delete!(genename_mm2words_counts, genenames_mm[i]);
                end
            end
        end
        saveMapFile("~/genename_mm2words_" * string(l) * ".tsv", genename_mm2words, genename_mm2words_counts);
        ncbiid2words_counts = Dict{Int64, Array{Int64, 1}}(i => Int64[] for i in allNCBIIDs);
        for i in 1:length(allNCBIIDs)
            if haskey(ncbiid2words, allNCBIIDs[i])
                if length(ncbiid2words[allNCBIIDs[i]])>0
                    ws, cs = rle(sort(ncbiid2words[allNCBIIDs[i]]));
                    ncbiid2words[allNCBIIDs[i]] = ws;
                    ncbiid2words_counts[allNCBIIDs[i]] = cs;
                else
                    delete!(ncbiid2words, allNCBIIDs[i]);
                    delete!(ncbiid2words_counts, allNCBIIDs[i]);
                end
            end
        end
        saveMapFile("~/NCBIID2words_" * string(l) * ".tsv", ncbiid2words, ncbiid2words_counts);

    return genename_hs2words, genename_hs2words_counts, genename_mm2words, genename_mm2words_counts, ncbiid2words, ncbiid2words_counts;
end