# Use this method to process the GWAS catalog by converting file names and 

function ConvertGWASCatalog(saveFileName="~/gwas_catalog_scfind.tsv", missense::Bool=false, hs::Bool=false)
    # This function is quite slow, but it does the job...
    gwas = readdlm("../data/gwas_catalog_short.tsv", '\t', header=true)[1];
    if missense
        gwas = readdlm("../data/gwas_catalog_missense.tsv", '\t', header=false);
    end
    traits = unique(gwas[:,1]);
    mm2hs = readdlm("../data/mousehuman_martexport.csv", ',', header=true)[1][:,[1,3]];
    writeFileHandle = open(saveFileName, "w");
    for i in 1:length(traits)
        rowInds = findall(gwas[:,1].==traits[i]);

        reportedGenesList = map(x->split(x, ","), map(x->string(x), gwas[rowInds,2]));
        reportedGenesArray = String[];
        ors = Float64[];
        foundGene = false;
        for j in 1:length(reportedGenesList)
            geneNamesHs = map(x->strip(x), reportedGenesList[j]); # remove white spaces
            for k in 1:length(geneNamesHs)
                geneInds = findall(mm2hs[:,2].==geneNamesHs[k]); # find the mouse gene 
                if hs
                    geneInds = collect(1:length(geneNamesHs));
                end
                if length(geneInds)>0
                    if !foundGene
                        foundGene = true;
                    end
                    for l in 1:length(geneInds)
                        geneName = geneNamesHs[l];
                        if !hs
                            geneName = convert(String, mm2hs[geneInds[l], 1]);
                        end
                        geneInd = findall(reportedGenesArray.==geneName);
                        if !missense
                            or = gwas[rowInds[j],4]; # Now add the odds ratio
                            if isempty(geneInd) # check if gene already found
                                push!(reportedGenesArray, geneName);
                                if or==""
                                    push!(ors, 1.0);
                                else
                                    push!(ors, convert(Float64, or));
                                end
                            elseif or!=""
                                if or>ors[geneInd[1]]
                                    ors[geneInd[1]] = or;
                                end
                            end
                        elseif isempty(geneInd) # check if gene already found
                            push!(reportedGenesArray, geneName);
                        end
                    end
                end
            end
        end
        if foundGene
            write(writeFileHandle, traits[i] * "\t");
            for j in 1:(length(reportedGenesArray)-1)
                write(writeFileHandle, reportedGenesArray[j] * ",");
            end
            if !missense
                write(writeFileHandle, reportedGenesArray[end] * "\t");
                for j in 1:(length(reportedGenesArray)-1)
                    write(writeFileHandle, string(ors[j]) * ",", );
                end
                write(writeFileHandle, string(ors[end]) * "\n");
            else
                write(writeFileHandle, reportedGenesArray[end] * "\n");
            end
        else
            println("Could not find any genes for " * traits[i]);
        end
    end
    close(writeFileHandle);
end

function ConvertGO(saveFileName="~/go_scfind.tsv")
    go = readdlm("../data/genenames_mm_go.tsv", '\t', header=true)[1];
    terms = unique(go[:,2]);
    writeFileHandle = open(saveFileName, "w");
    for i in 1:length(terms)
        if terms[i]!=""
            rowInds = findall(go[:,2].==terms[i]);
            reportedGenesArray = go[rowInds,1];
            if length(reportedGenesArray)>0
                write(writeFileHandle, terms[i] * "\t");
                for j in 1:(length(reportedGenesArray)-1)
                    write(writeFileHandle, reportedGenesArray[j] * ",");
                end
                write(writeFileHandle, reportedGenesArray[end] * "\n");
            end
        end
    end
    close(writeFileHandle);
end

function buildGWAS2variantMap(saveFileName="~/gwastraits2genename_hs.tsv", missense::Bool=false)
    # Use this method to parse the GWAS catalog and build a map for getting from trait to a list of variants. This can be combined with the map of variants to gene names to find associated genes based on literature rather than on proximity.
    # get the map from snps to gene names for hs and mm
    variants2genenames_hs, variants2genenames_counts_hs = readMapFile("../data/variant2genename_hs.tsv")
    traitCol = 8;
    snpCol = 22;
    pValCol = 29;
    contextCol = 25;
    maxpval = 7.3; # 5e-8;
    gwas = readdlm("../data/gwas_catalog_v1.0-associations_e93_r2018-08-04.tsv", '\t', header=true)[1][:,[traitCol, snpCol, contextCol, pValCol]];
    traits = unique(sort(gwas[:,1]));
    pvals = gwas[:,4]; 
    fileHandle = open(saveFileName, "w");

                        
    for i in 1:length(traits)
        inds = findall(gwas[:,1].==traits[i]);
        inds = try 
            inds[findall(pvals[inds].>maxpval)];
        catch e
            println(traits[i])
            println(string(inds))
            println(string(e))
            Int64[]
        end
        genenames_hs = String[];
        genenames_hs_counts = Int64[];
        for j in 1:length(inds)
            key = gwas[inds[j], 2];
            if haskey(variants2genenames_hs, key)
                genes = variants2genenames_hs[key];
                for k in 1:length(genes)
                    ind = findall(genenames_hs.==genes[k]);
                    if length(ind)==1
                        genenames_hs_counts[ind[1]] += variants2genenames_counts_hs[key][k];
                    else
                        push!(genenames_hs, genes[k]);
                        push!(genenames_hs_counts, variants2genenames_counts_hs[key][k]);
                    end
                end
            end
        end
        if length(genenames_hs)>0
            write(fileHandle, traits[i] * "\t")
            for j in 1:(length(genenames_hs)-1)
                write(fileHandle, genenames_hs[j] * ",");
            end
            write(fileHandle, genenames_hs[end] * "\t");
            for j in 1:(length(genenames_hs_counts)-1)
                write(fileHandle, string(genenames_hs_counts[j]) * ",");
            end
            write(fileHandle, string(genenames_hs_counts[end]) * "\n");
        end
    end
    close(fileHandle);
end