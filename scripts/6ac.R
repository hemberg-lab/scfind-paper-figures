generateUCSCBedFileForEnhancers <- function(object, query.gene = "Nr4a1", tissue = "Heart", window = c(1e5, 0), target.window = c("chr15", 101097277, 101105225 ), annotation) {
    #Create a .bed file containing tracks showing the "number of cells"
    t.atac <- object
    max.dist <- as.integer(max(window[1:2]))
    min.dist <- as.integer(min(window[1:2]))
    
    peak.locs <- t.atac@index$genes()
    peak.pos.list <- strsplit((peak.locs), '_')
    tmp <- unlist((peak.pos.list))
    peak.chr <- tmp[(0:(length(tmp)/3-1))*3+1]
    peak.start <- as.integer(tmp[(0:(length(tmp)/3-1))*3+2])
    peak.end <- as.integer(tmp[(0:(length(tmp)/3-1))*3+3])
    
    peaks <- c()
    
    genes.inds <- which(annotation$gene.names==query.gene)
    
    for( i in genes.inds  )
    {
        gene.inds <- i
        peaks <- c(peaks, peak.locs[which(annotation$chr[gene.inds]==peak.chr & abs(annotation$start[gene.inds]-peak.start)<max.dist & abs(annotation$start[gene.inds]-peak.end)<max.dist & abs(annotation$start[gene.inds]-peak.start)>min.dist & abs(annotation$start[gene.inds]-peak.end)>min.dist)])
}
    peaks <- unique(peaks)
    
    #Keep track of the number of peaks for each cell type in this list
    cell.type.peak.count <- list()
    for (i in 1:length(peaks)) {
        cell.type.peak.count[[peaks[i]]] <- list()
        res <- hyperQueryCellTypes(t.atac, peaks[i], tissue)
        for (j in 1:nrow(res)) {
            
            # Percentage of cells
            cell.type.peak.count[[peaks[i]]][[as.character(res$cell_type[j])]] <- res$cell_hits[j]/res$total_cells[j]
        }
    }
    
    peak.pos.list <- strsplit((peaks), '_')
    tmp <- unlist((peak.pos.list))
    peak.chr <- tmp[(0:(length(tmp)/3-1))*3+1]
    peak.start <- as.integer(tmp[(0:(length(tmp)/3-1))*3+2])
    peak.end <- as.integer(tmp[(0:(length(tmp)/3-1))*3+3])
    
    # open bed file for each cell type for writing
    cell.types <- unique(sub(paste0(".*\\.", tissue), tissue, names(unlist(cell.type.peak.count)))) 
    save.file.name <- list()
    save.file.handle <- list()
    
    # writing content for each bed file
    dir.create(paste0("~/", query.gene, "_", as.character(max.dist/1e3),"k"))
    
    for (i in 1:length(cell.types)) {
        save.file.name[[cell.types[i]]] <- paste0("~/scfind/EnhancerAnalysis/", query.gene, "_", as.character(max.dist/1e3),"k/", as.character(max.dist),"_",query.gene, "_" , gsub("/", "-", cell.types[i]), "_peaks.bed")
        save.file.handle[[cell.types[i]]] <- file(save.file.name[[i]], "w")
        writeLines(paste0( "track  type=broadPeak visibility=1 db=mm9 name=\"",cell.types[i] ,"\" ", "description=\"",  
                           cell.types[i],
                           "\"  color=255,0,0"), con=save.file.handle[[cell.types[i]]])
        writeLines(paste0("browser position ", target.window[1],":",as.integer(target.window[2]) - max.dist, "-", as.integer(target.window[3]) + max.dist), con=save.file.handle[[cell.types[i]]])
    }
    
    # writing peak positions for each opened file
    for (i in 1:length(peaks)) {
        for (j in 1:length(peaks[[i]])) {
            
            ind <- which(names(cell.type.peak.count[[i]])[j]==cell.types)
            writeLines(paste(peak.chr[i], peak.start[i], peak.end[i], paste0(round(cell.type.peak.count[[i]][[j]], 4)*100,"%"), "0", ".", cell.type.peak.count[[i]][[j]], -1, -1, sep="\t"), con=save.file.handle[[ind]])
        }
    }
    
    # close writing
    for (i in 1:length(cell.types)) {
        close(save.file.handle[[i]])
    }
    
}

loadAnnotation <- function(path) {
    #ensure that start corresponds to the TSS
    annotation.mm9 <- read.table(path, header=F)
    chr.col <- 3
    strand.col <- 4
    start.col <- 5
    end.col <- 6
    genename.col <- 13
    plus.inds <- which(annotation.mm9[,strand.col]=='+')
    minus.inds <- which(annotation.mm9[,strand.col]=='-')
    return( list(chr=annotation.mm9[c(plus.inds, minus.inds),chr.col], start=c(annotation.mm9[plus.inds,start.col], annotation.mm9[minus.inds,end.col]), end=c(annotation.mm9[plus.inds,end.col], annotation.mm9[minus.inds,end.col]), gene.names=as.character(annotation.mm9[c(plus.inds, minus.inds),genename.col]) ) )
}