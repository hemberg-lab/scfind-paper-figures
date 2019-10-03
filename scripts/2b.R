#!/usr/bin/env Rscript

insilico.gating <- function (index, cell.type, dataset, reduced.dim, gating = list(), ignore.na = T)
{
    if(length(index@metadata) == 0) 
    {
        message("Please make sure your index contains metadata")
        return(data.frame())
    }
    if(!any(index@datasets %in% dataset)) return(c())
    if(!any(cellTypeNames(index, dataset) %in% cell.type) )
    {
        message(toString(cellTypeNames(index, dataset)))
        return(data.frame())       
    }
    cell.type <- sub(paste0(dataset, '.'), "", cell.type)
    tissue.ind <- grep(paste0("^", dataset, "$"), index@metadata[[1]][[1]])
    map.ind <- index@metadata[[1]][[tolower(reduced.dim)]][tissue.ind]
    t.map <- index@metadata[[map.ind]]
    if(reduced.dim  == "pca") percentVar <- attr(t.map, "percentVar")
    t.map <- subset(t.map, rownames(t.map) == cell.type)
    
    if(length(gating) == 0) return(t.map)
    
    rownames(t.map) <- rep("NA", nrow(t.map))
    
    for( i in names(gating) )
    {
        gate <- findCellTypes(index, gating[[i]], dataset)[["Thymus.T cell"]]
        if(length(gate) != 0) rownames(t.map)[gate] <- as.character(i)
    }
    
    t.map <- if(ignore.na == T) subset(t.map, rownames(t.map) != "NA") 
    attr(t.map, "percentVar") <- if(reduced.dim == "pca") percentVar
    return(t.map)
}


plot.gating <- function(index, cell.type, dataset, reduced.dim, gating = list(), ignore.na = T)
{
    df <- insilico.gating(index = index, cell.type = cell.type, dataset = dataset,
                          reduced.dim = reduced.dim, gating = gating, ignore.na = ignore.na)
    percentVar <- if(reduced.dim == "pca" && !is.null(attr(df, "percentVar"))) attr(df, "percentVar")
    row.df <- rownames(df)
    df <- data.frame(df[,c(1,2)])
    
    gate <- as.character(row.df)
    colnames(df) <- c("dim1", "dim2")
    df$gate <- gate
    if(reduced.dim == "pca")
    {
        to_plot <- 1:2
        if ( is.null(percentVar) ) {
            labs <- sprintf("PC%i", to_plot)
        } else {
            labs <- sprintf("PC%i (%i%%)", to_plot, round(percentVar[to_plot] * 100))
        }
        gg <- ggplot(df, aes(x = dim1, y = dim2, colour = gate)) + geom_point(alpha = .6, size = 2) + 
            xlab(labs[1]) + ylab(labs[2]) +
            main + labs(colour = element_blank()) #+ facet_grid(.~gate) 
    }
    else
    {
        gg <- ggplot(df, aes(x = dim1, y = dim2, colour = gate)) + geom_point(alpha = .6, size = 2) + 
            xlab(paste0(toupper(reduced.dim), "1" )) + ylab(paste0(toupper(reduced.dim), "2" )) +
            main + labs(colour = element_blank()) #+ facet_grid(.~gate) 
    }
    return(gg)
}
