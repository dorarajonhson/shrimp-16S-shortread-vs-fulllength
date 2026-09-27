
# Filter taxa at least 1% relative abundance

filter.feature.table <- function(tab, cutoff) {
    indices = c()
    for(i in 1:nrow(tab)) {
        indices = c(indices, which(tab[i,] > cutoff))
    }
    indices = unique(indices)
    return(tab[,indices])
}

filter.feature.table.mean <- function(tab, cutoff) {
    indices = which(colMeans(tab) > cutoff)
    return(tab[,indices])
}


