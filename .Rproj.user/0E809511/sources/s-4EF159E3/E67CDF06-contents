load(url("https://github.com/iloveless1/methylXprs-Data/raw/main/COAD.rda"))

tmp <- do.call(rbind,COAD[,2])
names(COAD) <- tmp$ensembl_gene
out <- COAD[,1]

for(i in 1:length(out)){
    out[[i]][[2]] <- tmp[i,]
}

for(i in 1:length(out)){
    rownames(out[[i]][[1]]) <- out[[i]][[1]][,1]
}

for(i in 1:length(out)){
    out[[i]][[1]] <- out[[i]][[1]][,-c(1)]

}

for(i in 1:length(out)){
    out[[i]][[1]] <- data.frame(Weights = as.numeric(out[[i]][[1]]), row.names = names(out[[i]][[1]]))
}

COAD <- out
save(COAD,file = 'C:/Users/ilovele1/Desktop/EpiXprS Models/COAD.rda')
