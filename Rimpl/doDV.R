args <- commandArgs(trailingOnly = TRUE)
exprsMat_path <- args[1]
cellTypes_path <- args[2]
topN <- 50
pSig <- 0.001

exprsMat <- read.csv(exprsMat_path, row.names=1)
cellTypes <- readLines(cellTypes_path)

## doDV
cellTypes <- droplevels(as.factor(cellTypes))
tt <- list()
for (i in seq_len(nlevels(cellTypes))) {
    tmp_celltype <- (ifelse(cellTypes == levels(cellTypes)[i], 1, 0))

    meanPct <- do.call(cbind, lapply(c(0,1), function(i){
        Matrix::rowSums(exprsMat[, tmp_celltype == i, drop = FALSE] > 0) / sum(tmp_celltype == i)
    }))
    posNeg <- (meanPct[,2] - meanPct[,1]) > 0.05
    exprsMat_filt <- exprsMat[posNeg,]
    tt[[i]] <- apply(exprsMat_filt, 1, function(x) {
        df <- data.frame(gene = x, cellTypes = as.factor(tmp_celltype))
        stats::bartlett.test(gene~cellTypes, df)$p.value
    })

    tt[[i]] <- stats::p.adjust(tt[[i]], method = "BH")
}

tt <- lapply(tt, function(x)
    sort(x))
res <- Reduce(union, lapply(tt, function(t)
    names(t)[seq_len(min(topN, sum(t < pSig)))]))
write(res, file.path(dirname(exprsMat_path), "DV_features.txt"))