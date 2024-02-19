args <- commandArgs(trailingOnly = TRUE)
exprsMat_path <- args[1]
cellTypes_path <- args[2]
topN <- 50
pSig <- 0.001

exprsMat <- read.csv(exprsMat_path, row.names=1)
cellTypes <- readLines(cellTypes_path)

# Select genes by bimodal index
cellTypes <- droplevels(as.factor(cellTypes))
tt <- list()
for (i in seq_len(nlevels(cellTypes))) {
    tmp_celltype <- (ifelse(cellTypes == levels(cellTypes)[i], 1, 0))

    pi <- table(tmp_celltype)/length(tmp_celltype)

    agg_mean <- do.call(cbind, lapply(c(0,1), function(i){
        Matrix::rowMeans(exprsMat[, tmp_celltype == i, drop = FALSE])
    }))

    agg_sd2 <- do.call(cbind, lapply(c(0,1), function(i){
        apply(exprsMat[, tmp_celltype == i, drop = FALSE], 1, stats::var)
    }))

    bi <- abs(agg_mean[,2] - agg_mean[,1])/sqrt(pi[1]*agg_sd2[,1] +
                                                    pi[2]*agg_sd2[,2])

    bi <- unlist(bi)
    names(bi) <- rownames(exprsMat)
    bi <- bi[order(bi, decreasing = TRUE)]
    tt[[i]] <- bi
}

tt <- lapply(tt, function(x)
    x)
res <- Reduce(union, lapply(tt, function(t)
    names(t)[seq_len(topN)]))
write(res, file.path(dirname(exprsMat_path), "BI_features.txt"))
