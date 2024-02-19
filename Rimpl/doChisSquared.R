args <- commandArgs(trailingOnly = TRUE)
exprsMat_path <- args[1]
cellTypes_path <- args[2]
topN <- 50
pSig <- 0.001
threshold <- 1

exprsMat <- read.csv(exprsMat_path, row.names=1)
cellTypes <- readLines(cellTypes_path)

## doChisSquared
cellTypes <- droplevels(as.factor(cellTypes))
tt <- list()
for (i in seq_len(nlevels(cellTypes))) {
    tmp_celltype <- (ifelse(cellTypes == levels(cellTypes)[i], 1, 0))

    zerosMat <- ifelse(exprsMat > threshold, 1, 0)

    tt[[i]] <- apply(zerosMat,1,  function(x){
        tab <- c()
        for (i in c(0,1)) {
            tmp <- factor(x[tmp_celltype == i], levels = c(0, 1))
            tab <- rbind(tab, table(tmp))
        }
        suppressWarnings(stats::chisq.test(tab)$p.value)
    })
    tt[[i]] <- stats::p.adjust(tt[[i]], method = "BH")
}

tt <- lapply(tt, function(x)
    sort(x))
res <- Reduce(union, lapply(tt, function(t)
    names(t)[seq_len(min(topN, sum(t < pSig)))]))
write(res, file.path(dirname(exprsMat_path), "chisq_features.txt"))