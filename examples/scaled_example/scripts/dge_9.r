suppressMessages(library(edgeR))
options(scipen=999)

x <- read.delim('/data/./edgeR_matfile.csv', sep=',', stringsAsFactors=TRUE)

x$strain <- factor(x$strain)

x$cuminic_acid <- factor(x$cuminic_acid)

x$iptg <- factor(x$iptg)

x$timepoint <- factor(x$timepoint)

x$vanillic_acid <- factor(x$vanillic_acid)

x$xylose <- factor(x$xylose)

drops <- c('strain', 'cuminic_acid', 'iptg', 'timepoint', 'vanillic_acid', 'xylose')
counts <- x[, !(names(x) %in% drops)]
t_cts <- t(counts)
t_cts[is.na(t_cts)] <- 0

#colnames(t_cts) <- x$filename


group <- factor(c(1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 8, 7, 1, 2, 3, 4, 5, 6, 8, 7, 1, 2, 3, 4, 5, 6, 8, 7, 9, 10, 11, 12, 12, 12, 12, 9, 13, 14, 15, 15, 15, 15, 9, 13, 14, 11))

y <- DGEList(counts=t_cts, group=group)
y <- calcNormFactors(y)
design <- model.matrix(~0 + group)
keep <- rowSums(cpm(y[, c(7, 15, 23, 31, 1, 8, 16, 24)]) >1) >= 8
y <- y[keep, ,]
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, contrast=c(-1,0,0,0,0,0,1,0,0,0,0,0,0,0,0))
tab <- topTags(qlf, n=Inf)
write.table(tab, file="Bacillussubtilis168Marburg_False_False_0_False_False-vs-Bacillussubtilis168Marburg_False_True_5_False_False.txt")