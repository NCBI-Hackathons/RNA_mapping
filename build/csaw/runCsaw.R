--------------------------------------------
#Run csaw on user supplied bam files
--------------------------------------------

#quietly load the two packages
suppressMessages(require(csaw))
suppressMessages(require(edgeR))
suppressMessages(require(TxDb.Hsapiens.UCSC.hg38.knownGene))
suppressMessages(require(org.Hs.eg.db))
suppressMessages(require(locfit))
suppressMessages(require(statmod))

#gather user supplied bam files as object "args"
args <- commandArgs(trailingOnly = TRUE)

#create counts using csaw, save as data
data = windowCounts(args,width=1000,bin=TRUE, filter=10)
abundances <- aveLogCPM(asDGEList(data))
keep <- abundances > aveLogCPM(10, lib.size=mean(data$totals))
summary(keep)

#this section makes sure that the TxDb names match the names in the bam files
#if you're using the whole genome you would have to change the names of the 
#sequences prior to mapping

#for whole genome uncomment line 27 and comment out line 28
##row.good.names <- rowRanges(data)
row.good.names <- renameSeqlevels(rowRanges(data), "chr20")
hexons <- exons(TxDb.Hsapiens.UCSC.hg38.knownGene)
suppressWarnings(keep <- overlapsAny(row.good.names, hexons))
summary(keep)
data <- data[keep,]
ranges <- rowRanges(data)
out.counts <- assay(data)
mcols(ranges) = out.counts
colnames(mcols(ranges))=args
summary(out.counts)
r.df <- as(ranges, "data.frame")

#set the design
normfacs <- rep(1, each = length(args))
y <- asDGEList(data, norm.factors=normfacs)
x = data.frame(aligner=c(args[1], rep("other", each= (length(args)-1))))
design <- model.matrix(~aligner, x)
y <- estimateDisp(y,design = design)
fit <- glmQLFit(y, design, robust=TRUE)
results <- glmQLFTest(fit, contrast=c(0, 1))

#convert pvalues
pvalue.table <- results$table
p.fdr.adjust <- p.adjust(pvalue.table$PValue, method = "fdr")
r.df$Pvalue <- p.fdr.adjust

#sort pvalues
p.sorted <- order(r.df$Pvalue)
results.p.sorted <- r.df[p.sorted,]

#output table
write.table(results.p.sorted, file = "csaw.bam.results.csv", sep = ",")
