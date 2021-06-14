# import annotation.bed,
# here:
# X1 = chr X2 = strand, num, X3 = start, X4 = end, X6 = gene name


# this script merges the transcripts belonging to the same gene, on the same chromosome and same strand

res <- annotation %>% group_by(X1, X2, X6, X7) %>% summarise(n = n(), min = min(X3), max = max(X4))

out <- res %>% select(X1, X2, min, max, X6, X7)

out_geneid <- out %>% mutate(gene_id = paste(X1, X2, "_", X6, sep = ""))

write.table(out_geneid, file="out.bed", quote = F, row.names = F, col.names = F, sep = "\t")
