guide.seq = "UUAUAGAGCAAGAACACUGUUUU"

# Generate DESEQ2 results test data

get_example_data("sirna")

sirna.data = load_example_data("sirna")

anno.db = load_annotations(reference.name = "rnor6")

res <- sirna.data$Schlegel_2022_Ttr_D1_30mkg

res = filter_deseq(res, fdr.cutoff=1, fc.cutoff=0, rm.na.log2fc = T)

sub.genes = intersect(res$gene_id[1:200], names(anno.db$seqs))

anno.db$seqs[names(anno.db$seqs) %in% sub.genes]

# Generate seqs file as fasta

res = SeedMatchR(res = res,
                 gtf = anno.db$gtf,
                 seqs = anno.db$seqs,
                 sequence = guide.seq,
                 seed.name = "mer7m8")

match.df = SeedMatchR(gtf = anno.db$gtf,
                      seqs = anno.db$seqs[sub.genes],
                      sequence = guide.seq,
                      seed.name = "mer7m8")


write.table(res[res$gene_id %in% sub.genes,], file = "tests/testthat/testdata/DESEQ2_results.tsv", sep="\t", quote=F, row.names = F)

Biostrings::writeXStringSet(anno.db$seqs[sub.genes], format = "fasta",filepath = "tests/testthat/testdata/Test_seqs.fa")

gtf.df = as.data.frame(anno.db$gtf[anno.db$gtf$gene_id %in% sub.genes])

write.table(gtf.df, file = "tests/testthat/testdata/Test_gtf.tsv", sep="\t", quote=F, row.names = F)

gtf.granges = makeGRangesFromDataFrame(gtf.df, keep.extra.columns=TRUE)
