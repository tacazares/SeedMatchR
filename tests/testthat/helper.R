####### Variables for SeedMatchR tests ##########
#################################################
gtf.df = read.table(file.path("testdata", "Test_gtf.tsv"), sep="\t", header = 1)
gtf.df$tx_id = gtf.df$tx_id

gtf = GenomicRanges::makeGRangesFromDataFrame(gtf.df, keep.extra.columns=TRUE)

res = utils::read.table(file.path("testdata", "DESEQ2_results.tsv"), sep="\t", header = 1)

seqs = Biostrings::readDNAStringSet(file.path("testdata", "Test_seqs.fa"))

guide.seq = "UUAUAGAGCAAGAACACUGUUUU"

dna.seq = "TTAUAGAGCAAGAACACUGUUUU"

####### Variables for deseq_fc_ecdf tests #######
#################################################
# Gene set 1
mer7m8.list = res$gene_id[res$mer7m8 >= 1]

# Gene set 2
background.list = res$gene_id[res$mer7m8 == 0]

####### Variables for plotAlignment tests #######
#################################################
# Seed information for full sequence
seed.inf = SeedMatchR::get_seed(guide.seq, "Full", start = 1, stop = 23)

# Convert the target sequence to a character string
guide.rev.comp = as.character(seed.inf$Target.Seq)

# This version of the target string has missing characters
guide.rev.comp.dot = "...............CTCTATA."

# This version of the target string has missing characters
guide.rev.comp.space = "               CTCTATA "

# This version of the target string has missing characters
guide.rev.comp.x = "xxxxxxxxxxxxxxxCTCTATAx"

# This sequence has 3 mismatches to the guide sequence at pos 2, 9, 20
query.seq.dna = "ATAACAGTCTTCTTGCTCTGTAA"

# This version of the target string is longer than the match
long.query.seq.dna = "AATAACAGTCTTCTTGCTCTGTAA"

SampleTranscript = seqs[names(seqs) %in% c("ENSRNOG00000016275")]

seedmatchr.hit = SeedMatchR::SeedMatchR(sequence = guide.seq, seqs = SampleTranscript, res.format = "data.frame", seed.name = "Partial", start.pos = 2, stop.pos = 19, max.mismatch = 3)

seedmatchr.hit = SeedMatchR::SeedMatchR(sequence = guide.seq, seqs = SampleTranscript, res.format = "data.frame", seed.name = "mer7m8", max.mismatch = 1, match.df = seedmatchr.hit)

seedmatchr.hit.gr = SeedMatchR::SeedMatchR(sequence = guide.seq, seqs = SampleTranscript, res.format = "granges", seed.name = "mer7m8", max.mismatch = 1)

test.ranges = SeedMatchR::full_search(guide.seq, seqs, res_format = "iranges")

