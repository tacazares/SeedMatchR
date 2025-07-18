utils::globalVariables(c("gene_id", "utr3_len", "tx_id"))

#' Load species-specific annotations needed to run [SeedMatchR::SeedMatchR()]
#'
#' @description `load_annotations` uses [AnnotationHub::AnnotationHub()] to load
#' species-specific annotations. This function currently works for human, rat,
#' and mouse.
#'
#' @details
#' This function is designed to perform all annotation loading, parsing, and
#' sequence extractions for the reference of interest.
#'
#' `load_annotations` uses the `ensembldb` package for handling annotations.
#' This allows for the use of [AnnotationFilter::AnnotationFilter]s for
#' easily selecting transcripts of interest.
#'
#' Most of the arguments to `load_annotations` are geared towards customizing
#' the output features and sequence extracted from the [ensembldb::EnsDb]
#' object. Filtering parameters are inherited from
#' [SeedMatchR::build_annotation_filter()]`.
#'
#' @inheritParams load_species_anno_db
#' @inheritParams build_annotation_filter
#' @inheritParams get_3UTR_granges
#' @param feature.type The transcript feature type to extract. Options: "3UTR", "5UTR", "exons", "cds"
#' @param min.feature.width The minimum length for features in nucleotides.
#' @param return_gene_name If true, the gene_id will be used to name the sequences.
#'
#' @return The following are returned in the object:
#' * `$gtf`: A [GenomicRanges::GRanges] object containing the transcript information
#' * `$txdb`: A [ensembldb::EnsDb] object of transcript information
#' * `$dna`: The [rtracklayer::TwoBitFile] DNA sequence
#' * `$features`: A [GenomicRanges::GRangesList] object representing the features of interest
#' * `$seqs`: A [Biostrings::DNAStringSet] for the features in `$features`
#' * `$filter`: A [AnnotationFilter::AnnotationFilterList] used for filtering the [ensembldb::EnsDb] object
#' * `$tx.names`: A list of transcript names in `$seqs`
#'
#' @export
#'
#' @examplesIf interactive()
#' anno.db = load_annotations("hg38")
load_annotations <- function(reference.name = c("hg38", "hg38-old", "mm39", "mm10", "rnor6", "rnor7"),
                             feature.type = c("3UTR", "5UTR", "exons", "cds"),
                             standard.chroms = TRUE,
                             remove.na.rows = TRUE,
                             protein.coding = TRUE,
                             transcript.support = NULL,
                             canonical = TRUE,
                             gene_id = NULL,
                             tx_id = NULL,
                             symbol = NULL,
                             entrez_id = NULL,
                             longest.utr = FALSE,
                             reduce.features = FALSE,
                             add.filter = NULL,
                             min.feature.width = 8,
                             return_gene_name = TRUE
                             ){
  feature.type <- match.arg(feature.type)
  reference.name <- match.arg(reference.name)

  # Filter GTF
  afl <- build_annotation_filter(standard.chroms = standard.chroms,
                                  remove.na.rows = remove.na.rows,
                                  protein.coding = protein.coding,
                                  transcript.support = transcript.support,
                                  canonical = canonical,
                                  gene_id = gene_id,
                                  tx_id = tx_id,
                                  symbol = symbol,
                                  entrez_id = entrez_id,
                                  add.filter = add.filter)

  # Load species specific annotations from AnnotationHub: $gtf $dna $txdb
  annodb = load_species_anno_db(reference.name = reference.name, filter = afl)

  message("Extracting ", feature.type, " from ensembldb object.\n",
          if(longest.utr){"Keeping the longest UTR per gene.\n"},
          if(reduce.features){"Keep the longest UTR per gene.\n"})

  # Extract feature based on feature.type argument
  features <- switch(
    feature.type,
    "3UTR" = get_3UTR_granges(annodb$txdb,
                               filter = afl,
                               reduce.features = reduce.features,
                               longest.utr = longest.utr),
    "5UTR" = ensembldb::fiveUTRsByTranscript(annodb$txdb, filter = afl),
    "exons" = ensembldb::exonsBy(annodb$txdb, by = "tx", filter = afl),
    "cds" = ensembldb::cdsBy(annodb$txdb, by = "tx", filter = afl)
  )

  message("Extracting sequences for each feature.\n",
          "Keeping sequences that are >= ", min.feature.width, "\n"
  )

  # Extract the feature sequences
  features.seqs = GenomicFeatures::extractTranscriptSeqs(annodb$dna, features)

  # Remove sequences for which there is not 3' UTR or the UTR is not > 0
  features.seqs = features.seqs[BiocGenerics::width(features.seqs) >= min.feature.width]

  tx.gtf = annodb$gtf[annodb$gtf$tx_id %in% names(features.seqs)]

  if(return_gene_name == TRUE){
    names(features.seqs) <- tx.gtf$gene_id[match(names(features.seqs), tx.gtf$tx_id)]
  }

  # Create a list of results to return
  results.list = list(gtf = tx.gtf,
                      txdb = annodb$txdb,
                      features = features,
                      seqs = features.seqs,
                      dna = annodb$dna,
                      filter = afl,
                      tx.names = features.seqs$names)

  return(results.list)
}


#' Load species-specific annotations from `AnnotationHub`
#'
#' @description Use [AnnotationHub::AnnotationHub()] to load species-specific
#' [rtracklayer::TwoBitFile] DNA sequence and [ensembldb::EnsDb]
#' object. This function currently works for human, rat, and mouse.
#'
#' The specific versions and builds used for each reference are further defined
#' in `articles(SeedMatchR_Annotation_Databases)`.
#'
#' @param reference.name Reference build name. Options: hg38, rnor6, rnor7,
#' mm39, mm10
#' @param filter An [AnnotationFilter::AnnotationFilterList] used to refine annotations.

#' @return The function will return:
#' * `$gtf`: A [GenomicRanges::GRanges] object containing the GTF information
#' * `$dna`: A [rtracklayer::TwoBitFile] with reference DNA sequence
#' * `$txdb`: A [ensembldb::EnsDb] object with transcriptome information
#'
#' @export
#'
#' @examplesIf interactive()
#' anno.db = SeedMatchR::load_species_anno_db("hg38")
load_species_anno_db <- function(reference.name = c("hg38", "hg38-old", "mm39", "mm10", "rnor6", "rnor7", "rnor7.113"), filter = NULL){
  reference.name <- match.arg(reference.name)

  message("Loading annotations from AnnotationHub for ", reference.name)

  # Load annotation hub
  ah <- AnnotationHub::AnnotationHub()

  ah.df = data.frame(reference = c("hg38", "mm39", "mm10", "rnor6", "rnor7", "rnor7.113"),
                     dna = c("AH106283", "AH106440", "AH88475", "AH93578", "AH106786", "AH106786"),
                     txdb = c("AH109606", "AH109655", "AH89211", "AH95846", "AH109732", "AH119437"))

  dna <- ah[[ah.df$dna[ah.df$reference == reference.name]]]

  txdb <- ah[[ah.df$txdb[ah.df$reference == reference.name]]]

  gtf <- ensembldb::transcripts(txdb, filter = filter)

  return(list(gtf = gtf,
              dna = dna,
              txdb = txdb))
}

#' Reduce 3' UTRs to the gene level across all transcripts
#'
#' @param gr.list Input [GenomicRanges::GRangesList] file with features grouped by transcript
#' @param txdb An [ensembldb::EnsDb] object
#' @param filter An [AnnotationFilter::AnnotationFilterList] for selecting transcripts
#'
#' @return A gtf of unique feature intervals per gene
#' @export
#'
#' @examplesIf interactive()
#' anno.db = SeedMatchR::load_species_anno_db("hg38")
#'
#' features = SeedMatchR::get_features_gtf(anno.db$txdb)
#'
#' gtf = SeedMatchR::.reduce_features(features, anno.db$txdb)
#'
.reduce_features <- function(gr.list, txdb, filter){
  # Create a df to map tx_id, tx_name, and gene_id
  tx2gene <- S4Vectors::mcols(GenomicFeatures::transcripts(txdb,
                                                           columns=c("tx_id",
                                                                     "tx_name",
                                                                     "gene_id"),
                                                           filter = filter
  ))

  tx2gene$gene_id <- as.character(tx2gene$gene_id)

  # Unlist the GTF so that they are not longer grouped by transcript
  all_features <- BiocGenerics::unlist(gr.list)

  ## Add "tx_name" and "gene_id" metadata columns to feature specific gtf
  m <- match(names(all_features), tx2gene$tx_name)

  mcol.df <- S4Vectors::mcols(all_features)

  S4Vectors::mcols(all_features) <- cbind(mcol.df, tx2gene[m, -1L, drop=FALSE])

  names(all_features) <- NULL

  # Split the GTF into a list by gene_id
  features_by_genes <- GenomicRanges::split(all_features, all_features$gene_id)

  # Reduce the overlapping genomic ranges to only have unique intervals
  gtf = GenomicRanges::reduce(features_by_genes)

  return(gtf)
}

#' Build an [AnnotationFilter::AnnotationFilter] for selecting transcripts
#'
#' @description
#' This function helps create an [AnnotationFilter::AnnotationFilterList] that will combine
#' multiple standard filters based on the flags.
#'
#' @param standard.chroms Keep only standard chromosomes. Boolean.
#' @param remove.na.rows Remove transcripts with NA in the ID column. Boolean.
#' @param protein.coding Keep only protein coding genes. Boolean.
#' @param transcript.support Filter by transcript support level. Integer 1-4
#' @param canonical Keep only the ENSEMBL canonical transcript. Boolean.
#' @param gene_id Keep specific ENSEMBL ID(s). Input is a vector
#' @param tx_id Keep specific tx_id(s). Input is a vector
#' @param entrez_id Keep specific ENTREZ ID(s). Input is a vector
#' @param symbol Keep specific gene symbol(s). Input is a vector
#' @param add.filter Optional [AnnotationFilter::AnnotationFilterList] to include
#'
#' @return An [AnnotationFilter::AnnotationFilterList] object
#' @export
#'
#' @examples
#' # Default: standard chromosomes, removes na tx_ids, and keep protein coding
#' afl = SeedMatchR::build_annotation_filter()
#'
#' # Keep only transcripts with support level 1
#' afl = SeedMatchR::build_annotation_filter(transcript.support = 1)
#'
build_annotation_filter <- function(standard.chroms = TRUE,
                                     remove.na.rows = TRUE,
                                     protein.coding = TRUE,
                                     transcript.support = NULL,
                                     canonical = FALSE,
                                     gene_id = NULL,
                                     tx_id = NULL,
                                     symbol = NULL,
                                     entrez_id = NULL,
                                     add.filter = NULL){
  message("Build AnnotationFilter for transcript features based on the following parameters: \n",
          "Keep only standard chroms: ", standard.chroms ,"\n",
          "Remove rows with NA in transcript ID: ", remove.na.rows, "\n",
          "Keep only protein coding genes and transcripts: ", protein.coding, "\n",
          "Filtering for transcripts with support level: ", if(is.null(transcript.support) > 0){FALSE} else {transcript.support}, "\n",
          "Keep only the ENSEMBL canonical transcript: ", canonical, "\n",
          "Filtering for specific genes: ", if(length(gene_id) > 0){TRUE} else {FALSE}, "\n",
          "Filtering for specific transcripts: ", if(length(tx_id) > 0){TRUE} else {FALSE}, "\n",
          "Filtering for specific gene symbols: ", if(length(symbol) > 0){TRUE} else {FALSE}, "\n",
          "Filtering for specific entrez id: ", if(length(entrez_id) > 0){TRUE} else {FALSE}, "\n"
  )

  afl <- AnnotationFilter::AnnotationFilterList()

  # Keep only standard chromosomes
  if (standard.chroms) {
    # filter the gtf for only standard chromosomes
    std_filter <-
      AnnotationFilter::SeqNameFilter(c(paste0(1:22), "X", "Y"))

    afl = AnnotationFilter::AnnotationFilterList(afl, std_filter)
  }

  # Remove transcripts without a name
  if (remove.na.rows) {
    txna_filter = AnnotationFilter::TxIdFilter("", "!=")

    afl = AnnotationFilter::AnnotationFilterList(afl, txna_filter)
  }

  # Include only protein coding genes
  if (protein.coding) {
    txbio_filter = AnnotationFilter::TxBiotypeFilter("protein_coding")
    genebio_filter = AnnotationFilter::GeneBiotypeFilter("protein_coding")

    afl = AnnotationFilter::AnnotationFilterList(afl, genebio_filter, txbio_filter)
  }

  # Include only ENSEMBL canonical transcripts
  if (canonical) {
    txcanon_filter = ensembldb::TxIsCanonicalFilter(1)

    afl = AnnotationFilter::AnnotationFilterList(afl, txcanon_filter)
  }

  # Filter by transcript support level
  if (!is.null(transcript.support)) {
    stopifnot(transcript.support <= 5)

    txsupport_filter = ensembldb::TxSupportLevelFilter(transcript.support)

    afl = AnnotationFilter::AnnotationFilterList(afl, txsupport_filter)
  }

  # Include only target genes in input list
  if (!is.null(gene_id)) {
    geneid_filter = AnnotationFilter::GeneIdFilter(gene_id)

    afl = AnnotationFilter::AnnotationFilterList(afl, geneid_filter)
  }

  # Include only target transcripts in input list
  if (!is.null(tx_id)) {
    txid_filter = AnnotationFilter::TxIdFilter(tx_id)

    afl = AnnotationFilter::AnnotationFilterList(afl, txid_filter)
  }

  # Include only target gene symbols in input list
  if (!is.null(symbol)) {
    symbol_filter = AnnotationFilter::SymbolFilter(symbol)

    afl = AnnotationFilter::AnnotationFilterList(afl, symbol_filter)
  }

  # Include only target gene symbols in input list
  if (!is.null(entrez_id)) {
    entrez_filter = AnnotationFilter::EntrezFilter(entrez_id)

    afl = AnnotationFilter::AnnotationFilterList(afl, entrez_filter)
  }

  if (!is.null(add.filter)){
    afl = AnnotationFilter::AnnotationFilterList(afl, add.filter)
  }

  return(afl)
}

#' Extract [GenomicRanges::GRangesList] of 3' UTRs
#'
#' @description
#' `get_3UTR_granges` leverages the [GenomicFeatures::threeUTRsByTranscript()]
#' function to extract 3' UTR regions from a given [ensembldb::EnsDb] object.
#'
#' @inheritParams .reduce_features
#' @param longest.utr Reduce annotations to the longest UTR per gene
#' @param reduce.features Reduce 3' UTRs to non-overlapping intervals across all
#' transcripts per gene.
#'
#' @return A [GenomicRanges::GRangesList] object representing the genomic
#' features of interest. Results are returned grouped by transcript, unless
#' `longest.utr` is TRUE, then results are grouped by gene.
#' @export
#'
#' @examplesIf interactive()
#' annodb <- SeedMatchR::load_species_anno_db("hg38")
#'
#' features <- SeedMatchR::get_3UTR_granges(annodb$txdb)
#'
get_3UTR_granges <- function(txdb,
                             filter = NULL,
                             longest.utr = FALSE,
                             reduce.features = FALSE){
  utr.grlist = ensembldb::threeUTRsByTranscript(txdb, filter = filter)

  # If longest UTR == T, create a set of features that
  # reports the longest tx per gene.
  if (longest.utr){
    utr.grlist = .get_longest_utr(utr.grlist, txdb,  filter = filter)
  }

  # If reduce.features == T, create a set of features
  # that takes the union of all 3' UTRs per gene.
  if (reduce.features){
    utr.grlist = suppressWarnings(.reduce_features(utr.grlist, txdb,  filter = filter))
  }

  return(utr.grlist)

}

#' Get the longest transcript per gene
#'
#' @description
#' `.get_longest_utr` will extract the longest transcript per gene in a [GenomicRanges::GRangesList]
#' using an [ensembldb::EnsDb] object. The [GenomicFeatures::transcriptLengths()] function
#' is used to calculate the 3' UTR lengths for each transcript in the [ensembldb::EnsDb] object.
#'
#' @param utr.grlist A [GenomicRanges::GRangesList] of UTR exons grouped by `tx_id`
#' @param txdb An [ensembldb::EnsDb] object
#' @param filter Optional [AnnotationFilter::AnnotationFilterList] object
#'
#' @return A [GenomicRanges::GRangesList] filtered containing the longest 3' UTR
#'  per gene in the [ensembldb::EnsDb] object that pass filters in the optional
#'  [AnnotationFilter::AnnotationFilterList]
#' @export
#'
#' @examplesIf interactive()
#' annodb <- SeedMatchR::load_species_anno_db("hg38")
#'
#' features = GenomicFeatures::threeUTRsByTranscript(annodb$txdb)
#'
#' longest.tx = SeedMatchR:::.get_longest_tx(features, annodb$txdb)
.get_longest_utr <- function(utr.grlist, txdb, filter = NULL){
  # Extract transcript lengths for each gene
  tx_lens <- GenomicFeatures::transcriptLengths(txdb, with.utr3_len=TRUE, filter = filter)

  # Find the longest transcript per gene
  longest_utr_per_gene = tx_lens %>% dplyr::group_by(gene_id) %>% dplyr::summarize(longest_UTR = max(utr3_len),
                                                                            tx_id = tx_id[which.max(utr3_len)])

  # Subset annotations for the longest transcript
  utr.grlist = utr.grlist[longest_utr_per_gene$tx_id]

  return(utr.grlist)
}

