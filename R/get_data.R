#' Download example DESeq2 data from GEO
#' @description This function will download data that can be used  for SeedMatchR.
#' Choosing 'sirna' will download 3 DESeq2 results files from GSE184929. Choosing
#' 'mirna' will download the miRDB database as a tsv.
#'
#' @param example.type Name of the example to load. Options: sirna, mirna
#' @return None?
#' @export
#'
#' @examplesIf interactive()
#' get_example_data()
get_example_data <- function(example.type = c("sirna", "mirna")){
  options(timeout=180)

  example.type <- match.arg(example.type)
  data.path = tools::R_user_dir("SeedMatchR", "data")

  message("Example data directory being created at: ", data.path)

  dir.create(data.path, recursive = TRUE)

  if (example.type == "sirna"){
    Schlegel_2022_Ttr_D1_30mkg = download_parse_file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184929/suppl/GSE184929_Fig3-TTR-2_groupD1_30mkg.vs.Mock.DESeq2.txt.gz",
               file.path(data.path, "GSE184929_Fig3-TTR-2_groupD1_30mkg.vs.Mock.DESeq2.txt.gz"))

    # Download D4 data 30 mg/kg
    Schlegel_2022_Ttr_D4_30mkg = download_parse_file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184929/suppl/GSE184929_Fig3-TTR-2_groupD4_30mkg.vs.Mock.DESeq2.txt.gz",
                file.path(data.path, "GSE184929_Fig3-TTR-2_groupD4_30mkg.vs.Mock.DESeq2.txt.gz"))

    # Download D1 data 10 mg/kg
    Schlegel_2022_Ttr_D1_10mkg = download_parse_file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184929/suppl/GSE184929_Fig3-TTR-1_groupD1_10mkg.vs.Mock.DESeq2.txt.gz",
                file.path(data.path, "GSE184929_Fig3-TTR-1_groupD1_10mkg.vs.Mock.DESeq2.txt.gz"))

    save(Schlegel_2022_Ttr_D1_30mkg, Schlegel_2022_Ttr_D4_30mkg, Schlegel_2022_Ttr_D1_10mkg, file = file.path(data.path, "Schlegel_2022.Rdata"))

  } else if (example.type == "mirna") {
    mirdb.path = file.path(data.path, "miRDB_v6_prediction_results.txt.gz")

    utils::download.file("https://mirdb.org/download/miRDB_v6.0_prediction_result.txt.gz",
                  mirdb.path)

    # Load table
    mirdb = utils::read.table(mirdb.path,
                              sep="\t",
                              header=F,
                              stringsAsFactors=F)

    # Split the names into species and miRNA name
    mirdb[ , c("species", "miRNA.name")] <- stringr::str_split_fixed(mirdb[,1], "-", 2)

    # Update column names
    colnames(mirdb) = c("miRDB.ID", "target.REFSEQ.ID", "target.score", "species", "miRNA.name")

    save(mirdb, file = file.path(data.path, "mirdb.Rdata"))

  }
}


#' Download and parse DESeq2 output from GSE184929
#'
#' @param download.path File url to be downloaded
#' @param output.path Filename used for saving downloaded file
#'
#' @return DESeq2 results as a data.frame.
#' @export
#'
#' @examplesIf interactive()
#' download_parse_file()
download_parse_file <- function(download.path, output.path){
  utils::download.file(download.path,
                output.path)

  res = utils::read.csv(output.path, sep="\t")

  colnames(res) <- c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "symbol")

  return(res)
}

#' Load example DESeq2 data into the environment
#'
#' @param example.type Name of the example to load. Options: sirna, mirna
#' @return Loads either the Schlegel 2022 RNAseq data or miRDB into the environment.
#' @export
#'
#' @examplesIf interactive()
#' load_example_data()
load_example_data <- function(example.type = c("sirna", "mirna")){
  example.type <- match.arg(example.type)
  data.path = tools::R_user_dir("SeedMatchR", "data")

  sirna.path = file.path(data.path, "Schlegel_2022.Rdata")
  mirna.path = file.path(data.path, "mirdb.Rdata")

  if (example.type == "sirna") {
    ex.data = load(sirna.path)

    return(mget(ex.data))

  } else if (example.type == "mirna") {
    ex.data = load(mirna.path)

    return(get(ex.data))

  }
}
