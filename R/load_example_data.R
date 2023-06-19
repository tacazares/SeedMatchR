#' Load example DESeq2 data into the environment
#'
#' @param example.type Name of the example to load. Options: sirna, mirna
#' @return Loads either the Schlegel 2022 RNAseq data or miRDB into the environment.
#' @export
#'
#' @examplesIf interactive()
#' load_example_data()
load_example_data <- function(example.type){
  data.path = tempdir()

  message(paste0("Example data directory being created at: ", data.path))

  sirna.path = paste0(data.path, "/Schlegel_2022.Rdata")
  mirna.path = paste0(data.path, "/mirdb.Rdata")

  if (example.type == "sirna") {
         ex.data = load(sirna.path)

         return(mget(ex.data))

    } else if (example.type == "mirna") {
         ex.data = load(mirna.path)

         return(get(ex.data))

    } else {
      stop("Invalid option. Options: sirna, mirna")
    }
}
