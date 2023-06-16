#' Load example DESeq2 data into the environment
#'
#' @param example.type Name of the example to load. Options: sirna, mirna
#' @export
#'
#' @examplesIf interactive()
#' load_example_data()
load_example_data <- function(example.type){
  data.path = tools::R_user_dir("SeedMatchR", "data")

  sirna.path = paste0(data.path, "/Schlegel_2022.Rdata")
  mirna.path = paste0(data.path, "/mirdb.Rdata")

  if (example.type == "sirna") {
         load(sirna.path, envir = .GlobalEnv)
    } else if (example.type == "mirna") {
         load(mirna.path, envir = .GlobalEnv)
    } else {
      stop("Invalid option. Options: sirna, mirna")
    }
}
