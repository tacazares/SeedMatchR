x <- NULL
y.seed <- NULL
y.non.seed <- NULL
val <- NULL
group <- NULL

#' Function to compute overall, positive, and negative Wasserstein differences between ECDFs
#'
#' @param seed a vector of log2 FC values for the pre-specified seed matching genes.
#' @param non.seed a vector of log2FC values for the pre-specified non-seed matching genes.
#' @param n.interp number of interpolation points between (0,1) for integrating eCDF difference.
#'
#' @return A list containing:
#' list(abs.auc      = Absolute Wasserstein statistic,
#'      neg.cens.auc = Negative Wasserstein statistic,
#'      pos.cens.auc = Positive Wasserstein statistic,
#'      dens.plot    = Density plot,
#'      cdf.plot     = CDF plot,
#'      diff.plot    = Difference plot)
#'
#' @export
#'
#' @examplesIf interactive()
#' library(dplyr)
#'
#' guide.seq = "UUAUAGAGCAAGAACACUGUUUU"
#'
#' anno.db = load_annotations("rnor7")
#'
#' # Load test data
#' get_example_data("sirna")
#'
#' sirna.data = load_example_data("sirna")
#'
#' res <- sirna.data$Schlegel_2022_Ttr_D1_30mkg
#'
#' # Filter DESeq2 results for SeedMatchR
#' res = filter_res(res, fdr.cutoff=1, fc.cutoff=0, rm.na.log2fc = TRUE)
#'
#' res = SeedMatchR(res = res, gtf = anno.db$gtf, seqs = anno.db$seqs,
#' sequence = guide.seq, seed.name = "mer7m8", tx.id.col= FALSE)
#'
#' # Gene set 1
#' mer7m8.list = res$log2FoldChange[res$mer7m8 >= 1]
#'
#' # Gene set 2
#' background.list = res$log2FoldChange[!(res$mer7m8 %in% mer7m8.list)]
#'
#' ecdf.res = SeedMatchR::wass_dist(res, mer7m8.list, background.list)
wass_dist <- function(seed, non.seed, n.interp=10000){
  if (length(seed) < 10){
    stop("wass_dist: Error: length of the seed values is < 10")
  }
  if (length(non.seed) < 10){
    stop("wass_dist: Error: length of the non-seed values is < 10")
  }

  # ----------------------------------------------------
  # Compute one-sided area (assumes non-seed N > seed N)
  # ----------------------------------------------------
  non.seed.cdf <- stats::ecdf(non.seed)
  d.non.seed <- data.frame(x=non.seed, y.non.seed=non.seed.cdf(non.seed)) %>%
                dplyr::arrange(x)

  seed.cdf <- stats::ecdf(seed)
  d.seed <- data.frame(x=seed, y.seed=seed.cdf(seed)) %>%
            dplyr::arrange(x)

  d.approx <- data.frame(x=seq(from=min(seed, non.seed), to=max(seed, non.seed), length=n.interp))
  d.approx$y.seed <- stats::approx(x=d.seed$x, y=d.seed$y.seed, xout=d.approx$x, yleft=0, yright=1)$y
  d.approx$y.non.seed <- stats::approx(x=d.non.seed$x, y=d.non.seed$y.non.seed, xout=d.approx$x, yleft=0, yright=1)$y

  d.approx <- d.approx %>%
              dplyr::mutate(diff = y.seed - y.non.seed) %>%
              dplyr::mutate(abs.diff = abs(diff)) %>%
              dplyr::mutate(diff.neg.cens = dplyr::case_when(diff < 0 ~ 0, TRUE ~ diff)) %>%
              dplyr::mutate(diff.pos.cens = dplyr::case_when(diff > 0 ~ 0, TRUE ~ abs(diff)))

  neg.cens.auc <- round(caTools::trapz(d.approx$x, d.approx$diff.neg.cens), 3) # Negative values of the (seed - non.seed eCDF) set to 0 prior to integration
  pos.cens.auc <- round(caTools::trapz(d.approx$x, d.approx$diff.pos.cens), 3) # Positive values of the (seed - non.seed eCDF) set to 0 prior to integration
  abs.auc <- round(caTools::trapz(d.approx$x, d.approx$abs.diff), 3)           # Integral of the |seed - non.seed eCDF|

  # ----------------
  # Overlay CDF Plot
  # ----------------
  d.x <- data.frame(group=rep("non-seed", length(non.seed)), val=non.seed)
  d.y <- data.frame(group=rep("seed", length(seed)), val=seed)
  d <- rbind(d.x, d.y)
  cdf.plot <- ggplot2::ggplot(d, ggplot2::aes(x = val, color = group)) +
    ggplot2::scale_color_manual(values = c("black", "red")) +
    ggplot2::stat_ecdf(lwd = 1.25, geom = "line") +
    ggplot2::scale_x_continuous("log2(FC)") +
    ggplot2::scale_y_continuous("eCDF")

  # --------------------
  # Overlay density plot
  # --------------------
  dens.plot <- ggplot2::ggplot(d, ggplot2::aes(x=val, fill=group)) +
    ggplot2::geom_density(alpha=.25) +
    ggplot2::scale_x_continuous("log2(FC)") +
    ggplot2::scale_y_continuous("Density")

  # -------------------
  # CDF difference plot
  # -------------------
  main.title <- paste("Wass(- censored) = ", neg.cens.auc, "   Wass(+ censored) = ", pos.cens.auc, "   Wass (total) = ", abs.auc, sep="")
  diff.plot <- ggplot2::ggplot(d.approx, ggplot2::aes(x=x, y=diff)) +
    ggplot2::geom_line() +
    ggplot2::scale_x_continuous("log2(FC)") +
    ggplot2::scale_y_continuous("Seed - NonSeed eCDF") +
    ggplot2::ggtitle(main.title)

  # ------
  # Output
  # ------
  out <- list(abs.auc      = abs.auc,
              neg.cens.auc = neg.cens.auc,
              pos.cens.auc = pos.cens.auc,
              dens.plot    = dens.plot,
              cdf.plot     = cdf.plot,
              diff.plot    = diff.plot)
  return(out)
}
