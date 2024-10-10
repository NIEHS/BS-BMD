library(splines)
library(Matrix)
library(invgamma)
library(readr)
library(lattice)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)
library(reshape2)
library(data.table)
library(stringr)
library(tables)
require(gridExtra)




#' Creates box or ridge plots for the posterior BMD given the BMD MCMC chains.
#'
#' The plots will display hallmark sets on the Y axis and the log-scale dose on
#' the X axis.  Each individual plot for a chemical shows two sets of boxes, one
#' for liver and one for kidney
#'
#' @param full_df data frame with each column an MCMC iteration and each row a
#'   hallmark set, along with columns for the tissue and chemical and hallmark
#'   set name
#' @param BMD_MCMC_iter integer number of MCMC iterations per tissue/chemical
#' @param burnin number of initial MCMC iterations to discard
#' @param quantile_cut the quantile used to compute the BMD,  used for reference
#' @param sort boolean to sort the hallmark sets by their mean value (aggregates
#'   across tissues), default is True
#' @param use_ridgelets boolean to show ridge plot instead of box plots. Default
#'   is False
#' @param file_label a string descriptor to be pre-pended to the file name
#'   default is "biplots"
make_fancy_ggbox <- function(full_df, BMD_MCMC_iter, burnin = 0,
                             quantile_cut = 0.78, sort = TRUE,
                             use_ridgelets = FALSE,
                             file_label = "biplots") {
  effective_iter <- BMD_MCMC_iter - burnin
  full_df <- full_df[, c(1:3, (1 + burnin):BMD_MCMC_iter + 3)]
  significance <-
    apply(full_df[, 1:effective_iter + 3],
      MARGIN = 1,
      FUN = function(rx) ifelse(FALSE, "Not Significant", "Significant")
    )

  plotdata <- full_df %>%
    mutate(Hallmark = sapply(Hallmark, function(x) {
      unlist(str_split(x, pattern = "HALLMARK_"))[2]
    })) %>%
    reshape2::melt(., id.vars = c("Chemical", "Tissue", "Hallmark")) %>%
    rename(., c(BMD = value)) %>%
    mutate(BMD = as.numeric(BMD)) %>%
    dplyr::select(Chemical, Tissue, Hallmark, BMD)
  chemicals <- unique(plotdata$Chemical)

  jsh_hallmark_biplot <- function(mychem = chemicals[1]) {
    mydata <- plotdata %>% dplyr::filter(Chemical == mychem)
    plt_title <- mychem # paste0(mychem,  " - N=", BMD_MCMC_iter)
    if (length(grep("thujone", mychem))) {
      plt_title <- expression(paste(alpha, beta, " Thujone"))
    }
    if (nrow(mydata) > 5) {
      if (sort) {
        sort_function <- function(x) {
          quantile(x, 0.05) + 1e4 * (quantile(x, .90) == max(x))
        }
        mydata$Hallmark <- with(mydata, reorder(Hallmark, BMD, mean))
      }
      ggplot(data = mydata, aes(x = Hallmark, y = BMD)) +
        ylab("Pathway BMD - mg/kg") +
        xlab("Hallmark Pathway") +
        ggtitle(plt_title) +
        geom_boxplot(
          outlier.size = .5, outlier.alpha = .1,
          outlier.shape = 1, outlier.color = "red",
          color = "black"
        ) +
        stat_summary(
          fun = mean, geom = "point", shape = 23,
          size = 1.2, color = "purple", fill = NA, aes(colour = "meanLines")
        ) +
        stat_summary(
          fun = function(x) quantile(x, 0.05), geom = "point", shape = 3,
          size = 1.2, color = "blue", fill = "blue", aes(colour = "meanLines")
        ) +
        scale_x_discrete(limits = rev(levels(mydata$Hallmark))) +
        theme_bw() +
        theme(legend.position = "top") +
        theme(axis.text.y = element_text(size = 8)) +
        coord_flip() +
        scale_y_continuous(trans = "log10") +
        annotation_logticks(sides = "b") +
        facet_wrap(Chemical ~ Tissue) -> biplot
    }

    return(biplot)
  }


  hallmark_ridgeplot <- function(mychem = chemicals[1]) {
    mydata <- plotdata %>% dplyr::filter(Chemical == mychem)
    plt_title <- mychem
    # pull EPA values
    if (length(grep("thujone", mychem))) {
      plt_title <- expression(paste(alpha, beta, " Thujone"))
    }
    x_label <- "Pathway BMD - mg/kg"
    if (mychem == "EE") {
      x_label <- "Pathway BMD - ug/kg"
    }

    if (nrow(mydata) > 5) {
      if (sort) {
        sort_function <-
          function(x) quantile(x, 0.05) + 1e4 * (quantile(x, .90) == max(x))
        mydata$Hallmark <- with(mydata, reorder(Hallmark, BMD, mean))
      }
      ggplot(data = mydata, aes(y = Hallmark, x = BMD, fill = stat(x))) +
        geom_density_ridges_gradient(
          scale = 3, rel_min_height = 0.01) +
        xlab(x_label) +
        ylab("Hallmark Pathway") +
        scale_y_discrete(limits = rev(levels(mydata$Hallmark))) +
        theme_bw() +
        theme(legend.position = "None") +
        theme(axis.text.y = element_text(size = 8)) +
        scale_x_continuous(trans = "log10") +
        stat_summary(
          fun = function(x) log(mean(10^x), base = 10), geom = "point", shape = 23, stroke = .4,
          size = 1.2, color = "black", fill = "green", aes(colour = "meanLines")
        ) +
        stat_summary(
          fun = function(x) quantile(x, 0.05), geom = "point", shape = 24,
          size = 1.2, stroke = .4, color = "black", fill = "cyan",
          aes(colour = "meanLines")
        ) +
        #coord_trans(x = "log10", xlim = c(1, 1000)) +
        annotation_logticks(sides = "b") +
        scale_fill_viridis_c(
          name = "Log BMD", option = "C",
          direction = -1, trans = "log"
        ) +
        facet_wrap(Chemical ~ Tissue) -> biplot_ridge
    }
    return(biplot_ridge)
  }

  biplots_chems <- list()
  for (i in seq_along(chemicals)) {
    if (i == 4 || i == 20) next # skip BPAF, Triclosan
    if (use_ridgelets) {
      biplots_chems[[i]] <- try(hallmark_ridgeplot(mychem = chemicals[i]))
    } else {
      biplots_chems[[i]] <- try(jsh_hallmark_biplot(mychem = chemicals[i]))
    }
  }

  filename <- paste("output/figs/",
    file_label, "_",
    quantile_cut, "qtile_",
    BMD_MCMC_iter, "-", burnin, ".pdf",
    sep = ""
  )
  print(filename)
  pdf(file = filename, width = 8, height = 7)
  for (i in seq_along(biplots_chems)) {
    print(biplots_chems[i])
  }
  dev.off()
}


# plotting GO results:  may need scatte rplot?



#' Plot the gene set test statistic and individual gene contributions
#'
#' Running this script requires the MCMC parameter estimates to be in memory.
#' Ie: run the MCMC, stop at some iteration, and then run this function to
#' create a plot of the hallmark set test statistic broken down by gene.
#'
#' @param hallmark_idx Integer for index of a specific hallmark set to plot
#' @param log_plot boolean to apply log scaling to the x-axis, default is True
#' @param saveplot_pdf boolean to save the plot as a PDF, default is False
#' @param return_cut
#' @param zero_center
#' @param plt_title
set_bmd_breakdown <- function(hallmark_idx, log_plot = TRUE,
                              saveplot_pdf = FALSE,
                              return_cut = FALSE,
                              zero_center = TRUE,
                              plt_title = "Stacked Gene Contributions") {
  clust_idx <- clust_groups[[hallmark_idx]]
  gene_ids <- curr_data$Probe_ID[clust_idx]
  clust_size <- length(clust_idx)
  B <- solve(mat[clust_idx, clust_idx])
  R <- t(ju %*% betas[, clust_idx]) # allow for interpolator
  if (zero_center) R <- R[, ] - R[, 1]
  # Form my Test Statistic
  TS <- t(R) %*% B
  TS <- t(TS) * R
  G <- colSums(TS)
  # Null Response
  R_null <- t(ju %*% betas_null[, clust_idx])
  R_null <- R_null[, ] - R_null[, 1]
  # Form my Test Statistic
  TS_null <- t(R_null) %*% tau_mat[clust_idx, clust_idx]
  TS_null <- t(TS_null) * R_null
  G_null <- colSums(TS_null)
  # interpolated curves
  dose_seq <- seq(0, log(max(doses) + 1), 0.0525)
  hallmark_bubs <- splinefun(dose_seq, G)
  null_bubs <- splinefun(dose_seq, G_null)
  test_bubs <- splinefun(dose_seq, G - G_null)
  hallmark_curve <- hallmark_bubs(dose_seq)
  null_curve <- null_bubs(dose_seq)
  test_curve <- test_bubs(dose_seq)
  # collect various statistics to plot
  hallmark_df <- data.frame(
    "Dose" = exp(dose_seq[1:length(dose_seq)]) - 1,
    "Response" = hallmark_curve,
    "Noise" = null_curve,
    "Test" = test_curve
  )
  # individual gene curves
  gene_curves <- matrix(0, nrow = length(clust_idx), ncol = nrow(ju))
  decor <- t(chol(mat[clust_idx, clust_idx]))
  # compute decorrelated responses for genes in set
  R_decor <- t(ju %*% betas[, clust_idx] %*% solve(t(decor)))
  R_decor_cent <- R_decor
  if (zero_center) R_decor_cent <- R_decor - R_decor[, 1]
  # square to get the effective test statistic for the individual genes
  gene_contrib <- R_decor_cent^2
  # estimate percent contribute by comparing maximum responses for each gene
  max_gene_contrib <- apply(gene_contrib, MARGIN = 1, max)
  gene_ordering <- order(max_gene_contrib, decreasing = TRUE)
  rownames(gene_contrib) <- gene_ids
  # compute precent contribution by summing up the maximum contributions
  gene_percent_contr <- max_gene_contrib[gene_ordering] / sum(max_gene_contrib)
  # get names of top genes for annotation
  top_genes <- gene_ids[gene_ordering[1:5]]
  top_gene_labels <- paste(gene_ids[gene_ordering[1:5]], " ",
    round(100 * gene_percent_contr[1:5], 1),
    "%",
    sep = ""
  )
  df <- reshape2::melt(gene_contrib[gene_ordering, ])
  names(df) <- c("Gene", "Dose", "Response")
  df$Dose <- exp(dose_seq[df$Dose]) - 1
  hallmark_id <- unique(gene_probe_map$MSigDBName)[hallmark_idx]
  hallmark_id <- strsplit(hallmark_id, split = "HALLMARK_")[[1]][2]
  # use a regex to get the chemical and tissue from the data file name
  chem_tiss <- sub(".*/([^/]+)\\..*", "\\1", curr_file) # thanks chatGPT
  file_name_plot <- paste("Gene_breakdown_", hallmark_id, "_",
    chem_tiss, ".pdf",
    sep = ""
  )
  plot_title <- plt_title
  test_threshold <- qchisq(0.98, df = clust_size) -
    qchisq(0.78, df = clust_size)
  # the initial plot is not log-scaled and difficult to analyze
  my_gg <- ggplot(data = df, aes(x = Dose, y = Response, col = Gene)) +
    geom_line(
      data = hallmark_df,
      aes(x = Dose, y = Response, col = "Hallmark Set"),
      col = "black", linewidth = 1.5
    ) +
    geom_line(position = "stack") +
    geom_hline(
      yintercept = test_threshold,
      linetype = 2, color = "black",
      linewidth = 1.25
    ) +
    geom_line(
      data = hallmark_df, aes(x = Dose, y = Noise),
      col = "gray", linewidth = 1.25
    ) +
    geom_line(
      data = hallmark_df, aes(x = Dose, y = Test),
      col = "red", linewidth = 1.25
    ) +
    scale_color_discrete(
      breaks = top_genes,
      name = "Top Contributing Genes",
      labels = top_gene_labels
    ) +
    theme(legend.position = "bottom") +
    ggtitle(plot_title)
  if (log_plot) {
    # reference color scale that is safe for colorblind
    okabe <- c(
      "#E69F00", "#56B4E9", "#009E73", "#F0E442",
      "#0072B2", "#D55E00", "#CC79A7"
    )
    plotcolors <- c(
      "Aggregated Response" = "#009E73",
      "Null Response" = "#E69F00",
      "Test Statistic" = "#56B4E9",
      "Test Cutoff" = "black"
    )
    df$Dose <- df$Dose + 0.01
    hallmark_df$Dose <- hallmark_df$Dose + 0.01
    my_gg <- ggplot(data = df, aes(x = Dose, y = Response)) +
      geom_line(
        data = hallmark_df,
        aes(x = Dose, y = Response),
        col = "black", linewidth = 2
      ) +
      geom_line(
        data = hallmark_df,
        aes(x = Dose, y = Response, col = "Aggregated Response"),
        linewidth = 1.5
      ) +
      geom_line(position = "stack", aes(group = Gene), col = "gray") +
      geom_hline(aes(yintercept = test_threshold, col = "Test Cutoff"),
        linetype = 2, linewidth = 1.25
      ) +
      geom_line(
        data = hallmark_df, aes(x = Dose, y = Noise),
        col = "black", linewidth = 1.75
      ) +
      geom_line(data = hallmark_df, aes(
        x = Dose, y = Noise,
        col = "Null Response"
      ), linewidth = 1) +
      geom_line(
        data = hallmark_df, aes(x = Dose, y = Test),
        col = "black", linewidth = 1.75
      ) +
      geom_line(
        data = hallmark_df,
        aes(x = Dose, y = Test, col = "Test Statistic"),
        linewidth = 1.25
      ) +
      scale_color_manual(values = plotcolors) +
      theme(legend.position = "bottom") +
      ggtitle(plot_title) +
      xlab("Dose (Axis Log-Scaled)")
    my_gg <- my_gg + scale_x_log10() + annotation_logticks(sides = "b")

    # annotate the top genes
    y_vals <- rev(cumsum(gene_contrib[rev(gene_ordering), nrow(ju) - 2]))
    for (i in 1:3) {
      my_gg <- my_gg + annotate("segment",
        x = exp(dose_seq[15]),
        xend = exp(dose_seq[nrow(ju) - 3]),
        y = y_vals[i],
        yend = y_vals[i],
        colour = "black"
      ) +
        annotate("text",
          label = paste("Gene", top_genes[i]),
          y = y_vals[i], x = exp(dose_seq[1]) - .5
        )
    }
    my_gg
  }
  if (saveplot_pdf) {
    ggsave(file_name_plot, plot = my_gg, device = "pdf", width = 7, height = 5)
  }
  # other plots use the test
  if (return_cut) {
    return(list(my_gg, test_threshold, test_statistics))
  }
  my_gg
}



plot_corrmat <- function(lambda, taus, thinning = 10,
                         cuts = 10, q_order = FALSE) {
  require(ggplot2)
  require(tidyr)
  require(tibble)
  require(viridis)

  cormat <- 0
  n_genes <- 0
  if (typeof(lambda) == "list") {
    cor_mat_avg <- 0
    for (i in 1:length(lambda)) {
      print(paste("reading", i, "th", "lambda"))
      covmat <- lambda[[i]] %*% t(lambda[[i]]) +
        as.numeric(1 / taus[[i]]) * diag(length(taus[[i]]))
      cormat_i <- cov2cor(covmat)
      diag(cormat_i) <- 0
      cor_mat_avg <- cor_mat_avg + cormat_i
    }
    n_genes <- nrow(lambda[[1]])
    cormat <- cor_mat_avg / length(lambda)
  } else if (typeof(lambda) == "double" && length(dim(lambda)) == 3) {
    cor_mat_avg <- 0
    for (i in 1:dim(lambda)[1]) {
      print(paste("reading", i, "th", "lambda"))
      lambda_i <- lambda[i, , ]
      covmat <- lambda_i %*% t(lambda_i) +
        as.numeric(1 / taus[i, ]) * diag(length(taus[i, ]))
      cormat_i <- cov2cor(covmat)
      diag(cormat_i) <- 0
      cor_mat_avg <- cor_mat_avg + cormat_i
    }
    n_genes <- nrow(lambda[1, , ])
    cormat <- cor_mat_avg / dim(lambda)[1]
  } else {
    n_genes <- nrow(lambda)
    covmat <- lambda %*% t(lambda) + as.numeric(1 / taus) * diag(length(taus))
    cormat <- cov2cor(covmat)
    diag(cormat) <- 0
  }

  coul <- colorRampPalette(brewer.pal(8, "Spectral"))(24)
  clust <- hclust(as.dist(1 - as.matrix(cormat), diag = TRUE))
  vals <- cutree(clust, cuts)
  q <- order(vals)
  half_dens <- seq(1, n_genes, by = thinning)

  if (typeof(q_order) == "integer") {
    if (length(q_order) != length(vals)) {
      print("Length Mismatch!")
    } else {
      q <- q_order
    }
  }

  gplot <- cormat[q, q][half_dens, half_dens] %>%
    # Data wrangling
    as_tibble() %>%
    rowid_to_column(var = "X") %>%
    gather(key = "Y", value = "Z", -1) %>%
    # Change Y to numeric
    mutate(Y = as.numeric(gsub("V", "", Y))) %>%
    ggplot(aes(X, Y, fill = Z)) +
    geom_tile() +
    scale_fill_viridis(discrete = FALSE) +
    theme(
      legend.position = "right",
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(), # remove x axis labels
      axis.ticks.x = element_blank(), # remove x axis ticks
      axis.text.y = element_blank(), # remove y axis labels
      axis.ticks.y = element_blank() # remove y axis ticks
    )

  myplot <- fields::imagePlot(cormat[q, q][half_dens, half_dens],
    xaxt = "n", yaxt = "n"
  )
  myplot
  return(gplot)
}



get_corr_ordering <- function(lambda, taus, clust_idx, cuts = 10) {
  lambda_sub <- lambda[clust_idx, ]
  tau_sub <- taus[clust_idx]
  n_genes <- nrow(lambda_sub)
  covmat <- lambda_sub %*% t(lambda_sub) + as.numeric(1 / tau_sub) *
    diag(length(tau_sub))
  cormat <- cov2cor(covmat)
  diag(cormat) <- 0
  clust <- hclust(as.dist(1 - as.matrix(cormat), diag = TRUE))
  vals <- cutree(clust, cuts)
  q <- order(vals)
  return(q)
}


#' Collect BMD estimates into a single data frame
#'
#' Reads in BSBMD, EPA ETAP and/or Apical bmd estimates and creates an ordered
#' table for presenting and making comparisons between the methods
#'
#' @param hallmark_bigset_df data frame with chemical, tissue, hallmark, and BMD
#'   samples.  Number of rows equals (hallmark sets) x (tissues) x (chemicals),
#'   number of columns is number of MCMC iterations + 3 (chem, tiss, hallmark)
#' @param take_top_n boolean or integer to subset the collected data.  Defaults
#'   to False, ie take all data.
#' @param apical_top boolean to compare BSBMD to apical or EPA, default is False
#' @param epa_vs_ap boolean to compare EPA and apical, default is False
#' @param burnin integer MCMC burnin interations to discard, default is 500
#'
#' @return a data frame
#'
compare_bmd_EPA <- function(hallmark_bigset_df, take_top_n = FALSE,
                            apical_top = FALSE, epa_vs_ap = FALSE,
                            burnin = 500) {
  # Compare BMDL = 5 percentile, BMD Express 5%, and Apical endpoint
  # Correlation?

  # hallmark_df <- hallmark_bigset_df %>%
  #  mutate(BMD = apply(.data[[as.char(1:BMD_MCMC_iter+3]],
  #                     MARGIN=1,
  #                     FUN = function(x) quantile(x,0.05))) %>%
  #  dplyr::select(Chemical,Tissue, Hallmark, BMD)

  hallmark_df <- hallmark_bigset_df[, 1:3]
  hallmark_BMDLs <- apply(hallmark_bigset_df[, burnin:BMD_MCMC_iter],
    MARGIN = 1,
    FUN = function(x) quantile(x, 0.05)
  )
  hallmark_BMDs <- apply(hallmark_bigset_df[, burnin:BMD_MCMC_iter],
    MARGIN = 1,
    FUN = function(x) median(x)
  )
  significant_BMDs <- apply(hallmark_bigset_df[, burnin:BMD_MCMC_iter],
    MARGIN = 1,
    FUN = function(rx) {
      ifelse(quantile(rx, .95) == max(rx),
        "Not Significant",
        "Significant"
      )
    }
  )

  hallmark_df$significance <- significant_BMDs
  hallmark_df$BMD <- hallmark_BMDs
  hallmark_df$BMDL <- hallmark_BMDLs
  hallmark_df <- hallmark_df %>%
    group_by(Chemical, Tissue) %>%
    mutate(Rank = min_rank(BMDL))

  # Get teh top 3 bmd for each chemical/tissue
  top3_BMDL <- hallmark_df %>%
    # filter(significance == "Significant") %>%
    group_by(Chemical, Tissue) %>%
    slice_min(BMDL, n = 3, with_ties = FALSE) %>%
    select(Chemical, Tissue, Hallmark, BMD, BMDL)

  top3_BMDL$Hallmark <- unlist(lapply(
    top3_BMDL$Hallmark,
    FUN = function(x) {
      istr <- strsplit(x, split = "HALLMARK_")[[1]][2]
      nistr <- str_replace_all(string = istr, pattern = "_", replacement = " ")
      stringr::str_to_title(nistr)
    }
  ))


  # get the BMD for 5-day study
  etap_go_hlmk <- readxl::read_xlsx("../../scott_epa/ETAP Final GO + Hallmark.xlsx")
  ETAP_hallmark_df <- etap_go_hlmk %>% filter(`Gene Set` == "Hallmark")%>%
    transform("Active_Genes" = paste0(`Genes That Passed All Filters`,
                                    "/",
                                    `All Genes (Expression Data)`)) %>%
    select(Chemical, Organ,
           `GO/Pathway/Gene Set/Gene Name`,
           `BMD Median`,
           `BMDL Median`,
          `Active_Genes`) %>%
    rename( Tissue = `Organ`,
            Hallmark = `GO/Pathway/Gene Set/Gene Name`,
            BMDL_Median_epa = `BMDL Median`,
            BMD_Median_epa = `BMD Median`) 
  

  epa_hallmark_df <- ETAP_hallmark_df %>%
    group_by(Chemical, Tissue) %>%
    mutate(Rank_EPA = min_rank(BMDL_Median_epa))

  # join EPA and our data
  comp_df <- left_join(hallmark_df,
    epa_hallmark_df,
    by = c("Chemical", "Tissue", "Hallmark")
  )
  # format Hallmark col
  comp_df$Hallmark <- unlist(lapply(
    comp_df$Hallmark,
    FUN = function(x) {
      istr <- strsplit(x, split = "HALLMARK_")[[1]][2]
      nistr <- str_replace_all(string = istr, pattern = "_", replacement = " ")
      stringr::str_to_title(nistr)
    }
  ))
  
  

  if (take_top_n > 0) {
    # return a table with the lowest n BMDs from our method with the
    # corresponding BMD (ie matched by hallmark) from EPA
    sub_comp_df <- comp_df %>%
      # filter(significance == "Significant") %>%
      group_by(Chemical, Tissue) %>%
      slice_min(BMD, n = take_top_n) %>%
      select(
        Chemical, Tissue, Hallmark, BMD, BMDL,
        Rank, BMD_Median_epa, BMDL_Median_epa,
        Rank_EPA, Active_Genes
      )
    return(sub_comp_df)

    # create table: for each chemical
    primary_indices <- c()
    for (chem in unique(sub_comp_df$Chemical)) {
      for (tiss in unique(sub_comp_df$Tissue)) {
        chem_idx <- which(sub_comp_df$Chemical == chem)
        tiss_idx <- which(sub_comp_df$Tissue == tiss)
        my_idx <- intersect(chem_idx, tiss_idx)
        if (length(my_idx) == 0) next
        primary_indices <- c(primary_indices, my_idx[1])
      }
    }
    sub_chem_pt1 <- sub_comp_df[primary_indices, ]
    secondary_indices <- setdiff(1:nrow(sub_comp_df), primary_indices)
    sub_chem_pt2 <- sub_comp_df[secondary_indices, ]
    sub_comp_df <- sub_chem_pt1 %>%
      left_join(sub_chem_pt2,
        by = c("Chemical", "Tissue"), suffix = c(".1", ".2")
      ) %>%
      select(
        Chemical, Tissue, Hallmark.1,
        BMD.1, BMDL.1, Hallmark.2, BMD.2, BMDL.2
      )
  }

  if (apical_top > 0) {
    chem_names <- unique(hallmark_df$Chemical)
    # lookup apical data
    apical_df <- read_apical_input(chem_names, take_top_n = apical_top)
    # sub comp only takes BSBMD values
    sub_comp_df <- comp_df %>%
      # filter(significance == "Significant") %>%
      group_by(Chemical) %>%
      slice_min(BMD, n = min(apical_top, 3), with_ties = FALSE) %>%
      transform(
        BSBMD = BMD,
        BSBMDL = BMDL,
        Chemical.1 = Chemical
      ) %>%
      select(Chemical.1, Tissue, Hallmark, BSBMD, BSBMDL) %>%
      dplyr::mutate_if(is.numeric, round, 2)
    if (epa_vs_ap) {
      sub_comp_df <- comp_df %>%
        # filter(significance == "Significant") %>%
        group_by(Chemical) %>%
        slice_min(BMD_Median_epa, n = min(apical_top, 3)) %>%
        transform(
          ETAP.BMD = BMD_Median_epa,
          ETAP.BMDL = BMDL_Median_epa,
          Chemical.1 = Chemical
        ) %>%
        select(Chemical.1, Tissue, Hallmark, ETAP.BMD, ETAP.BMDL) %>%
        dplyr::mutate_if(is.numeric, round, 2)
    }
    # cant do inner join: duplicates the apical entries
    # apical_comp_df <- apical_df %>% inner_join(sub_comp_df, by = "Chemical")
    n_apical <- table(apical_df$Chemical)
    matched_hallmark_idx <-
      unlist(lapply(unique(apical_df$Chemical), FUN = function(x) {
        which(sub_comp_df$Chemical == x)[1:n_apical[x]]
      }))
    comp_apical_df <- cbind(
      apical_df,
      sub_comp_df[matched_hallmark_idx, ]
    )
    # write_csv(comp_apical_df, "apical_bmd.csv")
    return(comp_apical_df)
  }
  return(comp_df)
}

read_apical_input <- function(chem_names, tableid = 5, take_top_n = 1) {
  apical_df <- readxl::read_xlsx("../../scott_epa/apical_data.xlsx") %>%
    filter(Table_id == tableid) %>%
    select(!Table_id)
  # correct the abbreviations
  abbr_indx <- stringr::str_to_upper(chem_names) %>%
    stringr::str_sub(., end = 4)
  abbr_indx[1] <- "THU"
  abbr_indx[8] <- "EE2"
  abbr_indx[15] <- "TBBPA"


  map_chemID <-
    sapply(unique(apical_df$Chemical), FUN = function(x) grep(x, abbr_indx))
  chems_present <- which(sapply(map_chemID, length) > 0)
  # filter out missing chemicals
  apical_df <- apical_df %>%
    filter(Chemical %in% names(chems_present)) %>%
    mutate(
      Chemical = chem_names[unlist(map_chemID[apical_df$Chemical])],
      BMDL = as.numeric(BMDL),
      BMD = as.numeric(BMD)
    ) %>%
    group_by(Chemical) %>%
    slice_min(BMDL, n = take_top_n)

  if (tableid == 5 && take_top_n == 1) {
    # Fenofibrate BMDL 0.269
    apical_df[7, 4] <- 0.269 # from TGGATEs
    # DEHP: 2year Pancreas adenoma carcinoma 31.2  20.3
    apical_df[which(apical_df$Chemical == "DEHP"), 3:4] <- list(31.2, 20.3)
    apical_df[which(apical_df$Chemical == "DEHP"), 2] <-
      "Pancreas adenoma carcinoma"

    # HCB Chronic nephrosis severe  0.59 0.35
    apical_df[which(apical_df$Chemical == "HCB"), 3:4] <- list(0.59, 0.35)
    apical_df[which(apical_df$Chemical == "HCB"), 2] <-
      "Chronic nephrosis severe"

    # PFOA  Liver hepatocyte hypertrophy 0.5 0.41
    apical_df[which(apical_df$Chemical == "PFOA"), 3:4] <- list(0.5, 0.41)
    apical_df[which(apical_df$Chemical == "PFOA"), 2] <-
      "Liver hepatocyte hypertrophy"

    # Pulegon Nose olfactory epithelial degeneration 14.02 9.47
    apical_df[which(apical_df$Chemical == "Pulegone"), 3:4] <- list(14.02, 9.47)
    apical_df[which(apical_df$Chemical == "Pulegone"), 2] <-
      "Nose olfactory epithelial degeneration"

    # TBBPA 90day Total thyroid 56.5 46.1
    apical_df[which(apical_df$Chemical == "TBBPA"), 3:4] <- list(56.5, 46.1)
    apical_df[which(apical_df$Chemical == "TBBPA"), 2] <-
      "Total thyroid"

    # TCPP lung granulomatous focal inflammation 223.6 58.8
    apical_df[which(apical_df$Chemical == "TCPP"), 3:4] <- list(223.6, 58.8)
    apical_df[which(apical_df$Chemical == "TCPP"), 2] <-
      "Lung granulomatous focal inflammation"

    # Triclosan zorrila 2009 thyroid hormone 14.5 7.23 (BMR: 20%)
    apical_df[which(apical_df$Chemical == "TCPP"), 3:4] <- list(14.5, 7.23)
    apical_df[which(apical_df$Chemical == "TCPP"), 2] <-
      "Lung granulomatous focal inflammation"

    # MTE/ Ginseng
    # Bisphenol a f


    apical_df[which(apical_df$Chemical == "EE"), ] <-
      apical_df %>%
      filter(Chemical == "EE") %>%
      dplyr::mutate_if(is.numeric, function(x) x / 1000)
  }


  return(apical_df)
}

compare_GO_apical <- function(GO_df, n_terms = 3, use_ETAP = FALSE) {
  #  GO_df = GO_lowest5_df

  chem_names <- unique(GO_df$Chemical)
  # lookup apical data
  apical_df <- read_apical_input(chem_names)
  sub_comp_df <- GO_df %>%
    group_by(Chemical) %>%
    slice_min(BMDL, n = n_terms) %>%
    rename(
      GO.BSBMD = BMD,
      GO.BSBMDL = BMDL
    ) %>%
    select(Chemical, Tissue, GO.Term, GO.BSBMD, GO.BSBMDL) %>%
    dplyr::mutate_if(is.numeric, round, 2)

  if (use_ETAP) {
    etap_go_hlmk <- readxl::read_xlsx("../../scott_epa/ETAP Final GO + Hallmark.xlsx")
    ETAP_apical_GO_df_tot <- etap_go_hlmk %>% filter(`Gene Set` == "GOBP")%>%
      select(Chemical, Organ,
             `GO/Pathway/Gene Set/Gene ID`,
             `GO Level`,
             `GO/Pathway/Gene Set/Gene Name`,
             `BMD Median`,
             `BMDL Median`) %>%
      rename( Tissue = `Organ`,
              BMDL = `BMDL Median`,
              BMD = `BMD Median`) %>%
      filter(Chemical %in% unique(GO_df$Chemical))

    sub_comp_df <- ETAP_apical_GO_df_tot %>%
      group_by(Chemical) %>%
      slice_min(BMD, n = n_terms, with_ties = F) %>%
      select(1, 2, 3, 6, 7) %>%
      rename(
        ETAP.GO.BMD = BMD,
        ETAP.GO.BMDL = BMDL
      ) %>%
      mutate_if(is.numeric, round, 3)
  }


  sub_comp_df[which(sub_comp_df$Chemical == "EE"), ] <-
    sub_comp_df %>%
    filter(Chemical == "EE") %>%
    dplyr::mutate_if(is.numeric, function(x) x / 1000)


  if (n_terms == 1) {
    # cant do inner join: duplicates the apical entries
    # apical_comp_df <- apical_df %>% inner_join(sub_comp_df, by = "Chemical")
    n_apical <- table(apical_df$Chemical)
    matched_GO_idx <-
      unlist(lapply(unique(apical_df$Chemical), FUN = function(x) {
        which(sub_comp_df$Chemical == x)[1:n_apical[x]]
      }))
    comp_apical_df <- cbind(
      apical_df,
      sub_comp_df[matched_GO_idx, ]
    )
    return(comp_apical_df)
  }
  return(sub_comp_df)
}

#' Compare BMD ordering
#'
#' Creates a bipartite graph with the orderings of the EPA ETAP and BSBMD
#' method.  Hallmark sets are paired across the two methods and the change in
#' position is indicated by number and color.
#'
#' @param comp_df data frame, output of the function [compare_bmd_EPA]
#' @param chem string, the specific chemical to create a plot for
#'
create_ordering_vis <- function(comp_df, chem = "Fenofibrate") {
  BSBMD_ordered <- comp_df %>%
    # filter(significance == "Significant") %>%
    filter(Chemical == chem) %>%
    filter(Tissue == "Liver") %>%
    group_by(Hallmark) %>%
    slice_min(BMD, n = min(1)) %>%
    dplyr::mutate_if(is.numeric, round, 2) %>%
    arrange(BMDL) %>%
    mutate(Hallmark_BMD = paste0(
      Hallmark, " (",
      # BMD, ", ",
      BMDL, ")"
    )) %>%
    select(Chemical, Hallmark, BMDL, Hallmark_BMD)
  BSBMD_ordered$BSBMD_rank <- 1:nrow(BSBMD_ordered)

  ETAP_ordered <- comp_df %>%
    filter(Chemical == chem) %>%
    # filter(significance == "Significant") %>%
    filter(Tissue == "Liver") %>%
    filter(complete.cases(.)) %>%
    group_by(Hallmark) %>%
    slice_min(BMD_Median_epa, n = min(1)) %>%
    dplyr::mutate_if(is.numeric, round, 2) %>%
    arrange(BMDL_Median_epa) %>%
    mutate(Hallmark_BMD = paste0(
      Hallmark, " (",
      # BMD_Median_epa, ", ",
      BMDL_Median_epa, ")"
    )) %>%
    select(Chemical, Hallmark, BMDL_Median_epa, Hallmark_BMD)
  ETAP_ordered$ETAP_rank <- 1:nrow(ETAP_ordered)

  ordered_etap <- ETAP_ordered %>%
    arrange(Hallmark) %>%
    select(ETAP_rank)
  ordered_bsbmd <- BSBMD_ordered %>%
    arrange(Hallmark) %>%
    select(Chemical, Hallmark, BSBMD_rank)

  d <- left_join(ordered_bsbmd, ordered_etap, by = "Hallmark")

  # diff_vec <- (d$BSBMD_rank-d$ETAP_rank)[order(d$ETAP_rank)]
  diff_vec <- (d$ETAP_rank - d$BSBMD_rank)[order(d$BSBMD_rank)]
  colnames(d)[3:4] <- c("Bayesian Set BMD", "EPA ETAP")
  d <- d %>% pivot_longer(cols = c("Bayesian Set BMD", "EPA ETAP"))
  d$name <- factor(d$name, levels = c("EPA ETAP", "Bayesian Set BMD"))
  # Pad the names with spaces so they don’t overlap the points when plotted
  d <- d %>% mutate(Hallmark = str_c("        ", Hallmark, "  "))


  diff_col <- ifelse(diff_vec < 0, "red", "green")
  diff_shape <- rep(16, 100)
  diff_shape[which(is.na(d$value)) - 1] <- 1 # ifelse(is.na(diff_vec), 1, 16)

  plotd <- d %>% ggplot(aes(x = name, y = value)) +
    geom_point(size = 2, shape = diff_shape) +
    geom_line(aes(group = Hallmark)) +
    geom_text(aes(label = Hallmark), hjust = "outward") +
    annotate("text", x = 2.1, y = 1:50, label = diff_vec, col = diff_col) +
    # ylab("Rank")+xlab("Method")+
    ggtitle(chem) +
    theme(
      axis.line = element_blank(),
      #  axis.text.x=element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "none",
      panel.background = element_blank(),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_blank()
    )

  plotd + scale_y_reverse() + expand_limits(x = c(-.25, 3.5))

  file_name <- paste0("hallmark_ordering_vis_BMDL_", chem, ".pdf")
  ggsave(file_name, width = 8, height = 8)
}

prepare_MCMC_output_for_plot <- function(MCMC_output_file,
                                         BMD_MCMC_iter = 7000,
                                         includes_full_genome = F,
                                         use_GO = F) {
  col_names <- unique(gene_probe_map$MSigDBName)
  n_hlmk_sets <- length(col_names)
  if (includes_full_genome) {
    col_names <- c(col_names, "HALLMARK_FULL_GENOME")
    n_hlmk_sets <- n_hlmk_sets + 1
  }
  if (use_GO) {
    col_names <- unique(clust_groups_GO_df$GO)
    n_hlmk_sets <- length(col_names)
  }
  plt_title <- str_replace_all(MCMC_output_file,
    pattern = "_", replacement = " "
  )
  print(plt_title)
  curr_file_dir <- paste(output_dir, MCMC_output_file, sep = "")
  bmd_mat <- read.table(curr_file_dir)
  # original files need to be transposed:
  bmd_corrected <- t(matrix(array(t(bmd_mat)),
    nrow = n_hlmk_sets,
    ncol = BMD_MCMC_iter, byrow = TRUE
  ))
  chem <- strsplit(plt_title, split = " ")[[1]][1]
  tissue <- strsplit(plt_title, split = " ")[[1]][2]
  # for repeats:
  # chem <- paste(strsplit(plt_title, split = " ")[[1]][1],
  #               strsplit(plt_title, split = " ")[[1]][2])
  # tissue <- strsplit(plt_title, split = " ")[[1]][4]

  mydf <- data.frame(
    "Chemical" = chem,
    "Tissue" = tissue,
    "Hallmark" = col_names
  )
  if (use_GO) {
    colnames(mydf)[3] <- "GO Term"
  }
  widedf <- cbind(mydf, t(bmd_corrected))
  return(widedf)
}

#' Hill function
#'
#' A simple function to compute the Hill function response given parameters and
#' a concentration
#'
#' @param a the maximum effect or sill
#' @param b the EC50
#' @param c the slope
#' @param conc an input concentration or dose
#'
#' @return a real number prediction of the dose response
#' @export
#'
#' @examples
#' hill_function(1, 1.5, 2, 3)
hill_function <- function(a, b, c, conc) {
  a / (1 + (b / conc)^c)
}


get_comparison_dfs = function(){
  
}


scatter_plot_3way <- function(plot_name = "Compare_3way.pdf"){
  # create scatter plot for apical and epa comparison
  BMD_MCMC_iter <- 7000
  
  BSBMD_vs_apical <- compare_bmd_EPA(hallmark_bigset_df,
                                     take_top_n = FALSE,
                                     apical_top = 1
  )
  ETAP_vs_apical <- compare_bmd_EPA(hallmark_bigset_df,
                                    take_top_n = FALSE,
                                    apical_top = 1,
                                    epa_vs_ap = T
  )
    compare_alt <- data.frame(
      "x" = c(
        BSBMD_vs_apical$BSBMDL,
        BSBMD_vs_apical$BSBMDL,
        ETAP_vs_apical$ETAP.BMDL
      ),
      "y" = c(
        ETAP_vs_apical$ETAP.BMDL,
        BSBMD_vs_apical$BMDL,
        ETAP_vs_apical$BMDL
      ),
      "Comparison" = c(
        rep("BS-BMD vs ETAP", nrow(BSBMD_vs_apical)),
        rep("BS-BMD vs Apical", nrow(BSBMD_vs_apical)),
        rep("ETAP vs Apical", nrow(ETAP_vs_apical))
      )
    )
    BMDL_GO_ggplot_alt <- ggplot(compare_alt, aes(x = x, y = y,
                                                  color = Comparison,
                                                  shape = Comparison)) +
      geom_point(size = 3) +
      scale_x_log10() +
      scale_y_log10() +
      geom_abline(slope = 1, intercept = 0) +
      xlab("Approximate Method BMDL (mg/kg)") +
      ylab("Comparison BMDL (mg/kg)") +
      theme(legend.position = "bottom") +
      guides(color = guide_legend(nrow = 2))
    
    ggsave(plot_name, device = "pdf", width = 7, height = 5)
  }
  
  
  
  
  
scatter_plot_4way <- function(plot_name = "Compare_all_table5.pdf"){
  
  BSBMD_vs_apical <- compare_bmd_EPA(hallmark_bigset_df,
                                     take_top_n = FALSE,
                                     apical_top = 1)
  BSBMD_vs_apical<- BSBMD_vs_apical[-16,]
  BSBMD_vs_apical[BSBMD_vs_apical$Chemical == "EE",c(8,9)] <- 
    BSBMD_vs_apical[BSBMD_vs_apical$Chemical == "EE",c(8,9)]/1000
  ETAP_vs_apical <- compare_bmd_EPA(hallmark_bigset_df,
                                    take_top_n = FALSE,
                                    apical_top = 1,
                                    epa_vs_ap = T)
  ETAP_vs_apical<- ETAP_vs_apical[-16,]
  ETAP_vs_apical[ETAP_vs_apical$Chemical == "EE",c(8,9)] <- 
    ETAP_vs_apical[ETAP_vs_apical$Chemical == "EE",c(8,9)]/1000
  BSBMD_vs_apical_GO <- compare_GO_apical(GO_lowest5_df,
                                          n_terms = 1)
  ETAP_vs_apical_GO <- compare_GO_apical(GO_lowest5_df,
                                         n_terms = 1,
                                         use_ETAP = T)
  # compute RMSE to show on legend
  get_RMSE <- function(pt1, pt2) round(sqrt(mean((pt1 - pt2)^2, na.rm = T)), 2)
  
  RMSE_BSBMD_GO <- get_RMSE(
    BSBMD_vs_apical_GO$GO.BSBMDL,
    BSBMD_vs_apical_GO$BMDL
  )
  RMSE_ETAP_GO <- get_RMSE(
    ETAP_vs_apical_GO$ETAP.GO.BMDL,
    BSBMD_vs_apical_GO$BMDL
  )
  RMSE_BSBMD_Hall <- get_RMSE(
    BSBMD_vs_apical$BSBMDL,
    BSBMD_vs_apical$BMDL
  )
  RMSE_ETAP_Hall <- get_RMSE(
    ETAP_vs_apical$ETAP.BMDL,
    ETAP_vs_apical$BMDL
  )
  # legend labels
  cid_1 <- paste("BS-BMD (GO) vs Apical [RMSE: ", RMSE_BSBMD_GO, "]", sep = "")
  cid_2 <- paste("ETAP (GO) vs Apical [RMSE: ", RMSE_ETAP_GO, "]", sep = "")
  cid_3 <- paste("BS-BMD (Hallmark) vs Apical [RMSE: ", RMSE_BSBMD_Hall, "]", sep = "")
  cid_4 <- paste("ETAP (Hallmark) vs Apical [RMSE: ", RMSE_ETAP_Hall, "]", sep = "")
  # construct data frame for ggplot
  compare_all <- data.frame(
    "x" = c(
      BSBMD_vs_apical_GO$GO.BSBMDL,
      ETAP_vs_apical_GO$ETAP.GO.BMDL,
      BSBMD_vs_apical$BSBMDL,
      ETAP_vs_apical$ETAP.BMDL),
    "y" = c(
      BSBMD_vs_apical_GO$BMDL,
      ETAP_vs_apical_GO$BMDL,
      BSBMD_vs_apical$BMDL,
      ETAP_vs_apical$BMDL),
    "Comparison" = c(
      rep(cid_1, nrow(BSBMD_vs_apical_GO)),
      rep(cid_2, nrow(ETAP_vs_apical_GO)),
      rep(cid_3, nrow(BSBMD_vs_apical)),
      rep(cid_4, nrow(ETAP_vs_apical)))
  )
  compare_all <- compare_all[!is.na(compare_all$y), ]
  BMDL_GO_ggplot <- ggplot(compare_all, aes(x = x, y = y,
                                            color = Comparison,
                                            shape = Comparison)) +
    geom_point(size = 3) +
    scale_x_log10() +
    scale_y_log10() +
    geom_abline(slope = 1, intercept = 0) +
    xlab("Approximate Method BMDL (mg/kg)") +
    ylab("Apical BMDL (mg/kg)") +
    theme(legend.position = "bottom") +
    guides(color = guide_legend(nrow = 2))
  ggsave(plot_name, device = "pdf", width = 7.5, height = 5.5)
}


# Calls: ####
# here is where most of the functions are called
if (FALSE) {
  # GG Boxplot ####
  require(stringr)
  output_dir <- "output/"
  file_list_out <- list.files(path = output_dir)
  # use the output_id string or manually select only relevant files
  hallmark_list <- file_list_out[grep("7k_full_genom_clust", file_list_out)]
  hallmark_list_filenames <-
    sapply(hallmark_list,
           FUN = function(x) strsplit(x, split = "[.]")[[1]][1],
           USE.NAMES = FALSE
    )
  col_names <- unique(gene_probe_map$MSigDBName)
  # single data.frame:
  hallmark_bigset_df <- data.frame()
  for (idx in 1:length(hallmark_list_filenames)) {
    plt_title <- str_replace_all(hallmark_list_filenames[idx],
      pattern = "_", replacement = " "
    )
    print(plt_title)
    widedf <- prepare_MCMC_output_for_plot(hallmark_list[idx],
      BMD_MCMC_iter = 7000,
      includes_full_genome = T,
      use_GO = F
    )
    hallmark_bigset_df <- rbind(hallmark_bigset_df, widedf)
  }
  # trace plot test: suggests that variance is decreasing
  # plot(bmd_corrected[, 12])

  # save(hallmark_bigset_df, file= "output/chi2_final_hallmark_BMDs.RData")
  # load("../../chi2_final_hallmark_BMDs.RData")
  sub_idx <- c(103:204)
  make_fancy_ggbox(hallmark_bigset_df[sub_idx, ],
    BMD_MCMC_iter = BMD_MCMC_iter,
    burnin = 2000,
    quantile_cut = 0.78,
    sort = TRUE,
    use_ridgelets = TRUE,
    file_label = "Furan_7k_Full_genome"
  )

  # simple plot to show BMD vs hallmark size, in case of trend
  BMDs <- compute_BMD(hallmark_bigset_df, BMD_MCMC_iter, burnin = 2000)
  hallmark_sizes <- unlist(lapply(clust_groups, length))
  plot(rep(hallmark_sizes, each = 40), BMDs)
  cbind(hallmark_bigset_df[, 1:3], BMDs)

  # compute effective sample sizes
  MCMC_eff_ss <- function(rowvals) {
    1 / (2 + 2 * sum(pacf(rowvals, plot = FALSE)$acf))
  }
  ess_values <- apply(hallmark_bigset_df[, 2000:BMD_MCMC_iter + 3],
    MARGIN = 1, MCMC_eff_ss
  )
  summary(ess_values)
  # 0.2356  0.4469  0.4866  0.4775  0.5119  0.6461      74
  # minimum:  1178 samples
  plot(1:BMD_MCMC_iter, hallmark_bigset_df[226, 1:BMD_MCMC_iter + 3])



  # GO term calculations ####
  # similar to above loop but for GO terms instead of Hallmark sets
  output_dir <- "output/"
  output_dir <- "/Volumes/zilberds/bmd-contain/output/" # server:  
  input_dir <- "input/"
  input_files_dir <- "input/organ_data/"
  probe_map <- readr::read_table(file = "input/Probe File_Rat S1500+.txt")
  gene_set <- 
    readr::read_table(file = "input/Rat_Genome_Whole_H_msigdb_v5.0.txt")
  gene_probe_map <- dplyr::inner_join(gene_set, probe_map, "Entrez_Gene_ID")

  file_list_out <- list.files(path = output_dir)
  # for each data set, read in the corresponding probe sets
  file_list_in <- list.files(path = input_files_dir)

  # use the output_id string or manually select only relevant files
  GO_list <- file_list_out[grep("7k_GO_with_3", file_list_out)]
  GO_list_filenames <-
    sapply(GO_list,
      FUN = function(x) strsplit(x, split = "[.]")[[1]][1],
      USE.NAMES = FALSE
    )
  col_names <- unique(gene_probe_map$MSigDBName)
  # single data.frame:
  GO_lowest5_df <- data.frame()
  for (idx in 1:length(GO_list_filenames)) {
    # make sure probe are correct
    plt_title <- str_replace_all(GO_list_filenames[idx],
      pattern = "_", replacement = " "
    )
    print(plt_title)
    chem <- strsplit(plt_title, split = " ")[[1]][1]
    tissue <- strsplit(plt_title, split = " ")[[1]][2]
    input_index <- 
      intersect(grep(chem, file_list_in), grep(tissue, file_list_in))
    curr_file <- paste(input_files_dir, file_list_in[input_index], sep = "")
    curr_data <- readr::read_table(curr_file, col_names = FALSE, skip = 2)
    names(curr_data)[1] <- "Probe_ID"
    probe_idx <- 
      data.frame("Probe_ID" = curr_data$Probe_ID, "indx" = 1:nrow(curr_data))
    # use the GO to probe dataframe to collect relevant indices into lists
    clust_groups_GO_df <- GO2probe_long %>%
      left_join(probe_idx, by = "Probe_ID") %>%
      na.omit() %>% # some probes not in the data
      group_by(GO) %>%
      summarise(set_indx = list(indx),
                n_genes = length(indx)) %>%
      filter(n_genes > 2 & n_genes < 501)
    go_names <- clust_groups_GO_df$GO # may help for labeling

    GO_BMD_df <- prepare_MCMC_output_for_plot(GO_list[idx],
      BMD_MCMC_iter = 7000,
      includes_full_genome = F,
      use_GO = T
    )
    
    GO_BMDL <- apply(GO_BMD_df[, 2004:7003], MARGIN = 1,
                     FUN = function(rw) quantile(rw, 0.05))
    GO_BMD <- apply(GO_BMD_df[, 2004:7003], MARGIN = 1, mean)
    BMDL_seq <- order(GO_BMDL)[1:5]

    GO_BMD_top5 <- data.frame(
      "Chemical" = chem, "Tissue" = tissue,
      "GO Term" = go_names[BMDL_seq],
      "BMDL" = GO_BMDL[BMDL_seq],
      "BMD" = GO_BMD[BMDL_seq], # order the BMD with BMDL
      "Gene Count" = clust_groups_GO_df$n_genes[BMDL_seq]
    )
    GO_lowest5_df <- rbind(GO_lowest5_df, GO_BMD_top5)
  }

  #save(GO_lowest5_df, file = "output/final_GO_BMDs_with_3small_top5.RData")
  #load("output/final_GO_BMDs_with_3small_top5.RData")
  table_for_supplement <- compare_GO_apical(GO_lowest5_df, n_terms = 3)
  BSBMD_vs_apical_GO <- compare_GO_apical(GO_lowest5_df, n_terms = 1)
  ETAP_vs_apical_GO <- compare_GO_apical(GO_lowest5_df, n_terms = 1, use_ETAP = T)

  ggplot(BSBMD_vs_apical_GO, aes(x = GO.BSBMDL, y = BMDL)) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    geom_abline(slope = 1, intercept = 0) +
    geom_point(data = ETAP_vs_apical_GO, aes(x = ETAP.GO.BMDL, y = BMDL, color = "red"))







  # Gene breakdown ####
  title_6_TGF <- "Stacked gene responses in fenofibrate treated liver for \n the Hallmark TGF-Beta Signaling pathway"
  title_33_Fatty <- "Stacked gene responses in fenofibrate treated liver for \n the Hallmark Fatty Acid Metabolism pathway"

  set_bmd_breakdown(33,
    saveplot_pdf = TRUE, log = TRUE, zero_center = TRUE,
    plt_title = title_33_Fatty
  )

  # CovMat plot ####
  plot_corrmat(lambda, taus)
  cov_ggplot <- plot_corrmat(lambda_record, tau_record[2001:3000, ])
  ggsave(
    filename = "fenofibrate_correlation_thin10.pdf",
    plot = cov_ggplot,
    units = "in",
    width = 5,
    height = 4,
    device = "pdf"
  )


  # lambda list: 1,2,3,4,19,20, 20_noise, 3_noise, 1 noise
  lambda_array <- simplify2array(lambda_list)
  tau_array <- simplify2array(tau_list)
  lambda_mean <- apply(lambda_array, MARGIN = c(1, 2), mean)
  tau_mean <- apply(tau_array, MARGIN = c(1, 2), mean)
  plot1 <- plot_corrmat(lambda_mean, tau_mean)
  plot2 <- plot_corrmat(lambda_null) +
    theme(legend.position = "right") +
    guides(fill = guide_colorbar(title = "Corr"))
  plot3 <- plot_corrmat(lambda_list[[3]], thinning = 5)

  plot_corrmat(lambda_list[[3]][, c(1, 5, 6, 7)])
  ggsave(
    filename = "fenofibrate_correlation_thin4.pdf",
    plot = plot3,
    units = "in",
    width = 6,
    height = 5,
    device = "pdf"
  )

  load("Desktop/Wheeler/lambda_tau_listlist.RData")
  # compare different chemicals, run a bit of everthing
  Lambda_3ds <- lambda_list_list %>%
    lapply(., FUN = simplify2array) %>%
    lapply(., FUN = function(lmat) apply(lmat, MARGIN = c(1, 2), mean))

  tau_3ds <- tau_list_list %>%
    lapply(., FUN = simplify2array) %>%
    lapply(., FUN = function(lmat) apply(lmat, MARGIN = c(1, 2), mean))

  cov_ggplot <- plot_corrmat(Lambda_3ds[[4]], tau_3ds[[4]], thinning = 10)
  ggsave(
    filename = "fenofibrate_correlation_thin10.pdf",
    plot = cov_ggplot,
    units = "in",
    width = 5,
    height = 4,
    device = "pdf"
  )

  # unique(gene_probe_map$MSigDBName)
  clust_to_plot <- 33
  q_new <- get_corr_ordering(
    Lambda_3ds[[3]],
    tau_3ds[[3]],
    clust_groups[[clust_to_plot]]
  )
  # q_new = unlist(clust_groups)
  # hallmark gene set only: clust_groups[[33]]
  gplot_list <- lapply(1:length(tau_3ds),
    FUN = function(idx) {
      plot_corrmat(Lambda_3ds[[idx]][q_new, ],
        tau_3ds[[idx]][q_new],
        thinning = 1,
        q_order = 1:length(q_new)
      )
    }
  )
  cov_pair <- grid.arrange(gplot_list[[1]],
    gplot_list[[3]],
    gplot_list[[5]],
    gplot_list[[2]],
    gplot_list[[4]],
    gplot_list[[6]],
    ncol = 3,
    nrow = 2
  )

  # widths = c(.45,.55))
  ggsave(
    filename = "hallmark_fatty_DE71_Feno_PFOA.png",
    plot = cov_pair,
    units = "in",
    width = 7,
    height = 3,
    device = "png"
  )

  # QQ Plots ####
  qq_df <- data.frame(organ_data - t(X %*% betas) - lambda %*% eta)
  qq_mdf <- reshape2::melt(qq_df)
  rqq_mdf <- qq_mdf # [sample(nrow(qq_mdf), 400),]
  rqq_mdf$value <- (rqq_mdf$value - mean(rqq_mdf$value)) / sd(rqq_mdf$value)
  ggplot(rqq_mdf, aes(sample = value)) +
    stat_qq(shape = 1) +
    geom_abline(slope = 1, intercept = 0) +
    ggtitle("QQ Norm, Data - XB - LF")
  plot(density(rqq_mdf$value))
  lines(seq(-5, 5, by = 0.1),
    dnorm(seq(-5, 5, by = 0.1),
      mean = mean(rqq_mdf$value),
      sd = sd(rqq_mdf$value)
    ),
    col = 3
  )
  # ggsave("BMD_gene_express_qq.pdf", width = 7, height = 5)
  qq_df <- data.frame(organ_data - t(X %*% betas))
  ggplot(melt(qq_df), aes(sample = value)) +
    stat_qq() +
    geom_abline(slope = 1, intercept = 0) +
    ggtitle("QQ Norm, Data - XB")


  # Ordering visualization ####
  comp_df = compare_bmd_EPA(hallmark_bigset_df,
                            take_top_n = 50)
  create_ordering_vis(comp_df, chem = "Fenofibrate")
  create_ordering_vis(comp_df, chem = "Furan")
  create_ordering_vis(comp_df, chem = "Methyleugenol")
  
  # 3-way and 4-way Scatter ####
  scatter_plot_3way() 
  scatter_plot_4way()
  

  # Latex tables ####
  # GO + HLMK Tables ####
  BSBMD_vs_apical <- compare_bmd_EPA(hallmark_bigset_df,
                                     take_top_n = FALSE,
                                     apical_top = 1)
  BSBMD_vs_apical[BSBMD_vs_apical$Chemical == "EE",c(8,9)] <- 
    BSBMD_vs_apical[BSBMD_vs_apical$Chemical == "EE",c(8,9)]/1000
  
  ETAP_vs_apical <- compare_bmd_EPA(hallmark_bigset_df,
                                    take_top_n = FALSE,
                                    apical_top = 1,
                                    epa_vs_ap = T)
  ETAP_vs_apical[ETAP_vs_apical$Chemical == "EE",c(8,9)] <- 
    ETAP_vs_apical[ETAP_vs_apical$Chemical == "EE",c(8,9)]/1000
  
  BSBMD_vs_apical_GO <- compare_GO_apical(GO_lowest5_df,
                                          n_terms = 1)
  ETAP_vs_apical_GO <- compare_GO_apical(GO_lowest5_df,
                                         n_terms = 1,
                                         use_ETAP = T)
  
  
  # these are in the supplement
  colnames(ETAP_vs_apical_GO)[1] <- "Chemical"
  colnames(BSBMD_vs_apical_GO)[1] <- "Chemical"
  wide_comparison_BMDL <- BSBMD_vs_apical %>%
    select(-c(2, 3, 5:7, 8)) %>%
    full_join(ETAP_vs_apical[, c(1, 9)], by = "Chemical") %>%
    full_join(BSBMD_vs_apical_GO[, c(1, 9)], by = "Chemical") %>%
    full_join(ETAP_vs_apical_GO[, c(1, 9)], by = "Chemical") 


    
    exclude_chems = c(9, 16)
  bmdmat = sweep(wide_comparison_BMDL[-exclude_chems, 3:6],
        MARGIN = 1,
        STATS = wide_comparison_BMDL$BMDL[-exclude_chems],
        FUN = "/") %>%
    mutate_if(is.numeric, log) %>%
    dplyr::mutate_if(is.numeric, signif, 2)
  

  
  latex.tabular(as.tabular(
    rbind(cbind(wide_comparison_BMDL$Chemical[-exclude_chems], bmdmat),
          c("Mean" ,round(colMeans(bmdmat, na.rm = T), 2)))
    ))
  
  
  wide_comparison_BMDL <- wide_comparison_BMDL%>%
    dplyr::mutate_if(is.numeric, round, 2)
  
  wide_comparison_BMD <- BSBMD_vs_apical %>%
    select(-c(2, 4:7, 9)) %>%
    full_join(ETAP_vs_apical[, c(1, 8)], by = "Chemical") %>%
    full_join(BSBMD_vs_apical_GO[, c(1, 8)], by = "Chemical") %>%
    full_join(ETAP_vs_apical_GO[, c(1, 8)], by = "Chemical") %>%
    dplyr::mutate_if(is.numeric, round, 2)

  xtable::xtable(wide_comparison_BMDL, digits = 3, display = c("s", "s", rep("g", 5)))
  xtable::xtable(wide_comparison_BMD, digits = 3, display = c("s", "s", rep("g", 5)))
  
  latex.tabular(as.tabular(as.data.frame(wide_comparison_BMDL)))
  latex.tabular(as.tabular(as.data.frame(wide_comparison_BMD)))



  # Table 1, EPA_BMD comparison to apical and BSBMD
  EPA_BMD <- compare_bmd_EPA(hallmark_bigset_df, take_top_n = 5)
  sub_EPA <- EPA_BMD %>%
    filter(Chemical %in% c("Furan", "Fenofibrate", "Methyleugenol")) %>%
    filter(Tissue == "Liver") %>%
    select(!c(Tissue)) %>%
    dplyr::mutate_if(is.numeric, round, 2) %>%
    ungroup() %>%
    select(
      Chemical, Hallmark, BMD, BMDL, Rank,
      BMD_Median_epa, BMDL_Median_epa, Rank_EPA,
      Active_Genes)

  # for large tables, just save to excel
  # openxlsx::write.xlsx(sub_EPA, "BMD_EPA_ranked_comparison.xlsx")
  latex.tabular(as.tabular(as.data.frame(sub_EPA)))

  # get lowest EPA liver per 3 chems of interest
  compare_bmd_EPA(hallmark_bigset_df, take_top_n = 50) %>%
    filter(Chemical %in% c("Furan", "Fenofibrate", "Methyleugenol")) %>%
    filter(Tissue == "Liver") %>%
    group_by(Chemical) %>%
    slice_min(BMDL_Median_epa, n = 1, with_ties = FALSE) %>%
    dplyr::mutate_if(is.numeric, round, 2)

  # comparison of BSBMD and apical endpoints from EPA data
  EPA_apical_BMD <- compare_bmd_EPA(hallmark_bigset_df,
    take_top_n = FALSE,
    apical_top = 3)
  EPA_apical_BMD <- EPA_apical_BMD[, -5] %>%
    dplyr::mutate_if(is.numeric, round, 2)
  latex.tabular(as.tabular(as.data.frame(EPA_apical_BMD)))

  latex.tabular(as.tabular(as.data.frame(sub_comp_df)))
  top3_noBPAF_tric <- top3_BMDL %>%
    filter(!(Chemical %in% c("BPAF", "Triclosan")))
  latex.tabular(as.tabular(as.data.frame(top3_noBPAF_tric)))


    top3_BMDL <- top3_BMDL %>% mutate_if(is.numeric, round, 2)  
  latex.tabular(as.tabular(as.data.frame(top3_BMDL)))
  }





if(FALSE){
  # Old ETAP data filtering code
  # Old GO ####
  ETAP_apical_GO_df_pt1 <- read_delim(file = "../../scott_epa/ETAP GOBP v5 EPA Direct.txt") %>% # ,skip = 6
    filter(`Input Genes` > 2 & `Input Genes` < 501) %>%
    select(c(1, 2, 3, 4, 17, 18)) %>% # 23,24 for v3
    mutate(
      Chemical = unlist(lapply(Analysis, FUN = function(r) strsplit(r, "_")[[1]][2])),
      Tissue = unlist(lapply(Analysis, FUN = function(r) {
        tiss <- strsplit(r, "_")[[1]][4]
        if (tiss == "AA") tiss <- "liver" # TBBPA has "__" instead of "_"
        return(tiss)
      })),
      BMDL = `BMDL Median`,
      BMD = `BMD Median`
    ) %>%
    mutate(Chemical = unlist(lapply(Chemical,
                                    FUN = function(r) str_replace(r, "-", "")
    ))) %>%
    filter(Chemical %in% unique(GO_df$Chemical))
  
  
  ETAP_apical_go_miss <- "../../scott_epa/ETAP GO Final - Analysis in BMDE2 GO BP 3 Genes Extra chemicals not used in ETAP Combined Analyses_filtered.txt"
  ETAP_apical_GO_df_pt2 <- read_delim(ETAP_apical_go_miss, skip = 4) %>% # ,skip = 6
    filter(`Input Genes` > 2 & `Input Genes` < 501) %>%
    select(c(1, 2, 3, 4, 25, 31)) %>% # 23,24 for v3 #17, 18 for v4
    mutate(
      Chemical = unlist(lapply(Analysis, FUN = function(r) strsplit(r, "_")[[1]][1])),
      Tissue = unlist(lapply(Analysis, FUN = function(r) {
        tiss <- strsplit(r, "_")[[1]][2]
        if (tiss == "AA") tiss <- "liver" # TBBPA has "__" instead of "_"
        return(tiss)
      })),
      BMDL = `BMDL Median`,
      BMD = `BMD Median`
    ) %>%
    mutate(Chemical = unlist(lapply(Chemical,
                                    FUN = function(r) str_replace(r, "-", "")
    ))) %>%
    filter(Chemical %in% unique(GO_df$Chemical)) 
  
  ETAP_apical_GO_df_tot <- rbind(ETAP_apical_GO_df_pt1, ETAP_apical_GO_df_pt2)
  
  sub_comp_df <- ETAP_apical_GO_df_tot %>%
    group_by(Chemical) %>%
    slice_min(BMD, n = n_terms, with_ties = F) %>%
    select(7, 8, 2, 10, 9) %>%
    rename(
      ETAP.GO.BMD = BMD,
      ETAP.GO.BMDL = BMDL
    ) %>%
    mutate_if(is.numeric, round, 3)
  
  # Old HALLMARK ####
  epa_df_pt1 <- read_delim("../../scott_epa/EPA_hallmark.txt", delim = "\t")
  epa_hallmark_df_pt1 <- epa_df_pt1 %>%
    mutate(
      Chemical = unlist(lapply(
        Analysis,
        function(x) unlist(strsplit(x, split = "_"))[2]
      )),
      Tissue = unlist(lapply(
        Analysis,
        function(x) str_to_title(unlist(strsplit(x, split = "_"))[4])
      ))
    ) %>%
    transform(
      # BMD_epa = `BMD at 5th Percentile of Total Genes`,
      # BMD_Mean_epa = `BMD Mean`,
      BMD_Median_epa = `BMD Median`,
      # BMDL_Mean_epa = `BMDL Mean`,
      BMDL_Median_epa = `BMDL Median`,
      Hallmark = `GO/Pathway/Gene Set/Gene Name`,
      Active_Genes = paste0(
        `Genes That Passed All Filters`,
        "/", `All Genes (Expression Data)`
      )
    ) %>%
    dplyr::select(
      Chemical, Tissue, Hallmark, # BMD_Mean_epa, BMD_epa,BMDL_Mean_epa,
      BMD_Median_epa, BMDL_Median_epa, Active_Genes
    )
  
  
  epa_file = "../../scott_epa/ETAP Hallmark Final - Analysis in BMDE2 Hallmark 3 Genes Extra chemicals not used in ETAP Combined Analyses_filtered.txt"
  epa_df_pt2 <- read_delim(epa_file, skip = 4, delim = "\t")
  epa_hallmark_df_pt2 <- epa_df_pt2 %>%
    mutate(
      Chemical = unlist(lapply(
        Analysis,
        function(x) unlist(strsplit(x, split = "_"))[1]
      )),
      Tissue = unlist(lapply(
        Analysis,
        function(x) str_to_title(unlist(strsplit(x, split = "_"))[2])
      ))
    ) %>%
    transform(
      # BMD_epa = `BMD at 5th Percentile of Total Genes`,
      # BMD_Mean_epa = `BMD Mean`,
      BMD_Median_epa = `BMD Median`,
      # BMDL_Mean_epa = `BMDL Mean`,
      BMDL_Median_epa = `BMDL Median`,
      Hallmark = `GO/Pathway/Gene Set/Gene Name`,
      Active_Genes = paste0(
        `Genes That Passed All Filters`,
        "/", `All Genes (Expression Data)`
      )
    ) %>%
    dplyr::select(
      Chemical, Tissue, Hallmark, # BMD_Mean_epa, BMD_epa,BMDL_Mean_epa,
      BMD_Median_epa, BMDL_Median_epa, Active_Genes
    )
  
  epa_hallmark_df<- rbind(epa_hallmark_df_pt1, epa_hallmark_df_pt2)
  
  epa_hallmark_df$Chemical[grep(
    "thujone",
    epa_hallmark_df$Chemical
  )] <- "abthujone"
  
}