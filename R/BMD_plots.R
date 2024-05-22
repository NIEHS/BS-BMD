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
    apply(full_df[, 1:effective_iter + 3], MARGIN = 1,
          FUN = function(rx) ifelse(FALSE, "Not Significant", "Significant"))
  
  plotdata <- full_df %>%
    mutate(Hallmark = sapply(Hallmark, function(x) {
      unlist(str_split(x, pattern = "HALLMARK_"))[2]})) %>%
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
          scale = 3, rel_min_height = 0.01,
        ) +
        xlab(x_label) +
        ylab("Hallmark Pathway") +
        scale_y_discrete(limits = rev(levels(mydata$Hallmark))) +
        theme_bw() +
        scale_fill_viridis_c(
          name = "Log BMD", option = "C",
          direction = -1, trans = "log"
        ) +
        theme(legend.position = "None") +
        theme(axis.text.y = element_text(size = 8)) +
        scale_x_continuous(trans = "log10") +
        annotation_logticks(sides = "b") +
        stat_summary(
          fun = mean, geom = "point", shape = 23, stroke = .4,
          size = 1.2, color = "black", fill = "green", aes(colour = "meanLines")
        ) +
        stat_summary(
          fun = function(x) quantile(x, 0.05), geom = "point", shape = 24,
          size = 1.2, stroke = .4, color = "black", fill = "cyan",
          aes(colour = "meanLines")
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
                    BMD_MCMC_iter, "-", burnin, ".pdf", sep = "")
  print(filename)
  pdf(file = filename, width = 8, height = 7)
  for (i in seq_along(biplots_chems)) {
    print(biplots_chems[i])
  }
  dev.off()
}


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
    "Test" = test_curve)
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
                           sep = "")
  df <- reshape2::melt(gene_contrib[gene_ordering, ])
  names(df) <- c("Gene", "Dose", "Response")
  df$Dose <- exp(dose_seq[df$Dose]) - 1
  hallmark_id <- unique(gene_probe_map$MSigDBName)[hallmark_idx]
  hallmark_id <- strsplit(hallmark_id, split = "HALLMARK_")[[1]][2]
  # use a regex to get the chemical and tissue from the data file name
  chem_tiss <- sub(".*/([^/]+)\\..*", "\\1", curr_file) # thanks chatGPT
  file_name_plot <- paste("Gene_breakdown_", hallmark_id, "_",
                          chem_tiss, ".pdf",
                          sep = "")
  plot_title <- plt_title
  test_threshold <- qchisq(0.98, df = clust_size) -
    qchisq(0.78, df = clust_size)
  # the initial plot is not log-scaled and difficult to analyze
  my_gg <- ggplot(data = df, aes(x = Dose, y = Response, col = Gene)) +
    geom_line(data = hallmark_df,
              aes(x = Dose, y = Response, col = "Hallmark Set"),
              col = "black", linewidth = 1.5) +
    geom_line(position = "stack") +
    geom_hline(yintercept = test_threshold,
               linetype = 2, color = "black",
               linewidth = 1.25) +
    geom_line(data = hallmark_df, aes(x = Dose, y = Noise),
              col = "gray", linewidth = 1.25) +
    geom_line(data = hallmark_df, aes(x = Dose, y = Test),
              col = "red", linewidth = 1.25) +
    scale_color_discrete(
      breaks = top_genes,
      name = "Top Contributing Genes",
      labels = top_gene_labels) +
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
      "Test Cutoff" = "black")
    df$Dose <- df$Dose + 0.01
    hallmark_df$Dose <- hallmark_df$Dose + 0.01
    my_gg <- ggplot(data = df, aes(x = Dose, y = Response)) +
      geom_line(data = hallmark_df,
                aes(x = Dose, y = Response),
                col = "black", linewidth = 2) +
      geom_line(data = hallmark_df,
                aes(x = Dose, y = Response,col = "Aggregated Response"),
                linewidth = 1.5) +
      geom_line(position = "stack", aes(group = Gene), col = "gray") +
      geom_hline(aes(yintercept = test_threshold, col = "Test Cutoff"),
                 linetype = 2, linewidth = 1.25) +
      geom_line(data = hallmark_df, aes(x = Dose, y = Noise),
                col = "black", linewidth = 1.75) +
      geom_line(data = hallmark_df, aes(x = Dose, y = Noise,
                                        col = "Null Response"), linewidth = 1) +
      geom_line(data = hallmark_df, aes(x = Dose, y = Test),
                col = "black", linewidth = 1.75) +
      geom_line(data = hallmark_df,
                aes(x = Dose, y = Test, col = "Test Statistic"),
                linewidth = 1.25) +
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
                                colour = "black") +
        annotate("text",
                 label = paste("Gene", top_genes[i]),
                 y = y_vals[i], x = exp(dose_seq[1]) - .5)
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
                              xaxt = "n", yaxt = "n")
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

get_epa_bmd <- function(epa_input, chem, tissue, hallmark) {
  # read epa bmd
  epa_input <-
    read_delim("/Users/zilberds/Desktop/Wheeler/scott_epa/EPA_hallmark.txt",
                          delim = "\t"
  )
  
  
  # create filter
  my_regex <- paste("\b", chem, "\b.*\b", tissue, "\b", sep = "")
  acr_epa <- epa_input %>%
    filter(grepl(my_regex, Analysis, ignore.case = TRUE)) %>%
    filter(grepl(hallmark, `GO/Pathway/Gene Set/Gene Name`)) %>%
    select(`BMDU Mean`)
  
  if (length(acr_epa) > 0) {
    return(acr_epa)
  }
  return(NA)
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
                          FUN = function(x) quantile(x, 0.05))
  hallmark_BMDs <- apply(hallmark_bigset_df[, burnin:BMD_MCMC_iter],
                         MARGIN = 1,
                         FUN = function(x) median(x))
  significant_BMDs <- apply(hallmark_bigset_df[, burnin:BMD_MCMC_iter],
                            MARGIN = 1,
                            FUN = function(rx) {
                              ifelse(quantile(rx, .95) == max(rx),
                                     "Not Significant",
                                     "Significant")})

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
  epa_df <- read_delim("scott_epa/EPA_hallmark.txt", delim = "\t")
  epa_hallmark_df <- epa_df %>%
    mutate(
      Chemical = unlist(lapply(
        Analysis,
        function(x) unlist(strsplit(x, split = "_"))[1]
      )),
      Tissue = unlist(lapply(
        Analysis,
        function(x) str_to_title(unlist(strsplit(x, split = "_"))[4])
      ))
    ) %>%
    transform(
      BMD_epa = `BMD at 5th Percentile of Total Genes`,
      BMD_Mean_epa = `BMD Mean`,
      BMD_Median_epa = `BMD Median`,
      BMDL_Mean_epa = `BMDL Mean`,
      BMDL_Median_epa = `BMDL Median`,
      Hallmark = `GO/Pathway/Gene Set/Gene Name`,
      Active_Genes = paste0(`Genes That Passed All Filters`,
                            "/", `All Genes (Expression Data)`)
    ) %>%
    dplyr::select(
      Chemical, Tissue, Hallmark, BMD_epa, BMD_Mean_epa,
      BMDL_Mean_epa, BMD_Median_epa, BMDL_Median_epa, Active_Genes
    )
  
  epa_hallmark_df$Chemical[grep("thujone",
                                epa_hallmark_df$Chemical)] <- "abthujone"
  
  epa_hallmark_df <- epa_hallmark_df %>%
    group_by(Chemical, Tissue) %>%
    mutate(Rank_EPA = min_rank(BMDL_Median_epa))
  
  # join EPA and our data
  comp_df <- left_join(hallmark_df,
                       epa_hallmark_df,
                       by = c("Chemical", "Tissue", "Hallmark"))
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
    sub_comp_df <- comp_df %>%
      filter(significance == "Significant") %>%
      group_by(Chemical, Tissue) %>%
      slice_min(BMD, n = take_top_n) %>%
      select(Chemical, Tissue, Hallmark, BMD, BMDL,
             Rank, BMD_Median_epa, BMDL_Median_epa,
             Rank_EPA, Active_Genes)
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
                by = c("Chemical", "Tissue"), suffix = c(".1", ".2")) %>%
      select(Chemical, Tissue, Hallmark.1,
             BMD.1, BMDL.1, Hallmark.2, BMD.2, BMDL.2)
  }

  if (apical_top > 0) {
    chem_names <- unique(hallmark_df$Chemical)
    # lookup apical data
    apical_df <- readxl::read_xlsx("scott_epa/apical_data.xlsx") %>%
      filter(Table_id == 9) %>%
      select(!Table_id)
    # correct the abbreviations
    abbr_indx <- stringr::str_to_upper(chem_names) %>%
      stringr::str_sub(., end = 4)
    abbr_indx[1] <- "THU"
    abbr_indx[8] <- "EE2"
    map_chemID <-
      sapply(unique(apical_df$Chemical), FUN = function(x) grep(x, abbr_indx))
    chems_present <- which(sapply(map_chemID, length) > 0)
    # filter out missing chemicals
    apical_df <- apical_df %>% filter(Chemical %in% names(chems_present))
    apical_df$Chemical <- chem_names[unlist(map_chemID[apical_df$Chemical])]
    apical_df <- apical_df %>%
      group_by(Chemical) %>%
      slice_min(BMDL, n = min(apical_top, 3))
    sub_comp_df <- comp_df %>%
      # filter(significance == "Significant") %>%
      group_by(Chemical) %>%
      slice_min(BMD, n = min(apical_top, 3)) %>%
      transform(
        BMD.1 = BMD,
        BMDL.1 = BMDL,
        Chemical.1 = Chemical
      ) %>%
      select(Chemical.1, Tissue, Hallmark, BMD.1, BMDL.1) %>%
      dplyr::mutate_if(is.numeric, round, 2)
    if (epa_vs_ap) {
      sub_comp_df <- comp_df %>%
        # filter(significance == "Significant") %>%
        group_by(Chemical) %>%
        slice_min(BMD_Median_epa, n = min(apical_top, 3)) %>%
        transform(
          BMD.1 = BMD_Median_epa,
          BMDL.1 = BMDL_Median_epa,
          Chemical.1 = Chemical
        ) %>%
        select(Chemical.1, Tissue, Hallmark, BMD.1, BMDL.1) %>%
        dplyr::mutate_if(is.numeric, round, 2)
    }
    # cant do inner join: duplicates the apical entries
    # apical_comp_df <- apical_df %>% inner_join(sub_comp_df, by = "Chemical")
    n_apical <- table(apical_df$Chemical)
    matched_hallmark_idx <- 
      unlist(lapply(unique(apical_df$Chemical), FUN = function(x) {
        which(sub_comp_df$Chemical == x)[1:n_apical[x]]
      }))
    comp_apical_df <- cbind(apical_df, sub_comp_df[matched_hallmark_idx, ])
    # write_csv(comp_apical_df, "apical_bmd.csv")
    return(comp_apical_df)
  }
  return(comp_df)
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
    slice_min(BMD_epa, n = min(1)) %>%
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

# Calls: ####
# here is where most of the functions are called
if (FALSE) {
  # GG Boxplot ####
  require(stringr)
  file_list_out <- list.files(path = output_dir)
  # use the output_id string or manually select only relevant files
  hallmark_list <- file_list_out[grep("7k_correct", file_list_out)]
  hallmark_list_filenames <-
    sapply(hallmark_list, FUN = function(x) strsplit(x, split = "[.]")[[1]][1],
           USE.NAMES = FALSE)
  col_names <- unique(gene_probe_map$MSigDBName)
  # single data.frame:
  hallmark_bigset_df <- data.frame()
  for (idx in 1:length(hallmark_list_filenames)) {
    plt_title <- str_replace_all(hallmark_list_filenames[idx],
                                 pattern = "_", replacement = " ")
    print(plt_title)
    curr_file_dir <- paste(output_dir, hallmark_list[idx], sep = "")
    bmd_mat <- read.table(curr_file_dir)
    # original files need to be transposed:
    bmd_corrected <- t(matrix(array(t(bmd_mat)), nrow = 50,
                              ncol = BMD_MCMC_iter, byrow = TRUE))
    chem <- strsplit(plt_title, split = " ")[[1]][1]
    tissue <- strsplit(plt_title, split = " ")[[1]][2]
    # for repeats:
    # chem <- paste(strsplit(plt_title, split = " ")[[1]][1],
    #               strsplit(plt_title, split = " ")[[1]][2])
    # tissue <- strsplit(plt_title, split = " ")[[1]][4]

    mydf <- data.frame("Chemical" = chem,
                       "Tissue" = tissue,
                       "Hallmark" = col_names)
    widedf <- cbind(mydf, t(bmd_corrected))
    hallmark_bigset_df <- rbind(hallmark_bigset_df, widedf)

    rm_idx <- which(apply(bmd_corrected,
                          MARGIN = 1,
                          FUN = function(rx) all(rx == 0)))
    if (length(rm_idx) > 0) bmd_corrected <- bmd_corrected[-rm_idx, ]
    bmd_corrected <- bmd_corrected # [500:5000, ]
    # print(make_gg_boxplot(bmd_corrected, col_names, plt_title))
  }
  # trace plot test: suggests that variance is decreasing
  plot(bmd_corrected[, 12])

  # save(hallmark_bigset_df, file= "output/chi2_final_hallmark_BMDs.RData")
  # load("output/chi2_final_hallmark_BMDs.RData")
  sub_idx <- c(901:1000)
  make_fancy_ggbox(hallmark_bigset_df,
                   BMD_MCMC_iter = BMD_MCMC_iter,
                   burnin = 2000,
                   quantile_cut = 0.78,
                   sort = TRUE,
                   use_ridgelets = TRUE,
                   file_label = "chi2_final")

  # simple plot to show BMD vs hallmark size, in case of trend
  BMDs <- compute_BMD(hallmark_bigset_df, BMD_MCMC_iter, burnin = 2000)
  hallmark_sizes <- unlist(lapply(clust_groups, length))
  plot(rep(hallmark_sizes, each = 40), BMDs)
  cbind(hallmark_bigset_df[, 1:3], BMDs)

  # compute effective sample sizes
  MCMC_eff_ss <- function(rowvals) {
    1 / (2 + 2 * sum(pacf(rowvals, plot = FALSE)$acf))}
  ess_values <- apply(hallmark_bigset_df[, 2000:BMD_MCMC_iter + 3],
                      MARGIN = 1, MCMC_eff_ss)
  summary(ess_values)
  # 0.2356  0.4469  0.4866  0.4775  0.5119  0.6461      74
  # minimum:  1178 samples
  plot(1:BMD_MCMC_iter, hallmark_bigset_df[226, 1:BMD_MCMC_iter + 3])


  # Gene breakdown ####
  title_6_TGF <- "Stacked gene responses in fenofibrate treated liver for \n the Hallmark TGF-Beta Signaling pathway"
  title_33_Fatty <- "Stacked gene responses in fenofibrate treated liver for \n the Hallmark Fatty Acid Metabolism pathway"

  set_bmd_breakdown(33,
                    saveplot_pdf = TRUE, log = TRUE, zero_center = TRUE,
                    plt_title = title_33_Fatty)

  # CovMat plot ####
  plot_corrmat(lambda, taus)

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
    device = "pdf")

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
  q_new <- get_corr_ordering(Lambda_3ds[[3]],
                             tau_3ds[[3]],
                             clust_groups[[clust_to_plot]])
  # q_new = unlist(clust_groups)
  # hallmark gene set only: clust_groups[[33]]
  gplot_list <- lapply(1:length(tau_3ds),
                       FUN = function(idx) {
                         plot_corrmat(Lambda_3ds[[idx]][q_new, ],
                                      tau_3ds[[idx]][q_new],
                                      thinning = 1,
                                      q_order = 1:length(q_new))})
  cov_pair <- grid.arrange(gplot_list[[1]],
                           gplot_list[[3]],
                           gplot_list[[5]],
                           gplot_list[[2]],
                           gplot_list[[4]],
                           gplot_list[[6]],
                           ncol = 3,
                           nrow = 2)

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
              sd = sd(rqq_mdf$value)),
        col = 3)
  ggsave("BMD_gene_express_qq.pdf", width = 7, height = 5)
  qq_df <- data.frame(organ_data - t(X %*% betas))
  ggplot(melt(qq_df), aes(sample = value)) +
    stat_qq() +
    geom_abline(slope = 1, intercept = 0) +
    ggtitle("QQ Norm, Data - XB")

  # Latex tables ####
  EPA_BMD <- compare_bmd_EPA(hallmark_bigset_df, take_top_n = 5)
  sub_EPA <- EPA_BMD %>%
    # filter(significance=="Significant")%>%
    filter(Chemical %in% c("Furan", "Fenofibrate", "Methyleugenol")) %>%
    filter(Tissue == "Liver") %>%
    # select(!c(Tissue))%>%
    dplyr::mutate_if(is.numeric, round, 2) %>%
    ungroup() %>%
    select(Chemical, Hallmark, BMD, BMDL, Rank,
           BMD_Median_epa, BMDL_Median_epa, Rank_EPA,
           Active_Genes)

  latex.tabular(as.tabular(as.data.frame(sub_EPA)))

  is.na(EPA_BMD$BMD_epa) & (EPA_BMD$significance == "Not Significant")

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

  # 3-way Scatter ####
  # create scatter plot for apical and epa comparison
  BMD_MCMC_iter <- 7000
  EPA_BMD <- compare_bmd_EPA(hallmark_bigset_df, take_top_n = 1)
  EPA_apical_BMD <-
    compare_bmd_EPA(hallmark_bigset_df, take_top_n = FALSE, apical_top = 1)
  EPA_apical_BMD <- EPA_apical_BMD[, -5] %>%
    dplyr::mutate_if(is.numeric, round, 2)
  EPA_apical_BMD_compl <-
    data.frame(EPA_apical_BMD[(EPA_apical_BMD$BMDL != "NA"), ])

  EPA_apical <- compare_bmd_EPA(hallmark_bigset_df,
                                take_top_n = FALSE,
                                apical_top = 1, epa_vs_ap = TRUE)
  EPA_apical <- EPA_apical[, -5] %>%
    dplyr::mutate_if(is.numeric, round, 2)
  EPA_apical_compl <- data.frame(EPA_apical[(EPA_apical$BMDL != "NA"), ])

  compare_df_BMD <- data.frame(
    "x" = c(
      EPA_BMD$BMD,
      EPA_apical_BMD_compl[, 7],
      EPA_apical_compl[, 7]
    ),
    "y" = c(
      EPA_BMD$BMD_Median_epa,
      as.numeric(EPA_apical_BMD_compl[, 3]),
      as.numeric(EPA_apical_compl[, 3])
    ),
    "Comparison" = c(
      rep("BS-BMD vs ETAP", nrow(EPA_BMD)),
      rep("BS-BMD vs Apical", nrow(EPA_apical_BMD_compl)),
      rep("ETAP vs Apical", nrow(EPA_apical_compl))
    )
  )
  compare_df <- data.frame(
    "x" = c(
      EPA_BMD$BMDL,
      EPA_apical_BMD_compl[, 8],
      EPA_apical_compl[, 8]
    ),
    "y" = c(
      EPA_BMD$BMDL_Median_epa,
      as.numeric(EPA_apical_BMD_compl[, 4]),
      as.numeric(EPA_apical_compl[, 4])
    ),
    "Comparison" = c(
      rep("BS-BMD vs ETAP", nrow(EPA_BMD)),
      rep("BS-BMD vs Apical", nrow(EPA_apical_BMD_compl)),
      rep("ETAP vs Apical", nrow(EPA_apical_compl))
    )
  )
  BMD_ggplot <- ggplot(compare_df_BMD,
                       aes(x = x, y = y,
                           color = Comparison,
                           shape = Comparison)) +
    geom_point(size = 3) +
    scale_x_log10() +
    scale_y_log10() +
    geom_abline(slope = 1, intercept = 0) +
    xlab("Method 1 BMD (mg/kg)") +
    ylab("Method 2 BMD (mg/kg)")
  BMDL_ggplot <- ggplot(compare_df, aes(x = x, y = y,
                                        color = Comparison,
                                        shape = Comparison)) +
    geom_point(size = 3) +
    scale_x_log10() +
    scale_y_log10() +
    geom_abline(slope = 1, intercept = 0) +
    xlab("Method 1 BMDL (mg/kg)") +
    ylab("Method 2 BMDL (mg/kg)")
  ggsave("Compare_BMDL_points_top1_3way.pdf",
         device = "pdf", width = 7, height = 5)
}
