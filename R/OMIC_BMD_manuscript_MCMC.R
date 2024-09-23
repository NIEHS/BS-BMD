# for command line runs: pass in settings
myargs <- commandArgs(trailingOnly = TRUE)
print(myargs)
ext_file_idx <- 18
quantile_cutoff <- .78
use_GO <- FALSE
save_params <- FALSE
simulation_study <- FALSE
study_type <- 2
use_full_genome_as_cluster <- FALSE
output_id <- date()
if (length(myargs) > 0) {
  ext_file_idx <- as.numeric(myargs[1])
  if (length(myargs) > 1) {
    quantile_cutoff <- as.numeric(myargs[2])
    if (length(myargs) > 2) {
      output_id <- myargs[3]
      if (length(myargs) > 3) {
        sim_sets <- myargs[4]
        if (sim_sets == "use_GO") use_GO <- TRUE
        if (sim_sets == "simulate") simulation_study <- TRUE
      }
    }
  }
}

print(ext_file_idx)
print(quantile_cutoff)
print(output_id)

library(splines)
library(Matrix)
library(invgamma)
library(readr)
library(lattice)
library(RColorBrewer)
library(pheatmap)
library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(data.table)
library(stringr)
library(Rcpp)
sourceCpp("src/code.cpp")



print("reading files")


# setup GO
probe2gene <- scan(
  file = "input/GO_data/probe2gene.gz",
  what = "list", sep = "\n"
)
probe2gene_long <- probe2gene %>%
  lapply(., FUN = function(lent) unlist(strsplit(lent, ";|\t"))) %>%
  lapply(., FUN = function(list_input) {
    data.frame(
      "Probe_ID" = cbind(rep(
        list_input[1],
        length(list_input) - 1
      )),
      "gene" = list_input[-c(1)]
    )
  }) %>%
  bind_rows()

GO2gene <- scan(
  file = "input/GO_data/genes2gos.gz",
  what = "list", sep = c("\n")
)

GO2probe_long <- GO2gene[-1] %>%
  lapply(., FUN = function(lent) unlist(strsplit(lent, ";|\t"))) %>%
  lapply(., FUN = function(list_input) {
    data.frame(
      "GO" = list_input[-c(1:2)], # the 2nd entry is unnecessary text?
      "gene" = cbind(rep(list_input[1], length(list_input) - 2))
    )
  }) %>%
  bind_rows() %>%
  left_join(probe2gene_long, by = "gene", relationship = "many-to-many") %>%
  distinct(GO, Probe_ID, .keep_all = TRUE) %>%
  arrange(GO)


input_dir <- "input/"
input_files_dir <- "input/null_data/" # organ_data
#  null data can be used instead, input/null_data/
output_dir <- "output/"
file_list <- list.files(path = input_files_dir)
# read hallmark gene sets and link gene to probes
probe_map <- readr::read_table(file = "input/Probe File_Rat S1500+.txt")
gene_set <- readr::read_table(file = "input/Rat_Genome_Whole_H_msigdb_v5.0.txt")
gene_probe_map <- dplyr::inner_join(gene_set, probe_map, "Entrez_Gene_ID")
BMD_MCMC_iter <- 7000
set.seed(BMD_MCMC_iter)
# list to store empirical test statistic from permutation tests
gene_cluster_statistic_table <- list()
print("Beginning outer loop over files")
for (file_idx in 6) { # c(17, 18, 19, 20,25, 26 )) {
  # Load Data ####
  print(paste("Started:", file_list[file_idx]))
  curr_file <- paste(input_files_dir, file_list[file_idx], sep = "")
  curr_data <- readr::read_table(curr_file, col_names = FALSE, skip = 2)
  names(curr_data)[1] <- "Probe_ID"
  clust_groups <- list()
  # for each gene set, find the indices of the corresponding probes
  # this is done for each file in case there are different probes
  for (hlmk in unique(gene_probe_map$MSigDBName)) {
    probset <- gene_probe_map$Probe_ID[gene_probe_map$MSigDBName == hlmk]
    probe_group <- which(unlist(lapply(curr_data$Probe_ID,
                                       FUN = function(x) (x %in% probset))))
    clust_groups <- c(clust_groups, list(probe_group))
  }

  if (use_GO) {
    # create a simple map between probe and data index
    probe_idx <- data.frame("Probe_ID" = curr_data$Probe_ID,
                            "indx" = 1:nrow(curr_data))
    # use the GO to probe dataframe to collect relevant indices into lists
    clust_groups_GO_df <- GO2probe_long %>%
      left_join(probe_idx, by = "Probe_ID") %>%
      na.omit() %>% # some probes not in the data
      group_by(GO) %>%
      summarise(
        set_indx = list(indx),
        n_genes = length(indx)
      ) %>%
      filter(n_genes > 19 & n_genes < 501)
    go_names <- clust_groups_GO_df$GO # may help for labeling
    clust_groups <- clust_groups_GO_df$set_indx
  }


  # create obs data matrix
  organ_data <- as.matrix(curr_data[, -1])

  # Simulation section ####
  # the fixed effect should be added before the centering
  if (simulation_study) {
    genes_to_add_signal <- c()
    # apply fixed effects (hill function?) to 50% of the genes in some sets
    set.seed(1298)
    p_add_signal <- 1
    if (study_type == 1) {
      # choose 3 gene sets:
      true_signal_genesets <- sample(1:n_grps, 3) # 28, 11, 37
      genes_to_add_signal <- unique(unlist(clust_groups[true_signal_genesets]))
      p_add_signal <- 0.5
    }
    if (study_type == 2) {
      # choose one gene from each gene set
      # sample from set diff
      for (gp in 1:length(clust_groups)) {
        sample_set <- setdiff(clust_groups[[gp]], unlist(clust_groups[-gp]))
        sampled_gene <- sample(sample_set, 1)
        if (length(sample_set) == 1) sampled_gene <- sample_set
        genes_to_add_signal <- c(genes_to_add_signal, sampled_gene)
      }
    }
    # for each gene, add signals with probability 50%
    genes_with_added_signal <- c()
    for (gs in genes_to_add_signal) {
      if (runif(1) < p_add_signal) {
        amax <- 2 # runif(1, 1, 3)
        ec50 <- 100 # runif(1, min = 0, max =  max(doses)/2)
        effect <- hill_function(a = amax, b = ec50, c = 1, conc = doses)
        organ_data[gs, ] <- organ_data[gs, ] + effect
        genes_with_added_signal <- c(genes_with_added_signal, gs)
      }
    }
  }

  if (use_full_genome_as_cluster) {
    clust_groups <- c(clust_groups, list(c(1:nrow(curr_data))))
  }
  n_grps <- length(clust_groups)

  # center data to match centered prior
  r_means <- rowMeans(organ_data)
  organ_data <- organ_data - matrix(r_means, nrow(organ_data), ncol(organ_data))
  # extract doses
  dose_data <- readr::read_table(curr_file, col_names = FALSE, n_max = 2)
  dose_vec <- (as.numeric(dose_data[2, ]))
  # doses included as unlabeled row, needs correction
  doses <- dose_vec[-which(is.na(dose_vec))]
  # null data has some dose format issue
  if (any(grep("null", input_files_dir))) {
    #  the N0 null data sets are missing a dose
    if (file_idx < 6) {
      doses <- c(doses, 20)
    }
  }
  if (length(doses) < ncol(organ_data)) {
    warning("Dose Missing, added maximum dose to list")
    doses <- c(doses, max(doses))
  }




  X <- splines::bs(log(doses + 1), df = 5, intercept = TRUE)
  interp_knots <- attributes(X)$knots
  bdd_knots <- attributes(X)$Boundary.knots
  # for interpolation
  ju <- splines::bs(seq(0, log(max(doses) + 1), 0.0525),
    knots = interp_knots,
    Boundary.knots = bdd_knots,
    intercept = TRUE
  )
  # penalty matrix for spline coefficients beta:  ie prior covariance is pC*phi
  pC <- diag(ncol(X))
  taus <- matrix(1, nrow = nrow(organ_data), ncol = 1)
  betas <- matrix(0, nrow = ncol(X), ncol = nrow(organ_data))
  h_betas <- matrix(0, nrow = ncol(X))
  #####################################
  # Set up latent Factor over responses
  n_l_factors <- 15
  eta <- matrix(
    rnorm(n_l_factors * ncol(organ_data)),
    n_l_factors,
    ncol(organ_data)
  )
  lambda <- matrix(0, nrow(organ_data), n_l_factors)
  #####################################
  # Sparse matrices for our Latent factor model
  J <- Matrix::Matrix(1, nrow = ncol(organ_data), 1)
  Big_D <- Matrix::Diagonal(nrow(organ_data))
  Big_D <- kronecker(J, Big_D)
  BMD <- matrix(0, BMD_MCMC_iter, ncol = n_grps)
  BMD_individ <- matrix(0, BMD_MCMC_iter, ncol = nrow(organ_data))
  ####################################
  # for checking stationary, record parameters
  beta_record <- array(0, dim = c(BMD_MCMC_iter, nrow(betas), ncol(betas)))
  tau_record <- matrix(0, nrow = BMD_MCMC_iter, ncol = length(taus))
  eta_record <- array(0, dim = c(BMD_MCMC_iter, nrow(eta), ncol(eta)))
  if (save_params) {
    lambda_record <- array(0, dim = c(
      BMD_MCMC_iter - 2000,
      nrow(lambda),
      ncol(lambda)
    ))
  }
  ####################################
  # local scale variance
  theta <- matrix(1, nrow(organ_data), n_l_factors)
  # global scale variance
  gammas_prior <- rgamma(n_l_factors, 2, 1)
  gam_prior <- cumprod(gammas_prior) * 0 + 1
  h_phi <- rep(1, nrow(organ_data))
  bm <- rep(0, nrow(organ_data))
  # vector for sampling the test statistic
  test_threshold <- rep(0, n_grps)
  ##########################################################################
  # Begin MCMC ################################################################
  for (nn in 1:BMD_MCMC_iter) {
    if (nn %% 10 == 0) {
      print(nn)
    }
    ##### Beta Regression #####################################
    t_organ_data <- organ_data - lambda %*% eta
    betas <- sample_beta(t_organ_data, X, taus, pC, matrix(h_phi))
    ##### Hierarchical Variance (simple at first)
    t_betas <- betas
    b <- (t(t_betas) - bm)
    # need to account for prior variance, pC instead if I: BKB
    b <- 0.5 * rowSums((b %*% pC) * b) + .1
    a <- 0.5 * rep(nrow(betas) - 1, length(bm)) + .1
    h_phi <- rgamma(length(b), a, b)

    ##### Latent Factor Analysis (simple at first)############################
    # first the residual
    t_organ_data <- organ_data - t(X %*% betas)
    #####  Lambda ####
    lambda <- compute_lambda(
      t_organ_data,
      eta,
      taus,
      lambda,
      matrix(gam_prior),
      theta
    )

    ##### Eta ####
    eta <- compute_eta(t_organ_data, lambda, taus)
    ##### Multiplicative Gamma Process ####
    mult_lambda <- lambda * lambda * theta
    for (kk in 1:n_l_factors) {
      temp_g <- gam_prior / gammas_prior[kk]
      mult_A <- matrix(temp_g, nrow(lambda), ncol(lambda), byrow = TRUE)
      b_t <- 0.5 * sum(mult_A[, kk:n_l_factors] * 
                         mult_lambda[, kk:n_l_factors]) + 1
      a_t <- 0.5 * length(mult_A[, kk:n_l_factors]) + 3.3
      gammas_prior[kk] <- rgamma(1, a_t, b_t)
      gam_prior <- cumprod(gammas_prior)
    }

    #### Local Scale update ####
    t_gam_mat <- matrix(rep(gam_prior, nrow(mult_lambda)),
      ncol = length(gam_prior), byrow = TRUE
    )
    theta <- matrix(
      rgamma(
        length(mult_lambda),
        1 / 2 + 1 / 2,
        1 / 2 + array(mult_lambda * t_gam_mat) / 2
      ),
      ncol = length(gam_prior)
    )
    ##### update residual ####
    t_organ_data <- organ_data - t(X %*% betas) - lambda %*% eta
    #####  Variance update (conjugate) #################################
    bs <- rowSums(0.5 * t_organ_data * t_organ_data) + 2
    as <- 0.5 * rep(ncol(t_organ_data), nrow(t_organ_data)) + 2
    taus <- matrix(rgamma(length(bs), as, bs))
    ##### BMD Start ##################################
    # Estimate the BMD
    BMD_vec <- rep(0, n_grps)
    sig_mat <- as.numeric(1 / taus) * diag(length(taus))
    tau_mat <- as.numeric(taus) * diag(length(taus))
    mat <- lambda %*% t(lambda) + sig_mat
    permutation_test_list <- list()
    # fit spline to null data
    betas_null <- sample_beta(t_organ_data, X, taus, pC, matrix(h_phi))
    for (uv in 1:length(clust_groups)) {
      clust_idx <- clust_groups[[uv]]
      clust_size <- length(clust_idx)
      B <- 1
      if (clust_size == nrow(curr_data)) { # Woodbury identity
        A_i <- tau_mat # Matrix(as.numeric(taus)*Diagonal(length(taus)))
        B <- solve(diag(n_l_factors) + t(lambda) %*% A_i %*% lambda)
        B <- A_i - A_i %*% lambda %*% B %*% t(lambda) %*% A_i
      } else {
        B <- solve(mat[clust_idx, clust_idx])
      }

      # Pred Response
      R <- t(ju %*% betas[, clust_idx]) # allow for interpolator
      R <- R[, ] - R[, 1]
      # Form my Test Statistic
      TS <- t(R) %*% B
      TS <- t(TS) * R
      G <- colSums(TS)
      # Null Response
      R_null <- t(ju %*% betas_null[, clust_idx]) # allow for interpolator
      R_null <- R_null[, ] - R_null[, 1]
      # compute null test statistic
      TS_null <- t(R_null) %*% tau_mat[clust_idx, clust_idx]
      TS_null <- t(TS_null) * R_null
      G_null <- colSums(TS_null)
      bubs <- splinefun(seq(0, log(max(doses) + 1), 0.0525), G - G_null)
      chi2_threshold <- qchisq(0.98, df = clust_size) - qchisq(0.78,
        df = clust_size
      )

      solve_fun <- function(d) bubs(d) - chi2_threshold
      # compute the BMD
      if (sign(solve_fun(0)) != sign(solve_fun(log(max(doses) + 1)))) {
        BMD[nn, uv] <- exp(uniroot(
          solve_fun,
          c(log(1), log(max(doses) + 1))
        )$root) - 1
      } else {
        # else use max dose + flag
        BMD[nn, uv] <- max(doses)
      }
    }
    # for plotting traces: record additional parameters
    if (save_params) {
      beta_record[nn, , ] <- betas
      tau_record[nn, ] <- taus
      eta_record[nn, , ] <- eta
      if (nn > 2000) lambda_record[nn - 2000, , ] <- lambda
    }
  }
  print(paste("Finished:", file_list[file_idx]))
  file_id <- strsplit(file_list[file_idx], split = "[.]")[[1]][1]
  print(paste("Writing BMD chain for", file_list[file_idx]))
  file_name_bmds <- paste(output_dir, file_id, "_", output_id, ".txt", sep = "")
  write(BMD, file = file_name_bmds, ncolumns = ncol(BMD))
  if (save_params) {
    file_name_params <- paste(output_dir, file_id, "_",
      output_id, "_BMDparams[run1].Rdata",
      sep = ""
    )
    if (file.exists(file_name_params)) {
      for (i in 2:20) {
        file_name_params <- paste(output_dir, file_id, "_",
          output_id,
          "_BMDparams[run", i, "].Rdata",
          sep = ""
        )
        if (!file.exists(file_name_params)) break
      }
    }
    print(paste("Parameters saved to", file_name_params))
    save(beta_record, tau_record, eta_record, lambda_record,
      file = file_name_params
    )
    # if checking BMD across runs for stationarity, save the following instead
    # save(BMD, file = file_name_params)
  }
}

#### End MCMC ###
