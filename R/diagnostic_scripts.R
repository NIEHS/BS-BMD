# This code requires the MCMC to have been run to populate env vars
library(tensor)
library(coda)
library(R.utils)


# diagnostic plot:  check if splines fit well
check_gene_set_fits <- function(hallmark_idx,
                                file_save_name = "chem.pdf",
                                plot_out = FALSE,
                                decorrelate = TRUE,
                                use_coexpress = FALSE) {
  if (use_coexpress) {
    hallmark_name <- present_genesets[hallmark_idx]
  } else {
    hallmark_name <- unique(gene_probe_map$MSigDBName)[hallmark_idx]
    hallmark_name <- unlist(str_split(hallmark_name, pattern = "HALLMARK_"))[2]
  }
  # organ_data - t(X%*%betas) - lambda%*%eta
  clust_idx <- clust_groups[[hallmark_idx]]
  gene_ids <- curr_data$Probe_ID[clust_idx]
  clust_size <- length(clust_idx)
  dose_seq <- seq(0, log(max(doses) + 1), 0.0525)
  R <- t(ju %*% betas[, clust_idx]) # allow for interpolator
  
  # individual decorrelated gene curves
  decor <- t(chol(mat[clust_idx, clust_idx]))
  R_decor <- t(ju %*% betas[, clust_idx] %*% solve(t(decor)))
  R_decor_cent <- R_decor - R_decor[, 1]
  gene_contrib <- R_decor_cent^2
  max_gene_contrib <- apply(gene_contrib, MARGIN = 1, max)
  gene_order <- order(max_gene_contrib, decreasing = T)
  
  file_save_name <- paste(hallmark_name, "_", file_save_name, sep = "")
  
  defpar <- par(mfrow = c(4, 3), mar = c(1.5, 2.5, .5, .5))
  organ_data_t <- organ_data - lambda %*% eta
  if (decorrelate) { # plot original data after removing latent factor
    for (i in gene_order) {
      plot(log(doses + 1), organ_data_t[clust_idx[i], ])
      lines(dose_seq, R[i, ], lty = 3, col = "gray")
      text(paste(gene_ids[i]), x = 2.5,
           y = mean(organ_data_t[clust_idx[i], ]), col = 2)
    }
  } else {
    for (i in gene_order) { # plot original data directly
      #+ (lambda %*% eta)[clust_idx[i],]
      # lambdas have different dimension that ju beta, cant combine
      plot(log(doses + 1), organ_data[clust_idx[i], ])
      # points(log(doses + 1), organ_data_t[clust_idx[i], ], col=2)
      lines(dose_seq, R[i, ], lty = 3, col = "gray")
      text(paste(gene_ids[i]), x = 2.5,
           y = mean(organ_data[clust_idx[i], ]), col = 4)
    }
  }
  par(defpar)
}




extract_avg_rank <- function(full_df,
                             BMD_MCMC_iter,
                             burnin = 0,
                             quantile_cut = 0.78) {
  effective_iter <- BMD_MCMC_iter - burnin
  full_df <- full_df[, c(1:3, (1 + burnin):BMD_MCMC_iter + 3)]
  significance <- apply(full_df[, 1:effective_iter + 3],
                        MARGIN = 1,
                        FUN = function(rx) {
                          ifelse(quantile(rx, .95) == max(rx),
                                 "Not Significant",
                                 "Significant"
                          )
                        }
  )
  full_df <- cbind(significance, full_df)
  plotdata <- full_df %>%
    mutate(Hallmark = sapply(Hallmark, function(x) unlist(str_split(x, pattern = "HALLMARK_"))[2])) %>%
    reshape2::melt(.,
                   id.vars = c("Chemical", "Tissue", "Hallmark", "significance")
    ) %>%
    rename(., c(BMD = value)) %>%
    mutate(BMD = as.numeric(BMD)) %>%
    dplyr::select(Chemical, Tissue, Hallmark, BMD, significance)
  chemicals <- unique(plotdata$Chemical)
  
  ranking_list <- list()
  for (i in seq_along(chemicals)) {
    mydata <- plotdata %>% dplyr::filter(Chemical == chemicals[i])
    # reorder the factor, not the data?
    mydata$Hallmark <- with(mydata, reorder(Hallmark, BMD, mean))
    ranking_list[[i]] <- levels(mydata$Hallmark)
  }
  
  unique_hlmk <- unlist(lapply(unique(gene_probe_map$MSigDBName),
                               function(x) unlist(str_split(x, pattern = "HALLMARK_"))[2]))
  avg_rank <- matrix(0, nrow = length(unique_hlmk), ncol = 1)
  for (hmi in seq_along(unique_hlmk)) {
    my_hlmk <- unique_hlmk[hmi]
    avg_rank[hmi] <- mean(unlist(lapply(ranking_list,
                                        FUN = function(id) {
                                          which(id == my_hlmk)
                                        }
    )))
  }
  hallmark_sizes <- unlist(lapply(clust_groups, length))
  rank_df <- data.frame(
    "Hallmark" = unique_hlmk,
    "Size" = hallmark_sizes,
    "Rank" = avg_rank
  )
  plot(rank_df$Size, rank_df$Rank,
       main = "Rank vs Size, Null Hallmark Data",
       xlab = "Gene Set Size",
       ylab = "Avg Rank"
  )
}



create_null_data <- function() {
  # full residuals
  # organ_data - t(X %*% betas) - lambda %*% eta
  # correlated residuals: just fit beta
  XtX <- t(X) %*% X
  for (i in 1:nrow(organ_data)) {
    ##### Beta Regression ####
    tV <- solve(XtX * taus[i] + h_phi[i] * pC)
    tM <- tV %*% ((t(X) * taus[i]) %*% t(t_organ_data[i, , drop = F]) +
                    h_phi[i] * pC %*% matrix(bm[i], ncol(X), 1))
    betas[, i] <- as.numeric(t(chol(tV)) %*% as.matrix(rnorm(nrow(betas))) + tM)
  }
  # optionally: do not refit beta, use from current iteration
  corr_error <- organ_data - t(X %*% betas)
  R <- t(ju %*% betas) # allow for interpolator
  response_intercept <- R[, 1]
  null_data <- corr_error + response_intercept
  null_df <- data.frame("Probe_ID" = curr_data$Probe_ID, null_data)
  dose_df <- data.frame(t(dose_vec))
  names(dose_df) <- names(null_df)
  null_df <- rbind(dose_df, null_df)
  null_src <- str_split(file_list[file_idx], "[.]")[[1]][1]
  save_name <- paste("Null_from_", null_src, ".txt", sep = "")
  write_tsv(null_df, file = save_name)
}





plot_trace_bmd <- function() {
  pdf("trace_plot_bmd_feno_oxidative_phos.pdf", width = 5, height = 4)
  # 884: oxidative phos, fenofib liver
  # TGF Beta 856
  plot(array(t(hallmark_bigset_df[850 + 33, seq(2000, 7000, by = 2)])),
       pch = 1, ylab = "BMD Estimate",
       xlab = "Iteration",
       main = "Sample Trace Plot, Fenofibrate Liver \n TGF Beta Signaling"
  )
  dev.off()
}




coda_test = function(){
  
  coda_results <- array(dim = 0)
  for (i_idx in 0) {
    for (j_idx in 0) {
      chain_list <- list()
      mcmc_records <- list.files("MCMC_chains/")
      
      for (i in 13:16) {
        load(paste0("MCMC_chains/", mcmc_records[i]))
        record_mat <- BMD
        chain_list <- c(chain_list, list(record_mat))
      }
      coda_test <- matrix(0, nrow = ncol(record_mat), ncol = 2)
      for (i in 1:ncol(record_mat)) {
        mcmc_list <- as.mcmc.list(lapply(chain_list,
                                         FUN = function(mat) {
                                           coda::mcmc(mat[, i])
                                         }
        ))
        coda_test[i, ] <- coda::gelman.diag(mcmc_list)$psrf
      }
      coda_results <- rbind(coda_results, coda_test)
    }
  }
  
  mcmc_list <- as.mcmc.list(lapply(chain_list,
                                   FUN = function(mat) coda::mcmc(mat)))
  coda_test <- coda::gelman.diag(mcmc_list)$psrf
  
  summary(coda_test)
  beepr::beep()
  plot(chain_list[[2]][, 22], col = 2)
  for (i in c(1, 3, 4)) points(chain_list[[i]][, 22], col = i)
  
}



get_reduced_covmat <- function(lambda_record, tau_record, i_idx, j_idx) {
  if (i_idx > 5 || j_idx > 5) {
    warning("index too large")
    return()
  }
  n_iters <- dim(lambda_record)[1]
  n_genes <- dim(lambda_record)[2]
  # can't store full cov record in memory, compute each part separately
  cov_mat_record <- array(0, dim = c(n_iters, 50, 50))
  for (i in 1:n_iters) {
    if (i %% 500 == 0) print(i)
    lam_sqr <- lambda_record[i, , ] %*% t(lambda_record[i, , ])
    cov_mat_record[i, , ] <-
      (lam_sqr + diag(1 / tau_record[i, ]))[
        1:50 + i_idx * 50,
        1:50 + j_idx * 50
      ]
  }
  return(wrap.array(cov_mat_record, map = list(1, c(3, 2))))
}

coda_covmat = function(){
  chain_list <- list()
  mcmc_records <- list.files("MCMC_chains/")
  
  for (i in 5:8) {
    load(paste0("MCMC_chains/", mcmc_records[i]))
    # lambda record and tau record are part of the mcmc_records memory
    record_mat <- get_reduced_covmat(lambda_record, tau_record, 0, 0)
    chain_list <- c(chain_list, list(record_mat))
  }
  coda_test <- matrix(0, nrow = ncol(record_mat), ncol = 2)
  for (i in 1:ncol(record_mat)) {
    mcmc_list <- as.mcmc.list(lapply(chain_list,
                                     FUN = function(mat) coda::mcmc(mat[, i])))
    coda_test[i, ] <- coda::gelman.diag(mcmc_list)$psrf
  }
  
  summary(coda_test)
}


