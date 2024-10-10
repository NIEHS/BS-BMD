# standalone code to produce the posterior BMD ridge plots
library(dplyr)
library(stringr)
library(ggplot2)
library(ggridges)

make_fancy_ggbox <- function(full_df, BMD_MCMC_iter, burnin = 0,
                             quantile_cut = 0.78, sort = T,
                             use_ridgelets = F,
                             file_label = "biplots") {
  effective_iter <- BMD_MCMC_iter - burnin
  full_df <- full_df[, c(1:3, (1 + burnin):BMD_MCMC_iter + 3)]
  significance <- apply(full_df[, 1:effective_iter + 3],
                        MARGIN = 1,
                        #FUN = function(rx) ifelse(quantile(rx, .95) == max(rx), "Not Significant", "Significant")
                        FUN = function(rx) ifelse(FALSE, "Not Significant", "Significant")
  )
  
  #full_df <- cbind(significance, full_df)
  plotdata <- full_df %>%
    mutate(Hallmark = sapply(Hallmark, function(x) unlist(str_split(x, pattern = "HALLMARK_"))[2])) %>%
    reshape2::melt(.,
                   id.vars = c("Chemical", "Tissue", "Hallmark")#, "significance")
    ) %>%
    rename(
      .,
      c(BMD = value)
    ) %>%
    mutate(BMD = as.numeric(BMD)) %>%
    dplyr::select(Chemical, Tissue, Hallmark, BMD)#, significance)
  chemicals <- unique(plotdata$Chemical)
  
  jsh_hallmark_biplot <- function(mychem = chemicals[1]) {
    mydata <- plotdata %>% dplyr::filter(Chemical == mychem)
    plt_title <- mychem#paste0(mychem,  " - N=", BMD_MCMC_iter)
    if(length(grep("thujone",mychem))) plt_title <-
      expression(paste(alpha,beta, " Thujone"))
    if (nrow(mydata) > 5) {
      if (sort) {
        sort_function <- function(x) quantile(x, 0.05) +
          1e4 * (quantile(x, .90) == max(x))
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
          #aes(color = significance)
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

    plt_title <- mychem#paste0(mychem,  " - N=", BMD_MCMC_iter)

    if(length(grep("thujone",mychem))) plt_title <-
      expression(paste(alpha,beta, " Thujone"))
    
    if (nrow(mydata) > 5) {
      if (sort) {
        sort_function <- function(x) quantile(x, 0.05) +
          1e4 * (quantile(x, .90) == max(x))
        mydata$Hallmark <- with(mydata, reorder(Hallmark, BMD, mean))
      }
      ggplot(data = mydata, aes(y = Hallmark, x = BMD, fill = stat(x))) +
        geom_density_ridges_gradient(scale = 3,
                                     rel_min_height = 0.01)+
        xlab("Pathway BMD - mg/kg") +
        ylab("Hallmark Pathway") +
        #ggtitle(plt_title) +
        scale_y_discrete(limits = rev(levels(mydata$Hallmark))) +
        theme_bw() +
        scale_fill_viridis_c(name = "Log BMD", option = "C",
                             direction = -1, trans = "log")+
        theme(legend.position = "None") +
        theme(axis.text.y = element_text(size = 8)) +
        scale_x_continuous(trans = "log10") +
        annotation_logticks(sides = "b") +
        stat_summary(
          fun = function(x) log(mean(10^(x)), base=10), geom = "point", shape = 23,stroke=.4,
          size = 1.2, color = "black", fill = "green", aes(colour = "meanLines")
        ) +
        stat_summary(
          fun = function(x) quantile(x, 0.05), geom = "point", shape = 24,
          size = 1.2, stroke=.4, color = "black",
          fill = "cyan", aes(colour = "meanLines")
        ) +
        facet_wrap(Chemical ~ Tissue) -> biplot_ridge
    }
    
    return(biplot_ridge)
  }
  
  
  biplots_chems <- list()
  for (i in seq_along(chemicals)) {
    if(use_ridgelets){
      biplots_chems[[i]] <- try(hallmark_ridgeplot(mychem = chemicals[i]))
    }else{
      biplots_chems[[i]] <- try(jsh_hallmark_biplot(mychem = chemicals[i]))
    }
  }
  
  filename <- paste("output/", file_label, "_",
                    quantile_cut, "qtile_",
                    BMD_MCMC_iter, "-", burnin, ".pdf", sep = "")
  print(filename)
  pdf(file = filename, width = 8, height = 7)
  for (i in seq_along(biplots_chems)) {
    print(biplots_chems[i])
  }
  dev.off()
}



load("output/final_BMD_posteriors_fenofibrate_furan.Rdata")
fenofibrate_data = fenofib_furan_bmd_post[1:100,]
furan_data = fenofib_furan_bmd_post[101:200,]
make_fancy_ggbox(furan_data,
                 BMD_MCMC_iter = 7000,
                 burnin = 2000,
                 quantile_cut = 0.78,
                 sort = T,
                 use_ridgelets = T,
                 file_label = "chi2_final_furan"
)
