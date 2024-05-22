#' some comments:
#' - an extreme response can be positive or negative, so we take both tails
#' - squaring or abs is needed to avoid effects canceling out
#' - take a set of dimension (gene set size) n, we have a multivariate distribution
#' - the extreme response for a unimodal represents points outside of some sphere or boundary
#' - for the iid gaussian case, the boundary is a circle.
#' - For arbitrary case, can be smallest area that gives the desired mass (ie 0.89, etc)
#' - However, we use magnitude, so we can just find the quantile for that, accounting for two tails

utail_x <- seq(qnorm(.99), 4, by = 0.01)
ltail_x <- seq(-4, qnorm(.01), by = 0.01)
utail78_x <- seq(qnorm(.89), 4, by = 0.01)
ltail78_x <- seq(-4, qnorm(.11), by = 0.01)
gauss_df <- data.frame(x, "y" = dnorm(x), "lab" = "Gauss")
utail_df <- data.frame("x" = utail_x, "y" = dnorm(utail_x), "lab" = "Tail")
ltail_df <- data.frame("x" = ltail_x, "y" = dnorm(ltail_x), "lab" = "Tail")
utail78_df <- data.frame("x" = utail78_x, "y" = dnorm(utail78_x), "lab" = "Tail")
ltail78_df <- data.frame("x" = ltail78_x, "y" = dnorm(ltail78_x), "lab" = "Tail")
dses <- 10^(seq(-1, log(20, base = 10), by = 0.05))
error_terms <- rnorm(length(dses), sd = .5)
dose_response_df <- data.frame("x" = dses, "y" = 3 / (1 + (1.5 / dses)^2))
dose_response_noise_df <- data.frame(
  "x" = dses, "y" = 3 / (1 + (1.5 / dses)^2),
  "noisy" = 3 / (1 + (1.5 / dses)^2) + error_terms,
  "er" = error_terms
)
okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")




gp_noisy_exampe <- ggplot(dose_response_noise_df, aes(x = x, y = y)) +
  geom_line() +
  scale_x_log10() +
  annotation_logticks(sides = "b") +
  ylim(-1, 4) +
  geom_point(aes(x = x, y = noisy)) +
  xlab("Dose") +
  ylab("Noisy Response") +
  geom_segment(aes(
    x = x[6],
    y = y[6],
    xend = x[6],
    yend = noisy[6]
  ), data = dose_response_noise_df, col = "red")
gp_null_exampe <- ggplot(dose_response_noise_df, aes(x = x, y = er)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  scale_x_log10() +
  annotation_logticks(sides = "b") +
  ylim(-1.5, 1.5) +
  xlab("Dose") +
  ylab("Error") +
  geom_segment(aes(
    x = x[6],
    y = 0,
    xend = x[6],
    yend = er[6]
  ), data = dose_response_noise_df, col = "red")

polar_df <- data.frame(x = c(.1, .2, .3), y = c(
  quantile(test_stat_df$s, .78),
  quantile(test_stat_df$s, .98) - quantile(test_stat_df$s, .78), 2
), fill = c(NA, "green", "red"))
rand_stats <- rnorm(1000)
polar_df <- data.frame(x = c(0, 0, 0), y = c(
  quantile(rand_stats, .78),
  quantile(rand_stats, .98) - quantile(rand_stats, .78), 1
), fill = c(NA, "green", "red"))

gp2 <- ggplot(polar_df, aes(x = factor(1), y = y)) +
  geom_col(width = 1, fill = c(NA, "green", "red")) +
  coord_polar()


scatter_df <- data.frame(x = rnorm(10000), y = rnorm(10000))
full_rad <- seq(0, 3 * pi, by = .1)
circle_df <- data.frame(x = cos(full_rad), y = sin(full_rad))
width_green <- qchisq(.98, df = 2) - qchisq(.78, df = 2)
radius1 <- sqrt(qchisq(.78, df = 2))
radius2 <- sqrt(qchisq(.98, df = 2))


okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
pt_color <- rep("grey", nrow(scatter_df))
pt_mag <- sqrt(scatter_df$x^2 + scatter_df$y^2)
pt_color[which(pt_mag > radius1)] <- "#56B4E9"
pt_color[which(pt_mag > radius2)] <- "#E69F00"
gp2_alt <- ggplot(data = scatter_df, aes(x = x, y = y)) +
  geom_point(alpha = .75, col = pt_color) +
  geom_point(pch = ".") +
  geom_path(data = circle_df * (radius1), aes(x = x, y = y), linewidth = 2, alpha = 1, color = "#56B4E9") +
  geom_path(data = circle_df * (radius2), aes(x = x, y = y), linewidth = 2, alpha = 1, color = "#E69F00") +
  xlab("x1") +
  ylab("x2") +
  annotate("segment",
    x = 0, y = (radius1), xend = 0, yend = (radius2),
    arrow = arrow(type = "closed", length = unit(0.02, "npc")), color = "black"
  ) +
  theme_bw()


arrow1a <- annotate("segment",
  x = qnorm(.89), y = dnorm(qnorm(.89)), xend = qnorm(.99), yend = dnorm(qnorm(.89)),
  arrow = arrow(type = "closed", length = unit(0.02, "npc"))
)
arrow1b <- annotate("segment",
  xend = qnorm(.01), y = dnorm(qnorm(.11)), x = qnorm(.11), yend = dnorm(qnorm(.11)),
  arrow = arrow(type = "closed", length = unit(0.02, "npc"))
)


# pdf("explanation_sauron.pdf", width = 9, height = 4.25)
gridExtra::grid.arrange(gp1 + arrow1a + arrow1b,
  gp2_alt,
  nrow = 1,
  layout_matrix = matrix(c(1, 2), nrow = 1)
)

# dev.off()


# Crump95: explainer plot ####

gp1 <- ggplot(gauss_df, aes(x = x, y = y)) +
  geom_ribbon(data = utail78_df, aes(x = x, ymin = 0, ymax = y), fill = "#56B4E9") +
  geom_ribbon(data = ltail78_df, aes(x = x, ymin = 0, ymax = y), fill = "#56B4E9") +
  geom_ribbon(data = utail_df, aes(x = x, ymin = 0, ymax = y), fill = "#E69F00") +
  geom_ribbon(data = ltail_df, aes(x = x, ymin = 0, ymax = y), fill = "#E69F00") +
  ylab("p(x)") +
  theme_bw() +
  geom_line()


gp_crump1 <- ggplot(gauss_df, aes(x = x, y = y)) +
  geom_ribbon(data = utail_df, aes(x = x, ymin = 0, ymax = y), fill = "red") +
  ylab("p(x)") +
  geom_line()


gp_crump2 <- ggplot(gauss_df, aes(x = x, y = y)) +
  geom_ribbon(data = utail78_df, aes(x = x, ymin = 0, ymax = y), fill = "green") +
  geom_ribbon(data = utail_df, aes(x = x, ymin = 0, ymax = y), fill = "red") +
  ylab("p(x)") +
  geom_line()

gp_crump3 <- ggplot(dose_response_df, aes(x = x, y = y)) +
  geom_line() +
  scale_x_log10() +
  annotation_logticks(sides = "b") +
  ylim(-2, 4)


# pdf("crump_bits.pdf", width = 9, height = 9)
gridExtra::grid.arrange(gp_crump1,
  gp_crump2,
  gp_crump3,
  nrow = 3,
  layout_matrix = rbind(c(1, 2), c(3), c(3))
)
# dev.off()
