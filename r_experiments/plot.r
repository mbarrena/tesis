library(tibble) 

plot_nl_patched <- function (results_nl) 
{
    col_regime_1 <- "#21618C"
    col_regime_2 <- "#D68910"
    specs <- results_nl$specs
    if (specs$model_type == 0) {
        irf_s1_mean <- results_nl[[1]]
        irf_s1_low <- results_nl[[2]]
        irf_s1_up <- results_nl[[3]]
        irf_s2_mean <- results_nl[[4]]
        irf_s2_low <- results_nl[[5]]
        irf_s2_up <- results_nl[[6]]
        gg_s1 <- rep(list(NaN), specs$endog * specs$endog)
        gg_s2 <- rep(list(NaN), specs$endog * specs$endog)
        plot_num <- 1
        for (rr in 1:(specs$endog)) {
            for (ss in 1:(specs$endog)) {
                tbl_s1_mean <- as.matrix(t(irf_s1_mean[, 0:specs$hor+1, 
                  ss]))[, rr]
                tbl_s1_low <- as.matrix(t(irf_s1_low[, 0:specs$hor+1, 
                  ss]))[, rr]
                tbl_s1_up <- as.matrix(t(irf_s1_up[, 0:specs$hor+1, 
                  ss]))[, rr]
                tbl_s1 <- data.frame(x = 0:specs$hor, mean = tbl_s1_mean, 
                  low = tbl_s1_low, up = tbl_s1_up)
                tbl_s1_mean <- as.matrix(t(irf_s2_mean[, 0:specs$hor+1, 
                  ss]))[, rr]
                tbl_s2_low <- as.matrix(t(irf_s2_low[, 0:specs$hor+1, 
                  ss]))[, rr]
                tbl_s2_up <- as.matrix(t(irf_s2_up[, 0:specs$hor+1, 
                  ss]))[, rr]
                tbl_s2 <- data.frame(x = 0:specs$hor, mean = tbl_s1_mean, 
                  low = tbl_s2_low, up = tbl_s2_up)
                gg_s1[[plot_num]] <- ggplot() + geom_line(data = tbl_s1, 
                  aes(y = mean, x = x), col = col_regime_1) + 
                  geom_ribbon(data = tbl_s1, aes(x = x, ymin = low, 
                    ymax = up), col = "grey", fill = "grey", 
                    alpha = 0.3) + theme_classic() + ggtitle(paste(specs$column_names[ss], 
                  "on", specs$column_names[rr], sep = " ")) + 
                  xlab("") + ylab("") + theme(title = element_text(size = 6), 
                  plot.title = element_text(hjust = 0.5)) + scale_y_continuous(expand = c(0, 
                  0)) + scale_x_continuous(expand = c(0, 0), 
                  breaks = seq(0, specs$hor, 2)) + geom_hline(yintercept = 0, 
                  col = "black", linewidth = 0.25, linetype = "dashed")
                gg_s2[[plot_num]] <- ggplot() + geom_line(data = tbl_s2, 
                  aes(y = mean, x = x), col = col_regime_2) + 
                  geom_ribbon(data = tbl_s2, aes(x = x, ymin = low, 
                    ymax = up), col = "grey", fill = "grey", 
                    alpha = 0.3) + theme_classic() + ggtitle(paste(specs$column_names[ss], 
                  "on", specs$column_names[rr], sep = " ")) + 
                  xlab("") + ylab("") + theme(title = element_text(size = 6), 
                  plot.title = element_text(hjust = 0.5)) + scale_y_continuous(expand = c(0, 
                  0)) + scale_x_continuous(expand = c(0, 0), 
                  breaks = seq(0, specs$hor, 2)) + geom_hline(yintercept = 0, 
                  col = "black", linewidth = 0.25, linetype = "dashed")
                plot_num <- plot_num + 1
            }
        }
    }
    else if (specs$model_type == 1 | specs$model_type == 2) {
        gg_s1 <- rep(list(NaN), specs$endog)
        gg_s2 <- rep(list(NaN), specs$endog)
        plot_num <- 1
        for (rr in 1:(specs$endog)) {
            tbl_s1_mean <- results_nl$irf_s1_mean[rr, ]
            tbl_s1_low <- results_nl$irf_s1_low[rr, ]
            tbl_s1_up <- results_nl$irf_s1_up[rr, ]
            tbl_s1 <- data.frame(x = 0:specs$hor, mean = tbl_s1_mean, 
                low = tbl_s1_low, up = tbl_s1_up)
            tbl_s2_mean <- results_nl$irf_s2_mean[rr, ]
            tbl_s2_low <- results_nl$irf_s2_low[rr, ]
            tbl_s2_up <- results_nl$irf_s2_up[rr, ]
            tbl_s2 <- data.frame(x = 0:specs$hor, mean = tbl_s2_mean, 
                low = tbl_s2_low, up = tbl_s2_up)
            gg_s1[[rr]] <- ggplot() + geom_line(data = tbl_s1, 
                aes(y = mean, x = x), col = col_regime_1) + geom_ribbon(data = tbl_s1, 
                aes(x = x, ymin = low, ymax = up), col = "grey", 
                fill = "grey", alpha = 0.3) + theme_classic() + 
                ggtitle(paste("Shock", "on", specs$column_names[rr], 
                  sep = " ")) + xlab("") + ylab("") + theme(title = element_text(size = 6), 
                plot.title = element_text(hjust = 0.5)) + scale_y_continuous(expand = c(0, 
                0)) + scale_x_continuous(expand = c(0, 0), breaks = seq(0, 
                specs$hor, 2)) + geom_hline(yintercept = 0, col = "black", 
                linewidth = 0.25, linetype = "dashed")
            gg_s2[[rr]] <- ggplot() + geom_line(data = tbl_s2, 
                aes(y = mean, x = x), col = col_regime_2) + geom_ribbon(data = tbl_s2, 
                aes(x = x, ymin = low, ymax = up), col = "grey", 
                fill = "grey", alpha = 0.3) + theme_classic() + 
                ggtitle(paste("Shock", "on", specs$column_names[rr], 
                  sep = " ")) + xlab("") + ylab("") + theme(title = element_text(size = 6), 
                plot.title = element_text(hjust = 0.5)) + scale_y_continuous(expand = c(0, 
                0)) + scale_x_continuous(expand = c(0, 0), breaks = seq(0, 
                specs$hor, 2)) + geom_hline(yintercept = 0, col = "black", 
                linewidth = 0.25, linetype = "dashed")
        }
    }
    list(gg_s1 = gg_s1, gg_s2 = gg_s2)
}


plot_lin_patched <- function (results_lin) 
{
    irf_lin_mean <- results_lin[[1]]
    irf_lin_low <- results_lin[[2]]
    irf_lin_up <- results_lin[[3]]
    specs <- results_lin$specs
    if (specs$model_type == 0) {
        plot_num <- 1
        gg_lin <- rep(list(NaN), specs$endog * specs$endog)
        for (rr in 1:(specs$endog)) {
            for (ss in 1:(specs$endog)) {
                tbl_lin_mean <- as.matrix(t(irf_lin_mean[, 0:specs$hor+1, 
                  ss]))[, rr]
                tbl_lin_low <- as.matrix(t(irf_lin_low[, 0:specs$hor+1, 
                  ss]))[, rr]
                tbl_lin_up <- as.matrix(t(irf_lin_up[, 0:specs$hor+1, 
                  ss]))[, rr]
                tbl_lin <- tibble(x = 0:specs$hor, 
                  mean = tbl_lin_mean, low = tbl_lin_low, up = tbl_lin_up)
                gg_lin[[plot_num]] <- ggplot() + geom_line(data = tbl_lin, 
                  aes(y = mean, x = x)) + geom_ribbon(data = tbl_lin, 
                  aes(x = x, ymin = low, ymax = up), col = "grey", 
                  fill = "grey", alpha = 0.3) + theme_classic() + 
                  ggtitle(paste(specs$column_names[ss], "on", 
                    specs$column_names[rr], sep = " ")) + xlab("") + 
                  ylab("") + theme(title = element_text(size = 6), 
                  plot.title = element_text(hjust = 0.5)) + scale_y_continuous(expand = c(0, 
                  0)) + scale_x_continuous(expand = c(0, 0), 
                  breaks = seq(0, specs$hor, 2)) + geom_hline(yintercept = 0, 
                  col = "black", linewidth = 0.25, linetype = "dashed")
                plot_num <- plot_num + 1
            }
        }
    }
    else if (specs$model_type == 1 | specs$model_type == 2) {
        gg_lin <- rep(list(NaN), specs$endog)
        for (rr in 1:(specs$endog)) {
            tbl_lin_mean <- irf_lin_mean[rr, ]
            tbl_lin_low <- irf_lin_low[rr, ]
            tbl_lin_up <- irf_lin_up[rr, ]
            tbl_lin <- tibble(x = 0:specs$hor, mean = tbl_lin_mean, 
                low = tbl_lin_low, up = tbl_lin_up)
            gg_lin[[rr]] <- ggplot() + geom_line(data = tbl_lin, 
                aes(y = mean, x = x)) + geom_ribbon(data = tbl_lin, 
                aes(x = x, ymin = low, ymax = up), col = "grey", 
                fill = "grey", alpha = 0.3) + theme_classic() + 
                ggtitle(paste("Shock", "on", specs$column_names[rr], 
                  sep = " ")) + xlab("") + ylab("") + theme(title = element_text(size = 6), 
                plot.title = element_text(hjust = 0.5)) + scale_y_continuous(expand = c(0, 
                0)) + scale_x_continuous(expand = c(0, 0), breaks = seq(0, 
                specs$hor, 2)) + geom_hline(yintercept = 0, col = "black", 
                linewidth = 0.25, linetype = "dashed")
        }
    }
    return(gg_lin)
}
