combine_results_3 <- function(sim3_results_all,
                              sim3_results_reps,
                              sim3_results_all_wt,
                              sim3_results_reps_wt,
                              sim3_results_generalized_sim,
                              sim3_results_generalized_sim_reps,
                              sim3_results_generalized_sim_combine_first,
                              sim3_results_generalized_sim_combine_first_reps,
                              sim3_results_simex,
                              sim3_results_empirical_simex) {
    
    results_list <- list(sim3_results_all,
                         sim3_results_reps,
                         sim3_results_all_wt,
                         sim3_results_reps_wt,
                         sim3_results_generalized_sim,
                         sim3_results_generalized_sim_reps,
                         sim3_results_generalized_sim_combine_first,
                         sim3_results_generalized_sim_combine_first_reps,
                         sim3_results_simex,
                         sim3_results_empirical_simex)

    do.call(rbind, lapply(results_list, function(results){
        results %>% pivot_longer(cols = -c("scenario", "tar_batch", "tar_rep"))
    }))
}

# 3 BP
# Plotting Helpers
# Generate the plot for each of the different variables
# This will be passed data either with scenarios that will (eventually) be
# faceted, *or* it will be filtered down to the scenario to begin.
X_maker_3 <- function(data, true_vals=NA, nmethod=1){ 
    data %>% ggplot(aes(y = value), outlier.size=NA) +
            geom_boxplot(aes(x = Method)) + 
            geom_segment(aes(x = 0, y = true_vals["X"],  xend = (nmethod)+.5, yend = true_vals["X"]), lty=3) + 
            labs(x="", y = expression(beta ~ " estimates"), title=expression(X)) + 
            theme_minimal() +
            scale_y_continuous(guide = guide_axis(check.overlap = TRUE))+
            theme(axis.text.y = element_text(size=10), 
                    axis.text.x = element_text(size=10),
                    panel.spacing.x = unit(4, "mm"))
}
int_maker_3 <- function(data, true_vals=NA, nmethod=1){
    data %>% ggplot(aes(y = value), outlier.size=NA) +
            geom_boxplot(aes(x = Method)) + 
            geom_segment(aes(x = 0, y = true_vals["Intercept"],  xend = (nmethod)+.5, yend = true_vals["Intercept"]), lty=3) + 
            labs(x="", y = expression(beta ~ " estimates"), title="Intercept") + 
            theme_minimal() +
            scale_y_continuous(guide = guide_axis(check.overlap = TRUE))+
            theme(axis.text.y = element_text(size=10), 
                    axis.text.x = element_text(size=10),
                    panel.spacing.x = unit(4, "mm"))
}
## Produce a Boxplot from the Results of SIMULATION 1
generate_plot_3_boxplot <- function(sim3) {
    # Specify the values here
    true_vals <- c(Intercept = 0.5, X = -0.5)
    file_start<- "sim3_new_boxplot"
    folder <- "media"
    combined_plot <- paste0(folder, "/", file_start, ".png")

    # Actual Plot Code
    sim3 <- sim3 %>% 
                separate(name, c("Variable", "Method"), sep="\\s", extra="merge") %>% 
                mutate(Method = str_sub(Method, 2, -2))

    nmethod <- length(unique(sim3$Method))

    X <- sim3 %>% filter(Variable == "X") %>% arrange(scenario) %>% X_maker_3(true_vals=true_vals, nmethod=nmethod) + facet_wrap(~scenario,  nrow = 1, scales="free_x") + coord_flip()
    int <- sim3 %>% filter(Variable == "Intercept") %>% arrange(scenario) %>% int_maker_3(true_vals=true_vals, nmethod=nmethod) + facet_wrap(~scenario,  nrow = 1, scales="free_x") + coord_flip()

    
    png(filename=combined_plot, width=14, height=9, units="in", res=300)
    grid.arrange(arrangeGrob(int, X, nrow=2, left = textGrob("Estimation technique", rot = 90, vjust = 1)))
    dev.off()

    plots <- c(combined_plot)

    for (sc in unique(sim3$scenario)) {
        X <- sim3 %>% filter(Variable == "X" && scenario == sc) %>% X_maker_3(true_vals=true_vals, nmethod=nmethod) + coord_flip()
        int <- sim3 %>% filter(Variable == "Intercept" && scenario == sc) %>% int_maker_3(true_vals=true_vals, nmethod=nmethod) + coord_flip()

        pname <- paste0(folder, "/", file_start, "_", sc, ".png")
        png(filename=pname, width=14, height=9, units="in", res=300)
        grid.arrange(arrangeGrob(int, X, nrow=2, left = textGrob("Estimation technique", rot = 90, vjust = 1)))
        dev.off()
        plots <- c(plots, pname)
    }

    plots
}


generate_plot_3_probplot <- function(sim3) {
    test_x <- seq(-2, 8, length.out=100)
    true_probs <- expit(0.5 - 0.5*test_x)

    est_probs_upper <- list()
    est_probs_lower <- list()
    which_q <- c(0.025, 0.975)

    sim3 <- sim3 %>% 
                separate(name, c("Variable", "Method"), sep="\\s", extra="merge") %>% 
                mutate(Method = str_sub(Method, 2, -2))
    
    for(sc in unique(sim3$scenario)) {
        est_probs_lower[[sc]] <- data.frame(X = test_x)
        est_probs_upper[[sc]] <- data.frame(X = test_x)

        for (me in unique(sim3$Method)) {
            est_probs_lower[[sc]][me] <- 0 
            est_probs_upper[[sc]][me] <- 0 
            
            filtered_df <- sim3 %>% 
                            filter(scenario == sc && Method == me) %>%
                            pivot_wider(names_from = Variable, values_from = value, id_cols = c(tar_batch, tar_rep, Method,scenario))
            
            all_results  <- vapply(X = test_x, FUN=function(x) expit(filtered_df$Intercept + filtered_df$X*x), FUN.VALUE=numeric(nrow(filtered_df)))
            qs <- apply(all_results, FUN=quantile, MARGIN=2, na.rm=TRUE, probs = which_q)
            est_probs_lower[[sc]][me] <- qs[1,]
            est_probs_upper[[sc]][me] <- qs[2,]
        }
    }

    plots <- c()
    for(idx in 1:length(est_probs_lower)) {
        long_lower <- est_probs_lower[[idx]] %>% pivot_longer(cols=-1)
        long_upper <- est_probs_upper[[idx]] %>% pivot_longer(cols=-1)

        
        intervals <- long_lower %>% 
            inner_join(long_upper, by=c(X="X",name="name")) %>% 
            left_join(data.frame(X=test_x, true=true_probs), by=c(X="X"))

        p <- intervals %>% 
            ggplot(aes(x = X, group = name)) + 
            geom_rect(aes(xmin=3-qnorm(0.975), xmax=3+qnorm(0.975), ymin=0, ymax=Inf), fill = "grey70", alpha=0.01) +
            geom_line(aes(y = value.x), lty=2) + 
            geom_line(aes(y = value.y), lty=2) + 
            geom_line(aes(y = true), lwd=1.5) + 
            facet_wrap(~name, nrow=2) + 
            labs(y='Estimated Probability', x='X Value') + 
            scale_y_continuous(expand = c(0, 0)) + 
            theme_minimal() +
            theme(axis.text.y = element_text(size=10), axis.text.x = element_text(size=10))

        pname <- paste0("media/sim3_probs", names(est_probs_lower)[idx], ".png")
        ggsave(pname, p, dpi=300, width=15, height=7)
        plots <- c(plots, pname)
    }

    plots
}
