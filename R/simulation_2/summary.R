combine_results_2 <- function(sim2_results_truth,
                              sim2_results_RC,
                              sim2_results_generalized_RC,
                              sim2_results_generalized_RC_noZ,
                              sim2_results_generalized_RC_wtd,
                              sim2_results_generalized_RC_wtd_noZ,
                              sim2_results_standard_simex,
                              sim2_results_empirical_simex,
                              sim2_results_generalized_SIMEX,
                              sim2_results_generalized_SIMEX_combine_first) {
    
    results_list <- list(sim2_results_truth,
                         sim2_results_RC,
                         sim2_results_generalized_RC,
                         sim2_results_generalized_RC_noZ,
                         sim2_results_generalized_RC_wtd,
                         sim2_results_generalized_RC_wtd_noZ,
                         sim2_results_standard_simex,
                         sim2_results_empirical_simex,
                         sim2_results_generalized_SIMEX,
                         sim2_results_generalized_SIMEX_combine_first)

    do.call(rbind, lapply(results_list, function(results){
        results %>% pivot_longer(cols = -c("scenario", "tar_batch", "tar_rep"))
    }))
}



# Plotting Helpers
# Generate the plot for each of the different variables
# This will be passed data either with scenarios that will (eventually) be
# faceted, *or* it will be filtered down to the scenario to begin.
X_maker_2 <- function(data, true_vals=NA, nmethod=1){ 
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
Z_maker_2 <- function(data, true_vals=NA, nmethod=1){
    data %>% ggplot(aes(y = value), outlier.size=NA) +
            geom_boxplot(aes(x = Method)) + 
            geom_segment(aes(x = 0, y = true_vals["Z"],  xend = (nmethod)+.5, yend = true_vals["Z"]), lty=3) + 
            labs(x="", y = expression(beta ~ " estimates"), title=expression(Z)) + 
            theme_minimal() +
            scale_y_continuous(guide = guide_axis(check.overlap = TRUE))+
            theme(axis.text.y = element_text(size=10), 
                    axis.text.x = element_text(size=10),
                    panel.spacing.x = unit(4, "mm"))
}
int_maker_2 <- function(data, true_vals=NA, nmethod=1){
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
## Produce a Boxplot from the Results of SIMULATION 2
generate_plot_2 <- function(sim2) {
    # Specify the values here
    true_vals <- c(Intercept = 2, X = 2, Z = -3)
    file_start<- "sim2_new_boxplot"
    folder <- "media"
    combined_plot <- paste0(folder, "/", file_start, ".png")

    # Actual Plot Code
    sim2 <- sim2 %>% 
                separate(name, c("Variable", "Method"), sep="\\s", extra="merge") %>% 
                mutate(Method = str_sub(Method, 2, -2))

    nmethod <- length(unique(sim2$Method))

    Z <- sim2 %>% filter(Variable == "Z") %>% arrange(scenario) %>% Z_maker_2(true_vals=true_vals, nmethod=nmethod) + facet_wrap(~scenario,  nrow = 1, scales="free_x") + coord_flip()
    X <- sim2 %>% filter(Variable == "X") %>% arrange(scenario) %>% X_maker_2(true_vals=true_vals, nmethod=nmethod) + facet_wrap(~scenario,  nrow = 1, scales="free_x") + coord_flip()
    int <- sim2 %>% filter(Variable == "Intercept") %>% arrange(scenario) %>% int_maker_2(true_vals=true_vals, nmethod=nmethod) + facet_wrap(~scenario,  nrow = 1, scales="free_x") + coord_flip()

    
    png(filename=combined_plot, width=14, height=9, units="in", res=300)
    grid.arrange(arrangeGrob(int, X, Z, nrow=3, left = textGrob("Estimation technique", rot = 90, vjust = 1)))
    dev.off()

    plots <- c(combined_plot)

    for (sc in unique(sim2$scenario)) {
        Z <- sim2 %>% filter(Variable == "Z" && scenario == sc) %>% Z_maker_2(true_vals=true_vals, nmethod=nmethod) + coord_flip()
        X <- sim2 %>% filter(Variable == "X" && scenario == sc) %>% X_maker_2(true_vals=true_vals, nmethod=nmethod) + coord_flip()
        int <- sim2 %>% filter(Variable == "Intercept" && scenario == sc) %>% int_maker_2(true_vals=true_vals, nmethod=nmethod) + coord_flip()

        pname <- paste0(folder, "/", file_start, "_", sc, ".png")
        png(filename=pname, width=14, height=9, units="in", res=300)
        grid.arrange(arrangeGrob(int, X, Z, nrow=3, left = textGrob("Estimation technique", rot = 90, vjust = 1)))
        dev.off()
        plots <- c(plots, pname)
    }

    plots
}


generate_table_2 <- function(sim2) {
    # Specify the values here
    true_vals <- c(Intercept = 2, X = 2, Z = -3)

    # Rebasing Function:
    # Let's you re-scale the MSE in the table, based on the input.
    # Each 'x' will be a numeric vector corresponding to a single variable,
    # so the re-scaling can differ based on properties of that value. 
    # e.g.: rebase <- function(x) { 100*x } will multiply all MSEs by 100
    rebase <- function(x) { x }
    
    sim2 %>% 
        separate(name, c("Variable", "Method"), sep="\\s", extra="merge") %>% 
        mutate(Method = str_sub(Method, 2, -2)) %>% 
        mutate(square_error = (value - true_vals[Variable])^2) %>% 
        select(scenario, Variable, Method, value, square_error) %>%
        group_by(scenario, Variable, Method) %>% 
        summarize(MSE = mean(square_error)) %>% 
        pivot_wider(id_cols = c(scenario, Method), names_from = Variable, values_from = MSE) %>% 
        mutate_if(is.numeric, rebase)
}
