## Produce a Summarized Version of the Results of simulation 1
combine_results_1 <- function(sim1_results_RC,
                              sim1_results_generalized_RC,
                              sim1_results_weighted_generalized_RC,
                              sim1_results_SIMEX,
                              sim1_results_ESIMEX,
                              sim1_results_generalized_SIMEX,
                              sim1_results_combine_first_generalized_SIMEX) {
    results_list <- list(
                sim1_results_RC,
                sim1_results_generalized_RC,
                sim1_results_weighted_generalized_RC,
                sim1_results_generalized_SIMEX,
                sim1_results_combine_first_generalized_SIMEX,
                sim1_results_SIMEX,
                sim1_results_ESIMEX)

    do.call(rbind, lapply(results_list, function(results){
        results %>% pivot_longer(cols = -c("scenario", "tar_batch", "tar_rep"))
    }))
}

# Plotting Helpers
# Generate the plot for each of the different variables
# This will be passed data either with scenarios that will (eventually) be
# faceted, *or* it will be filtered down to the scenario to begin.
X3_maker_1 <- function(data, true_vals=NA, nmethod=1) {
    data %>% ggplot(aes(y = value)) +
    geom_boxplot(aes(x = Method)) + 
    geom_segment(aes(x = 0, y = true_vals["X3"],  xend = nmethod+0.5, yend = true_vals["X3"]), lty=3) +
    labs(x="", y = expression(beta ~ " estimates"), title=expression(X[3])) + 
    theme_minimal() +
    scale_y_continuous(guide = guide_axis(check.overlap = TRUE))+
            theme(axis.text.y = element_text(size=10), 
                    axis.text.x = element_text(size=10),
                    panel.spacing.x = unit(4, "mm"))
}
X2_maker_1 <- function(data, true_vals=NA, nmethod=1){ 
    data %>% ggplot(aes(y = value)) +
    geom_boxplot(aes(x = Method)) + 
    geom_segment(aes(x = 0, y = true_vals["X2"],  xend = nmethod+0.5, yend = true_vals["X2"]), lty=3) +
    labs(x="", y = expression(beta ~ " estimates"), title=expression(X[2])) + 
    theme_minimal() +
    scale_y_continuous(guide = guide_axis(check.overlap = TRUE))+
            theme(axis.text.y = element_text(size=10), 
                    axis.text.x = element_text(size=10),
                    panel.spacing.x = unit(4, "mm"))
}
X1_maker_1 <- function(data, true_vals=NA, nmethod=1){
    data %>% ggplot(aes(y = value)) +
            geom_boxplot(aes(x = Method)) + 
            geom_segment(aes(x = 0, y = true_vals["X1"],  xend = nmethod+0.5, yend = true_vals["X1"]), lty=3) +
            labs(x="", y = expression(beta ~ " estimates"), title=expression(X[1])) + 
            theme_minimal() +
            scale_y_continuous(guide = guide_axis(check.overlap = TRUE))+
            theme(axis.text.y = element_text(size=10), 
                    axis.text.x = element_text(size=10),
                    panel.spacing.x = unit(4, "mm"))
}
int_maker_1 <- function(data, true_vals=NA, nmethod=1){
    data %>% ggplot(aes(y = value)) +
            geom_boxplot(aes(x = Method)) + 
            geom_segment(aes(x = 0, y = true_vals["int"],  xend = nmethod+0.5, yend = true_vals["int"]), lty=3) +
            labs(x="", y = expression(beta ~ " estimates"), title="Intercept") + 
            theme_minimal() +
            scale_y_continuous(guide = guide_axis(check.overlap = TRUE))+
            theme(axis.text.y = element_text(size=10), 
                    axis.text.x = element_text(size=10),
                    panel.spacing.x = unit(4, "mm"))
}
## Produce a Boxplot from the Results of SIMULATION 1
generate_plot_1 <- function(sim1) {
    # Specify the values here
    true_vals <- c(int = 2, X1 = -1, X2 = 2, X3 = 0.5)
    file_start<- "sim1_new_boxplot"
    folder <- "media"
    combined_plot <- paste0(folder, "/", file_start, ".png")

    # Actual Plot Code
    sim1 <- sim1 %>% 
                separate(name, c("Variable", "Method"), sep="\\s", extra="merge") %>% 
                mutate(Method = str_sub(Method, 2, -2))

    nmethod <- length(unique(sim1$Method))

    X3 <- sim1 %>% 
            filter(Variable == "X3") %>% 
            arrange(scenario) %>% 
            X3_maker_1(true_vals=true_vals, nmethod=nmethod) + facet_wrap(~scenario,  nrow = 1, scales="free_x") + coord_flip()

    X2 <- sim1 %>% 
            filter(Variable == "X2") %>% 
            arrange(scenario) %>% 
            X2_maker_1(true_vals=true_vals, nmethod=nmethod) + facet_wrap(~scenario,  nrow = 1, scales="free_x") + coord_flip()

    X1 <- sim1 %>% 
            filter(Variable == "X1") %>% 
            arrange(scenario) %>% 
            X1_maker_1(true_vals=true_vals, nmethod=nmethod) + facet_wrap(~scenario,  nrow = 1, scales="free_x") + coord_flip()

    int <- sim1 %>% 
            filter(Variable == "Intercept") %>% 
            arrange(scenario) %>% 
            int_maker_1(true_vals=true_vals, nmethod=nmethod) + facet_wrap(~scenario,  nrow = 1, scales="free_x") +  coord_flip()

    png(filename=combined_plot, width=14, height=9, units="in", res=300)
    grid.arrange(arrangeGrob(int, X1, X2, X3, nrow=4, left = textGrob("Estimation technique", rot = 90, vjust = 1)))
    dev.off()

    plots <- c(combined_plot)

    for (sc in unique(sim1$scenario)) {
        X3 <- sim1 %>% filter(Variable == "X3" && scenario == sc) %>% arrange(scenario) %>% X3_maker_1(true_vals=true_vals, nmethod=nmethod) + coord_flip()
        X2 <- sim1 %>%  filter(Variable == "X2" && scenario == sc) %>% arrange(scenario) %>% X2_maker_1(true_vals=true_vals, nmethod=nmethod) + coord_flip()
        X1 <- sim1 %>% filter(Variable == "X1" && scenario == sc) %>% arrange(scenario) %>% X1_maker_1(true_vals=true_vals, nmethod=nmethod) + coord_flip()
        int <- sim1 %>% filter(Variable == "Intercept" && scenario == sc) %>% arrange(scenario) %>% int_maker_1(true_vals=true_vals, nmethod=nmethod) + coord_flip()

        pname <- paste0(folder, "/", file_start, "_", sc, ".png")
        png(filename=pname, width=14, height=9, units="in", res=300)
        grid.arrange(arrangeGrob(int, X1, X2, X3, nrow=4, left = textGrob("Estimation technique", rot = 90, vjust = 1)))
        dev.off()
        plots <- c(plots, pname)
    }

    plots
}


generate_table_1 <- function(sim1) {
    # Specify the values here
    true_vals <- c(Intercept = 2, X1 = -1, X2 = 2, X3 = 0.5)

    # Rebasing Function:
    # Let's you re-scale the MSE in the table, based on the input.
    # Each 'x' will be a numeric vector corresponding to a single variable,
    # so the re-scaling can differ based on properties of that value. 
    # e.g.: rebase <- function(x) { 100*x } will multiply all MSEs by 100
    rebase <- function(x) { x }
    
    sim1 %>% 
        separate(name, c("Variable", "Method"), sep="\\s", extra="merge") %>% 
        mutate(Method = str_sub(Method, 2, -2)) %>% 
        mutate(square_error = (value - true_vals[Variable])^2) %>% 
        select(scenario, Variable, Method, value, square_error) %>%
        group_by(scenario, Variable, Method) %>% 
        summarize(MSE = mean(square_error)) %>% 
        pivot_wider(id_cols = c(scenario, Method), names_from = Variable, values_from = MSE) %>% 
        mutate_if(is.numeric, rebase)
}