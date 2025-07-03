# Required packages
library(lavaan)
library(dplyr)

# to view the forward selection process and model reduction the code down to line 172 shows the process. 
# If you just want to see the final model, run code from model 4 and onward

# Inputs
x_vars <- c(paste0("dim.", 1:5)) # direct variables

mediators = c(test_list_combo[c(1,3,4,7,8,9,14,16,17)]) # forward selection here

response <- "Y"

# Build X → Y (direct effects)
x_to_y <- paste0(response, " ~ ", 
                 paste0("X", seq_along(x_vars), "_to_Y*", x_vars, collapse = " + "))

# Build X → M (paths from predictors to each mediator)
x_to_m <- unlist(lapply(mediators, function(m) {
  paste0(m, " ~ ", 
         paste0("X", seq_along(x_vars), "_to_", m, "*", x_vars, collapse = " + "))
}))

# Build M → Y (mediators to outcome)
m_to_y <- paste0(response, " ~ ", 
                 paste0(mediators, "_to_Y*", mediators, collapse = " + "))

# Indirect effect expressions
indirect_lines <- sapply(mediators, function(m) {
  paste0("(", paste0("X", seq_along(x_vars), "_to_", m, collapse = " + "), ")*", m, "_to_Y")
})
indirect_expr <- paste(indirect_lines, collapse = " + ")

# Define effect labels
effects <- c(
  "direct := "   %>% paste0(paste0("X", seq_along(x_vars), "_to_Y", collapse = " + ")),
  "indirect := " %>% paste0(indirect_expr),
  "total := direct + indirect",
  "prop_mediated := indirect / total"
)

# Combine full model
model <- paste(
  "# X to Y (direct)",
  x_to_y,
  "\n# X to M",
  paste(x_to_m, collapse = "\n"),
  "\n# M to Y",
  m_to_y,
  "\n# Effects",
  paste(effects, collapse = "\n"),
  sep = "\n"
)


model_core <- paste(
  x_to_y,
  paste(x_to_m, collapse = "\n"),
  m_to_y,
  sep = "\n"
)

model_effects <- paste(effects, collapse = "\n")

fit_core <- sem(model_core, data =data_sem, std.lv = T)
summary(fit_core,rsquare=T,fit.measures=T)

lavResiduals(fit_core, type="cor.bollen" ) # check residuals are not inheritly signifcant wrt Y


# manually build the model from model_core to asses indirect effects (all mediators)

model3 <- "
# Direct effect of X (abiotic) on Y
Y ~ X1_to_Y*dim.1 + X2_to_Y*dim.2 + X3_to_Y*dim.3 + X4_to_Y*dim.4 + X5_to_Y*dim.5

# Mediator regressions (X to M)
M9  ~ X1_to_M9*dim.1 + X2_to_M9*dim.2 + X3_to_M9*dim.3 + X4_to_M9*dim.4 + X5_to_M9*dim.5
M22 ~ X1_to_M22*dim.1 + X2_to_M22*dim.2 + X3_to_M22*dim.3 + X4_to_M22*dim.4 + X5_to_M22*dim.5
M35 ~ X1_to_M35*dim.1 + X2_to_M35*dim.2 + X3_to_M35*dim.3 + X4_to_M35*dim.4 + X5_to_M35*dim.5
M53 ~ X1_to_M53*dim.1 + X2_to_M53*dim.2 + X3_to_M53*dim.3 + X4_to_M53*dim.4 + X5_to_M53*dim.5
M87 ~ X1_to_M87*dim.1 + X2_to_M87*dim.2 + X3_to_M87*dim.3 + X4_to_M87*dim.4 + X5_to_M87*dim.5
M40 ~ X1_to_M40*dim.1 + X2_to_M40*dim.2 + X3_to_M40*dim.3 + X4_to_M40*dim.4 + X5_to_M40*dim.5



# Mediators to Y
Y ~ M9_to_Y*M9 + M22_to_Y*M22 + M35_to_Y*M35 + M53_to_Y*M53 + M87_to_Y*M87 + M40_to_Y*M40

# Total direct effect
direct := X1_to_Y + X2_to_Y + X3_to_Y + X4_to_Y + X5_to_Y

# Total mediator effect
direct_fungi := M9_to_Y + M22_to_Y + M35_to_Y + M53_to_Y + M87_to_Y + M40_to_Y

# Indirect effects: dim.1
indirect1_1 := X1_to_M9  * M9_to_Y
indirect1_2 := X1_to_M22 * M22_to_Y
indirect1_3 := X1_to_M35 * M35_to_Y
indirect1_4 := X1_to_M53 * M53_to_Y
indirect1_5 := X1_to_M87 * M87_to_Y
indirect1_6 := X1_to_M40 * M40_to_Y


ind_X1 := indirect1_1 + indirect1_2 + indirect1_3 + indirect1_4 + indirect1_5 + indirect1_6
total_X1 := direct + ind_X1
prop_mediated_X1 := ind_X1 / total_X1

# Indirect effects: dim.2
indirect2_1 := X2_to_M9  * M9_to_Y
indirect2_2 := X2_to_M22 * M22_to_Y
indirect2_3 := X2_to_M35 * M35_to_Y
indirect2_4 := X2_to_M53 * M53_to_Y
indirect2_5 := X2_to_M87 * M87_to_Y
indirect2_6 := X2_to_M40 * M40_to_Y


ind_X2 := indirect2_6 
total_X2 := direct  + ind_X2
prop_mediated_X2 := ind_X2 / total_X2

# Indirect effects: dim.3
indirect3_1 := X3_to_M9  * M9_to_Y
indirect3_2 := X3_to_M22 * M22_to_Y
indirect3_3 := X3_to_M35 * M35_to_Y
indirect3_4 := X3_to_M53 * M53_to_Y
indirect3_5 := X3_to_M87 * M87_to_Y
indirect3_6 := X3_to_M40 * M40_to_Y


# nothing significant
ind_X3 := indirect3_1 + indirect3_2 + indirect3_3 + indirect3_4 + indirect3_5 + indirect3_6
total_X3 := direct + ind_X3
prop_mediated_X3 := ind_X3 / total_X3

# Indirect effects: dim.4
indirect4_1 := X4_to_M9  * M9_to_Y
indirect4_2 := X4_to_M22 * M22_to_Y
indirect4_3 := X4_to_M35 * M35_to_Y
indirect4_4 := X4_to_M53 * M53_to_Y
indirect4_5 := X4_to_M87 * M87_to_Y
indirect4_6 := X4_to_M40 * M40_to_Y


# nothing signifciant
ind_X4 := indirect4_1 + indirect4_2 + indirect4_3 + indirect4_4 + indirect4_5 + indirect4_6
total_X4 := direct  + ind_X4
prop_mediated_X4 := ind_X4 / total_X4

# Indirect effects: dim.5
indirect5_1 := X5_to_M9  * M9_to_Y
indirect5_2 := X5_to_M22 * M22_to_Y
indirect5_3 := X5_to_M35 * M35_to_Y
indirect5_4 := X5_to_M53 * M53_to_Y
indirect5_5 := X5_to_M87 * M87_to_Y
indirect5_6 := X5_to_M40 * M40_to_Y


# nothing signifcant
ind_X5 := indirect5_1 + indirect5_2 + indirect5_3 + indirect5_4 + indirect5_5 + indirect5_6
total_X5 := direct  + ind_X5
prop_mediated_X5 := ind_X5 / total_X5

"

# view model
fit_ind <- sem(model3, data =data_sem, std.lv = T)
summary(fit_ind,rsquare=T,fit.measures=T)



# the reduced model 

model4 <- "
# Direct effect of X (abiotic) on Y
Y ~ X1_to_Y*dim.1 + X2_to_Y*dim.2 + X5_to_Y*dim.5

# Mediator regressions (X to M)
M9  ~ X1_to_M9*dim.1 
M22 ~ X1_to_M22*dim.1 
M35 ~ X1_to_M35*dim.1 
M53 ~ X1_to_M53*dim.1 
M87 ~ X1_to_M87*dim.1 
M40 ~ X1_to_M40*dim.1 + X2_to_M40*dim.2


# Mediators to Y
Y ~ M9_to_Y*M9 + M22_to_Y*M22 + M35_to_Y*M35 + M53_to_Y*M53 + M87_to_Y*M87 + M40_to_Y*M40
# Total direct effect
direct := X1_to_Y + X2_to_Y + X5_to_Y

# Total mediator effect
direct_fungi := M9_to_Y + M22_to_Y + M35_to_Y + M53_to_Y + M87_to_Y + M40_to_Y

# Indirect effects: dim.1
indirect1_1 := X1_to_M9  * M9_to_Y
indirect1_2 := X1_to_M22 * M22_to_Y
indirect1_3 := X1_to_M35 * M35_to_Y
indirect1_4 := X1_to_M53 * M53_to_Y
indirect1_5 := X1_to_M87 * M87_to_Y
indirect1_6 := X1_to_M40 * M40_to_Y


ind_X1 := indirect1_1 + indirect1_2 + indirect1_3# + indirect1_4 + indirect1_5 + indirect1_6
total_X1 := direct + ind_X1
prop_mediated_X1 := ind_X1 / total_X1

# Indirect effects: dim.2
indirect2_6 := X2_to_M40 * M40_to_Y



ind_X2 := indirect2_6 
total_X2 := direct  + ind_X2
prop_mediated_X2 := ind_X2 / total_X2

# Indirect effects: dim.3
#indirect3_1 := X3_to_M9  * M9_to_Y
#indirect3_2 := X3_to_M22 * M22_to_Y
#indirect3_3 := X3_to_M35 * M35_to_Y
#indirect3_4 := X3_to_M53 * M53_to_Y
#indirect3_5 := X3_to_M87 * M87_to_Y
#indirect3_6 := X3_to_M40 * M40_to_Y
#indirect3_7 := X3_to_M19 * M19_to_Y
#indirect3_8 := X3_to_M56 * M56_to_Y
#indirect3_9 := X3_to_M65 * M65_to_Y

# nothing significant
#ind_X3 := indirect3_1 + indirect3_2 + indirect3_3 + indirect3_4 + indirect3_5 + indirect3_6
#total_X3 := direct + ind_X3
#prop_mediated_X3 := ind_X3 / total_X3

# Indirect effects: dim.4
#indirect4_1 := X4_to_M9  * M9_to_Y
#indirect4_2 := X4_to_M22 * M22_to_Y
#indirect4_3 := X4_to_M35 * M35_to_Y
#indirect4_4 := X4_to_M53 * M53_to_Y
#indirect4_5 := X4_to_M87 * M87_to_Y
#indirect4_6 := X4_to_M40 * M40_to_Y
#indirect4_7 := X4_to_M19 * M19_to_Y
#indirect4_8 := X4_to_M56 * M56_to_Y
#indirect4_9 := X4_to_M65 * M65_to_Y

# nothing signifciant
#ind_X4 := indirect4_1 + indirect4_2 + indirect4_3 + indirect4_4 + indirect4_5 + indirect4_6
#total_X4 := direct  + ind_X4
#prop_mediated_X4 := ind_X4 / total_X4

# Indirect effects: dim.5
#indirect5_1 := X5_to_M9  * M9_to_Y
#indirect5_2 := X5_to_M22 * M22_to_Y
#indirect5_3 := X5_to_M35 * M35_to_Y
#indirect5_4 := X5_to_M53 * M53_to_Y
#indirect5_5 := X5_to_M87 * M87_to_Y
#indirect5_6 := X5_to_M40 * M40_to_Y
#indirect5_7 := X5_to_M19 * M19_to_Y
#indirect5_8 := X5_to_M56 * M56_to_Y
#indirect5_9 := X5_to_M65 * M65_to_Y

# nothing signifcant
#ind_X5 := indirect5_1 + indirect5_2 + indirect5_3 + indirect5_4 + indirect5_5 + indirect5_6
#total_X5 := direct  + ind_X5
#prop_mediated_X5 := ind_X5 / total_X5

"

# view model
fit_ind <- sem(model4, data =data_sem, std.lv = T)
summary(fit_ind,rsquare=T,fit.measures=T)

# get standardized residuals
resi_z=lavResiduals(fit_ind, type="cor.bollen")
resi_z$cov.z

# case-wise residuals (varepsilon1)
preds=lavPredictY(fit_ind, ynames = c("Y"),
            xnames =colnames(lavInspect(fit_ind,what="data"))[-1])

y_obs=lavInspect(fit_ind,what="data")[,1]
resid_y=preds-y_obs
resid_y=as.numeric(resid_y)
# test
1-rmse(preds,y_obs)^2/var(y_obs) # ok!

library(car)
group_field_type  <- merged_data$field_type
group_season  <- merged_data$season

df_test <- data.frame(resid_y = resid_y, season = group_season,field_type=group_field_type)

leveneTest(resid_y ~ season, data = df_test)
leveneTest(resid_y ~ field_type, data = df_test)
leveneTest(resid_y ~ field_type:season, data = df_test)


library(patchwork)

df_plot <- data.frame(
  preds = preds,
  y_obs = y_obs,
  resid_y = resid_y,
  season = as.factor(merged_data$season),
  field_type = as.factor(merged_data$field_type),
  x10 = lavInspect(fit_ind, what = "data")[, 10]
)

# Create each plot
p1 <- ggplot(df_plot, aes(x = preds, y = resid_y, color = season, shape = field_type)) +
  geom_point(alpha = 0.7, size = 2.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Residuals vs Fitted",
    x = "Fitted (Predicted)", 
    y = "Residuals", 
    color = "Season", 
    shape = "Field Type"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

p2 <- ggplot(df_plot, aes(x = preds, y = y_obs, color = season, shape = field_type)) +
  geom_point(alpha = 0.7,size = 2.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
  labs(
    title = "Predicted vs Observed",
    x = "Predicted", 
    y = "Observed", 
    color = "Season",
    shape = "Field Type"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

p3 <- ggplot(df_plot, aes(sample = resid_y)) +
  stat_qq(alpha = 0.7) +
  stat_qq_line() +
  labs(title = "QQ Plot of Residuals") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Combine plots into one figure
(p1 | p2) / p3+
  plot_annotation(tag_levels = "A")  # Automatically adds A, B, C



# get case-wise residuals for varepilon2
M_preds=lavPredictY(fit_ind, ynames = c(colnames(lavInspect(fit_ind,what="data"))[c(2:7)]),
            xnames =colnames(lavInspect(fit_ind,what="data"))[c(8:9)])

M_preds=data.frame(M_preds)

M_obs_data=lavInspect(fit_ind,what="data")[,c(2:7)]

M_obs_data=data.frame(M_obs_data)

# test!
1-rmse(M_obs_data$M9,M_preds$M9)^2/var(M_obs_data$M9) # ok!
1-rmse(M_obs_data$M22,M_preds$M22)^2/var(M_obs_data$M22) # ok!
1-rmse(M_obs_data$M35,M_preds$M35)^2/var(M_obs_data$M35) # ok!
1-rmse(M_obs_data$M53,M_preds$M53)^2/var(M_obs_data$M53) # ok!
1-rmse(M_obs_data$M87,M_preds$M87)^2/var(M_obs_data$M87) # ok!
1-rmse(M_obs_data$M40,M_preds$M40)^2/var(M_obs_data$M40) # ok!

resid_M=M_preds-M_obs_data

# plot residuals
library(tidyr)
# Assuming merged_data (or similar) has season and field_type in the same row order
df_long <- data.frame(
  fitted = as.vector(as.matrix(M_preds)),
  residuals = as.vector(as.matrix(resid_M)),
  variable = rep(colnames(M_preds), each = nrow(M_preds)),
  season = rep(merged_data$season, times = ncol(M_preds)),
  field_type = rep(merged_data$field_type, times = ncol(M_preds))
)


library(gridExtra)

resid_vs_fitted_plots <- lapply(seq_along(unique(df_long$variable)), function(i) {
  var <- unique(df_long$variable)[i]
  df_sub <- subset(df_long, variable == var)
  ggplot(df_sub, aes(x = fitted, y = residuals, color = season, shape = field_type)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
   # annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = 1.5, size = 6, fontface = "bold") +
    labs(
      title = paste("Residuals vs Fitted for", var),
      x = "Fitted values",
      y = "Residuals",
      color = "Season",
      shape = "Field Type"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
})

qq_plots <- lapply(seq_along(unique(df_long$variable)), function(i) {
  var <- unique(df_long$variable)[i]
  df_sub <- subset(df_long, variable == var)
  ggplot(df_sub, aes(sample = residuals)) +
    stat_qq(alpha = 0.7) +
    stat_qq_line(color = "black") +
 #   annotate("text", x = -Inf, y = Inf, label = "B", hjust = -0.5, vjust = 1.5, size = 6, fontface = "bold") +
    labs(
      title = paste("QQ Plot of Residuals for", var),
      color = "Season",
      shape = "Field Type"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
})
# Combine plots side by side for each variable (resid vs fitted + QQ)
combined_plots <- mapply(function(rvf, qq) {
  rvf + qq + plot_annotation(tag_levels = "A")
}, resid_vs_fitted_plots, qq_plots, SIMPLIFY = FALSE)

# Display all combined plots
for (p in combined_plots) {
  print(p)
}

library(MVN)
mvn(resid_M, mvn_test = "mardia")

mvn(resid_M, mvn_test="mardia", multivariatePlot="qq")

# 1. Compute Mahalanobis distances
mahal_dist <- mahalanobis(resid_M, 
                          center = colMeans(resid_M), 
                          cov = cov(resid_M))

# 2. Compute Chi-squared quantiles
n <- nrow(resid_M)
p <- ncol(resid_M)
chi_sq_quantiles <- qchisq(ppoints(n), df = p)

# 3. Q-Q Plot
qqplot(chi_sq_quantiles, mahal_dist,
       main = "Mahalanobis Distances Q-Q Plot",
       xlab = expression(paste(chi^2, " Quantiles")),
       ylab = "Ordered Mahalanobis Distances",
       pch = 19, col = "darkblue")
abline(0, 1, col = "red", lwd = 2)  # Reference line


# do a bootstrap (it takes a while!- note results might varry a little)
fit_boot4 <- sem(model4, data = data_sem, se = "bootstrap", bootstrap = 500,std.lv = T)
summary(fit_boot4)



# estimate p.value from FDR
params <- parameterEstimates(fit_boot4)

params$pvalue_fdr <- p.adjust(params$pvalue, method = "fdr",n=length(params$pvalue))

# check FDR vs label
cbind(params$label,params$pvalue, params$pvalue_fdr)


# estimate post hoc power
library(semPower)
fitted_cov <- fitted(fit_boot4)$cov
sample_cov <- lavInspect(fit_boot4, "sampstat")$cov

h <- semPower.postHoc(alpha = .05, N = 72, df =26, SigmaHat = fitted_cov,Sigma = sample_cov)
h$power


# Monte Carlo simulation to asses fit further
library(simsem)
Output=sim(model=fit_ind,nRep=5000, n=72, sequential = T, seed = 123, lavaanfun="sem")
summary(Output)


# plot the model results and check model results are identical
library(semPlot)
library(semptools)
library(qgraph)

pe <- parameterEstimates(fit_ind, standardized = TRUE)

# Create color vector: black for significant, transparent for non-significant
edge_colors <- ifelse(pe$pvalue < 0.05 & pe$op %in% c("~", "=~"), "black", "transparent")

# using sempaths
semPaths(fit_ind, "std", edge.label.cex = 0.5, curvePivot = TRUE,style = "ram",nCharNodes = 0, nCharEdges = 0,layout = "spring")


# Get standardized estimates
est_std <- standardizedSolution(fit_ind,type = "std.lv")

# Directed edges: regressions and loadings
edges_directed <- est_std %>%
  filter(op %in% c("~", "=~"), pvalue < 0.05) %>%
  select(from = rhs, to = lhs, weight = est.std) %>%
  mutate(across(c(from, to), ~ gsub("^dim\\.(\\d+)$", "X\\1", .)))

# Covariance edges: significant residual covariances

# Combine all variables involved
vars <- unique(c(edges_directed$from, edges_directed$to))

# Define node shapes
shapes <- ifelse(grepl("^X\\d+$", vars), "square", "triangle")
shapes[10]="circle"
# Create empty adjacency matrix
adj <- matrix(0, nrow = length(vars), ncol = length(vars),
              dimnames = list(vars, vars))

# Add directed edges
for (i in 1:nrow(edges_directed)) {
  adj[edges_directed$from[i], edges_directed$to[i]] <- edges_directed$weight[i]
}


# Plot with qgraph
par(mar = c(5, 4, 6, 10))  # bottom, left, top, right
qgraph(adj,
       layout = "spring",
       edge.labels = TRUE,
       curvePivot = TRUE,
       vsize = 5,
       label.prop = 0.8,
       posCol = "red",
       negCol = "blue",
       shape = shapes)

edge_legend <- matrix(c(
  "Mortierella",      "M9", 
  "Trichocladium",    "M22",
  "Leucothecium",     "M35",
  "Phialemonium",     "M40",
  "Candida",          "M53",
  "Thermomyces",      "M87"
), ncol = 2, byrow = TRUE)


edge_legend_PC <- matrix(c(
  "Component 1",    "X1", 
  "Component 2",    "X2", 
  "Component 5", "X5"
), ncol = 2, byrow = TRUE)


par(xpd = TRUE)

# Set up side-by-side layout: 1 row, 2 columns
layout(matrix(c(1, 2), 1, 2), widths = c(4, 1))  # Left = plot, Right = legend


# Plot scaled edges WITHOUT labels
qgraph(adj,
      # title = "SEM mediation analysis model result",
       layout = "spring",
       edge.labels = TRUE,
       edge.label.cex = 0.9,     # Smaller edge label text
       edge.label.position = 0.7, # Label midpoint
       curve = 2.0,              # Stronger curvature
       curvePivot = TRUE,
       vsize = 5,
       asize = 3,                # Smaller arrowhead size
       fade=F,
       label.prop = 0.8,
       posCol = "red",
       negCol = "blue",
       shape = shapes)

# Legend panel
par(mar = c(0, 0, 0, 0))  # No margins
plot.new()                # Open a blank plotting area

legend("center",
       legend = paste(edge_legend[,2], "=", edge_legend[,1]),
       title = "\u25B3   Mediators",   # Unicode followed by space and title
       cex = 1.2,
       bty = "n",
       ncol = 1,
       pch = 2,
       inset = c(0.02, 0.1)) 

legend("top",
       legend = paste(edge_legend_PC[,2], "=", edge_legend_PC[,1]),
       title = "\u25A1   Abiotic effects",    # Unicode open triangle
       pch = 0,                       # triangle symbol
       cex = 1.2,
       bty = "n",
       inset = c(0.02, 0.2))

legend("bottom",
       legend = "Y = NDVI",
       title = "\u25CB  Response variable",    # Unicode open triangle
       pch = 1,                       # triangle symbol
       cex = 1.2,
       bty = "n",
       inset = c(0.02, 0.25))

