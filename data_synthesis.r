# %% load necessary libs
library(phyloseq)
library(SPRING)
library(SpiecEasi)
library(JGL)
library(EstimateGroupNetwork)
library(SpiecEasi)
library(igraph)
library(caret)     # For confusionMatrix function
library(pROC)      # For ROC and AUC calculation
library(MLmetrics) # For F1 Score, MCC

data_soil_raw <- readRDS("data\\soil_sample_gr2.rds")
# Remove taxa not seen more than 3 times in at least 20% of the samples.
data_soil = filter_taxa(data_soil_raw, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
otu_tax <- tax_table(data_soil)[,1:5]
otu_Ab <- t(otu_table(data_soil)) #any(rowSums(otu_Ab) == 0)

# Split by soil conditions and filter 0 counts taxas
soil_info <- sample_data(data_soil)[, "Soil_Type"]
otu_Ab_naural <- otu_Ab[soil_info == "Forest",
                        colSums(otu_Ab[soil_info == "Forest",]) != 0]
otu_Ab_potting <- otu_Ab[soil_info == "Potting",
                        colSums(otu_Ab[soil_info == "Potting",]) != 0]

# Keep only shared taxa/cols between otu_Ab_naural and otu_Ab_potting
shared_taxa <- intersect(colnames(otu_Ab_naural), colnames(otu_Ab_potting))
otu_Ab_naural <- otu_Ab_naural[, shared_taxa]
otu_Ab_potting <- otu_Ab_potting[, shared_taxa]

preprocess_and_estimate_network <- function(data_list, labels = NULL,  nlambda1 = 10, nlambda2 = 10) {
  # Helper function for common scaling normalization
  norm_to_total <- function(x) x / sum(x)
  common_scaling <- function(data) {
    depths <- rowSums(data)
    data_normalized <- t(apply(data, 1, norm_to_total))
    common_depth <- min(depths)  # Calculate only once
    data_common_scaled <- round(data_normalized * common_depth) 
    return(data_common_scaled)
  }
  # Preprocessing steps for each data set in the list
  processed_data_list <- lapply(data_list, function(data) {
    scaled_data <- common_scaling(data)
    mclr_data <- mclr(scaled_data) # Using SPRING's mclr
    Kcor <- mixedCCA::estimateR(mclr_data, type = "trunc", method = "approx", tol = 1e-6, verbose = FALSE)$R
    return(Kcor)
  })

  # Sample sizes for each dataset.
  n_samples <- sapply(data_list, nrow)

  # Estimate the network using the preprocessed data
  Res <- EstimateGroupNetwork(
    processed_data_list, 
    inputType = "list.of.covariance.matrices",
    n = n_samples,
    labels = labels, # Pass labels directly 
    nlambda1 = nlambda1, 
    nlambda2 = nlambda2, 
    truncate = 1e-10, 
    criterion = 'aic'
  )

  return(Res)
}
data_list <- list(naural = otu_Ab_naural, potting = otu_Ab_potting)
network_results <- preprocess_and_estimate_network(data_list, labels = shared_taxa)
network_naural <- network_results$naural
network_potting <- network_results$potting

# %% evaluation
# Create the true metwork matrix obtaned from the data
true_adj_naural <- (network_naural !=0)*1
true_adj_potting <- (network_potting !=0)*1

# %%Create a list to store the simulated matrices
set.seed(10010)
synthesize_scaled_data <- function(dat, net)
  {
    graph <- (net != 0)*1     # SpiecEasi::make_graph('cluster', dat$n_OTUs, dat$n_edge)
    attr(graph, "class") <- "graph"
    Prec <- SpiecEasi::graph2prec(graph)
    Cor <- cov2cor(SpiecEasi::prec2cov(Prec))
    X <- SpiecEasi::synth_comm_from_counts(dat, mar = 2, distr = 'zinegbin', Sigma = Cor, n = nrow(dat))
    return(X)
  }

Sim_list <- list()
for (i in 1:10)
  {
    Sim_list[[i]] <- list(
      naural = synthesize_scaled_data(otu_Ab_naural, network_naural),
      potting = synthesize_scaled_data(otu_Ab_potting, network_potting)
    )
  }
# %%Simulation: estimated matrix as starting point close to true network.  TP/FP/Recall ... HV
# Create a list to store the simulated matrices network
Res_sim <- list()
for (i in 1:10)
  {
    Res_sim[[i]] <- preprocess_and_estimate_network(Sim_list[[i]], labels = shared_taxa)
  }

# %%Create a list to store the simulated adjacency matrices
Sim_adj <- list()
for (i in 1:10)
  {
    Sim_adj[[i]] <- list(
      naural = (Res_sim[[i]]$naural !=0)*1,
      potting = (Res_sim[[i]]$potting !=0)*1
    )
  }
 
# %%Define a function to calculate all performance metrics
calculate_mcc <- function(tp, tn, fp, fn)
  {
    numerator <- (tp * tn) - (fp * fn)
    denominator <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  
    # If the denominator is 0, return 0 to avoid NaN errors
    if (denominator == 0) {
      return(0)
    } else {
      return(numerator / denominator)
    }
  }

calculate_metrics <- function(true_adj, sim_adj)
  {
    # Flatten the matrices into vectors (upper triangular part, excluding diagonal)
    true_edges <- as.vector(true_adj[upper.tri(true_adj)])
    sim_edges <- as.vector(sim_adj[upper.tri(sim_adj)])
  
    # Confusion matrix
    cm <- confusionMatrix(as.factor(sim_edges), as.factor(true_edges), positive = "1")
    tn <- as.numeric(cm$table[1,1]) #true negatives
    fp <- as.numeric(cm$table[1,2]) #false positives
    fn <- as.numeric(cm$table[2,1]) #false negatives
    tp <- as.numeric(cm$table[2,2]) #true positives
  
    # Calculate TPR, FPR, Precision, Recall
    tpr <- tp / (tp + fn)  # Sensitivity / Recall
    fpr <- fp / (fp + tn)  # 1 - Specificity
    precision <- tp / (tp + fp)
    recall <- tpr
  
    # Calculate ROC and AUC
    roc_obj <- roc(true_edges, sim_edges)
    auc <- auc(roc_obj) # plot here
  
    # F1 Score
    f1 <- F1_Score(sim_edges, true_edges)
  
    # MCC (Matthews correlation coefficient)
    mcc <- calculate_mcc(tp, tn, fp, fn)
  
    # Return metrics as a list
    return(list(TPR = tpr, FPR = fpr, Precision = precision, Recall = recall,
                F1 = f1, AUC = auc, MCC = mcc))
  }

#%%Compare the true graph with each simulated graph
confusion_results <- lapply(1:10, function(i) {
  # Calculate metrics for both natural and potting
  natural_metrics <- calculate_metrics(true_adj_naural, Sim_adj[[i]]$naural)
  potting_metrics <- calculate_metrics(true_adj_potting, Sim_adj[[i]]$potting)
  # Return both results as a named list
  return(list(natural = natural_metrics, potting = potting_metrics))
})

# Convert the results list into a data frame for easy analysis
results_df <- do.call(rbind, lapply(confusion_results, as.data.frame))


# %% Plot edge weights distribution
library(ggplot2)

# Convert the correlation matrices to vectors
cor_values_naural <- as.vector(Kcor_naural)
cor_values_potting <- as.vector(Kcor_potting)

# Create a data frame for plotting
cor_df <- data.frame(
  correlation = c(cor_values_naural, cor_values_potting),
  group = factor(rep(c("Natural", "Potting"), each = length(cor_values_naural)))
)

# Plot histogram
hist_plot <- ggplot(cor_df, aes(x = correlation, fill = group)) +
  geom_histogram(position = "dodge", bins = 30, alpha = 0.7) +
  scale_fill_manual(values = c("Natural" = "blue", "Potting" = "red")) +
  theme_minimal() +
  labs(title = "Distribution of Correlation Values",
       x = "Correlation",
       y = "Frequency",
       fill = "Soil Type")

# Save the plot
ggsave("histogram_comparison.pdf", hist_plot, width = 10, height = 6)

# Return the plot
return(hist_plot)

# %%Plot Network
library(ggraph)
library(igraph)
library(tidyverse)
library(RColorBrewer)

# Extract taxonomic information
otu_tax_df <- otu_tax %>%
  as.data.frame() %>%
  rownames_to_column("OTU") %>%
  select(-OTU, everything(), OTU)

pairs <- list(
  c("Kingdom", "Phylum"),
  c("Phylum", "Class"),
  c("Class", "Order"),
  c("Order", "Family"),
  c("Family", "OTU"))
# Function to create edges for a single pair
create_edges <- function(pair, data)
  {
    from <- data[[pair[1]]]
    to <- data[[pair[2]]]
    data.frame(from = from, to = to)
  }

# Apply the function to all pairs and combine the results
edges <- unique(do.call(rbind, lapply(pairs, create_edges, data = otu_tax_df[otu_tax_df$OTU %in% shared_taxa, ])))

# Extract lower triangular part
lower_tri <- lower.tri(network_naural, diag = FALSE)
# Get non-zero elements and their indices
non_zero <- which(lower_tri & network_naural != 0, arr.ind = TRUE)
# Create the new table
connect <- data.frame(
  from = rownames(network_naural)[non_zero[, 1]],
  to = colnames(network_naural)[non_zero[, 2]],
  score = network_naural[non_zero])

# create a vertices data.frame. One line per object of our hierarchy
vertices  <-  data.frame(
  name = unique(c(as.character(edges$from), as.character(edges$to))) , 
  value = runif(length(unique(c(as.character(edges$from), as.character(edges$to))))))
vertices$group  <-  edges$from[ match( vertices$name, edges$to ) ]

vertices$id <- NA
myleaves <- which(is.na(match(vertices$name, edges$from)))
vertices$value[myleaves] <- colSums(network_naural) / sum(network_naural)
nleaves <- length(myleaves)
vertices$id[myleaves] <- seq(1:nleaves)
vertices$angle <- 90 - 360 * vertices$id / nleaves
vertices$angle <- ifelse(vertices$angle < -90, vertices$angle+180, vertices$angle)
vertices$hjust <- ifelse( vertices$angle < -90, 1, 0)
mygraph <- igraph::graph_from_data_frame( edges, vertices=vertices )
from  <-  match( connect$from, vertices$name)
to  <-  match( connect$to, vertices$name)

ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
  geom_conn_bundle(data = get_con(from = from, to = to), alpha=0.2, width=0.1, aes(colour = after_stat(index))) +
  scale_edge_colour_gradient(low = "red", high = "blue") +
  geom_node_text(aes(x = x*1.15, y=y*1.15, filter = leaf, label=name, angle = angle, hjust=hjust, colour=group), size=2, alpha=1) +
  geom_node_point(aes(filter = leaf, x = x*1.07, y=y*1.07, colour=group, size=value, alpha=0.2)) +
  scale_colour_manual(values= rep( brewer.pal(14,"Paired") , 30)) +
  scale_size_continuous(range = c(0.1,10) ) +

  theme_void() +
  theme(
    legend.position="none",
    plot.margin=unit(c(0,0,0,0),"cm"),
  ) +
  expand_limits(x = c(-1.3, 1.3), y = c(-1.3, 1.3))
ggsave("naural_plot_synthesized.pdf", width = 12, height = 12, units = "in")

# igraph package for network analysis
# packages for ML network
# extend simulated data to large data set
