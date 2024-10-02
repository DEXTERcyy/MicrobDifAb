# %% load necessary libs
library(phyloseq)
library(SPRING)
library(SpiecEasi)
library(JGL)
library(EstimateGroupNetwork)
library(SpiecEasi)

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

# %%scaling
common_scaling <- function(data)
  {
    depths <- rowSums(data)
    data_normalized <- t(apply(data, 1, norm_to_total))
    data_common_scaled <- round(data_normalized * min(depths))
    
    d <- ncol(data_common_scaled)
    n <- nrow(data_common_scaled)
    e <- d

    return(list(
      processed_data = data_common_scaled,
      n_OTUs = d,
      n_samples = n,
      n_edge = e
    ))
  }

otu_Ab_naural_scaled <- common_scaling(otu_Ab_naural)
otu_Ab_potting_scaled <- common_scaling(otu_Ab_potting)

# synthesize data
synthesize_scaled_data <- function(dat)
  {
    set.seed(10010)
    graph <- SpiecEasi::make_graph('cluster', dat$n_OTUs, dat$n_edge)
    Prec <- SpiecEasi::graph2prec(graph)
    Cor <- cov2cor(SpiecEasi::prec2cov(Prec))
    
    X <- SpiecEasi::synth_comm_from_counts(dat$processed_data, mar = 2, distr = 'zinegbin', Sigma = Cor, n = dat$n_samples)
    
    return(X)
  }
otu_Ab_naural_synthesized <- synthesize_scaled_data(otu_Ab_naural_scaled)
otu_Ab_potting_synthesized <- synthesize_scaled_data(otu_Ab_potting_scaled)

# mClr
mclr <- function(dat, base = exp(1), tol = 1e-16, eps = NULL, atleast = 1) #from SPRING.r
  {
    dat <- as.matrix(dat)
    nzero <- (dat >= tol)  # index for nonzero part
    LOG <- ifelse(nzero, log(dat, base), 0.0) # take log for only nonzero values. zeros stay as zeros.

    # centralize by the log of "geometric mean of only nonzero part" # it should be calculated by each row.
    if (nrow(dat) > 1){
      clrdat <- ifelse(nzero, LOG - rowMeans(LOG)/rowMeans(nzero), 0.0)
    } else if (nrow(dat) == 1){
      clrdat <- ifelse(nzero, LOG - mean(LOG)/mean(nzero), 0.0)
    }

    if (is.null(eps)){
      if(atleast < 0){
        warning("atleast should be positive. The functions uses default value 1 instead.")
        atleast = 1
      }
      if( min(clrdat) < 0 ){ # to find the smallest negative value and add 1 to shift all data larger than zero.
        positivecst <- abs(min(clrdat)) + atleast # "atleast" has default 1.
      }else{
        positivecst <- 0
      }
      # positive shift
      ADDpos <- ifelse(nzero, clrdat + positivecst, 0.0) ## make all non-zero values strictly positive.
      return(ADDpos)
    } else if(eps == 0){
      ## no shift. clr transform applied to non-zero proportions only. without pseudo count.
      return(clrdat)
    } else if(eps > 0){
      ## use user-defined eps for additional positive shift.
      ADDpos <- ifelse(nzero, clrdat + eps, 0.0)
      return(ADDpos)
    } else {
      stop("check your eps value for additional positive shift. Otherwise, leave it as NULL.")
    }
  }
otu_Ab_naural_mclr <- mclr(otu_Ab_naural_synthesized)
otu_Ab_potting_mclr <- mclr(otu_Ab_potting_synthesized)
Kcor_naural <- mixedCCA::estimateR(otu_Ab_naural_mclr, type = "trunc", method = "approx", tol = 1e-6, verbose = TRUE)$R
Kcor_potting <- mixedCCA::estimateR(otu_Ab_potting_mclr, type = "trunc", method = "approx", tol = 1e-6, verbose = TRUE)$R
# Graphical Lasso
Res <- EstimateGroupNetwork(list(naural = Kcor_naural, potting = Kcor_potting), inputType = "list.of.covariance.matrices",
                          n = c(dim(otu_Ab_naural)[1],dim(otu_Ab_potting)[1]), labels = shared_taxa,
                          nlambda1 = 10, nlambda2 = 10, truncate = 1e-10, criterion = 'aic')
network_naural <- Res$naural
network_potting <- Res$potting

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