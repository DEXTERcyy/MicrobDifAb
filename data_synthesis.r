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

# scaling
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

# %% mclr
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