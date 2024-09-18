# %% load necessary libs
library(phyloseq)
library(SPRING)
library(SpiecEasi)
library(JGL)
library("qgraph")
library("parallel")
library("psych")
library("mvtnorm")
library("EstimateGroupNetwork")
ncores <- detectCores() -1
options(mc.cores = ncores)

data_soil_raw <- readRDS("data\\soil_sample_gr2.rds")

# Remove taxa not seen more than 3 times in at least 20% of the samples.
data_soil = filter_taxa(data_soil_raw, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
otu_tax <- tax_table(data_soil) # Assign variables
otu_Ab <- t(otu_table(data_soil))
#any(rowSums(otu_Ab) == 0)
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
#any(colSums(otu_Ab_naural) == 0) || any(colSums(otu_Ab_potting) == 0)  

# %% SPRC
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
otu_Ab_naural_mclr <- mclr(as(otu_Ab_naural,"matrix"))
otu_Ab_potting_mclr <- mclr(as(otu_Ab_potting,"matrix"))
Kcor_naural <- mixedCCA::estimateR(otu_Ab_naural_mclr, type = "trunc", method = "approx", tol = 1e-6, verbose = TRUE)$R
Kcor_potting <- mixedCCA::estimateR(otu_Ab_potting_mclr, type = "trunc", method = "approx", tol = 1e-6, verbose = TRUE)$R

# %%
library("qgraph")
library("parallel")
library("psych")
library("mvtnorm")
Res <- EstimateGroupNetwork(list(naural = Kcor_naural, potting = Kcor_potting), inputType = "list.of.covariance.matrices",
                          n = c(dim(otu_Ab_naural)[1],dim(otu_Ab_potting)[1]),
                          nlambda1 = 2, nlambda2 = 2, truncate = 1e-10)
Layout <- averageLayout(Res$naural,Res$potting)
layout(t(1:2))
qgraph(Res$naural, layout = Layout, title = "Naural (JGL)")
qgraph(Res$potting, layout = Layout, title = "Potting (JGL)")

# pdf("data\\EGN_SPRC_JGL.pdf")
# layout(t(1:2))
# qgraph(Res$naural, layout = Layout, title = "Naural (JGL)")
# qgraph(Res$potting, layout = Layout, title = "Potting (JGL)")
# dev.off()
# %% example
# data(bfi)
# bfi2 <- bfi[rowSums(is.na(bfi[,1:26])) == 0,]
# CorMales <- cor_auto(bfi2[bfi2$gender == 1,1:25])
# CorFemales <- cor_auto(bfi2[bfi2$gender == 2,1:25])
# Res_exa <- EstimateGroupNetwork(list(males = CorMales, females = CorFemales),
# n = c(sum(bfi2$gender == 1),sum(bfi2$gender == 2)))
# qgraph(Res_exa$males, layout = Layout, title = "Males")
