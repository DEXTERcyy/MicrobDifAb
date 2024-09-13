# %% load necessary libs
library(phyloseq)
library(SPRING)
library(SpiecEasi)
library(JGL)

data_soil_raw <- readRDS("data\\soil_sample_gr2.rds")
# Remove taxa not seen more than 3 times in at least 20% of the samples.
data_soil = filter_taxa(data_soil_raw, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
otu_tax <- tax_table(data_soil) # Assign variables
otu_Ab <- t(otu_table(data_soil))
# any(rowSums(otu_Ab) == 0)

soil_info <- sample_data(data_soil)[, "Soil_Type"]
otu_Ab_natural <- otu_Ab[soil_info == "Forest",
                        colSums(otu_Ab[soil_info == "Forest",]) != 0]
otu_Ab_potting <- otu_Ab[soil_info == "Potting",
                        colSums(otu_Ab[soil_info == "Potting",]) != 0]
# any(colSums(otu_Ab_naural) == 0) || any(colSums(otu_Ab_potting) == 0)))

# %% SPRING
otu_Ab_naural_fit <- SPRING(otu_Ab_naural, Rmethod = "approx", quantitative = TRUE, 
              lambdaseq = "data-specific", nlambda = 20, rep.num = 20)
opt.K_naural <- otu_Ab_naural_fit$output$stars$opt.index
# Estimated adjacency matrix from sparse graphical modeling technique ("mb" method) (1 = edge, 0 = no edge)
adj.K_naural <- as.matrix(otu_Ab_naural_fit$fit$est$path[[opt.K_naural]])
# Estimated partial correlation coefficient, same as negative precision matrix.
pcor.K_naural <- as.matrix(SpiecEasi::symBeta(otu_Ab_naural_fit$output$est$beta[[opt.K_naural]], mode = 'maxabs'))

otu_Ab_potting_fit <- SPRING(otu_Ab_potting, Rmethod = "approx", quantitative = TRUE, 
              lambdaseq = "data-specific", nlambda = 20, rep.num = 20)
opt.K_potting <- otu_Ab_potting_fit$output$stars$opt.index
adj.K_potting <- as.matrix(otu_Ab_potting_fit$fit$est$path[[opt.K_potting]])
pcor.K_potting <- as.matrix(SpiecEasi::symBeta(otu_Ab_potting_fit$output$est$beta[[opt.K_potting]], mode = 'maxabs'))
# isSymmetric(pcor.K_naural) && isSymmetric(adj.K_potting)
# anyNA(pcor.K_naural) || anyNA(pcor.K_potting)

# %% Save progress
# write.csv(pcor.K_naural, "pcor_naural.csv")
# write.csv(pcor.K_potting, "pcor_potting.csv")
# save.image(file = "workspace.RData")

# %% Load progress
# load("workspace.RData")
# library(JGL)
# data(example.data)
# str(example.data)
# JGL(Y=list(example.data[[1]],example.data[[2]]),penalty="fused",lambda1=.25,lambda2=.1)

# %%
# JGL(Y=list(as(otu_table(otu_Ab_naural), "matrix"),as(otu_table(otu_Ab_potting), "matrix")),penalty="fused",lambda1=.25,lambda2=.1)
result <- JGL(Y = list(pcor.K_naural, pcor.K_potting), 
              penalty = "fused", # FGL with λ1 =0.2, λ2 = 0.1 and GGL with λ1 = 0.05, λ2 = 0.25
              lambda1 = 0.25,  # sparsity of the estimated precision matrices
              lambda2 = 0.1, # fused penalty parameter / similarity between the estimated precision matrices
              penalize.diagonal = FALSE,
              return.whole.theta = FALSE) #full estimated precision matrices
