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
# %% SPRING
otu_Ab_naural_fit <- SPRING(otu_Ab_naural, Rmethod = "approx", quantitative = TRUE, 
              lambdaseq = "data-specific", nlambda = 20, rep.num = 20)
opt.K_naural <- otu_Ab_naural_fit$output$stars$opt.index
# Estimated adjacency matrix from sparse graphical modeling technique ("mb" method) (1 = edge, 0 = no edge)
adj.K_naural <- as.matrix(otu_Ab_naural_fit$fit$est$path[[opt.K_naural]])
# Estimated partial correlation coefficient, same as negative precision matrix.
pcor.K_naural <- as.matrix(SpiecEasi::symBeta(otu_Ab_naural_fit$output$est$beta[[opt.K_naural]], mode = 'maxabs')) 
colnames(pcor.K_naural) <- colnames(otu_Ab_naural) # set pcor.K_naural colnames

otu_Ab_potting_fit <- SPRING(otu_Ab_potting, Rmethod = "approx", quantitative = TRUE, 
              lambdaseq = "data-specific", nlambda = 20, rep.num = 20)
opt.K_potting <- otu_Ab_potting_fit$output$stars$opt.index
adj.K_potting <- as.matrix(otu_Ab_potting_fit$fit$est$path[[opt.K_potting]])
pcor.K_potting <- as.matrix(SpiecEasi::symBeta(otu_Ab_potting_fit$output$est$beta[[opt.K_potting]], mode = 'maxabs'))
colnames(pcor.K_potting) <- colnames(otu_Ab_potting)
#isSymmetric(pcor.K_naural) && isSymmetric(adj.K_potting)
#anyNA(pcor.K_naural) || anyNA(pcor.K_potting)

# %% Save progress
#write.csv(pcor.K_naural, "pcor_naural.csv")
#write.csv(pcor.K_potting, "pcor_potting.csv")
save.image(file = "workspace_JGL_pcor.RData")
# %% Load progress
load("workspace_JGL_pcor.RData")
# Fused Graphical Lasso on Partial Correlation Matrices
soft <- function(a,lam,penalize.diagonal)
  { # if last argument is FALSE, soft-threshold a matrix but don't penalize the diagonal
    out <- sign(a)*pmax(0, abs(a)-lam)
    if(!penalize.diagonal) diag(out) <- diag(a)
    return(out)
  }


flsa2 <-function(A,L,lam1,lam2,penalize.diagonal)  #A is a list of 2 matrices from which we apply an L2 penalty to departures
  {
  # 1st apply fused penalty:
  S1 = abs(A[[1]]-A[[2]])<=2*lam2/L
  X1 = (A[[1]]+A[[2]])/2
  Y1 = X1

  S2 = (A[[1]] > A[[2]]+2*lam2/L)
  X2 = A[[1]] - lam2/L
  Y2 = A[[2]] + lam2/L

  S3 = (A[[2]] > A[[1]]+2*lam2/L)
  X3 = A[[1]] + lam2/L
  Y3 = A[[2]] - lam2/L

  X = soft(a = S1*X1 + S2*X2 + S3*X3, lam = lam1/L, penalize.diagonal=penalize.diagonal)
  Y = soft(a = S1*Y1 + S2*Y2 + S3*Y3, lam = lam1/L, penalize.diagonal=penalize.diagonal)

  return(list(X,Y))
  }
penalty.as.matrix <-function(lambda,p,penalize.diagonal) # from JGL.r
  {
    # for matrix penalties:  check dim and symmetry:
    if(is.matrix(lambda))
    {
      if(sum(lambda!= t(lambda))>0) {stop("error: penalty matrix is not symmetric")}
      if(sum(abs(dim(lambda)-p))!=0 ) {stop("error: penalty matrix has wrong dimension")}
    }
    # for scalar penalties: convert to matrix form:
    if(length(lambda)==1) {lambda=matrix(lambda,p,p)}
    # apply the penalize.diagonal argument:
    if(!penalize.diagonal) {diag(lambda)=0}
    return(lambda)
  }
JGL_modified <- function(corr_matrices, lambda1, lambda2, rho = 1, weights = "equal", penalize.diagonal = FALSE, maxiter = 500, tol = 1e-5, warm = NULL, return.whole.theta = FALSE, screening = "fast", truncate = 1e-5) 
  {
    # Initialize:
    K <- length(corr_matrices)
    p <- dim(corr_matrices[[1]])[1]
    n <- rep(1, K)  # Set sample sizes to 1 since we're using correlation matrices

    # Set weights:
    if (length(weights) == 1) {
      if (weights == "equal") {
        weights <- rep(1, K)
      }
    }

    # Define S as the correlation matrices
    S <- corr_matrices

    # Initialize theta:
    theta <- lapply(1:K, function(k) diag(p))

    # Initialize Z and W:
    Z <- W <- lapply(1:K, function(k) matrix(0, p, p))

    # Initialize lambdas:
    lam1 <- penalty.as.matrix(lambda1, p, penalize.diagonal = penalize.diagonal)
    lam2 <- penalty.as.matrix(lambda2, p, penalize.diagonal = TRUE)

    # ADMM iterations:
    iter <- 0
    diff_value <- 10
    while ((iter == 0) || (iter < maxiter && diff_value > tol)) {
      # Update theta:
      theta.prev <- theta
      for (k in 1:K) {
        edecomp <- eigen(S[[k]] - rho * Z[[k]] / weights[k] + rho * W[[k]] / weights[k])
        D <- edecomp$values
        V <- edecomp$vectors
        D2 <- weights[k] / (2 * rho) * (-D + sqrt(D^2 + 4 * rho / weights[k]))
        theta[[k]] <- V %*% diag(D2) %*% t(V)
      }

      # Update Z:
      A <- lapply(1:K, function(k) theta[[k]] + W[[k]])
      Z <- flsa2(A, rho, lam1, lam2, penalize.diagonal = TRUE)

      # Update W:
      for (k in 1:K) {
        W[[k]] <- W[[k]] + (theta[[k]] - Z[[k]])
      }

      # Bookkeeping:
      iter <- iter + 1
      diff_value <- sum(sapply(1:K, function(k) sum(abs(theta[[k]] - theta.prev[[k]])) / sum(abs(theta.prev[[k]]))))
      
      # Increment rho:
      rho <- rho * 1.01
    }

    # Prepare output:
    diff <- sum(sapply(1:K, function(k) sum(abs(theta[[k]] - Z[[k]]))))
    
    if (!return.whole.theta) {
      out <- list(theta = Z, diff = diff, iters = iter)
    } else {
      out <- list(theta = Z, diff = diff, iters = iter)
    }
    
    class(out) <- "jgl"
    return(out)
  }
jgl_result <- JGL_modified(corr_matrices = list(pcor.K_naural, pcor.K_potting), 
                       lambda1 = 0.1, 
                       lambda2 = 0.2)

# %% Visualizations
library(igraph)
library(RColorBrewer)
theta1 <- jgl_result$theta[[1]]
theta2 <- jgl_result$theta[[2]]

# Adjacency matrices
threshold <- 0.01
adj1 <- abs(theta1) > threshold
adj2 <- abs(theta2) > threshold
adj_combined <- adj1 | adj2

# Create igraph
g <- graph_from_adjacency_matrix(adj_combined, mode = "undirected", diag = FALSE)

# Set edge colors
edge_colors <- rep("grey", ecount(g))
edge_colors[which(adj1[lower.tri(adj1)] & adj2[lower.tri(adj2)])] <- "purple"  # Edges in both conditions
edge_colors[which(adj1[lower.tri(adj1)] & !adj2[lower.tri(adj2)])] <- "red"    # Edges only in naural
edge_colors[which(!adj1[lower.tri(adj1)] & adj2[lower.tri(adj2)])] <- "blue"   # Edges only in potting

# Set vertex colors
vertex_colors <- brewer.pal(9, "Set1")[1]

# Plot the graph
plot(g,
      vertex.color = vertex_colors,
      vertex.size = 5,
      vertex.label = V(g)$name,
      vertex.label.cex = 0.8,
      edge.color = edge_colors,
      main = "Network Graph")

# Add legend
# legend("bottomright", 
#       legend = c("Both conditions", "Naurual only", "Potting only"),
#       col = c("purple", "red", "blue"),
#       lwd = 2,
#       cex = 0.8)