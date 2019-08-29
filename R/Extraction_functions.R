### Extraction Code ###
#
# Set of functions used to extract information from INLA which will be used
# as arguments for the sampling call
#

# selection matrix
makeAMat <- function(inla_model, effect_name, contraint_point = 2){
     return(inla_model %>%
               reID(effect_name) %>%
               createTransform(inla_model, contraint_point))
}

# sample size manipulation
s.weights <- function(inla_model){
     weights <- c()
     for (i in 1:inla_model$misc$configs$nconfig) {
          weights[i] <- inla_model$misc$configs$config[[i]]$log.posterior
     }
     return(weights/sum(weights))
}

sampSizes <- function(inla_model, n = 1){
     return(as.vector(ceiling(s.weights(inla_model)*n) ) )
}

# extract the necessary components for sampling from INLA model
extractAllMeans <- function(inla_model){
     return(purrr::map(inla_model$misc$configs$config, function(xx) xx$mean))
}

extractAllCovMat <- function(inla_model){
     return(purrr::map(inla_model$misc$configs$config, function(xx) xx$Q))
}

# keep heavy multiplication steps seperate
extractEffectMeans <- function(mu, Amat){
     # takes in one mean vector and selects the bits we need based on the AMat
     return(as.vector(Matrix::crossprod(Amat,as.matrix(mu) ) ) )
}

extractEffectCovMat <- function(sigma, Amat){
     # takes in one precision matrix and selects the submatrix we need based on the AMat
     return(new("dsCMatrix",
          x = sigma@x,
          i = sigma@i,
          p = sigma@p,
          Dim = sigma@Dim) %>%
               Matrix::Cholesky(LDL = FALSE, perm = FALSE) %>%
               Matrix::solve(Amat) %>%
               Matrix::crossprod() )
}
