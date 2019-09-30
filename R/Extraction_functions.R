### Extraction Code ###
#
# Set of functions used to extract information from INLA which will be used
# as arguments for the sampling call
#

# selection matrix
makeAMat <- function(inla_model, effect_name,
                     constraint_point = which(day_constraint == 1)){
        return(inla_model %>%
                        reID(effect_name) %>%
                        createTransform(inla_model, constraint_point))
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

extractEffectCovMat <- function(prec, Amat){
        # takes in one precision matrix and selects the submatrix we need based on the AMat
        new_matrix <- new( "dsCMatrix",
                           x = prec@x,
                           i = prec@i,
                           p = prec@p,
                           Dim = prec@Dim )
        new_matrix_chol <- Matrix::Cholesky(new_matrix, LDL = TRUE, perm = TRUE)

        PA <- Matrix::solve(new_matrix_chol, Amat, system = "P")
        LinvPA <- Matrix::solve(new_matrix_chol, PA, system = "L")
        Dinvhalf <- Matrix::Diagonal( dim(new_matrix_chol)[1],
                                      1 / sqrt( new_matrix_chol@x[new_matrix_chol@p[1:nrow(new_matrix)]+1] ) )
        DinvhalfLinvPA <- Dinvhalf %*% LinvPA
        theVar <- Matrix::crossprod(DinvhalfLinvPA)

        return(theVar)
}



### A quick test
# prec_mat_list <- extractAllCovMat(real_rw_model)
# prec <- prec_mat_list[[1]]
# Amat <- makeAMat(inla_model = real_rw_model,
#                 effect_name = "day_count",
#                 constraint_point = which(day_constraint == 1))
# extractEffectCovMat(prec, Amat)
