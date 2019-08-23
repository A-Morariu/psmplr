# psmplr package
#
# Posterior distribution sampling framework created to take in an INLA model
# and output a sample to mimic all configurations of the model
#
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

library(MASS)
library(Matrix)
library(tidyverse)

reID <- function(inla_object, re_name){ ### FIX TO WHICH INSTEAD OF GREP
        ### return a vector of the indicies (row and column) of the sub-matrix
        ### corresponding to the named random effect (ONLY works for 1 effect)
        name <- grep(re_name,inla_object$misc$configs$contents$tag)
        indicies <- seq(inla_object$misc$configs$contents$start[name],
                len = inla_object$misc$configs$contents$length[name])

        return(indicies)
}

#meanOffset <- function(re_id, inla_object){
#  ### returns the mean of the latent field from the inla_object
#  return(as.vector(inla_object$misc$configs$config[[1]]$mean)[re_id])
#}

createTransform <- function(re_index, inla_object, constraint_point){
        ### return the A matrix which selects the sub-matrix corresponding to
        ### the named random effect while placing a constraint at a single
        ### point in the walk
        # Take the sub-matrix of Q that we want
        # number of columns - length of indices (for the random effect of interest)
        # number of rows - to match dimensions of Q (for multiplication)

        dimension <- nrow(inla_object$misc$configs$config[[1]]$Q)

        ### Top
        top_block <- Matrix::Matrix(0,
                nrow = (min(re_index)-1),
                ncol = length(re_index))

        ### Bottom
        bot_block <- Matrix::Matrix(0,
                nrow = (dimension - max(re_index)),
                ncol = length(re_index))

        ### Middle
        if(is.null(constraint_point)){
                mid_block <- Matrix::Diagonal(n = max(re_index) - min(re_index) + 1)
        } else{
                ### Middle
                differencing_mat <- Matrix(0,
                        nrow = max(re_index) - min(re_index) + 1,
                        ncol = max(re_index) - min(re_index) + 1)
                differencing_mat[,as.numeric(constraint_point)] <- 1

                # get +/- 1 entries on each row except in at entry
                # [fixed_point, constraint_point] where we get a row of zeros
                mid_block <- Matrix::Diagonal(n = max(re_index) - min(re_index) + 1) - differencing_mat
        }
        return(rbind(top_block, mid_block, bot_block))
}

s.weights <- function(inla_object){
        weights <- c()
        for (i in 1:inla_object$misc$configs$nconfig) {
                weights[i] <- inla_object$misc$configs$config[[i]]$log.posterior
        }
        return(exp(weights)/sum(exp(weights)))
}

PosteriorSampler <- function(inla_object, index = 1, effect_name, n=1, constraint_point=2){
        Amat <- inla_object %>%
                reID(effect_name) %>%
                createTransform(inla_object, constraint_point)

        mean_vec <- as.vector(Matrix::crossprod(Amat,
                as.matrix(inla_object$misc$configs$config[[index]]$mean)))

        # To bypass ForceSymmetric
        # Cholesky is "slow"
        cov_mat <- new("dsCMatrix",
                x = inla_object$misc$configs$config[[index]]$Q@x,
                i = inla_object$misc$configs$config[[index]]$Q@i,
                p = inla_object$misc$configs$config[[index]]$Q@p,
                Dim = inla_object$misc$configs$config[[index]]$Q@Dim) %>%
                Matrix::Cholesky(LDL = FALSE, perm = FALSE) %>%
                Matrix::solve(Amat) %>%
                Matrix::crossprod()

        as.data.frame(MASS::mvrnorm(n=n, mu = mean_vec, Sigma = cov_mat))
}

psmplr <- function(inla_object, effect_name, n = 1, constraint_point = 2){
        samp_size <- ceiling(s.weights(inla_object)*n)

        paths <- list()
        n.theta <- inla_object$misc$configs$nconfig

        for (ii in 1:inla_object$misc$configs$nconfig) {
                paths <- append(paths, PosteriorSampler(inla_object = inla_object,
                        index = ii,
                        effect_name = effect_name,
                        n = samp_size[ii],
                        constraint_point = constraint_point))
        }
        return(paths)
}

##### Helper functions #####

makeAMat <- function(inla_model, effect_name, contraint_point = 2){
        return(inla_model %>%
                reID(effect_name) %>%
                createTransform(inla_object, constraint_point))
}

s.weights <- function(inla_object){
        weights <- c()
        for (i in 1:inla_object$misc$configs$nconfig) {
                weights[i] <- inla_object$misc$configs$config[[i]]$log.posterior
        }
        return(exp(weights)/sum(exp(weights)))
}

sampSizes <- function(inla_model, n = 1){
        return(ceiling(s.weights(inla_model)*n))
}

extractAllMeans <- function(inla_model){
        return(purrr::map(inla_model$misc$configs$config, function(xx) xx$mean))
}

extractAllCovMat <- function(inla_model){
        return(purrr::map(inla_model$misc$configs$config, function(xx) xx$Q))
}

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
                Matrix::crossprod())
}
