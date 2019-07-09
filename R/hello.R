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

reID <- function(inla_object, re_name){
     ### return the indicies (row and column index) of the sub-matrix
     ### corresponding to the named random effect (ONLY works for 1 effect)
     name <- grep(re_name,inla_object$misc$configs$contents$tag)
     indicies <- seq(inla_object$misc$configs$contents$start[name],
                     len = inla_object$misc$configs$contents$length[name])

     return(indicies)
}

meanOffset <- function(re_id, inla_object){
     ### returns the mean of the latent field from the inla_object
     return(as.vector(inla_object$misc$configs$config[[1]]$mean)[re_id])
}

createTransform <- function(re_index, inla_object, constraint_point){
     ### return the A matrix which selects the sub-matrix corresponding to
     ### the named random effect while placing a constraint at a single
     ### point in the walk

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

PosteriorSampler <- function(inla_object, effect_name, n=1, constraint_point){
     Amat <- inla_object %>%
          reID("u1") %>%
          createTransform(inla_object, constraint_point)

     mean_vec <- as.vector(Matrix::crossprod(Amat,
                                             as.matrix(inla_object$misc$configs$config[[1]]$mean)))

     cov_mat <- new("dsCMatrix",
                    x = inla_object$misc$configs$config[[1]]$Q@x,
                    i = inla_object$misc$configs$config[[1]]$Q@i,
                    p = inla_object$misc$configs$config[[1]]$Q@p,
                    Dim = inla_object$misc$configs$config[[1]]$Q@Dim) %>%
          Matrix::Cholesky(LDL = FALSE, perm = FALSE) %>%
          Matrix::solve(Amat) %>%
          Matrix::crossprod()

     MASS::mvrnorm(n=n, mu = mean_vec, Sigma = cov_mat)
}
