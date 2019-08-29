### Create "AMat" ###
#
# Series of functions used to create the AMat - the matrix which contrains and
# selects the sub-matrix of the preicions from INLA
#

reID <- function(inla_model, re_name){ ### FIX TO WHICH INSTEAD OF GREP
     ### return a vector of the indicies (row and column) of the sub-matrix
     ### corresponding to the named random effect (ONLY works for 1 effect)
     name <- grep(re_name,inla_model$misc$configs$contents$tag)
     indicies <- seq(inla_model$misc$configs$contents$start[name],
          len = inla_model$misc$configs$contents$length[name])

     return(indicies)
}

createTransform <- function(re_index, inla_model, constraint_point){
     ### return the A matrix which selects the sub-matrix corresponding to
     ### the named random effect while placing a constraint at a single
     ### point in the walk
     # Take the sub-matrix of Q that we want
     # number of columns - length of indices (for the random effect of interest)
     # number of rows - to match dimensions of Q (for multiplication)

     dimension <- nrow(inla_model$misc$configs$config[[1]]$Q)

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
