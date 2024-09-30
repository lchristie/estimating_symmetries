
K_hat_sym_diff <- function ( TestRejections ) {
  ## This function tracks the error of the confidence region
  ## The true Test rejection vector is (1,1,1,-1,-1,-1)
  ## This can be sped up as a direct computation rather than conditionals, but this is fast enough for now
  
  mathcal_E <- 0
  if ( sum( TestRejections[1:3] ) < 3 ) {
    ## This case covers the error of $\langle R_h, R_v \rangle, $R_h$, $R_v$ and $R_\pi$
    mathcal_E <- mathcal_E + 1 + (3 - sum( TestRejections[1:3] ))
  }
  if ( TestRejections[6] * TestRejections[3] == 1) {
    ## This is the error of $\langle R_{\pi / 2} \rangle$
    mathcal_E <- mathcal_E + 1
  }
  if ( sum( TestRejections[3:5] ) == 3) {
    ## This is the error of $\langle R_/, R_\ \rangle$, $R_/$, $R_\$, and $R_\pi$
    mathcal_E <- mathcal_E + 1 + sum( TestRejections[4:5] )
  }
  if ( sum( TestRejections ) == 6 ) {
    ## This case covers the error of $D_4$
    mathcal_E <- mathcal_E + 1
  }
  return( mathcal_E )
}

G_hat_from_TR <- function ( TestRejections ) {
  ## This returns a vector indicating inclusion from the lattice base, G_hat is the join of these elements
  if ( sum(TestRejections) <= 1 ) {
    ## This case is the easiest, either returning $I$ or the singular accepted group (except R_270)
    output <- TestRejections
    output[6] <- 0
    return( output )
  }
  if ( TestRejections[3] == 0 ) {
    ## This case has only order two maximal subgroups
    possibilities = which( TestRejections[1:5] == 1, arr.ind = TRUE )
    output <- rep(0, 6)
    ind <- true_sample( possibilities )
    output[ind] <- 1
    return( output )
  } 
  if ( sum(TestRejections) == 2) {
    ## Recall that one of the accepted symmetries is G_3
    if ( TestRejections[6] == 1 ) {
      return( TestRejections )
    } else {
      possibilities = which( TestRejections == 1, arr.ind = TRUE )
      output <- rep(0, 6)
      ind <- true_sample( possibilities )
      output[ind] <- 1
      return( output )
    }
  }
  if ( (sum( TestRejections * c(1,1,1,0,0,0) ) == 3) * (sum(TestRejections) == 3) ) {
    return( TestRejections )
  }
  if ( (sum( TestRejections * c(0,0,1,1,1,0) ) == 3) * (sum(TestRejections) == 3) ) {
    return( TestRejections )
  }
  if ( (sum( TestRejections * c(1,1,1,1,1,0) ) == 5) * (sum(TestRejections) == 5) ) {
    if ( sample(2,1) == 1 ) {
      return( c(1,1,1,0,0,0) )
    } else {
      return( c(0,0,1,1,1,0) )
    }
  }
  if ( sum(TestRejections) == 6 ) {
    ## This is the case where $G_hat = D_4$
    return(TestRejections)
  }
  if ( sum(TestRejections) == 3 ) {
    ## Recall TestRejections[3] == 1, but cannot form a large maximal group
    possibilities = which( TestRejections == 1, arr.ind = TRUE )
    output <- rep(0, 6)
    ind <- true_sample( possibilities )
    output[ind] <- 1
    return( output )
  }
  if ( sum(TestRejections[1:3]) == 3) {
    possibilities = which( TestRejections[3:6] == 1, arr.ind = TRUE ) + 2
    output <- rep(0, 6)
    ind <- true_sample( possibilities )
    output[ind] <- 1
    if (output[3] == 1) {
      output[1] <- 1
      output[2] <- 1
    }
    return( output )
  }
  if ( sum(TestRejections[3:5]) == 3) {
    possibilities = which( TestRejections[c(1,2,3,6)] == 1, arr.ind = TRUE )
    output <- rep(0, 6)
    ind <- true_sample( possibilities )
    if (ind == 4) {ind <- 6}
    output[ind] <- 1
    if (output[3] == 1) {
      output[4] <- 1
      output[5] <- 1
    }
    return( output )
  }
  return( rep(-1, 6) ) ## This is an error code
}