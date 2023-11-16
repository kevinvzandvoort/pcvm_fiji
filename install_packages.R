#' The model depends on several packages, which all need to be installed
#' - I like to use the pacman package to install other packages, it automatically
#'   installs a package if not yet installed, loads it if installed, and updates
#'   it if a new version is available
install.packages("pacman")
#' - The pacman package can now be used to manage all other dependencies
pacman::p_load(data.table, Rcpp, RcppArmadillo, inline, deSolve, rootSolve,
               readxl, magrittr, binom, qs, units, BayesianTools, ggplot)

#' There is one other package that needs to be installed: rootsolve
#' - While this package is available on CRAN, for a technical reason, we are
#'   using a slightly adjusted version of the package, which I've put in a 
#'   Github repository
#' 1. install the remotes package, which allows to install packages from Github
pacman::p_load(remotes)
#' 2. delete the old version of the package, if it was installed
remove.packages("rootsolve")
#' 3. install the Github version
remotes::install_github("https://github.com/kevinvzandvoort/rootSolve")
#' 4. if all went well, it can now successfully be loaded
pacman::p_load(rootSolve)