class.ind <- function (cl) 
{
    n <- length(cl)
    cl <- as.factor(cl)
    x <- matrix(0, n, length(levels(cl)))
    x[(1L:n) + n * (unclass(cl) - 1L)] <- 1
    dimnames(x) <- list(names(cl), levels(cl))
    x
}
