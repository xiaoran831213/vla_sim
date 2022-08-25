library(car)
#' Levene's Robust Test
#'
#' @param y1 numeric vector of phenotype values
#' @param x1 integer vector of allale dosages in {0, 1, 2}
#' @param c1 numeric vector of covariates
#' @return p-value of putative GxE variant Test
LVT <- function(y1, x1, c1)
{
    m1 <- lm(y1 ~ c1 + 1)
    y2 <- resid(m1)
    x2 <- as.factor(x1)
    leveneTest(y2 ~ x2)[1, 3]
}

#' Deviation Regression Model
#'
#' @param y1 numeric vector of phenotype values
#' @param x1 integer vector of allale dosages in {0, 1, 2}
#' @param c1 numeric vector of covariates
#' @return p-value of putative GxE variant Test
DRM <- function(y1, x1, c1)
{
    m1 <- lm(y1 ~ c1 + 1)
    r1 <- resid(m1)
    y2 <- numeric(length(r1))
    for(.x in unique(x1))
    {
        i <- which(x1 == .x)
        y2[i] <- abs(r1[i] - median(r1[i]))
    }
    m2 <- lm(y2 ~ x1 + c1 + 1)
    z2 <- coef(m2)[2] / sqrt(vcov(m2)[2, 2])
    d2 <- length(y1) - 5
    t2p(z2, d2, 2)
}

#' Double Linear Model
#'
#' @param y1 numeric vector of phenotype values
#' @param x1 integer vector of allale dosages in {0, 1, 2}
#' @param c1 numeric vector of covariates
#' @return p-value of putative GxE variant Test
DLM <- function(y1, x1, c1)
{
    m1 <- lm(y1 ~ x1 + c1 + 1)
    r1 <- resid(m1)
    r2 <- r1^2
    m2 <- lm(r2 ~ x1 + c1 + 1)
    z2 <- coef(m2)[2] / sqrt(vcov(m2)[2, 2])
    d2 <- length(y1) - 5
    t2p(z2, d2, 2)
}

#' Double Generalized Linear Model
#'
#' DGLM with a linear dosage term.
#' 
#' @param y1 numeric vector of phenotype values
#' @param x1 integer vector of allale dosages in {0, 1, 2}
#' @param c1 numeric vector of covariates
#' @param fam family of the GLM (i.e., Gaussian, Bionomial, etc.)
#' @return p-value of putative GxE variant Test
DG1 <- function(y1, x1, c1, fam=gaussian())
{
    m1 <- try(suppressWarnings(dglm(y1 ~ x1 + c1, ~ x1 + c1, fam)))
    if('try-error' %in% class(m1))
    {
        NA
    }
    else
    {
        s1 <- summary(m1)$dispersion.summary
        wld(s1, 2:2)
    }
}

#' Double Generalized Linear Model
#'
#' DGLM with a linear and quadratic dosage term.
#' 
#' @param y1 numeric vector of phenotype values
#' @param x1 integer vector of allale dosages in {0, 1, 2}
#' @param c1 numeric vector of covariates
#' @param fam family of the GLM (i.e., Gaussian, Bionomial, etc.)
#' @return p-value of putative GxE variant Test
DG2 <- function(y1, x1, c1, fam=gaussian())
{
    x2 <- x1^2
    m1 <- try(suppressWarnings(dglm(y1 ~ x1 + c1, ~ x1 + x2 + c1, fam)))
    if('try-error' %in% class(m1))
    {
        NA
    }
    else
    {
        s1 <- summary(m1)$dispersion.summary
        if(nrow(coef(s1)) < 4)
            NA
        else
        {
            wld(s1, 2:3)
        }
    }
}

wld <- function(mdl, idx=NULL)
{
    x <- as.matrix(coef(mdl))[, 1] # coefficients
    V <- vcov(mdl)                 # variance covariance structure
    R <- diag(nrow(V))             # combination matrix

    ## part of the coefficients
    if(!is.null(idx))
        R <- R[idx, , drop=FALSE]
    ## chisqure statistics
    r <- R %*% x
    ksq <- crossprod(r, solve(R %*% tcrossprod(V, R), r))
    pvl <- pchisq(ksq, nrow(R), low=FALSE)
    pvl
}

## t-value to p-value
t2p <- function(tvl, dof, sds=2)
{
    if(sds > 1)
        unname(2 * pt(abs(tvl), dof, low=FALSE))
    else
        unname(pt(tvl, dof, low=sds < 0))
}


mtq <- function(Z, C, D=NULL, tol.egv=NULL, ...)
{
    tol.egv <- tol.egv %||%  sqrt(.Machine$double.eps)
    D <- D %||% 1
    dim(Z) <- NULL

    C <- eigen(C, TRUE, TRUE)$values
    D <- eigen(D, TRUE, TRUE)$values
    
    d <- kronecker(D, C)
    dim(d) <- NULL

    ## . <- d > max(d) * tol.egv
    ## if(!all(.))
    ##     d <- d[  .]
    L <- length(d)              # effective number of eigen

    ## P <- imhof(sum(Z^2), d, delta=rep(0, L))$Qq
    P <- davies(sum(Z^2), lambda = d)$Qq
    list(d=d, L=L, P=P)
}

pow <- function(rpt)
{
    rpt <- subset(rpt, se=-itr)
    grp <- subset(rpt, se=-c(pvl, egv))
    rpt <- by(rpt, grp, function(g)
    {
        cfg <- subset(g, se=-c(pvl, egv))[1, ]
        pow <- with(g, mean(pvl <= 0.05))
        egv <- with(g, mean(egv))
        cbind(cfg, pow=pow, egv=egv, rep=nrow(g))
    })
    rpt <- do.call(rbind, rpt)
    rpt
}


## balance data {x} grouped by {f}, scale the group size by {s}
blc <- function(f, x, s=1)
{
    x <- split(x, f)
    m <- max(sapply(x, nrow)) * s
    
    x <- lapply(x, function(g)
    {
        rbind(g, g[sample.int(nrow(g), m - nrow(g), TRUE), ])
    })
    do.call(rbind, x)
}

## mix data {x} grouped by {f}
mix <- function(x, f=1, r=1)
{
    f <- rep(f, len=NROW(x))
    x <- as.matrix(x)
    for(. in unique(f))
    {
        m <- f == .
        n <- sum(m)
        w <- matrix(sample(c(0, 1), n * n, TRUE), n, n)
        x[m, ] <- w %*% x[m, ] / rowSums(w)
    }
    x
}
