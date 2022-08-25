## test GxE models
library(simgen)
library(ggplot2)
library(car)
library(dglm)
library(simgen)
## source("lnk.R")
## source("mtd.R")

if(!exists("C17"))
    load("C17.RData")

psn <- function(x, ...) rpois(length(x), exp(x))             # exp link, possion
gma <- function(x, ...) rgamma(length(x), 1, scale=exp(x))   # exp link, gamma 
gsn <- function(x, ...) x
sgm <- function(x, ...) binomial()$linkinv(x)                # logit function (sigmoid)
bin <- function(x, ...) rbinom(length(x), 1, sgm(x))         # sgm link, binomial

#' t-value to p-value
t2p <- function(tvl, dof, sides=2)
{
    if(sides > 1)
        unname(2 * pt(abs(tvl), dof, low=FALSE))
    else
        unname(pt(tvl, dof, low=sides < 0))
}

#' Wald test
#'
#' @param mdl model object returned by R's glm(), lm(), or dglm()
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

#' simulate variance loci
#'
#' @param loc type of variance loci, "ADD" or "MLT".
#' @param rsp type of response variable, bin, psn, or gsu
sim <- function(N=5e2, a=0, b=1, c=0, d=1, m0=0, ve=0, loc="ADD", rsp='bin', rep=5e2, ...)
{
    
    cfg <- data.frame(N=N, a=a, b=b, c=c, d=d, m0=m0, ve=ve, loc=loc, rsp=rsp, rep=rep)
    dot <- list(...)
    set.seed(dot$rsd)
    ret <- replicate(rep,
    {
        ## genotype & environment
        x1 <- kgp(N, 1, MAF=.01); x0 <- mean(x1)
        u1 <- rnorm(N); u0 <- mean(u1); vu <- var(u1)
        u1 <- u1 - u0
        vx <- var(x1); # x1 <- x1 - x0
        xu <- x1 * u1
        c1 <- rnorm(N)
        ## response by loci type
        y1 <- m0 + rnorm(1, sd=sqrt(d)) * c1 + rnorm(1, sd=sqrt(a)) * x1
        if(loc == "ADD")
            y1 <- y1 + rnorm(1, sd=sqrt(b)) * xu +  rnorm(1, sd=sqrt(c)) * u1 + rnorm(N, sd=sqrt(ve))
        else
        {
            l0 <- log(ve + c * vu)
            ## l1 <- log(ve + c * vu + b * vu)
            l1 <- log(ve + c * vu + b * vx * vu)
            ld <- l1 - l0
            if(sample(c(TRUE, FALSE), 1))
                y1 <- y1 + rnorm(N, 0, exp((l1 - ld / 2 * x1) * 0.5))
            else
                y1 <- y1 + rnorm(N, 0, exp((l0 + ld / 2 * x1) * 0.5))
        }
        ## transforme the mean and sample
        y1 <- do.call(rsp, list(y1))
        fam <- c(gsn="gaussian", bin="binomial", psn="poisson")[rsp]
        pvl <- list()
        pvl$DL1 <- DLM(y1, x1, c1)      # double linear model
        pvl$VL0 <- VLA(y1, x1, c1)      # vla
        pvl$LV1 <- LVT(y1, x1, c1)      # Levene's Robust Test
        pvl$DR1 <- DRM(y1, x1, c1)      # deviation regression
        pvl$DG1 <- DG1(y1, x1, c1, fam) # double generalized linear model

        ## reject H_0: x1 is not a variance locus?
        unlist(pvl) < 0.05
    })
    set.seed(NULL)
    rep <- rowSums(!is.na(ret))      # actual repetition
    rej <- rowMeans(ret, na.rm=TRUE) # reject h0
    cbind(cfg, ..., mtd=names(rej), rej)
}
## gx2(1e3, m0=0, a=.5, gno=k05, env=b50, rsp=psn, rep=1e2, d=2.0, b=+2.0, c=0.5)[, -(1:9)]

VLA <- function(y1, x1, c1, ...)
{
    x2 <- x1^2
    m1 <- lm(y1 ~ x1 + c1 + 1)
    r1 <- resid(m1)
    r2 <- r1^2
    m2 <- lm(r2 ~ x2 + c1 + 1)

    z2 <- coef(m2)[2] / sqrt(vcov(m2)[2, 2])
    d2 <- length(y1) - 4
    t2p(z2, d2, 1)
}

## gamma regression with log link
LVT <- function(y1, x1, c1)
{
    m1 <- lm(y1 ~ c1 + 1)
    y2 <- resid(m1)
    x2 <- as.factor(x1)
    leveneTest(y2 ~ x2)[1, 3]
}

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

DG1 <- function(y1, x1, c1, fam=gaussian())
{
    m1 <- try(suppressWarnings(dglm(y1 ~ x1 + c1, ~ x1 + c1, fam)), silent=TRUE)
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
