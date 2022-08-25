## normalization
orq <- function(x)
{
    qnorm((rank(x) - 0.5) / length(x))
}

## GLM family object
sgm <- function(x, ...)
{
    if(length(x) > 0) binomial()$linkinv(x) else numeric()
}

## Softplut
sfp <- function(x, ...) log(1 + exp(x))

bin <- function(x, n=1, ...)
{
    n <- sample(n, length(x), TRUE)
    rbinom(length(x), n, sgm(x)) / n
}
bn1 <- function(x, ...) bin(x, 1, ...)
bn2 <- function(x, ...) bin(x, 2, ...)
bn3 <- function(x, ...) bin(x, 3, ...)
bn9 <- function(x, ...) bin(x, 9, ...)
dct <- function(x, ...)
{
    rbinom(length(x), 1, sgm(x)) - 1
}

psn <- function(x, ...)
{
    rpois(length(x), exp(x))
}
gma <- function(x, ...)
{
    rgamma(length(x), 1, scale=exp(x))
}
gsn <- function(x, ...) x



## data generating functions
bng <- function(N, L=1, MAF=NULL, drop=1, std=2, ...)
{
    MAF <- MAF %||% 0.05
    gmx <- rbinom(N * L, 2L, MAF)
    if(std > 0)
        gmx <- gmx - 2 * MAF
    if(std > 1)
        gmx <- gmx / sqrt(2 * MAF * (1 - MAF))
    if(L > 1 || !drop)
        dim(gmx) <- c(N, L)
    gmx
}

g01 <- function(N) bng(N, 1, MAF=.01, std=0)
g02 <- function(N) bng(N, 1, MAF=.02, std=0)
g05 <- function(N) bng(N, 1, MAF=.05, std=0)
g10 <- function(N) bng(N, 1, MAF=.10, std=0)
g20 <- function(N) bng(N, 1, MAF=.20, std=0)
g50 <- function(N) bng(N, 1, MAF=.50, std=0)
gu0 <- function(N) runif(N, -.5, .5)
k01 <- function(N) kgp(N, 1, MAF=.01)
k05 <- function(N) kgp(N, 1, MAF=.05)
k10 <- function(N) kgp(N, 1, MAF=.10)
k20 <- function(N) kgp(N, 1, MAF=.20)
gs0 <- function(N) rnorm(N, 0, 1.0) # normal
gs2 <- function(N) rnorm(N, 0, 0.2) # normal (almost) within (-1, +1)
gs5 <- function(N) rnorm(N, 0, 0.5) # normal (almost) within (-1, +1)

b01 <- function(N) rbinom(N, 1, .01) / sqrt(.01 * .99)
b02 <- function(N) rbinom(N, 1, .02) / sqrt(.02 * .98)
b05 <- function(N) rbinom(N, 1, .05) / sqrt(.05 * .95)
b10 <- function(N) rbinom(N, 1, .10) / sqrt(.10 * .90)
b20 <- function(N) rbinom(N, 1, .20) / sqrt(.20 * .80)
b50 <- function(N) rbinom(N, 1, .50) / sqrt(.50 * .50)
d01 <- function(N) (rbinom(N, 1, .01) - .01) / sqrt(.01 * .99)
d02 <- function(N) (rbinom(N, 1, .02) - .02) / sqrt(.02 * .98)
d05 <- function(N) (rbinom(N, 1, .05) - .05) / sqrt(.05 * .95)
d10 <- function(N) (rbinom(N, 1, .10) - .10) / sqrt(.10 * .90)
d20 <- function(N) (rbinom(N, 1, .20) - .20) / sqrt(.20 * .80)
d50 <- function(N) (rbinom(N, 1, .50) - .50) / sqrt(.50 * .50)
u01 <- function(N) runif(N, 0, 1)

dsp <- function(N) NULL
sn3 <- function(N) sample(0:2, N, TRUE) - 1
