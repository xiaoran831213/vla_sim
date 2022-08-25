library(ggplot2)
## Gaussian Phenotype
#' demonstrate Gausian phenotypes
#'
#' Simulate a total of 6 groups: 2 dosage {-1, 0, 1}, and 2 environment {-1, 1}.
#'
#' @param N group size
#' @param a genotype main effect
#' @param b gxe effect
#' @param c environment main effect
#' @param d covariate effect
#' @param ve variance of noise
#' 
dmo <- function(N=5e2, a=1, b=1, c=0, d=1, ve=1, rsp="gsn", m0=0, rs=NULL, ...)
{
    set.seed(rs)
    x1 <- sample(rep(c(-1, 0, 1), N * 2)); x2 <- x1 * x1 # genotype, and squared
    u1 <- sample(rep(c(-1,    1), N * 3)); xu <- x1 * u1 # enviroment and GxE
    c1 <- scale(rnorm(N * 6)) * sd(x1)                   # covariate 1
    qc <- qr(c1)                                         # QRD of covariates
    
    ## simulate outcomes homoscedastic noise
    e1 <- replicate(6, scale(rnorm(N))) * sqrt(ve)
    y1 <- m0 + a * x1 + b * xu + c * u1 + d * c1 + c(e1)
    ## transform
    if(rsp == "bin")
    {
        y1 <- rbinom(length(y1), 1, 1 / (1 + exp(-y1)))
    }
    else
    {
        y1 <- 4 * (y1 - min(y1))/(max(y1) - min(y1)) # move y1 to [0, 4]
    }

    ## residual squared as tier 2 "phenotype"
    if(d != 0)
        r2 <- lm(y1 ~ x1 + c1)$resid^2
    else
        r2 <- lm(y1 ~ x1)$resid^2
    ## if(rsp == "bin")
    ##     r2 <- (r2 - min(r2))/(max(r2) - min(r2)) # move y1 to [0, 4]
    ## r2 <- r2 * 2
    
    u1 <- ifelse(u1 < 1, "Env.A", "Env.B")
    d <- data.frame(env=u1)
    d <- rbind(
        data.frame(x=x1, y=y1, e=u1, mtd="LM1"),
        data.frame(x=x1, y=r2, e=u1, mtd="LM2"))
    set.seed(NULL)

    ## pack and return
    cbind(d, ...)
}

## dmo(5e2, a=1, b=1, c=0, d=0, ve=2)
plt <- function(d, x2=1, out=NULL)
{
    gsm <- geom_smooth
    thm <- theme_bw() +
        theme(strip.text.x = element_text(size = 17, face = "bold"),
              strip.text.y = element_text(size = 17, face = "bold"),
              legend.position = "none",
              legend.title = element_blank(),
              legend.margin = margin(0, 0, 0, 0, "cm"),
              axis.text.x = element_text(size = 14, face = "bold"),
              axis.text.y = element_text(size = 14, face = "bold"),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              plot.margin= margin(0, 0, 0, 0, "cm")) # unit(x=c(0,0,0,0), units="cm"))

    g <- ggplot(d, aes(x=x, y=y, color=e))
    ## scatter points
    g <- g + geom_point(position=position_jitter(.2, seed=1), alpha=.2)
    g <- g + scale_color_manual(values=c("magenta", "orange"))
    ## regression line for each environment, for LM1 only
    g <- g + gsm(data=subset(d, mtd=="LM1"), method=lm, formula=y ~ x, size=1, linetype=2)
    ## marginal regression line
    g <- g + gsm(data=subset(d, mtd=="LM1"), method=lm, color="gray", formula=y ~ x)
    g <- g + gsm(data=subset(d, mtd=="LM2"), method=lm, color="green", formula=y ~ x)
    ## VLA regression line
    g <- g + gsm(data=subset(d, mtd=="LM2"), method=lm, color="red", formula=y ~ I(x^2))

    g <- g + facet_grid(rows=vars(tag), cols=vars(key, dsc), scales="free_y")
    g <- g + scale_x_continuous(breaks=c(-1, 0, 1), labels=c("AA", "Aa", "aa"))
    g <- g + thm
    g
}

tmp0 <- function()
{
    rsd <- 4
    dat <- rbind(
        dmo(100, a=2.0, b=2, c=0.0, d=1, ve=1.0, m0=0.0, rsp="gsn", rs=rsd, key="GXE", tag="GSN"),
        dmo(100, a=2.0, b=2, c=1.0, d=1, ve=1.0, m0=0.0, rsp="gsn", rs=rsd, key="MIX", tag="GSN"),
        dmo(100, a=2.0, b=0, c=1.0, d=1, ve=1.0, m0=0.0, rsp="gsn", rs=rsd, key="NUL", tag="GSN"),
        dmo(2e2, a=1.5, b=4, c=0.0, d=4, ve=0.0, m0=0.0, rsp="bin", rs=rsd, key="GXE", tag="BIN"),
        dmo(100, a=1.5, b=4, c=2.0, d=4, ve=0.0, m0=0.0, rsp="bin", rs=rsd, key="MIX", tag="BIN"),
        dmo(200, a=2.0, b=0, c=1.0, d=1, ve=0.1, m0=-.5, rsp="bin", rs=rsd, key="NUL", tag="BIN"))
    KEY <- c(NUL="G+E", MIX="G+E + GxE", GXE="GxE")
    TAG <- c(GSN="Quantitative", BIN="Case/Control")
    DSC <- c(LM1="Phenotype", LM2="Variation")
    dat <- within(dat,
    {
        key <- factor(KEY[key], KEY)
        tag <- factor(TAG[tag], TAG)
        dsc <- factor(DSC[mtd], DSC)
    })

    p <- plt(dat)
    ggsave("dem_mix.png", p, width=24, height=8, scale=.6)
}
