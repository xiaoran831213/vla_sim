library(ggplot2)

#' collect simulation results
get.rpt <- function(rpt)
{
    rds <- paste0(rpt, ".rds")
    if(file.exists(rds))
        rpt <- readRDS(rds)
    else
    {
        ## gather reports
        rpt <- do.call(rbind, lapply(dir(rpt, "[.]rpt$", full.names = TRUE), readRDS))
        
        ## get average
        GRP <- c('N', 'rsp', 'key', 'tag', 'sim', "mtd")
        VAL <- c("rep", "rej")
        rpt <- by(rpt[, c(GRP, VAL)], rpt[, GRP], function(g)
        {
            ttl <- sum(g[, "rep"], na.rm=TRUE)
            avg <- sum(g[, "rep"] * g[, "rej"], na.rm=TRUE) / ttl
            within(g[1, ], {rep <- ttl; rej=avg})
        })
        rpt <- do.call(rbind, rpt)

        ## get net positive rate
        GRP <- c('N', 'rsp', 'key', 'tag', "mtd")
        rpt <- by(rpt, rpt[, GRP], function(g)
        {
            npr <- g[g[, "sim"] == "TPR", "rej"] - g[g[, "sim"] == "FPR", "rej"]
            rbind(g, within(g[1, ], {sim <- "NPR"; rej <- npr}))
        })
        rpt <- do.call(rbind, rpt)
        
        ## pack and return
        rownames(rpt) <- NULL
        saveRDS(rpt, rds, compress="xz")
    }
    invisible(rpt)
}

thm <- theme_bw() +
    theme(strip.text.x = element_text(size = 11, face = "bold"),
          strip.text.y = element_text(size = 11, face = "bold"),
          legend.text = element_text(size = 11, face = "bold"),
          legend.position = "bottom",
          legend.title = element_blank(),
          legend.margin = margin(0, 0, 0, 0, "cm"),
          legend.spacing.y = unit(0, "cm"),
          legend.box.margin = margin(0, 0, 0, 0, "cm"),
          legend.background = element_rect(fill = 'transparent'),
          legend.key = element_rect(fill = 'transparent'),
          axis.text.x = element_text(size = 11, face = "bold"),
          axis.title.x = element_text(size = 11, face = "bold"),
          axis.title.y = element_blank(),
          plot.margin= margin(0, 0, 0, 0, "cm"))

SIM <- c(TPR="True Positive Rate", FPR="False Positive Rate", NPR="Net Positive Rate")
KEY <- c(GSN="Gaussian", P00="36% zeros", P02="92% zeros", B50="50% cases rate", B07="7.6% case rate")
TAG <- c(NUL="No GxE", GXE="Pure GxE", FUL="Environment and GxE", MLT="Multiplicative")
MTD <- c(DL1="DLM", VL0="VLA", LV1="LVT", DR1="DRM", DG1="DGLM")
CLR <- c(DL1="green", VL0="red", LV1="blue", DR1="orange", DG1="cyan")
names(CLR) <- MTD[names(CLR)]

plt <- function(rpt, out=NULL, type=1, ...)
{
    ## read plotting data
    if(is.null(out))
        out <- paste0(rpt, ".png")
    rpt <- get.rpt(rpt)

    ## type of plot
    if(type == 2)
        rpt <- subset(rpt, sim %in% c("FPR", "NPR") & tag %in% c("FUL", "GXE"))
    dat <- within(subset(rpt, mtd %in% names(MTD)),
    {
        l2N <- log(N, 2)
        key <- factor(KEY[key], KEY)
        mtd <- MTD[mtd]
        tag <- factor(TAG[tag], TAG)
        sim <- factor(SIM[sim], SIM)
    })

    ## produce the figure
    g <- ggplot(dat, aes(x=l2N, y=rej))
    g <- g + geom_line(aes(colour=mtd), size=1.0)
    g <- g + geom_hline(yintercept=0.05, linetype="dashed", size=0.5, color="red")
    if(type == 2)
    {
        if(length(unique(rpt$tag)) > 1)
            g <- g + facet_grid(rows=vars(key), cols=vars(tag, sim))
        else
            g <- g + facet_grid(rows=vars(key), cols=vars(sim))
    }
    else
    {
        if(length(unique(rpt$key)) > 1)
            g <- g + facet_grid(rows=vars(key, sim), cols=vars(tag))
        else
            g <- g + facet_grid(rows=vars(sim), cols=vars(tag))
    }
    g <- g + coord_cartesian(ylim=c(0, 1), expand=TRUE)
    g <- g + xlab("Sample Size in power of 2")
    g <- g + scale_color_manual(values=CLR)
    g <- g + guides(color=guide_legend(nrow=1, byrow=TRUE))
    g <- g + thm
    ## print to file
    options(bitmapType = 'cairo') # , device = 'png')
    if(type == 2)
    {
        nc <- length(unique(rpt$tag)) * length(unique(rpt$sim))
        nr <- length(unique(rpt$key))
        eh <- 0.4
    }
    else
    {
        nc <- length(unique(rpt$tag))
        nr <- length(unique(rpt$key)) * length(unique(rpt$sim))
        eh <- 0.2
    }
    ggsave(out, g, width=4 * nc, height=3.3 * nr + eh, scale=.6)
    invisible(rpt)
}
