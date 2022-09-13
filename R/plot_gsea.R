##' plot_gsea.R is based on enrichplot: https://github.com/YuLab-SMU/enrichplot
##' devloped by Guangchuang YU
##'

##' Plot GSEA result
##' @param x GSEA result object
##' @param disease_name the tested disease
##' @param sub_type sub_type of tested disease
##' @param cluster_id trained immunity sub-population
##' @param reg_dir regulated direction of TIGs
##' @param geneSetID same as geneSetID in GSEA
##' @param color color of line segments
##' @param color.line color of running enrichment score line
##' @param color.vline color of vertical line which indicating the
##' maximum/minimal running enrichment score
##' @return ggplot2 object
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 geom_linerange
##' @importFrom ggplot2 geom_line
##' @importFrom ggplot2 geom_vline
##' @importFrom ggplot2 geom_hline
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 xlim
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 ggplotGrob
##' @importFrom ggplot2 geom_segment
##' @importFrom ggplot2 ggplot_gtable
##' @importFrom ggplot2 ggplot_build
##' @importFrom ggplot2 ggtitle
##' @importFrom ggplot2 element_text
##' @importFrom ggplot2 rel
##' @importFrom aplot plot_list
##' @importFrom DOSE theme_dose
##' @export
##'
##'
GSEAplot_Basic <- function (x, geneSetID, disease_name, sub_type, cluster_id, reg_dir,
                            by = "all",
                            title = "Enrichment Testing for ",
                                 color='black', color.line="green",
                                 color.vline="#FA5860", ...){

  title = paste0(title, toupper(reg_dir), "-regulated TIGs of ", cluster_id)
  if(sub_type == "all")
    x_label = paste0("Ranked genes in ", disease_name)
  else
    x_label = paste0("Ranked genes in ", sub_type, " of ", disease_name)
  num_decimal = function(x) trimws(format(round(x, 3), nsmall = 3))
  p.adj.add = paste0("Adjusted P-value: ", num_decimal(x@result$p.adjust))
  gene.rank = paste0("Gene rank: ", x@result$rank)
  enr.score = paste0("Enrichment score: ", num_decimal(x@result$enrichmentScore))

  annos <- data.frame(
    x = c(length(x@geneList), length(x@geneList), length(x@geneList)),
    y = c(x@result$enrichmentScore+0.1, (x@result$enrichmentScore+0.1)/3*2,
          (x@result$enrichmentScore+0.1)/3),
    text = c(p.adj.add, gene.rank, enr.score)
  )

  by <- match.arg(by, c("runningScore", "preranked", "all"))
  gsdata <- gsInfo(x, geneSetID)
  p <- ggplot(gsdata, aes_(x = ~x)) +
    theme_dose() + xlab(x_label)
  if (by == "runningScore" || by == "all") {
    p.res <- p + geom_linerange(aes_(ymin=~ymin, ymax=~ymax), color=color)
    p.res <- p.res + geom_line(aes_(y = ~runningScore), color=color.line,
                               size=1)
    enrichmentScore <- x@result[geneSetID, "enrichmentScore"]
    es.df <- data.frame(es = which.min(abs(p$data$runningScore - enrichmentScore)))
    p.res <- p.res + geom_vline(data = es.df, aes_(xintercept = ~es),
                                colour = color.vline, linetype = "dashed")
    p.res <- p.res + ylab("Enrichment Score")
    p.res <- p.res + geom_hline(yintercept = 0)+
      geom_text(data = annos, aes(x = x, y = y, label = text),
                vjust = "inward", hjust = "inward")
  }
  if (by == "preranked" || by == "all") {
    df2 <- data.frame(x = which(p$data$position == 1))
    df2$y <- p$data$geneList[df2$x]
    p.pos <- p + geom_segment(data=df2, aes_(x=~x, xend=~x, y=~y, yend=0),
                              color=color)
    p.pos <- p.pos + ylab("Ranked Values in Input Dataset") +
      xlim(0, length(p$data$geneList))
  }
  if (by == "runningScore")
    return(p.res + ggtitle(title))
  if (by == "preranked")
    return(p.pos + ggtitle(title))

  p.pos <- p.pos + xlab(NULL) + theme(axis.text.x = element_blank(),
                                      axis.ticks.x = element_blank())
  p.pos <- p.pos + ggtitle(title) +
    theme(plot.title=element_text(hjust=0.5, size=15))
  #plot_list(gglist =  list(p.pos, p.res), ncol=1)

  aplot::gglist(gglist = list(p.pos, p.res), ncol=1)
}


##' extract gsea result of selected geneSet
##'
##'
##' @title gsInfo
##' @param object gseaResult object
##' @param geneSetID gene set ID
##' @return data.frame
##' @export
gsInfo <- function(object, geneSetID) {
  geneList <- object@geneList

  if (is.numeric(geneSetID))
    geneSetID <- object@result[geneSetID, "ID"]

  geneSet <- object@geneSets[[geneSetID]]
  exponent <- object@params[["exponent"]]
  df <- gseaScores(geneList, geneSet, exponent, fortify=TRUE)
  df$ymin <- 0
  df$ymax <- 0
  pos <- df$position == 1
  h <- diff(range(df$runningScore))/20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$geneList <- geneList
  if (length(object@gene2Symbol) == 0) {
    df$gene <- names(geneList)
  } else {
    df$gene <- object@gene2Symbol[names(geneList)]
  }

  df$Description <- object@result[geneSetID, "Description"]
  return(df)
}


gseaScores <- getFromNamespace("gseaScores", "DOSE")


##' GSEA plot that mimic the plot generated by broad institute's GSEA software
##'
##'
##' @title GSEAplot
##' @param x gseaResult object
##' @param geneSetID gene set ID
##' @param sub_type tested sub_type of disease
##' @param disease_name tested disease
##' @param cluster_id trained immunity sub-population
##' @param reg_dir regulated direction of TIGs
##' @param title plot title
##' @param color color of running enrichment score line
##' @param base_size base font size
##' @param rel_heights relative heights of subplots
##' @param subplots which subplots to be displayed
##' @param pvalue_table whether add pvalue table
##' @param ES_geom geom for plotting running enrichment score,
##' one of 'line' or 'dot'
##' @return plot
##' @export
##' @importFrom ggplot2 theme_classic
##' @importFrom ggplot2 element_line
##' @importFrom ggplot2 element_text
##' @importFrom ggplot2 element_blank
##' @importFrom ggplot2 element_rect
##' @importFrom ggplot2 scale_x_continuous
##' @importFrom ggplot2 scale_y_continuous
##' @importFrom ggplot2 scale_color_manual
##' @importFrom ggplot2 theme_void
##' @importFrom ggplot2 geom_rect
##' @importFrom ggplot2 margin
##' @importFrom ggplot2 annotation_custom
##' @importFrom stats quantile
##' @importFrom RColorBrewer brewer.pal
##' @importFrom gridExtra tableGrob
GSEAplot <- function(x, geneSetID, sub_type, disease_name,
                      cluster_id, reg_dir,
                      title = "Enrichment Testing for ",
                      color="green", base_size = 11,
                      rel_heights=c(1.5, .5, 1), subplots = 1:3,
                      pvalue_table = FALSE, ES_geom="line") {
  ES_geom <- match.arg(ES_geom, c("line", "dot"))

  geneList <- position <- NULL ## to satisfy codetool

  if (length(geneSetID) == 1) {
    gsdata <- gsInfo(x, geneSetID)
  } else {
    gsdata <- do.call(rbind, lapply(geneSetID, gsInfo, object = x))
  }

  p <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) +
    theme_classic(base_size) +
    theme(panel.grid.major = element_line(colour = "grey92"),
          panel.grid.minor = element_line(colour = "grey92"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()) +
    scale_x_continuous(expand=c(0,0))

  if (ES_geom == "line") {
    es_layer <- geom_line(aes_(y = ~runningScore, color= ~Description),
                          size=1)
  } else {
    es_layer <- geom_point(aes_(y = ~runningScore, color= ~Description),
                           size=1, data = subset(gsdata, position == 1))
  }

  p.res <- p + es_layer +
    theme(legend.position = c(.8, .8), legend.title = element_blank(),
          legend.background = element_rect(fill = "transparent"))

  p.res <- p.res + ylab("Enrichment Score") +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.line.x=element_blank(),
          plot.margin=margin(t=.2, r = .2, b=0, l=.2, unit="cm"))

  i <- 0
  for (term in unique(gsdata$Description)) {
    idx <- which(gsdata$ymin != 0 & gsdata$Description == term)
    gsdata[idx, "ymin"] <- i
    gsdata[idx, "ymax"] <- i + 1
    i <- i + 1
  }
  p2 <- ggplot(gsdata, aes_(x = ~x)) +
    geom_linerange(aes_(ymin=~ymin, ymax=~ymax, color=~Description)) +
    xlab(NULL) + ylab(NULL) + theme_classic(base_size) +
    theme(legend.position = "none",
          plot.margin = margin(t=-.1, b=0,unit="cm"),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line.x = element_blank()) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0))

  if (length(geneSetID) == 1) {

    v <- seq(1, sum(gsdata$position), length.out=9)
    inv <- findInterval(rev(cumsum(gsdata$position)), v)
    if (min(inv) == 0) inv <- inv + 1

    col <- c(rev(brewer.pal(5, "Blues")), brewer.pal(5, "Reds"))

    ymin <- min(p2$data$ymin)
    yy <- max(p2$data$ymax - p2$data$ymin) * .3
    xmin <- which(!duplicated(inv))
    xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
    d <- data.frame(ymin = ymin, ymax = yy,
                    xmin = xmin,
                    xmax = xmax,
                    col = col[unique(inv)])
    p2 <- p2 + geom_rect(
      aes_(xmin=~xmin,
           xmax=~xmax,
           ymin=~ymin,
           ymax=~ymax,
           fill=~I(col)),
      data=d,
      alpha=.9,
      inherit.aes=FALSE)
  }


  if(sub_type == "all")
    x_label = paste0("Ranked genes in ", disease_name)
  else
    x_label = paste0("Ranked genes in ", sub_type, " of ", disease_name)


  df2 <- p$data #data.frame(x = which(p$data$position == 1))
  df2$y <- p$data$geneList[df2$x]
  p.pos <- p + geom_segment(data=df2, aes_(x=~x, xend=~x, y=~y, yend=0),
                            color="grey")
  p.pos <- p.pos + ylab("Ranked Values in Input Dataset") +
    xlab(x_label) +
    theme(plot.margin=margin(t = -.1, r = .2, b=.2, l=.2, unit="cm"))

  title = paste0(title, toupper(reg_dir), "-regulated TIGs of ", cluster_id)

  if (!is.null(title) && !is.na(title) && title != "")
    p.res <- p.res + ggtitle(title)

  if (length(color) == length(geneSetID)) {
    p.res <- p.res + scale_color_manual(values=color)
    if (length(color) == 1) {
      p.res <- p.res + theme(legend.position = "none")
      p2 <- p2 + scale_color_manual(values = "black")
    } else {
      p2 <- p2 + scale_color_manual(values = color)
    }
  }

  if (pvalue_table) {
    pd <- x[geneSetID, c("enrichmentScore", "rank", "p.adjust")]
    rownames(pd) <- ""

    for (i in seq_len(ncol(pd))) {
      pd[, i] <- format(pd[, i], digits = 3)
    }
    tp <- tableGrob(pd)

    p.res <- p.res + theme(legend.position = "none") +
      annotation_custom(tp,
                        xmin = quantile(p.res$data$x, .5),
                        xmax = quantile(p.res$data$x, .95),
                        ymin = quantile(p.res$data$runningScore, .75),
                        ymax = quantile(p.res$data$runningScore, .9))
  }


  plotlist <- list(p.res, p2, p.pos)[subplots]
  n <- length(plotlist)
  plotlist[[n]] <- plotlist[[n]] +
    theme(axis.line.x = element_line(),
          axis.ticks.x=element_line(),
          axis.text.x = element_text())

  if (length(subplots) == 1)
    return(plotlist[[1]] + theme(plot.margin=margin(t=.2, r = .2, b=.2,
                                                    l=.2, unit="cm")))


  if (length(rel_heights) > length(subplots))
    rel_heights <- rel_heights[subplots]

  aplot::gglist(gglist = plotlist, ncol=1, heights=rel_heights)
}

