pairs.bce <- function(bce,              # BCE object
                      sample=1,
                      gap=0,
                      upper.panel = NA,
                      diag.panel = NA,
                      ...)
{
panel.cor <- function(x, y) text(x=mean(range(x)),y=mean(range(y)),
                                 labels=format(cor(x,y),digits=2))

panel.hist <- function(x,...) {
        usr <- par("usr")
        on.exit(par(usr))
        par(usr = c(usr[1:2], 0, 2))
        h <- hist(x, plot = FALSE)
        breaks <- h$breaks
        nB <- length(breaks)
        y <- h$counts
        y <- y/max(y)
        rect(breaks[-nB], 0, breaks[-1], y, col = "cyan")

    }
if (!is.null(upper.panel) && is.na(upper.panel))upper.panel <- panel.cor
if (!is.null(diag.panel) && is.na(diag.panel))diag.panel <- panel.hist

ifelse(is.matrix(bce$X),x <- t(bce$X),x <- t(bce$X[sample,,]))      # if only one sample, X is a matrix

labels <- colnames(x)
pairs(x, diag.panel =diag.panel, labels = labels,
      gap = gap, upper.panel = upper.panel,...)

}

