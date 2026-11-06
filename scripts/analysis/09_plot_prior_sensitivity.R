#!/usr/bin/env Rscript

base <- file.path(
  "outputs",
  "bayesian",
  "prior_sensitivity"
)

input <- file.path(
  base, "tables",
  "Table_OCT2026_prior_sensitivity.tsv"
)

output <- file.path(base, "figures")
dir.create(output, recursive=TRUE, showWarnings=FALSE)

d <- read.delim(input, check.names=FALSE)

required <- c(
  "model", "prior_sd", "lineage",
  "OR_median", "CrI_low", "CrI_high"
)

stopifnot(all(required %in% names(d)))
stopifnot(nrow(d) == 60)
stopifnot(all(d$CrI_low > 0))
stopifnot(all(d$CrI_high >= d$CrI_low))

models <- c(
  "B1_primary",
  "B2_primary",
  "B3_primary",
  "B3_geo"
)

titles <- c(
  "A. Any genomic resistance",
  "B. MDR or higher",
  "C. Advanced resistance",
  "D. Advanced resistance + geography"
)

draw_plot <- function() {

  oldpar <- par(no.readonly=TRUE)
  on.exit(par(oldpar))

  par(
    mfrow=c(2,2),
    mar=c(4.5,5.5,3,1.5),
    oma=c(0,0,2,0)
  )

  for (i in seq_along(models)) {

    x <- d[d$model == models[i], ]

    lineages <- c("L6","L5","L3","L2","L1")
    priors <- sort(unique(x$prior_sd))

    positions <- seq_along(lineages)
    offsets <- seq(-0.22,0.22,length.out=length(priors))

    xrange <- range(
      c(x$CrI_low,x$CrI_high,1),
      finite=TRUE
    )

    plot(
      NA,
      xlim=xrange,
      ylim=c(0.5,5.5),
      log="x",
      yaxt="n",
      xlab="Posterior odds ratio (log scale)",
      ylab="",
      main=titles[i],
      bty="l"
    )

    axis(2,at=positions,labels=lineages,las=1)
    abline(v=1,lty=2,col="grey50")

    for (j in seq_along(priors)) {

      z <- x[x$prior_sd == priors[j], ]

      for (k in seq_along(lineages)) {

        r <- z[z$lineage == lineages[k], ]

        if (nrow(r) != 1) next

        yy <- positions[k] + offsets[j]

        segments(
          r$CrI_low, yy,
          r$CrI_high, yy,
          lwd=1.5,
          col=j
        )

        points(
          r$OR_median,
          yy,
          pch=16,
          cex=0.8,
          col=j
        )
      }
    }

    legend(
      "topright",
      legend=paste("Prior SD =",priors),
      col=seq_along(priors),
      pch=16,
      lty=1,
      cex=0.7,
      bty="n"
    )
  }

  mtext(
    "Bayesian prior sensitivity: lineage contrasts relative to L4",
    outer=TRUE,
    cex=1.15,
    font=2
  )
}

pdf(
  file.path(output,"Figure_OCT2026_prior_sensitivity.pdf"),
  width=12,
  height=9
)
draw_plot()
dev.off()

png(
  file.path(output,"Figure_OCT2026_prior_sensitivity.png"),
  width=3600,
  height=2700,
  res=300
)
draw_plot()
dev.off()

cat("Forest plots generated successfully.\n")
