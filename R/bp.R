# Arguments
# ---------
# files: character vector of (paths to) BP files
# legend: show a legend (or label systolic and diastolic directly)?
# smooth: show loess smooths?
# span: for loess smooths; if there is more than one, the span is
#       selected by 10-fold cross-validation
# confint: show pointwise confidence envelope(s) around smooths?
# levels: for confidence envelope(s)
# alpha: transparency for confidence envelopes (cumulative)
# vlines: list of lists for vertical lines to draw
#         ("at" component of each sub-list), labels for
#         the lines ("text" component), and (optionally, "at" component)
#         for vertical adjustment of text
#
bpPlot(c("2025-03.txt", "2025-04.txt", "2025-05.txt",
"2025-06.txt", "2025-07.txt", "2025-08.txt"),
span=0.2, legend=TRUE,
vlines=list(list(at="2025-07-08",
                 text="pause \u2192\nLenvima "),
            list(at="2025-07-15",
                 text="resume\u2192\nLenvima ",
                 adjust=15)
            #             list(at="2025-06-07",
            #                  text="stop\u2192\nRitalin "),
            #             list(at="2025-06-20",
            #                  text="restart\u2192\nRitalin ")
)
)


#' @export
bpPlot <- function(files, legend=FALSE,
                   smooth = nrow(BP) >= 15, span=seq(0.1, 0.9, by = 0.1),
                   confint=TRUE, level=c(0.5, 0.9), alpha=0.25,
                   vlines=NULL){

  BP <- read.table(files[1], header=TRUE, fill=TRUE)
  for (file in files[-1]){
    bp <- read.table(file, header=TRUE, fill=TRUE)
    BP <- rbind(BP, bp)
  }
  BP$date <- as.Date(BP$date)
  BP <- na.omit(BP)

  with(BP, {

    plot(range(date), range(c(sys, dia)), type="n",
         xlab="Date", ylab="Blood Pressure", xaxt="n", yaxt="n")
    lines(date, sys, lty=1, col="blue", type="b", pch=16)
    lines(date, dia, lty=2, col="magenta", type="b", pch=17)
    ticks <- axis.Date(1, date, format="%Y-%m-%d")
    axis(2, at=seq(60, 160, by=10), las=1)
    abline(h=seq(60, 160, by=10), col="gray30", lty="dotted")
    abline(v=ticks, col="gray30", lty="dotted")
    abline(h=120, col="blue", lty=3)
    abline(h=80, col="magenta", lty=3)

    if (smooth){

       if ((nmods <- length(span)) > 1){
         source("cv.loess.R")
         mods <- vector(nmods, mode="list")
         for (i in 1:nmods){
           cmd <- paste0('loess(sys ~ as.numeric(date), data=BP, span=', span[i],
                         ', degree=1, family="symmetric")')
           mods[[i]] <- eval(parse(text=cmd))
         }
         cvs <- sapply(mods, function (m) cv::cv(m, k="n")$"CV crit")
         cvs <-  cbind(span, cvs)
         colnames(cvs) <- c("span", "CV MSE")
         print(cvs)
         span <- span[which.min(cvs[, "CV MSE"])]
       }

      main <- paste0("robust local linear regression (span = ",
                     round(span, 2), ")")
      lo.sys <- loess(sys ~ as.numeric(date), data=BP, span=span,
                      degree=1, family="symmetric")
      lo.dia <- loess(dia ~ as.numeric(date), data=BP, span=span,
                      degree=1, family="symmetric")
      sys.fit <- predict(lo.sys, newdata=data.frame(date=BP$date), se=TRUE)
      dia.fit <- predict(lo.dia, newdata=data.frame(date=BP$date), se=TRUE)
      lines(date, sys.fit$fit, col="blue", lwd=2)
      lines(date, dia.fit$fit, col="magenta", lwd=2)

      if (confint){
        main <- paste0(main, "\nconfidence envelope",
                       if (length(level) > 1) "s",
                       " (level", if (length(level) > 1) "s", " = ",
                       paste(level, collapse=(", " )), ")")
        for (lev in level){
          z <- qnorm((1 - lev)/2, lower.tail=FALSE)
          sys.band <- c(sys.fit$fit + z*sys.fit$se,
                        rev(sys.fit$fit - z*sys.fit$se))
          polygon(c(date, rev(date)), sys.band,
                  col=rgb(0, 0, 1, alpha=alpha), border=NA)
          dia.band <- c(dia.fit$fit + z*dia.fit$se,
                        rev(dia.fit$fit - z*dia.fit$se))
          polygon(c(date, rev(date)), dia.band,
                  col=rgb(1, 0, 1, alpha=alpha), border=NA)
        }
      }

      title(main=main, cex.main=1, font.main=1)

    }

    if (legend){
      legend(as.Date("2025-03-20"), mean(par("usr")[3:4]),
             legend=c("Systolic", "Diastolic"), lty=1:2,
             col=c("blue", "magenta"), text.col=c("blue", "magenta"),
             pch=16:17, bty="n")
    } else{
      s.min <- which.min(sys)
      pos <- as.numeric(date[s.min]) >
        mean(as.numeric(date[1]), as.numeric(date[length(date)]))
      text(date[s.min], sys[s.min], labels="Systolic",
           col="blue", offset=1, xpd=TRUE,
           adj = if (pos) c(0, 1.5) else c(1, 1.5))
      d.max <- which.max(dia)
      pos <- as.numeric(date[d.max]) >
        mean(as.numeric(date[1]), as.numeric(date[length(date)]))
      text(date[d.max], dia[d.max], labels="Diastolic",
           col="magenta", offset=1, xpd=TRUE,
           adj = if (pos) c(1, -1) else c(0, -1))
    }

    if (!is.null(vlines)){
      for (vline in vlines){
        at <- as.numeric(as.Date(vline[["at"]]))
        pos <- if (at >
          mean(as.numeric(date[1]),
               as.numeric(date[length(date)]))) 2 else 4
        abline(v=at, col="gray30")
        adjust <- vline[["adjust"]]
        if (is.null(adjust)) adjust <- 0
        text(at, mean(range(sys, dia)) + adjust,
             labels=paste0(vline[["text"]]), pos=pos, col="gray30",
             offset=0.1)
      }
    }

  })

}

