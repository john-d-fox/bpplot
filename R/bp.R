#' Plot Daily Blood-Pressure Readings
#'
#' @param files character vector of (paths to) BP files.
#' @param legend (logical) if \code{TRUE} (the default) show a legend; if \code{FALSE} label
#'        systolic and diastolic values directly.
#' @param smooth (logical) show loess smooths (default \code{TRUE} if there are
#'        at least 15 dates)?
#' @param span for loess smooths; if there is more than one, the span is
#'        selected by 10-fold cross-validation; the default is 9 values
#'        between 0.1 and 0.9. Cross-validation doesn't work well if there
#'        are discontinuities. At present, the systolic data are used to
#'        pick the span.
#' @param confint (logical) show pointwise confidence envelope(s) around smooths?
#'        the default is the value of \code{smooth}. The confidence envelope
#'        assumes independent errors, which isn't generally the case.
#' @param level for confidence envelope(s); the default is 0.5 and 0.9.
#' @param alpha transparency for confidence envelopes (cumulative); the default
#'        is 0.25.
#' @param vlines list of lists for vertical lines to draw
#'        (the \code{"at"} component of each sub-list), labels for
#'        the lines (the \code{"text"} component), and (optionally, the
#'        \code{"at"} component) for vertical adjustment of text;
#'        see examples.
#' @param xlab,ylab horizontal axis and vertical axis labels.
#' @param show.main display main label (logical, default \code{TRUE}).
#'
#' @returns invisibly returns a data frame with columns for date, systolic
#'          blood pressure, and diastolic blood pressure.
#'          If cross-validation is used to pick the span, a table of spans
#'          and CV mean-squared error values for the systolic data is printed,
#'          along with residual autocorrelations.
#'          Otherwise, just the residual autocorrelations for the systolic
#'          data are printed. The latter can be distorted if there are
#'          discontinuities.
#'
#' @description
#' \code{bpplot()} graphs daily blood-pressure readings. For more detailed
#' analysis of blood-pressure data, see the \pkg{bp} package
#' (\url{https://github.com/johnschwenck/bp} and \url{https://cran.r-project.org/package=bp}).
#'
#' @examples
#' dir <- system.file("etc", package="bpplot")
#' (files <- list.files(dir)) # example data files provided
#' files <- paste0(dir, "/", files)
#' cat(paste(readLines(files[5], 10), "\n")) # part of one data file
#'
#' bpplot(files[3:5],
#'        span=0.3, # picked to smooth across discontinuity
#'        vlines=list(list(at="2025-07-08",
#'                         text="pause \u2192\nDrug      "),
#'                    list(at="2025-07-15",
#'                         text="resume\u2192\nDrug      ",
#'                         adjust=15))
#' )
#'
#' bpplot(files[1:3], legend=FALSE)


#' @export
bpplot <- function(files,
                   legend = TRUE,
                   smooth = nrow(BP) >= 15,
                   span = seq(0.1, 0.9, by = 0.1),
                   confint = smooth,
                   level = c(0.5, 0.9),
                   alpha = 0.25,
                   vlines = NULL,
                   xlab = "Date",
                   ylab = "Blood Pressure (mm Hg)",
                   show.main = TRUE
                     ){

  BP <- read.table(files[1], header=TRUE, fill=TRUE)
  for (file in files[-1]){
    bp <- read.table(file, header=TRUE, fill=TRUE)
    BP <- rbind(BP, bp)
  }
  BP$date <- as.Date(BP$date)
  BP <- na.omit(BP)

  with(BP, {
    plot(range(date), range(c(sys, dia)), type="n",
         xlab=xlab, ylab=ylab, xaxt="n", yaxt="n")
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
         mods <- vector(nmods, mode="list")
         autocor <- matrix(0, nmods, 10L)
         nr <- 10L
         for (i in 1L:nmods){
           cmd <- paste0('loess(sys ~ as.numeric(date), data=BP, span=', span[i],
                         ', degree=1, family="symmetric")')
           mods[[i]] <- eval(parse(text=cmd))
           r <- as.vector(acf(residuals(mods[[i]]), plot=FALSE)$acf)
           nr <- min(length(r) - 1, nr)
           autocor[i, 1L:nr] <- r[2L:(nr + 1)]
         }
         cvs <- sapply(mods, function (m) cv::cv(m, k="n")$"CV crit")
         cvs <-  cbind(span, cvs, round(autocor, 2))
         colnames(cvs) <- c("span", "CV MSE", paste0("r[", 1:nr, "]"))
         rownames(cvs) <- rep("", nrow(cvs))
         cat("\nFor smooths fit to the systolic data with",
             "residual autocorrelations:\n")
         print(cvs)
         span <- span[which.min(cvs[, "CV MSE"])]
       } else {
         mod <- loess(sys ~ as.numeric(date), data=BP, span=span,
                    degree=1, family="symmetric")
         r <- as.vector(acf(residuals(mod), plot=FALSE)$acf)
         nr <- min(length(r) - 1, 10L)
         r <- round(r[2L:(nr + 1)], 2)
         names(r) <- paste0("r[", 1L:nr, "]")
         cat("\nResidual autocorrelations for the systolic data:\n")
         print(r)
       }

      main <- if (show.main) {
       paste0("robust local linear regression (span = ",
                     round(span, 2), ")")
      } else {
        ""
      }
      lo.sys <- loess(sys ~ as.numeric(date), data=BP, span=span,
                      degree=1, family="symmetric")
      lo.dia <- loess(dia ~ as.numeric(date), data=BP, span=span,
                      degree=1, family="symmetric")
      sys.fit <- predict(lo.sys, newdata=data.frame(date=BP$date), se=TRUE)
      dia.fit <- predict(lo.dia, newdata=data.frame(date=BP$date), se=TRUE)
      lines(date, sys.fit$fit, col="blue", lwd=2)
      lines(date, dia.fit$fit, col="magenta", lwd=2)

      if (confint){
        if (show.main){
          main <- paste0(main, "\nconfidence envelope",
                         if (length(level) > 1) "s",
                         " (level", if (length(level) > 1) "s", " = ",
                         paste(level, collapse=(", " )), ")")
        }
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
      legend(min(date), mean(par("usr")[3:4]),
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

  return(invisible(BP))

}

#' @importFrom stats coef loess na.omit predict start update
#' @importFrom utils read.table globalVariables
#' @exportS3Method cv::cv
cv.loess <- function(model, data = insight::get_data(model), criterion = cv::mse,
                     k = 10L, reps = 1L, seed = NULL, details = k <= 10L, confint = n >=
                       400L, level = 0.95, ...) {

  f <- function(i_, ...) {
    # helper function to compute cv criterion for each fold
    indices.i <- fold(folds, i_)
    model.i <- if (start) {
      update(model, data = data[-indices.i, ], start = b)
    } else {
      update(model, data = data[-indices.i, ])
    }
    fit.all.i <- predict(model.i, newdata = data, type = type, ...)
    fit.i <- fit.all.i[indices.i]
    NAs <- is.na(fit.i)
    fit.i[NAs] <- y[indices.i][NAs]
    list(
      fit.i = fit.i,
      crit.all.i = criterion(y, fit.all.i),
      coef.i = coef(model.i)
    )
  }

  n <- nrow(data)

  cv::cvCompute(
    model = model,
    data = data,
    criterion = criterion,
    k = k,
    reps = reps,
    seed = seed,
    details = details,
    confint = confint,
    level = level,
    ncores = 1L,
    f = f,
    model.function = loess,
    model.function.name = "loess",
    ...
  )

}

globalVariables(c("fold", "folds", "b", "type", "y"))
