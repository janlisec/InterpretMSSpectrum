#'@export
plot.findMAIN <-
  function (x,
            rank = 1,
            correct_mass = NULL,
            ...)
  {
    if (length(rank) > 1) {
      opar <- graphics::par(mfrow = grDevices::n2mfrow(length(rank)))
      graphics::par(mar = c(2, 2, 2, 1))
      on.exit(graphics::par(opar))
    }
    idx <- rank
    if (any(idx > length(x))) { stop("object shorter than requested numer of ranks") }
    lidx <- length(idx)
    legend_text_col <- matrix(1, nrow = ncol(attr(x[[1]], "scores")), ncol = lidx)
    if (lidx > 1) {
      sm <- summary(x)[idx, ]
      legend_text_col[7, which(sm$mass_score == max(sm$mass_score))] <- 3
      legend_text_col[8, which(sm$int_score == max(sm$int_score))] <- 3
      legend_text_col[9, which(sm$supp_isos == max(sm$supp_isos))] <- 3
    }
    for (i in 1:lidx) {
      cols <- x[[idx[i]]][,"charge"]
      cols[is.na(cols)] <- 0
      cols <- cols+1
      cols[cols>=3] <- 6
      cols[!is.na(x[[idx[i]]][,"label"])] <- 4
      InterpretMSSpectrum::PlotSpec(x = x[[idx[i]]], cutoff = 0, cols = cols, txt = x[[idx[i]]][,c("mz","label")], ...)
      # browser()
      # the legend could be stablished based on he summary result...
      # leg <- summary(x)[idx[i],]
      # leg_lab <- colnames(leg)
      # leg_val <- leg
      leg_lab <- colnames(attr(x[[idx[i]]], "scores"))
      #leg_val <- as.character(round(x = as.numeric(attr(x[[idx[i]]], "scores")), digits = c(3, 0, 3, 0, 3, 2, 2, 2, 2, 2)))
      leg_val <- sapply(1:10, function(j) { formatC(x = as.numeric(attr(x[[idx[i]]], "scores"))[j], digits = c(4, 0, 4, 0, 2, 2, 2, 2, 2, 2)[j], format = "f") })
      leg_val[2] <- rownames(attr(x[[idx[i]]], "scores"))
      graphics::legend("topright", legend = paste(leg_lab, leg_val, sep=": "), bty = "n", cex = 0.75, text.col = legend_text_col[, i])
      mhyp <- attr(x[[idx[i]]], "scores")[, "neutral_mass"]
      color <- if (!is.null(correct_mass)) {
        if (abs(mhyp - correct_mass) < 0.01) 3 else 2
      } else {
        1
      }
      graphics::mtext(sprintf("[%d] %.4f", idx[i], round(mhyp, 4)), col = color, cex = 0.9, line = 0.05)
    }
  }