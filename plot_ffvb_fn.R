#plot simmr output for ffvb output
plot.simmr_output_ffvb <-
  function(x,
           type = c(
             "isospace",
             "histogram",
             "density",
             "matrix",
             "boxplot"
           ),
           group = 1,
           binwidth = 0.05,
           alpha = 0.5,
           title = if (x$group == 1) {
             "simmr output plot"
           } else {
             paste("simmr output plot: group", group)
           },
           ggargs = NULL,
           ...) {
    
    # Get the specified type
    type <- match.arg(type, several.ok = TRUE)
    
    # Iso-space plot is special as all groups go on one plot
    # Add in extra dots here as they can be sent to this plot function
    if ("isospace" %in% type) graphics::plot(x$input, group = group, title = title, ...)
    
    # Get group names
    group_names <- levels(x$groupnames)[group]
    
    for (i in 1:(x$group)) {
      
      # Prep data
      n <- 3600 
      K <- x$n_sources
      all_vb <- sim_thetacpp(n, x$lambda, x$n_sources, x$n_tracers)
      f_VB <- all_vb[,1:K]
      p_fun <- function(x) exp(x)/sum(exp(x))
      p_VB <- t(apply(all_vb[,1:K], 1, p_fun))

      
      out_all <- p_VB
      
      colnames(out_all) <- x$source_names
      df <- reshape2::melt(out_all)
      colnames(df) <- c("Num", "Source", "Proportion")
      if ("histogram" %in% type) {
        g <- ggplot(df, aes_string(
          x = "Proportion", y = "..density..",
          fill = "Source"
        )) +
          scale_fill_viridis(discrete = TRUE) +
          geom_histogram(binwidth = binwidth, alpha = alpha) +
          theme_bw() +
          ggtitle(title[i]) +
          facet_wrap("~ Source") +
          theme(legend.position = "none") +
          ggargs
        print(g)
      }
      
      if ("density" %in% type) {
        g <- ggplot(df, aes_string(
          x = "Proportion", y = "..density..",
          fill = "Source"
        )) +
          scale_fill_viridis(discrete = TRUE) +
          geom_density(alpha = alpha, linetype = 0) +
          theme_bw() +
          theme(legend.position = "none") +
          ggtitle(title[i]) +
          ylab("Density") +
          facet_wrap("~ Source") +
          ggargs
        print(g)
      }
      
      if ("boxplot" %in% type) {
        g <- ggplot(df, aes_string(
          y = "Proportion", x = "Source",
          fill = "Source", alpha = "alpha"
        )) +
          scale_fill_viridis(discrete = TRUE) +
          geom_boxplot(alpha = alpha, notch = TRUE, outlier.size = 0) +
          theme_bw() +
          ggtitle(title[i]) +
          theme(legend.position = "none") +
          coord_flip() +
          ggargs
        print(g)
      }
      
      # if ('convergence'%in%type) {
      #   coda::gelman.plot(x$output[[group[i]]],transform=TRUE)
      # }
      
      if ("matrix" %in% type) {
        # These taken from the help(pairs) file
        panel.hist <- function(x, ...) {
          usr <- graphics::par("usr")
          on.exit(graphics::par(usr))
          graphics::par(usr = c(usr[1:2], 0, 1.5))
          h <- graphics::hist(x, plot = FALSE)
          breaks <- h$breaks
          nB <- length(breaks)
          y <- h$counts
          y <- y / max(y)
          graphics::rect(breaks[-nB], 0, breaks[-1], y, col = "lightblue", ...)
        }
        panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
          usr <- graphics::par("usr")
          on.exit(graphics::par(usr))
          graphics::par(usr = c(0, 1, 0, 1))
          r <- stats::cor(x, y)
          txt <- format(c(r, 0.123456789), digits = digits)[1]
          txt <- paste0(prefix, txt)
          if (missing(cex.cor)) cex.cor <- 0.8 / graphics::strwidth(txt)
          graphics::text(0.5, 0.5, txt, cex = cex.cor * abs(r))
        }
        panel.contour <- function(x, y, ...) {
          usr <- graphics::par("usr")
          on.exit(graphics::par(usr))
          graphics::par(usr = c(usr[1:2], 0, 1.5))
          kd <- MASS::kde2d(x, y)
          kdmax <- max(kd$z)
          graphics::contour(kd, add = TRUE, drawlabels = FALSE, levels = c(kdmax * 0.1, kdmax * 0.25, kdmax * 0.5, kdmax * 0.75, kdmax * 0.9))
        }
        graphics::pairs(out_all, xlim = c(0, 1), ylim = c(0, 1), main = title[i], diag.panel = panel.hist, lower.panel = panel.cor, upper.panel = panel.contour)
      }
    }
    if (exists("g")) invisible(g)
  }


#test--------------------------
plot.simmr_output_ffvb(simmr_out1, type = "matrix")



