# convenience functions / wrappers

# Not In 
'%ni%' = Negate('%in%')

## save pheatmap output 
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# Wrapper around Seurat GenePlot function with flowjo-like density, no jitter, for Seurat V2
# https://slowkow.com/notes/ggplot2-color-by-density/
####### Map identity on the same plot: 
get_density = function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
require(viridis)
GenePlot4 = function(object, gene1, gene2, pch.use = 20, pt.size = 0.5, col.use = NULL,
                     plot.title = NULL, plot.ident = F, ident.plot = NULL, ... ) {
  if (plot.ident == T) {
    cell.ids = object@cell.names
    data.use = as.data.frame(FetchData(object = object, 
                                       vars.all = c(gene1, gene2), 
                                       cells.use = cell.ids))
    names(data.use) = c("x", "y")
    object = SetAllIdent(object, id = ident.plot)
    data.plot = tibble(data.use$x, data.use$y, as.factor(object@ident[cell.ids])) 
    names(data.plot) = c("x", "y", "ident.use")
    p = ggplot2::ggplot(data = data.plot, mapping = aes(x = x, y = y)) + 
      geom_point(size = pt.size, alpha = 0.9, aes(color = ident.use)) +
      theme_light() + 
      labs(title = plot.title, x = gene1, y = gene2) +
      guides(shape = guide_legend(override.aes = list(size = 4))) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) 
    return(p)
  } else { 
    cell.ids = object@cell.names
    data.use = as.data.frame(FetchData(object = object, 
                                       vars.all = c(gene1, gene2), 
                                       cells.use = cell.ids))
    data.plot = data.use[cell.ids, c(gene1, gene2)]
    names(data.plot) = c("x", "y")
    data.plot = data.plot %>% mutate(density = get_density(data.plot$x, data.plot$y, n = 100))
    p <- ggplot2::ggplot(data = data.plot, mapping = aes(x = x, y = y))
    p <- p + geom_point(aes(x,y, shape = pch.use, color = density), size = pt.size) + 
      scale_color_viridis(option = "inferno") +  
      scale_shape_identity() 
    p <- p + 
      labs(title = plot.title, x = gene1, y = gene2) + 
      scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + 
      theme_light()
    return(p)     
  }
}
