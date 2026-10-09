## Generates man/figures/logo.png (source of the pkgdown favicons).
## Motif: a skew-normal density whose shape estimate is pulled back from
## the half-normal limit (dashed curve) to a finite value (solid curve).
## Base grid graphics only (hexSticker and magick are not required).
library(grid)
hex <- function(r = 0.5) {
  a <- seq(30, 330, by = 60) * pi / 180
  list(x = 0.5 + r * cos(a), y = 0.5 + r * sin(a))
}
dsn <- function(x, a) 2 * dnorm(x) * pnorm(a * x)
png("man/figures/logo.png", width = 518, height = 600, res = 144, bg = "transparent")
grid.newpage()
pushViewport(viewport(width = unit(1, "snpc") * 0.98, height = unit(1, "snpc") * 0.98 * 600 / 518))
h <- hex(0.48)
grid.polygon(h$x, h$y, gp = gpar(fill = "#14365d", col = "#e8a33d", lwd = 5))
x <- seq(-1.2, 3.2, length.out = 300)
sx <- function(v) 0.20 + (v + 1.2) / 4.4 * 0.60
sy <- function(v) 0.38 + v / 0.9 * 0.36
grid.lines(sx(x), sy(ifelse(x < 0, 0, 2 * dnorm(x))), gp = gpar(col = "#e8a33d", lwd = 3, lty = 2))
grid.polygon(c(sx(x), sx(rev(x))[1]), c(sy(dsn(x, 4)), sy(0)), gp = gpar(fill = "#3b74ad", col = NA))
grid.lines(sx(x), sy(dsn(x, 4)), gp = gpar(col = "#f4f1ea", lwd = 4))
grid.lines(c(0.2, 0.8), c(sy(0), sy(0)), gp = gpar(col = "#f4f1ea", lwd = 2))
grid.text("snbr", x = 0.5, y = 0.24, gp = gpar(col = "#f4f1ea", fontsize = 30, fontface = "bold", fontfamily = "sans"))
popViewport()
invisible(dev.off())
