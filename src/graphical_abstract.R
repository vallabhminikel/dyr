# Graphical abstract: survival benefit of PrP lowering across in vivo models.
# Reads display_items/graphabs_summary_data.tsv; writes display_items/graphical_abstract.png (1200 x 1200 px, 300 dpi)
if (interactive()) setwd('~/d/sci/src/dyr')
suppressPackageStartupMessages({ library(tidyverse); library(png); library(magick) })

icon_dir <- '../native_function/data/iconography'
dat <- read_tsv('display_items/graphabs_summary_data.tsv', col_types = cols())

# ---- per-bar display spec (colors are the main colors used in the respective figures) ----
# rows are listed in the order they are drawn (top to bottom) within each experiment
spec <- tribble(
  ~panel,      ~condition_key,                                                              ~label,             ~color,
  'Figure 1F', 'Het KO',                                                                    'het KO',           '#501378',
  'Figure 1H', 'ASO 6, 300 ug, single dose at 127 dpi',                                     'low dose ASO',     '#57f781',
  'Figure 1H', 'ASO 6, 600 ug, single dose at 127 dpi',                                     'mid dose ASO',     '#1ac589',
  'Figure 1H', 'ASO 6, 1000 ug, single dose at 127 dpi',                                    'high dose ASO',    '#007a75',
  'Figure 2C', 'ASO 1, 300 ug, single dose at 64 dpi',                                      'late ASO',         '#1fa902',
  'Figure 2C', 'ASO 1, 300 ug, single dose at 43 dpi',                                      'early ASO',        '#93c47d',
  'Figure 2C', 'Het KO',                                                                    'het KO',           '#501378',
  'Figure 3F', 'ASO 6, 300 ug, 2 doses at 60 and 150 days of age',                          'ASO',              '#0077b6',
  'Figure 3G', 'ASO 6, 100 ug, 2 doses at 30 and 120 days of age',                          'early low ASO',    '#6baed6',
  'Figure 3G', 'ASO 6, 300 ug, 2 doses at 30 and 120 days of age',                          'early mid ASO',    '#3182bd',
  'Figure 3G', 'ASO 6, 500 ug, 2 doses at 30 and 120 days of age',                          'early high ASO',   '#08519c',
  'Figure 4C', 'Hemizygous transgene (half PrP dose; not a KO allele)',                     '10 copies',   '#3a6ea5',
  'Figure 4C', 'Hemizygous transgene, different line (Tg25109, 3 copies; not a KO allele)', '3 copies',    '#eb4d4b',
  'Figure 4F', 'Het KO (Ki817/KO)',                                                         'het KO',           '#3a6ea5'
)
ctl_label <- c('Figure 1F' = 'WT', 'Figure 1H' = 'vehicle', 'Figure 2C' = 'untreated WT',
               'Figure 3F' = 'PBS', 'Figure 3G' = 'PBS',
               'Figure 4C' = 'MM1, 20 copies', 'Figure 4F' = 'VV2, Ki817 hom')
ctl_color <- '#A9A9A9'

models <- tribble(
  ~model,             ~sub,                         ~icon,     ~bubble, ~panels,
  'Rat',              'RRML prions',                'rat',     'none',  list(c('Figure 1F', 'Figure 1H')),
  'Syrian hamster',   '263K prions',                'hamster', 'none',  list('Figure 2C'),
  'Tg(Vole) mouse',   'spontaneous disease',        'mouse',   'vole',  list(c('Figure 3F', 'Figure 3G')),
  'Humanized mouse',  'sCJD MM1 / VV2 prions',      'mouse',   'human', list(c('Figure 4C', 'Figure 4F'))
)

# ---- canvas: 4 x 4 in (1200 x 1200 px @ 300 dpi). y data units are INCHES (top = 0), x units are 100 per inch ----
W <- 4; H <- 4
mai <- c(0.50, 0.30, 0.04, 0.04)
plot_h <- H - mai[1] - mai[3]
xl <- c(-105, 100 * (W - mai[2] - mai[4]) - 105)       # left 1.05 in = icon column, 200% = 2.0 in
ix <- function(inch) inch * 100
iy <- function(inch) inch

# ---- text sizes (cex at pointsize 12) ----
cex_bar <- 0.56; cex_ctl <- 0.50; cex_model <- 0.62; cex_sub <- 0.50; cex_axis <- 0.58; cex_lab <- 0.62

# ---- bar lengths and percent increases use medians rounded to whole days before dividing ----
get_rows <- function(p) {
  d <- dat %>% filter(figure_panel == p)
  ctl_med <- round(d$control_median_survival[1])
  tr <- spec %>% filter(panel == p) %>%
    left_join(d, by = c('panel' = 'figure_panel', 'condition_key' = 'condition')) %>%
    mutate(len = 100 * round(test_median_survival) / ctl_med,
           pct = round(len - 100), is_ctl = FALSE) %>%
    select(label, color, len, pct, is_ctl)
  stopifnot(!anyNA(tr$len), nrow(tr) > 0)
  bind_rows(tibble(label = ctl_label[[p]], color = ctl_color, len = 100, pct = NA_real_, is_ctl = TRUE), tr)
}

# ---- icons: hi-res, smoothed edges; full square canvases are kept so relative sizes match the native_function repo ----
load_icon <- function(f, up = 3) {
  im <- image_read(file.path(icon_dir, paste0(f, '.png')))
  big <- image_resize(im, paste0(300 * up, 'x', 300 * up, '!'), filter = 'Lanczos')
  tf <- tempfile(fileext = '.png')
  image_write(big, tf, format = 'png'); arr <- readPNG(tf)               # h, w, 4
  # smooth the alpha edge: blur, then re-sharpen with a steep sigmoid, so upscaled outlines are smooth rather than blocky
  writePNG(arr[, , 4], tf)
  image_write(image_blur(image_read(tf), 0, up * 1.2), tf, format = 'png')
  sa <- readPNG(tf); if (length(dim(sa)) == 3) sa <- sa[, , 1]
  arr[, , 4] <- plogis((sa - 0.5) * 18)
  list(img = arr,
       rows = range(which(rowSums(arr[, , 4] > 0.5) > 0)) / nrow(arr),
       cols = range(which(colSums(arr[, , 4] > 0.5) > 0)) / ncol(arr))
}
icons <- map(unique(models$icon), load_icon) %>% set_names(unique(models$icon))
icon_canvas_in <- 1.0   # same canvas size for every rodent, as in the native_function figures

# ---- stack geometry (inches) ----
txt_h <- 0.24                 # model name + subtitle
bub_h <- 0.36; bub_gap <- 0.08
stack_h <- function(m) {
  ic <- icons[[models$icon[m]]]
  ch <- diff(ic$rows) * icon_canvas_in
  ch + txt_h + 0.04 + ifelse(models$bubble[m] != 'none', bub_h + bub_gap, 0)
}
rows_by_model <- map(seq_len(nrow(models)), function(m) map(unlist(models$panels[[m]]), get_rows))
n_rows <- map_int(rows_by_model, ~ sum(map_int(.x, nrow)))
n_gaps <- map_int(rows_by_model, ~ length(.x) - 1)
gap_exp_f <- 0.6; model_gap <- 0.15
blk_h <- function(p) pmap_dbl(list(seq_along(n_rows), n_rows, n_gaps), function(m, r, g) max(r * p + g * gap_exp_f * p, stack_h(m)))
tot <- function(p) sum(blk_h(p)) + model_gap * (nrow(models) - 1)
p <- uniroot(function(p) tot(p) - plot_h, c(0.03, 0.4))$root
cat(sprintf('row pitch %.3f in\n', p))

# ---- assign y positions ----
bars <- list(); blocks <- list(); y <- 0
bh_list <- blk_h(p)
for (m in seq_len(nrow(models))) {
  rows_h <- n_rows[m] * p + n_gaps[m] * gap_exp_f * p
  yb <- y + (bh_list[m] - rows_h) / 2          # center bars within block
  for (r in rows_by_model[[m]]) {
    r$y <- yb + (seq_len(nrow(r)) - 0.5) * p
    yb <- yb + nrow(r) * p + gap_exp_f * p
    bars[[length(bars) + 1]] <- r
  }
  blocks[[m]] <- c(top = y, bottom = y + bh_list[m])
  y <- y + bh_list[m] + model_gap
}
bars <- bind_rows(bars)

# ---- drawing helpers (all sizes in inches) ----
circle <- function(cx, cy, rin, n = 60) { t <- seq(0, 2 * pi, length.out = n); cbind(cx + ix(rin * cos(t)), cy + iy(rin * sin(t))) }
ell <- function(cx, cy, ax, ay, rot = 0, n = 60) {
  t <- seq(0, 2 * pi, length.out = n); xx <- ax * cos(t); yy <- ay * sin(t)
  cbind(cx + ix(xx * cos(rot) - yy * sin(rot)), cy + iy(xx * sin(rot) + yy * cos(rot)))
}
ink <- '#2b2d3a'
# smooth closed outline through control points (template coordinates, y up), scaled by s inches
smooth_poly <- function(cx, cy, pts, s = 1, shape = 1, col = ink) {
  xs <- cx + ix(pts[, 1] * s); ys <- cy - iy(pts[, 2] * s)
  xsp <- xspline(xs, ys, shape = shape, open = FALSE, draw = FALSE, repEnds = TRUE)
  polygon(xsp$x, xsp$y, col = col, border = NA)
}
# bank vole, facing right (blunt snout, small round ear, short tail), flat silhouette like the rodent icons
draw_vole <- function(cx, cy, s) {
  body <- rbind(c(-1.00, .10), c(-.85, .50), c(-.45, .72), c(.05, .80), c(.50, .74), c(.80, .58),
                c(1.02, .42), c(1.20, .24), c(1.20, .08), c(1.06, -.02), c(.85, -.05), c(.62, -.20),
                c(.60, -.38), c(.42, -.42), c(.30, -.28), c(-.10, -.28), c(-.40, -.40), c(-.62, -.42),
                c(-.78, -.30), c(-.95, -.12))
  smooth_poly(cx, cy, body, s)
  tt <- seq(0, 2 * pi, length.out = 13)[-13]
  smooth_poly(cx, cy, cbind(.50 + .11 * cos(tt), .82 + .12 * sin(tt)), s)
  xs <- cx + ix(s * c(-.97, -1.2, -1.42)); ys <- cy - iy(s * c(.08, .02, -.12))
  xsp <- xspline(xs, ys, shape = 1, open = TRUE, draw = FALSE)
  lines(xsp$x, xsp$y, lwd = 1.6, col = ink, lend = 'round')
}
# full-body person pictogram (round head, rounded shoulders/arms, straight legs); s = total height in inches
draw_human <- function(cx, cy, s) {
  lw <- function(w) w * 96                      # inches -> lwd
  seg <- function(x0, y0, x1, y1, w) segments(cx + ix(x0 * s), cy - iy(y0 * s), cx + ix(x1 * s), cy - iy(y1 * s), lwd = lw(w * s), col = ink, lend = 'round')
  polygon(circle(cx, cy - iy(.40 * s), .105 * s), col = ink, border = NA)              # head
  seg(0, .20, 0, -.02, .30)                                                            # torso
  seg(-.19, .22, -.19, .02, .075); seg(.19, .22, .19, .02, .075)                     # arms
  seg(-.055, -.05, -.055, -.42, .095); seg(.055, -.05, .055, -.42, .095)               # legs
}
thought_bubble <- function(cx, cy, w, h, what) {
  n <- 8; t <- seq(0, 2 * pi, length.out = n + 1)[-(n + 1)]
  for (k in 1:2) {
    col <- if (k == 1) 'black' else 'white'
    for (i in seq_along(t))
      polygon(ell(cx + ix(w * .36 * cos(t[i])), cy + iy(h * .36 * sin(t[i])), w * .24, h * .26), col = col, border = if (k == 1) 'black' else NA, lwd = 1.1)
    polygon(ell(cx, cy, w * .4, h * .4), col = col, border = NA)
  }
  if (what == 'vole') draw_vole(cx, cy + iy(h * .02), w * .21) else draw_human(cx, cy + iy(h * .02), h * .78)
}
trail <- function(x0, y0, x1, y1) {
  for (f in c(.30, .72)) polygon(circle(x0 + (x1 - x0) * f, y0 + (y1 - y0) * f, .02 + .018 * f), col = 'white', border = 'black', lwd = 1)
}
lum <- function(hex) { v <- col2rgb(hex) / 255; sum(c(.299, .587, .114) * v) }

# ---- render ----
png('display_items/graphical_abstract.png', width = W * 300, height = H * 300, res = 300)
par(mai = mai, family = 'sans', xpd = FALSE)
plot(NA, xlim = xl, ylim = c(plot_h, 0), axes = FALSE, ann = FALSE, xaxs = 'i', yaxs = 'i')
stopifnot(abs(par('pin')[2] - plot_h) < 1e-6, abs(par('pin')[1] * 100 - diff(xl)) < 1e-3)

segments(100, 0, 100, plot_h, lty = 3, col = 'gray50', lwd = .8)
segments(0, 0, 0, plot_h, col = 'black', lwd = 1)   # y axis line
bh <- 0.82 * p
for (i in seq_len(nrow(bars))) {
  b <- bars[i, ]
  rect(0, b$y - bh / 2, b$len, b$y + bh / 2, col = b$color, border = NA)
  if (b$is_ctl) {
    text(ix(.02), b$y, b$label, adj = c(0, .5), cex = cex_ctl, col = 'white')
  } else {
    text(b$len - ix(.02), b$y, sprintf('+%.0f%%', b$pct), adj = c(1, .5), cex = cex_bar, font = 2,
         col = ifelse(lum(b$color) > .55, 'black', 'white'))
    text(b$len + ix(.03), b$y, b$label, adj = c(0, .5), cex = cex_bar, col = 'black')
  }
}

axis(1, at = seq(0, 200, 50), labels = paste0(seq(0, 200, 50), '%'), pos = plot_h, cex.axis = cex_axis, tck = -.015, padj = -3.2)
mtext('survival (% of untreated median)', side = 1, line = 0.65, cex = cex_lab, at = 100)
par(xpd = NA)
text(xl[1] - ix(.19), plot_h / 2, 'in vivo model', srt = 90, cex = cex_lab)

# ---- model icons + names ----
for (m in seq_len(nrow(models))) {
  bl <- blocks[[m]]; cy <- mean(bl)
  ic <- icons[[models$icon[m]]]
  has_b <- models$bubble[m] != 'none'
  ch <- diff(ic$rows) * icon_canvas_in
  sh <- stack_h(m)
  top <- cy - sh / 2
  cx <- -52
  itop_c <- top + ifelse(has_b, bub_h + bub_gap, 0)           # top of icon content
  cbot <- itop_c + ch                                          # bottom of icon content
  canvas_bot <- cbot + (1 - ic$rows[2]) * icon_canvas_in
  canvas_top <- canvas_bot - icon_canvas_in
  rasterImage(as.raster(ic$img), cx - ix(icon_canvas_in / 2), canvas_bot, cx + ix(icon_canvas_in / 2), canvas_top, interpolate = TRUE)
  text(cx, cbot + .07, models$model[m], font = 2, cex = cex_model, adj = c(.5, .5))
  text(cx, cbot + .07 + .115, models$sub[m], cex = cex_sub, adj = c(.5, .5), col = 'gray30')
  if (has_b) {
    head_x <- cx + ix(icon_canvas_in * (mean(ic$cols) - .5) + .06)
    bx <- head_x - ix(.05)
    by <- top + bub_h / 2
    trail(head_x - ix(.04), itop_c - .005, bx, top + bub_h)
    thought_bubble(bx, by, .63, bub_h, models$bubble[m])
  }
  if (m < nrow(models)) {
    ys <- bl[['bottom']] + model_gap / 2
    segments(xl[1] + 2, ys, xl[2] - 2, ys, col = 'gray85', lwd = .6)
  }
}
dev.off()
cat('wrote display_items/graphical_abstract.png\n')
