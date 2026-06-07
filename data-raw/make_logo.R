## Generate BAMP hex-sticker logo — Lexis diagram
## Run once with: source("data-raw/make_logo.R")
##
## Required: install.packages(c("hexSticker", "showtext", "ggplot2", "scales"))

library(ggplot2)
library(hexSticker)
library(showtext)
library(scales)
library(bamp)

font_add_google("Montserrat", "montserrat")
showtext_auto()

# ── Data ──────────────────────────────────────────────────────────────────────
data(apc)
nop <- nrow(cases)   # 10 periods
noa <- ncol(cases)   # 8 age groups
vdb <- 2             # cohort diagonals at ~63°

rate <- cases / population
rate[is.nan(rate)] <- 0

df <- expand.grid(period = seq_len(nop), age = seq_len(noa))
df$rate <- as.vector(rate)

# ── Lexis cohort diagonals: period = vdb*(age-1) + k  (vdb=1 → 45°) ──────────
cohort_offsets <- seq(1 - vdb * (noa - 1), nop, by = vdb)

lex_lines <- do.call(rbind, lapply(cohort_offsets, function(k) {
  age_lo <- max(1, ceiling((1 - k) / vdb + 1))
  age_hi <- min(noa, floor((nop - k) / vdb + 1))
  if (age_lo >= age_hi) return(NULL)
  data.frame(
    age    = c(age_lo, age_hi),
    period = c(vdb * (age_lo - 1) + k, vdb * (age_hi - 1) + k),
    cohort = k
  )
}))

# ── Lexis diagram ─────────────────────────────────────────────────────────────
panel <- ggplot(df, aes(x = age, y = period, fill = rate)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_gradientn(
    colours = c("#05103a", "#1e3a8a", "#7c3aed", "#db2777", "#f97316", "#fde047"),
    values  = rescale(quantile(df$rate, c(0, .30, .55, .75, .90, 1))),
    guide   = "none"
  ) +
  # Period lines (vertical)
  geom_vline(
    xintercept  = seq(1.5, noa - 0.5, 1),
    colour      = "white", alpha = 0.30, linewidth = 0.25
  ) +
  # Age lines (horizontal)
  geom_hline(
    yintercept  = seq(1.5, nop - 0.5, 1),
    colour      = "white", alpha = 0.30, linewidth = 0.25
  ) +
  # Cohort diagonals
  geom_line(
    data        = lex_lines,
    aes(x = age, y = period, group = cohort),
    colour      = "white",
    alpha       = 0.55,
    linewidth   = 0.45,
    inherit.aes = FALSE
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_fixed(expand = FALSE) +
  theme_void() +
  theme(legend.position = "none")

# ── Hex sticker ───────────────────────────────────────────────────────────────
hexSticker::sticker(
  panel,
  package    = "bamp",
  p_size     = 22,
  p_color    = "#fde047",
  p_family   = "montserrat",
  p_fontface = "bold",
  p_x = 1.00, p_y = 0.36,
  s_x = 1.00, s_y = 1.08,
  s_width  = 1.60, s_height = 1.15,
  h_fill  = "#05103a",
  h_color = "#7c3aed",
  h_size  = 1.6,
  filename = "man/figures/logo.png",
  dpi = 320
)

message("Logo written to man/figures/logo.png")
