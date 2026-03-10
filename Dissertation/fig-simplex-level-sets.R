library(ggplot2)
library(dplyr)
library(tidyr)

# -----------------------------------------------------------------------
# Coordinate system
# -----------------------------------------------------------------------

u <- c(1, -1, 0) / sqrt(2) # horizontal: e1 left, e2 right
w <- c(-1, -1, 2) / sqrt(6) # vertical:   e3 at top
c_vec <- c(1 / 3, 1 / 3, 1 / 3)

project_2d <- function(x) {
  v <- x - c_vec
  c(s = sum(v * u), t = sum(v * w))
}

# -----------------------------------------------------------------------
# Simplex vertices in 2D
# -----------------------------------------------------------------------

vertices_3d <- list(e1 = c(1, 0, 0), e2 = c(0, 1, 0), e3 = c(0, 0, 1))
vertices_2d <- as.data.frame(t(sapply(vertices_3d, project_2d)))

simplex_polygon <- rbind(
  vertices_2d[1, c("s", "t")],
  vertices_2d[2, c("s", "t")],
  vertices_2d[3, c("s", "t")],
  vertices_2d[1, c("s", "t")]
)

simplex_edges <- bind_rows(
  data.frame(
    s = c(vertices_2d$s[1], vertices_2d$s[2]),
    t = c(vertices_2d$t[1], vertices_2d$t[2]),
    edge = "e1e2"
  ),
  data.frame(
    s = c(vertices_2d$s[2], vertices_2d$s[3]),
    t = c(vertices_2d$t[2], vertices_2d$t[3]),
    edge = "e2e3"
  ),
  data.frame(
    s = c(vertices_2d$s[3], vertices_2d$s[1]),
    t = c(vertices_2d$t[3], vertices_2d$t[1]),
    edge = "e3e1"
  )
)

# -----------------------------------------------------------------------
# Circles / arcs
# -----------------------------------------------------------------------

psi_low <- 0.38
psi_mid <- 0.50
psi_high <- 0.65

r_low <- sqrt(psi_low - 1 / 3)
r_mid <- sqrt(psi_mid - 1 / 3)
r_high <- sqrt(psi_high - 1 / 3)

circle_points <- function(psi_hat, n = 2000) {
  r <- sqrt(psi_hat - 1 / 3)
  phi <- seq(0, 2 * pi, length.out = n)
  pts <- lapply(phi, function(p) {
    x <- c_vec + r * (cos(p) * u + sin(p) * w)
    data.frame(
      s = cos(p) * r,
      t = sin(p) * r,
      x1 = x[1],
      x2 = x[2],
      x3 = x[3],
      phi = p,
      inside = all(x > 0)
    )
  })
  bind_rows(pts)
}

tag_segments <- function(df) {
  df <- df |> arrange(phi)
  seg <- cumsum(c(1, abs(diff(as.integer(df$inside))) > 0))
  df$seg <- seg
  df
}

pts_low <- circle_points(psi_low)
pts_mid <- circle_points(psi_mid)
pts_high <- circle_points(psi_high) |> tag_segments()

# -----------------------------------------------------------------------
# Tangency points
# -----------------------------------------------------------------------

tangency_3d <- list(c(1 / 2, 1 / 2, 0), c(1 / 2, 0, 1 / 2), c(0, 1 / 2, 1 / 2))
tangency_2d <- as.data.frame(t(sapply(tangency_3d, project_2d)))

# -----------------------------------------------------------------------
# Annotation positions
# psi labels on left, vertically at each circle's top
# arrow_len just past the small circle
# -----------------------------------------------------------------------

label_offset <- 0.05
label_x <- -0.09
arrow_len <- r_low + 0.12

# -----------------------------------------------------------------------
# Bounding box: driven purely by triangle vertices, circles, arrows,
# and psi labels — no vertex label text to account for
# -----------------------------------------------------------------------

pad <- 0.08

# psi label text extends left of label_x by approx 0.28 units
psi_label_left <- label_x - 0.28

extent_s_max <- max(vertices_2d$s, r_high, arrow_len + 0.08)
extent_s_min <- min(vertices_2d$s, -r_high, psi_label_left)
extent_t_max <- max(
  vertices_2d$t,
  r_high,
  arrow_len + 0.06,
  r_high + label_offset + 0.04
)
extent_t_min <- min(vertices_2d$t, -r_high)

plane_smin <- extent_s_min - pad
plane_smax <- extent_s_max + pad
plane_tmin <- extent_t_min - pad
plane_tmax <- extent_t_max + pad

plane_rect <- data.frame(
  s = c(plane_smin, plane_smax, plane_smax, plane_smin, plane_smin),
  t = c(plane_tmin, plane_tmin, plane_tmax, plane_tmax, plane_tmin)
)

# -----------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------

p <- ggplot() +

  # Affine plane Pi
  geom_polygon(
    data = plane_rect,
    aes(x = s, y = t),
    fill = "#eef2f7",
    colour = "#9bacc4",
    linewidth = 0.5
  ) +

  annotate(
    "text",
    x = plane_smax - 0.03,
    y = plane_tmax - 0.03,
    label = "Pi",
    parse = TRUE,
    size = 5.5,
    colour = "black",
    hjust = 1,
    vjust = 1,
    fontface = "italic"
  ) +

  # Simplex fill
  geom_polygon(
    data = simplex_polygon,
    aes(x = s, y = t),
    fill = "#f2f2f2",
    colour = NA
  ) +

  # Simplex edges
  geom_path(
    data = simplex_edges,
    aes(x = s, y = t, group = edge),
    colour = "#666666",
    linewidth = 0.9,
    linetype = "dashed"
  ) +

  # Full C_psi_hat for psi_high (faint dashed)
  geom_path(
    data = pts_high,
    aes(x = s, y = t),
    colour = "#b05a40",
    linewidth = 0.5,
    linetype = "dashed"
  ) +

  # psi = 1/2 circle
  geom_path(
    data = pts_mid,
    aes(x = s, y = t),
    colour = "#2a9d8f",
    linewidth = 1.2
  ) +

  # Omega_psi_hat: complete circle for psi_low
  geom_path(
    data = filter(pts_low, inside),
    aes(x = s, y = t),
    colour = "#2166ac",
    linewidth = 1.2
  ) +

  # Omega_psi_hat: disjoint arcs for psi_high
  geom_path(
    data = filter(pts_high, inside),
    aes(x = s, y = t, group = seg),
    colour = "#d6604d",
    linewidth = 1.2
  ) +

  # Tangency points
  geom_point(
    data = tangency_2d,
    aes(x = s, y = t),
    shape = 21,
    size = 2.5,
    colour = "#2a9d8f",
    fill = "white",
    stroke = 0.9
  ) +

  # Vertex labels: e1 and e2 below, e3 to the right
  annotate(
    "text",
    x = vertices_2d$s[1],
    y = vertices_2d$t[1] - 0.06,
    label = "bold(e)[2]",
    parse = TRUE,
    size = 5,
    colour = "#3a2a1a",
    hjust = 0.5,
    vjust = 1
  ) +
  annotate(
    "text",
    x = vertices_2d$s[2],
    y = vertices_2d$t[2] - 0.06,
    label = "bold(e)[1]",
    parse = TRUE,
    size = 5,
    colour = "#3a2a1a",
    hjust = 0.5,
    vjust = 1
  ) +
  annotate(
    "text",
    x = vertices_2d$s[3] + 0.06,
    y = vertices_2d$t[3],
    label = "bold(e)[3]",
    parse = TRUE,
    size = 5,
    colour = "#3a2a1a",
    hjust = 0,
    vjust = 0.5
  ) +

  # Centroid
  geom_point(aes(x = 0, y = 0), colour = "#222222", size = 1.8) +
  annotate(
    "text",
    x = -0.03,
    y = -0.05,
    label = "bold(c)",
    parse = TRUE,
    size = 4.5,
    colour = "#222222"
  ) +

  # Basis vector arrows
  annotate(
    "segment",
    x = 0,
    y = 0,
    xend = arrow_len,
    yend = 0,
    arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
    colour = "#888888",
    linewidth = 0.5
  ) +
  annotate(
    "text",
    x = arrow_len - 0.02,
    y = -0.04,
    label = "bold(u)",
    parse = TRUE,
    size = 3.8,
    colour = "#888888",
    hjust = 0,
    vjust = 0.5
  ) +
  annotate(
    "segment",
    x = 0,
    y = 0,
    xend = 0,
    yend = arrow_len,
    arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
    colour = "#888888",
    linewidth = 0.5
  ) +
  annotate(
    "text",
    x = -0.05,
    y = arrow_len - 0.02,
    label = "bold(w)",
    parse = TRUE,
    size = 3.8,
    colour = "#888888",
    hjust = 0,
    vjust = 0
  ) +

  # psi_hat labels
  annotate(
    "text",
    x = label_x,
    y = -r_low + label_offset + 0.02,
    label = "hat(psi) == 0.38",
    parse = TRUE,
    size = 4,
    colour = "#2166ac",
    hjust = 0
  ) +
  annotate(
    "text",
    x = label_x,
    y = -r_mid + label_offset + 0.01,
    label = "hat(psi) == 0.50",
    parse = TRUE,
    size = 4,
    colour = "#2a9d8f",
    hjust = 0
  ) +
  annotate(
    "text",
    x = label_x,
    y = -r_high + label_offset,
    label = "hat(psi) == 0.65",
    parse = TRUE,
    size = 4,
    colour = "#d6604d",
    hjust = 0
  ) +

  coord_equal(
    xlim = c(plane_smin, plane_smax),
    ylim = c(plane_tmin, plane_tmax),
    expand = FALSE
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(4, 4, 4, 4),
    plot.background = element_rect(fill = "white", colour = NA)
  )

ggsave(
  "Dissertation/files/Images/simplex-level-sets.png",
  plot = p,
  width = 4.5,
  height = 5.2,
  device = "png",
  dpi = 300,
  units = "in"
)
