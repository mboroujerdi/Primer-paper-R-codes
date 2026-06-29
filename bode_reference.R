###############################################################################
# Bode Plot Reference Script
# Boroujerdi M. — Systems Endocrinology / Control Theory Primer
#
# PURPOSE
#   Given only the open-loop transfer function (numerator & denominator
#   polynomial coefficients), compute and plot:
#     • Open-loop Bode plot  (magnitude + phase)
#     • Closed-loop Bode plot (unity negative feedback)
#     • Stability metrics    (phase margin, gain margin, bandwidth)
#
# TWO ENTRY POINTS
#   1. Direct polynomial entry  — supply num/den vectors directly
#   2. MCR / (omegan, zeta) entry — build the transfer function from
#      endocrine parameters using G-type or H-type architecture,
#      then pass to the same plotting pipeline
#
# SIGN CONVENTION
#   Polynomials are written in DESCENDING powers of s, matching R's
#   standard and the paper's notation.
#   Example: s^2 + 2s + 1  ->  c(1, 2, 1)
#
# UNITS
#   Use rad/min throughout (consistent with MCR values in Table 1).
#   Change freq_range if working in rad/s or rad/h.
###############################################################################

library(ggplot2)
library(patchwork)

# =============================================================================
# SECTION 1 — CORE UTILITIES
# =============================================================================

# Evaluate a polynomial p at complex values z
# p is a numeric vector of coefficients in DESCENDING power order.
polyval <- function(p, z) {
  n <- length(p) - 1
  sapply(z, function(zi) sum(p * zi^(n:0)))
}

# Phase unwrapping (prevents ±180° jumps in phase plot)
unwrap_phase <- function(phase_rad) {
  n   <- length(phase_rad)
  out <- numeric(n)
  out[1] <- phase_rad[1]
  for (i in 2:n) {
    d      <- phase_rad[i] - phase_rad[i - 1]
    d      <- d - 2 * pi * round(d / (2 * pi))
    out[i] <- out[i - 1] + d
  }
  out
}

# =============================================================================
# SECTION 2 — FREQUENCY RESPONSE
# =============================================================================
# Inputs:
#   num, den  : polynomial coefficient vectors (descending powers of s)
#   omega     : frequency vector (rad/time-unit)
#   dc_norm   : if TRUE, normalise so DC magnitude = 0 dB
#
# Returns a list:
#   $H        : complex frequency response vector
#   $mag_dB   : magnitude in dB
#   $phase_deg: unwrapped phase in degrees
#   $omega    : frequency vector (echoed back for convenience)

freq_response <- function(num, den, omega, dc_norm = FALSE) {
  s   <- 1i * omega
  H   <- polyval(num, s) / polyval(den, s)
  if (dc_norm) {
    dc <- Mod(H[1])
    if (dc > 0) H <- H / dc
  }
  ph_rad <- unwrap_phase(Arg(H))
  list(
    H         = H,
    mag_dB    = 20 * log10(Mod(H)),
    phase_deg = ph_rad * 180 / pi,
    omega     = omega
  )
}

# =============================================================================
# SECTION 3 — CLOSED-LOOP (unity negative feedback)
# =============================================================================
# Given the open-loop complex response vector H_ol,
# returns the closed-loop response:  H_cl = H_ol / (1 + H_ol)

closed_loop_from_ol <- function(fr_ol) {
  H_cl      <- fr_ol$H / (1 + fr_ol$H)
  ph_rad    <- unwrap_phase(Arg(H_cl))
  list(
    H         = H_cl,
    mag_dB    = 20 * log10(Mod(H_cl)),
    phase_deg = ph_rad * 180 / pi,
    omega     = fr_ol$omega
  )
}

# =============================================================================
# SECTION 4 — STABILITY METRICS
# =============================================================================
# Computed from the OPEN-LOOP frequency response.
#
# Phase margin (PM):
#   Phase at the gain crossover frequency (where |H| = 0 dB) plus 180°.
#   PM > 0 → stable.
#
# Gain margin (GM):
#   Negative of the gain (dB) at the phase crossover frequency
#   (where phase = –180°).
#   GM > 0 → stable.
#
# –3 dB bandwidth:
#   Frequency at which open-loop magnitude first drops 3 dB below DC.

stability_metrics <- function(fr) {
  mag   <- fr$mag_dB
  phase <- fr$phase_deg
  omega <- fr$omega

  # Gain crossover (0 dB crossing)
  gc_idx       <- which.min(abs(mag))
  gc_freq      <- omega[gc_idx]
  phase_at_gc  <- phase[gc_idx]
  PM           <- 180 + phase_at_gc

  # Phase crossover (–180° crossing)
  pc_idx       <- which.min(abs(phase + 180))
  pc_freq      <- omega[pc_idx]
  gain_at_pc   <- mag[pc_idx]
  GM           <- -gain_at_pc

  # –3 dB bandwidth relative to DC
  dc_gain <- mag[1]
  bw_idx  <- which(mag < dc_gain - 3)[1]
  BW      <- if (!is.na(bw_idx)) omega[bw_idx] else NA_real_

  # Slope at gain crossover (dB/decade) — Bode stability check: should be ~–20
  if (gc_idx > 1 && gc_idx < length(omega)) {
    slope <- (mag[gc_idx + 1] - mag[gc_idx - 1]) /
             (log10(omega[gc_idx + 1]) - log10(omega[gc_idx - 1]))
  } else {
    slope <- NA_real_
  }

  stable <- (PM > 0) && (GM > 0)

  list(
    phase_margin_deg      = round(PM,    1),
    gain_margin_dB        = round(GM,    2),
    gc_freq               = gc_freq,
    pc_freq               = pc_freq,
    bandwidth             = BW,
    crossover_slope_dB_dec= round(slope, 1),
    stable                = stable
  )
}

# =============================================================================
# SECTION 5 — BODE PLOT
# =============================================================================
# Produces a two-panel Bode plot (magnitude / phase) with stability annotations.
#
# Arguments:
#   fr_ol       : open-loop freq_response() list
#   fr_cl       : closed-loop freq_response() list (or NULL to omit)
#   title_text  : character string for plot title
#   x_label     : axis label for frequency (default "Frequency (rad/min)")

bode_plot <- function(fr_ol, fr_cl = NULL,
                      title_text = "Bode Plot",
                      x_label    = "Frequency (rad/min)") {

  m <- stability_metrics(fr_ol)

  # -- Build data frames --
  df_ol_mag <- data.frame(freq = fr_ol$omega, val = fr_ol$mag_dB,   loop = "Open-loop")
  df_ol_ph  <- data.frame(freq = fr_ol$omega, val = fr_ol$phase_deg, loop = "Open-loop")

  col_ol <- "#2E86AB"   # blue
  col_cl <- "#C84B31"   # red
  col_gc <- "#A23B72"   # crossover marker

  # -- Magnitude panel --
  p_mag <- ggplot(df_ol_mag, aes(x = freq, y = val)) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "black") +
    geom_line(aes(color = "Open-loop"), linewidth = 1) +
    geom_vline(xintercept = m$gc_freq, linetype = "dotted",
               color = col_gc, linewidth = 0.7) +
    scale_color_manual(name = NULL,
                       values = c("Open-loop" = col_ol, "Closed-loop" = col_cl)) +
    scale_x_log10() +
    labs(title = title_text, x = x_label, y = "Magnitude (dB)") +
    theme_minimal(base_size = 11) +
    theme(
      plot.title        = element_text(hjust = 0.5, face = "bold"),
      legend.position   = "top",
      panel.grid.minor  = element_line(color = "grey92", linewidth = 0.3)
    )

  if (!is.null(fr_cl)) {
    df_cl_mag <- data.frame(freq = fr_cl$omega, val = fr_cl$mag_dB, loop = "Closed-loop")
    p_mag <- p_mag +
      geom_line(data = df_cl_mag, aes(color = "Closed-loop"),
                linewidth = 1, linetype = "dashed")
  }

  # Stability annotation
  stab_col   <- if (m$stable) "#1A7A4A" else "#CC3311"
  stab_label <- if (m$stable) "[STABLE]"  else "[UNSTABLE]"
  p_mag <- p_mag +
    annotate("label",
             x     = min(fr_ol$omega) * 3,
             y     = min(fr_ol$mag_dB) * 0.85,
             label = paste0(stab_label, "\n",
                            "PM = ",   m$phase_margin_deg,       "\u00b0\n",
                            "GM = ",   m$gain_margin_dB,          " dB\n",
                            "Slope \u2248 ", m$crossover_slope_dB_dec, " dB/dec"),
             hjust = 0, vjust = 0, size = 3.2,
             color = stab_col, fill = "white",
             label.padding = unit(0.3, "lines"), fontface = "bold")

  # -- Phase panel --
  p_ph <- ggplot(df_ol_ph, aes(x = freq, y = val)) +
    geom_hline(yintercept = -180, linetype = "dotted",
               color = "grey50", linewidth = 0.5) +
    geom_line(aes(color = "Open-loop"), linewidth = 1) +
    geom_vline(xintercept = m$gc_freq, linetype = "dotted",
               color = col_gc, linewidth = 0.7) +
    scale_color_manual(name = NULL,
                       values = c("Open-loop" = "#F18F01", "Closed-loop" = col_cl)) +
    scale_x_log10() +
    labs(x = x_label, y = "Phase (degrees)") +
    theme_minimal(base_size = 11) +
    theme(
      legend.position  = "top",
      panel.grid.minor = element_line(color = "grey92", linewidth = 0.3)
    ) +
    annotate("text",
             x = m$gc_freq * 2, y = m$phase_deg_at_gc %||% (fr_ol$phase_deg[which.min(abs(fr_ol$mag_dB))]),
             label  = sprintf("PM = %.1f\u00b0", m$phase_margin_deg),
             hjust  = 0, size = 3.2, fontface = "bold",
             color  = stab_col)

  if (!is.null(fr_cl)) {
    df_cl_ph <- data.frame(freq = fr_cl$omega, val = fr_cl$phase_deg, loop = "Closed-loop")
    p_ph <- p_ph +
      geom_line(data = df_cl_ph, aes(color = "Closed-loop"),
                linewidth = 1, linetype = "dashed")
  }

  p_mag / p_ph + plot_layout(heights = c(1.5, 1))
}

# Null-coalescing helper (base R doesn't have %||%)
`%||%` <- function(a, b) if (!is.null(a)) a else b

# =============================================================================
# SECTION 6 — MCR / (omegan, zeta) HELPERS
# =============================================================================
# These convert endocrine parameters directly to transfer function polynomials.
# Use these if you prefer not to write num/den by hand.
#
# Second-order denominator (shared by both types):
#   den(s) = s^2 + k2*s + k3*k4   where k2 = 2*zeta*omegan, k3=k4=omegan
#
# G-type (feed-forward, phase lead — acute signalling hormones):
#   num(s) = k3*(s + k2)  =  k3*s + k3*k2
#
# H-type (feedback, pure lag — maintenance/integrating hormones):
#   num(s) = k2*k4

tf_G <- function(omegan, zeta) {
  k2  <- 2 * zeta * omegan
  k3  <- omegan
  k4  <- omegan
  list(
    num = c(k3, k3 * k2),          # k3*s + k3*k2
    den = c(1,  k2, k3 * k4),      # s^2 + k2*s + omegan^2
    k2  = k2, k3 = k3, k4 = k4
  )
}

tf_H <- function(omegan, zeta) {
  k2  <- 2 * zeta * omegan
  k3  <- omegan
  k4  <- omegan
  list(
    num = c(k2 * k4),              # constant numerator
    den = c(1,  k2, k3 * k4),      # s^2 + k2*s + omegan^2
    k2  = k2, k3 = k3, k4 = k4
  )
}

# Cascade two transfer functions (polynomial multiplication)
tf_cascade <- function(tf1, tf2) {
  list(
    num = convolve(tf1$num, rev(tf2$num), type = "open"),
    den = convolve(tf1$den, rev(tf2$den), type = "open")
  )
}

# =============================================================================
# SECTION 7 — FREQUENCY GRID
# =============================================================================
# Covers T4 (~7-day half-life) through epinephrine (~2-min half-life).
# Adjust freq_range to suit your axis of interest.

omega <- 10^seq(-5, 1, length.out = 3000)   # rad/min

# =============================================================================
# SECTION 8 — WORKED EXAMPLES
# =============================================================================

cat("=== Bode Plot Reference Script ===\n\n")

# -----------------------------------------------------------------------------
# EXAMPLE A — Direct polynomial entry
# -----------------------------------------------------------------------------
# Simple second-order system: H(s) = 1 / (s^2 + 1.4s + 1)
# omegan = 1 rad/min, zeta = 0.7  (generic reference)

cat("Example A: Direct polynomial entry (generic 2nd-order)\n")

num_A <- c(1)            # numerator:   1
den_A <- c(1, 1.4, 1)   # denominator: s^2 + 1.4s + 1

fr_ol_A <- freq_response(num_A, den_A, omega, dc_norm = TRUE)
fr_cl_A <- closed_loop_from_ol(fr_ol_A)
m_A     <- stability_metrics(fr_ol_A)

cat(sprintf("  Stable = %s | PM = %.1f deg | GM = %.1f dB\n\n",
            m_A$stable, m_A$phase_margin_deg, m_A$gain_margin_dB))

plot_A <- bode_plot(fr_ol_A, fr_cl_A,
                    title_text = "Example A — Direct Polynomial (2nd-order, zeta=0.7)",
                    x_label    = "Frequency (rad/min)")

# -----------------------------------------------------------------------------
# EXAMPLE B — HPA axis via MCR helpers (paper worked example)
# -----------------------------------------------------------------------------
# ACTH     : omegan = 0.0693 rad/min, zeta = 0.5, G-type
# Cortisol : omegan = 0.0077 rad/min, zeta = 0.7, H-type

cat("Example B: HPA axis — ACTH (G) -> Cortisol (H)\n")
cat("  Parameters from Boroujerdi Table 1 & 2 (omega_n = MCR in rad/min)\n")

acth     <- tf_G(omegan = 0.0693, zeta = 0.5)
cortisol <- tf_H(omegan = 0.0077, zeta = 0.7)
hpa_tf   <- tf_cascade(acth, cortisol)

cat(sprintf("  HPA open-loop num: [%s]\n",
            paste(round(hpa_tf$num, 6), collapse = ", ")))
cat(sprintf("  HPA open-loop den: [%s]\n",
            paste(round(hpa_tf$den, 6), collapse = ", ")))

fr_ol_B <- freq_response(hpa_tf$num, hpa_tf$den, omega, dc_norm = TRUE)
fr_cl_B <- closed_loop_from_ol(fr_ol_B)
m_B     <- stability_metrics(fr_ol_B)

cat(sprintf("  Stable = %s | PM = %.1f deg | GM = %.1f dB | BW = %.5f rad/min\n\n",
            m_B$stable, m_B$phase_margin_deg,
            m_B$gain_margin_dB, m_B$bandwidth))

plot_B <- bode_plot(fr_ol_B, fr_cl_B,
                    title_text = "Example B — HPA Axis: ACTH (G-type) \u2192 Cortisol (H-type)",
                    x_label    = "Frequency (rad/min)")

# -----------------------------------------------------------------------------
# EXAMPLE C — H-type only (single cortisol compartment, for reference)
# -----------------------------------------------------------------------------
# This is the note in the paper that H-type alone for cortisol shows
# pure integrating / lag behaviour.  Provided here so the user can see
# what H-type looks like in isolation.

cat("Example C: Cortisol H-type compartment alone\n")

fr_ol_C <- freq_response(cortisol$num, cortisol$den, omega, dc_norm = TRUE)
fr_cl_C <- closed_loop_from_ol(fr_ol_C)
m_C     <- stability_metrics(fr_ol_C)

cat(sprintf("  Stable = %s | PM = %.1f deg | GM = %.1f dB\n\n",
            m_C$stable, m_C$phase_margin_deg, m_C$gain_margin_dB))

plot_C <- bode_plot(fr_ol_C, fr_cl_C,
                    title_text = "Example C — Cortisol H-type alone (pure lag)",
                    x_label    = "Frequency (rad/min)")

# =============================================================================
# SECTION 9 — SAVE OUTPUTS
# =============================================================================

out_dir <- "bode_reference_outputs"
if (!dir.exists(out_dir)) dir.create(out_dir)

ggsave(file.path(out_dir, "ExA_generic_2ndorder.png"),  plot_A, width=9, height=7, dpi=300)
ggsave(file.path(out_dir, "ExB_HPA_ACTH_Cortisol.png"), plot_B, width=9, height=7, dpi=300)
ggsave(file.path(out_dir, "ExC_Cortisol_H_alone.png"),  plot_C, width=9, height=7, dpi=300)

# Combined reference panel
fig_combined <- (plot_A | plot_B | plot_C) +
  plot_annotation(
    title    = "Bode Plot Reference — Open-loop (blue) & Closed-loop (red dashed)",
    subtitle = "A: direct polynomial  |  B: HPA cascade via MCR helpers  |  C: H-type alone",
    theme    = theme(
      plot.title    = element_text(hjust=0.5, face="bold", size=13),
      plot.subtitle = element_text(hjust=0.5, size=10, color="grey40")
    )
  )

ggsave(file.path(out_dir, "Reference_panel.png"),
       fig_combined, width=22, height=7, dpi=300)

cat(sprintf("Outputs saved to: %s/\n", out_dir))

# =============================================================================
# QUICK-REFERENCE SUMMARY
# =============================================================================
cat("\n")
cat("=================================================================\n")
cat(" QUICK REFERENCE\n")
cat("=================================================================\n")
cat(" STEP 1 — Define your transfer function\n")
cat("   Option (a) Direct polynomials:\n")
cat("     num <- c(...)   # descending powers of s\n")
cat("     den <- c(...)\n")
cat("   Option (b) From MCR / endocrine params:\n")
cat("     tf  <- tf_G(omegan, zeta)   # G-type: acute signalling\n")
cat("     tf  <- tf_H(omegan, zeta)   # H-type: maintenance\n")
cat("     num <- tf$num ; den <- tf$den\n")
cat("   Option (c) Cascade two elements:\n")
cat("     tf12 <- tf_cascade(tf1, tf2)\n")
cat("\n")
cat(" STEP 2 — Compute frequency response\n")
cat("   fr_ol <- freq_response(num, den, omega, dc_norm=TRUE)\n")
cat("\n")
cat(" STEP 3 — Closed-loop\n")
cat("   fr_cl <- closed_loop_from_ol(fr_ol)\n")
cat("\n")
cat(" STEP 4 — Stability metrics\n")
cat("   m <- stability_metrics(fr_ol)\n")
cat("   # m$phase_margin_deg, m$gain_margin_dB, m$stable, m$bandwidth\n")
cat("\n")
cat(" STEP 5 — Plot\n")
cat("   bode_plot(fr_ol, fr_cl, title_text='My axis')\n")
cat("=================================================================\n")
