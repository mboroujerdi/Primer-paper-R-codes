###############################################################################
# Endocrine Cascade Bode Plot — Open-Loop & Closed-Loop Stability Analysis
#
# Each element in the cascade is defined directly by:
#   omegan  — natural frequency (rad/s)
#   zeta    — damping ratio
#
# Two transfer-function types per element:
#   G-type (band-pass / lead):   G(s) = k3*(s + k2) / [s*(s+k2) + k3*k4]
#   H-type (low-pass / lag):     H(s) = k2*k4       / [s*(s+k2) + k3*k4]
#
# where  k2 = 2*zeta*omegan,  k3 = k4 = omegan  (second-order factoring)
#
# Optional zeros can be appended to the open-loop transfer function to ensure
# the Bode gain plot crosses 0 dB at –20 dB/decade (Bode stability condition).
#
# ENDOCRINE EXAMPLES
#   • HPT  axis : Hypothalamus → Pituitary → Thyroid  (TRH → TSH → T4/T3)
#   • HPA  axis : Hypothalamus → Pituitary → Adrenal  (CRH → ACTH → Cortisol)
#   • HPG  axis : Hypothalamus → Pituitary → Gonad    (GnRH → LH/FSH → E2/T)
#   • Mixed demo: H + G + H + zero correction
#
# Closed-loop transfer function:  CL(s) = OL(s) / (1 + OL(s))
# Stability check: gain margin > 0 dB  AND  phase margin > 0°
###############################################################################

library(ggplot2)
library(patchwork)
library(pracma)
library(scales)

# ============================================================
# 1.  LOW-LEVEL TRANSFER FUNCTIONS
# ============================================================

# Convert (omegan, zeta) to rate constants
params_from_omegaz <- function(omegan, zeta) {
  k2 <- 2 * zeta * omegan
  k3 <- omegan
  k4 <- omegan
  list(k2 = k2, k3 = k3, k4 = k4)
}

# G-type element  (band-pass / lead character)
tf_G <- function(k2, k3, k4, s_vec) {
  k3 * (s_vec + k2) / (s_vec * (s_vec + k2) + k3 * k4)
}

# H-type element  (low-pass / lag character)
tf_H <- function(k2, k3, k4, s_vec) {
  (k2 * k4) / (s_vec * (s_vec + k2) + k3 * k4)
}

# A real zero at -z_val adds factor (s + z_val)
tf_zero <- function(z_val, s_vec) {
  (s_vec + z_val)
}

# ============================================================
# 2.  CASCADE OPEN-LOOP RESPONSE
# ============================================================
# elements  : list of lists, each with fields:
#               type   = "G" | "H"
#               omegan = natural frequency (rad/s)
#               zeta   = damping ratio
# zeros     : numeric vector of zero locations (rad/s); NULL = none
# dc_norm   : if TRUE, normalise so |H(j*0)| = 1  (0 dB at DC)

cascade_open_loop <- function(elements, zeros = NULL, omega_vec, dc_norm = FALSE) {

  if (length(elements) < 1 || length(elements) > 4)
    stop("Cascade must have 1-4 elements")

  s_vec <- 1i * omega_vec
  H_total <- rep(1 + 0i, length(omega_vec))

  for (el in elements) {
    p  <- params_from_omegaz(el$omegan, el$zeta)
    H_total <- H_total * switch(el$type,
      "G" = tf_G(p$k2, p$k3, p$k4, s_vec),
      "H" = tf_H(p$k2, p$k3, p$k4, s_vec),
      stop(paste("Unknown element type:", el$type))
    )
  }

  # Append zeros
  if (!is.null(zeros)) {
    for (z in zeros) {
      H_total <- H_total * tf_zero(z, s_vec)
    }
  }

  # Optional DC normalisation
  if (dc_norm) {
    dc_val <- Mod(H_total[1])
    if (dc_val > 0) H_total <- H_total / dc_val
  }

  mag_db       <- 20 * log10(Mod(H_total))
  phase_rad    <- Arg(H_total)
  phase_unwrap <- unwrap_phase_vec(phase_rad)
  phase_deg    <- phase_unwrap * 180 / pi

  list(H = H_total, magnitude = mag_db, phase = phase_deg)
}

# ============================================================
# 3.  CLOSED-LOOP RESPONSE
# ============================================================
cascade_closed_loop <- function(H_ol) {
  # Unity negative feedback:  CL = OL / (1 + OL)
  H_cl  <- H_ol / (1 + H_ol)
  mag_db   <- 20 * log10(Mod(H_cl))
  phase_rad <- Arg(H_cl)
  phase_deg <- unwrap_phase_vec(phase_rad) * 180 / pi
  list(H = H_cl, magnitude = mag_db, phase = phase_deg)
}

# ============================================================
# 4.  STABILITY METRICS
# ============================================================
stability_metrics <- function(mag_db, phase_deg, omega_vec) {

  # --- Gain crossover (0 dB) ---
  gc_idx  <- which.min(abs(mag_db))
  gc_freq <- omega_vec[gc_idx]
  phase_at_gc <- phase_deg[gc_idx]
  phase_margin <- 180 + phase_at_gc          # PM > 0 ⟹ stable

  # --- Phase crossover (–180°) ---
  pc_idx  <- which.min(abs(phase_deg + 180))
  pc_freq <- omega_vec[pc_idx]
  gain_at_pc <- mag_db[pc_idx]
  gain_margin <- -gain_at_pc                  # GM > 0 ⟹ stable

  # -3 dB bandwidth
  dc_gain <- mag_db[1]
  bw_idx  <- which(mag_db < dc_gain - 3)[1]
  bandwidth <- if (!is.na(bw_idx)) omega_vec[bw_idx] else NA

  # Slope at gain crossover  (dB / decade)
  if (gc_idx > 1 && gc_idx < length(omega_vec)) {
    slope <- (mag_db[gc_idx + 1] - mag_db[gc_idx - 1]) /
             (log10(omega_vec[gc_idx + 1]) - log10(omega_vec[gc_idx - 1]))
  } else {
    slope <- NA
  }

  stable    <- (gain_margin > 0) && (phase_margin > 0)
  slope_ok  <- !is.na(slope) && (slope > -25) && (slope < -15)  # near –20 dB/dec

  list(
    gc_freq        = gc_freq,
    phase_at_gc    = round(phase_at_gc, 1),
    phase_margin   = round(phase_margin, 1),
    pc_freq        = pc_freq,
    gain_at_pc     = round(gain_at_pc, 2),
    gain_margin    = round(gain_margin, 2),
    bandwidth      = bandwidth,
    crossover_slope_dB_dec = round(slope, 1),
    stable         = stable,
    slope_ok       = slope_ok
  )
}

# ============================================================
# 5.  UTILITY: phase unwrapping
# ============================================================
unwrap_phase_vec <- function(phase_rad) {
  n <- length(phase_rad)
  if (n <= 1) return(phase_rad)
  out <- numeric(n)
  out[1] <- phase_rad[1]
  for (i in 2:n) {
    d <- phase_rad[i] - phase_rad[i - 1]
    d <- d - 2 * pi * round(d / (2 * pi))
    out[i] <- out[i - 1] + d
  }
  out
}

# ============================================================
# 6.  PLOTTING HELPERS
# ============================================================

theme_endocrine <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40"),
    axis.title   = element_text(face = "bold", size = 11),
    axis.text    = element_text(face = "bold", size = 9),
    legend.text  = element_text(size = 9),
    panel.grid.minor = element_line(color = "gray92", linewidth = 0.3),
    panel.grid.major = element_line(color = "gray82", linewidth = 0.5)
  )
}

plot_magnitude <- function(df_ol, df_cl = NULL, metrics, title_text, x_label,
                           show_stability = TRUE) {

  p <- ggplot(df_ol, aes(x = frequency, y = magnitude)) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.4) +
    geom_line(aes(color = "Open-loop"), linewidth = 1) +
    scale_color_manual(name = NULL,
                       values = c("Open-loop" = "#2E86AB", "Closed-loop" = "#C84B31")) +
    scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
    labs(x = x_label, y = "Magnitude (dB)", title = title_text) +
    theme_endocrine()

  if (!is.null(df_cl)) {
    p <- p + geom_line(data = df_cl, aes(color = "Closed-loop"), linewidth = 1,
                       linetype = "dashed")
  }

  # Gain crossover marker
  p <- p + geom_vline(xintercept = metrics$gc_freq,
                      linetype = "dotted", color = "#A23B72", linewidth = 0.7)

  # Stability annotation box
  if (show_stability) {
    stab_label  <- if (metrics$stable) "[STABLE]" else "[UNSTABLE]"
    stab_colour <- if (metrics$stable) "#1A7A4A" else "#CC3311"
    slope_label <- sprintf("Slope ~ %.0f dB/dec", metrics$crossover_slope_dB_dec)

    p <- p +
      annotate("label",
               x = min(df_ol$frequency) * 3,
               y = min(df_ol$magnitude) * 0.85,
               label = paste0(stab_label, "\n",
                              "PM = ", metrics$phase_margin, "°\n",
                              "GM = ", metrics$gain_margin, " dB\n",
                              slope_label),
               hjust = 0, vjust = 0, size = 3, color = stab_colour,
               fill = "white", label.padding = unit(0.3, "lines"), fontface = "bold")
  }

  p
}

plot_phase <- function(df_ol, df_cl = NULL, metrics, x_label) {
  p <- ggplot(df_ol, aes(x = frequency, y = phase)) +
    geom_hline(yintercept = -180, linetype = "dotted", color = "gray50", linewidth = 0.5) +
    geom_line(aes(color = "Open-loop"), linewidth = 1) +
    scale_color_manual(name = NULL,
                       values = c("Open-loop" = "#F18F01", "Closed-loop" = "#C84B31")) +
    scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
    labs(x = x_label, y = "Phase (degrees)") +
    theme_endocrine()

  if (!is.null(df_cl)) {
    p <- p + geom_line(data = df_cl, aes(color = "Closed-loop"), linewidth = 1,
                       linetype = "dashed")
  }

  p <- p + geom_vline(xintercept = metrics$gc_freq,
                      linetype = "dotted", color = "#A23B72", linewidth = 0.7)

  # Phase margin annotation
  pm_colour <- if (metrics$phase_margin > 0) "#1A7A4A" else "#CC3311"
  p <- p +
    annotate("text",
             x = metrics$gc_freq * 2,
             y = metrics$phase_at_gc,
             label = sprintf("PM = %.1f°", metrics$phase_margin),
             hjust = 0, size = 3.2, color = pm_colour, fontface = "bold")

  p
}

# ============================================================
# 7.  FULL BODE PANEL (open-loop + closed-loop, mag + phase)
# ============================================================
bode_panel <- function(elements, zeros = NULL, omega_vec,
                       title_text, subtitle_text = "",
                       x_label = "Frequency (rad/s)") {

  # Open-loop  (DC-normalised so 0 dB crossing is always in-band)
  ol <- cascade_open_loop(elements, zeros, omega_vec, dc_norm = TRUE)
  df_ol_mag <- data.frame(frequency = omega_vec, magnitude = ol$magnitude)
  df_ol_ph  <- data.frame(frequency = omega_vec, phase    = ol$phase)

  # Closed-loop
  cl <- cascade_closed_loop(ol$H)
  df_cl_mag <- data.frame(frequency = omega_vec, magnitude = cl$magnitude)
  df_cl_ph  <- data.frame(frequency = omega_vec, phase    = cl$phase)

  # Metrics from open-loop
  m <- stability_metrics(ol$magnitude, ol$phase, omega_vec)

  # Compose title with subtitle
  full_title <- if (nchar(subtitle_text) > 0)
    bquote(bold(.(title_text)) ~ "\n" ~ .(subtitle_text))
  else
    title_text

  p_mag <- plot_magnitude(df_ol_mag, df_cl_mag, m,
                          title_text, x_label, show_stability = TRUE)
  p_ph  <- plot_phase   (df_ol_ph,  df_cl_ph,  m, x_label)

  list(panel  = p_mag / p_ph + plot_layout(heights = c(1.5, 1)),
       metrics = m)
}

# ============================================================
# 8.  ENDOCRINE CONFIGURATIONS
# ============================================================
# All parameters derived from Tables 1 & 2 of:
#   Boroujerdi M. "A control theoretic primer for systems endocrinology"
#   Frontiers in Endocrinology (submitted).
#
# Key mapping (first approximation, scaling factor = 1):
#   omega_n  = MCR  (rad/min)
#   k2       = 2 * zeta * omega_n
#   k3 = k4  = omega_n
#
# Architecture assignment (Table 2):
#   G-type: acute signalling hormones (phase lead; rapid secretion into circulation)
#   H-type: maintenance/integrating hormones (pure lag; endocytosis-mediated)
#
# All frequencies are in rad/min to match physiological half-life units.
# -------------------------------------------------------------

# --- HPA axis (PRIMARY WORKED EXAMPLE from paper, Section: Worked Example) ---
#
# Two-element cascade: ACTH (G-type) --> Cortisol (H-type)
# with long-loop negative feedback from cortisol.
#
# ACTH:    MCR = 0.0693 min^-1, omega_n = 0.0693 rad/min, zeta = 0.5, G-type
#          (Table 1: MCR_min = 0.0693; Table 2: zeta = 0.5, suggested_type = G)
#          half-life ~10 min; rapid pituitary signalling
#
# Cortisol: MCR = 0.0077 min^-1, omega_n = 0.0077 rad/min, zeta = 0.7, H-type
#           (Table 1: MCR_min = 0.0077; Table 2: zeta = 0.7, suggested_type = H)
#           half-life ~90 min; stable maintenance, strong negative feedback
#
# Predicted ultradian period ~60-90 min; dexamethasone suppression ~2-4 h.

HPA_elements <- list(
  list(type = "G", omegan = 0.0693, zeta = 0.5, label = "Pituitary (ACTH)"),
  list(type = "H", omegan = 0.0077, zeta = 0.7, label = "Adrenal (Cortisol)")
)
# No external zero required: G-type ACTH already provides phase lead.
# The G-type zero is intrinsic at s = -k2 = -2*0.5*0.0693 = -0.0693 rad/min.
HPA_zeros <- NULL

# --- HPT axis (Hypothalamus–Pituitary–Thyroid) ---
#
# Three-element cascade: TSH (H-type) --> T4/Thyroxine (H-type) --> T3 (H-type)
# (CRH/TRH not separately listed in Table 1; cascade modelled from pituitary onward)
#
# TSH:    MCR = 0.01160 min^-1, omega_n = 0.01160 rad/min, zeta = 0.7, H-type
#         (Table 1: MCR_min = 0.0116; Table 2: zeta = 0.7, type = H)
#         half-life ~60 min; pituitary integrator
#
# T4:     MCR = 0.00007 min^-1, omega_n = 0.00007 rad/min, zeta = 0.9, H-type
#         (Table 1: MCR_min = 0.00007; Table 2: zeta = 0.9, type = H)
#         half-life ~7 days (~10080 min); slowest element in endocrine table
#
# T3:     MCR = 0.00048 min^-1, omega_n = 0.00048 rad/min, zeta = 0.8, H-type
#         (Table 1: MCR_min = 0.00048; Table 2: zeta = 0.8, type = H)
#         half-life ~1 day (~1440 min); active form, faster than T4
#
# All three are H-type (pure lag), reflecting slow maintenance regulation.
# A stabilising zero is added near the gain crossover frequency to achieve
# -20 dB/dec crossing and adequate phase margin (see paper, Stability Analysis).

HPT_elements <- list(
  list(type = "H", omegan = 0.01160, zeta = 0.7, label = "Pituitary (TSH)"),
  list(type = "H", omegan = 0.00007, zeta = 0.9, label = "Thyroid (T4/Thyroxine)"),
  list(type = "H", omegan = 0.00048, zeta = 0.8, label = "Peripheral (T3/Triiodothyronine)")
)
# Stabilising zero placed near the geometric mean of TSH and T3 omega_n values
# to recover phase lead and ensure -20 dB/dec gain crossover.
HPT_zeros <- c(sqrt(0.01160 * 0.00048))   # ~ 0.00236 rad/min

# --- HPG axis (Hypothalamus–Pituitary–Gonad, female) ---
#
# Three-element cascade: LH (G-type) --> Estradiol (H-type) --> Progesterone (H-type)
# (GnRH not separately listed; cascade modelled from gonadotropin onward)
#
# LH:           MCR = 0.02310 min^-1, omega_n = 0.02310 rad/min, zeta = 0.5, G-type
#               (Table 1: MCR_min = 0.0231; Table 2: zeta = 0.5, type = G)
#               half-life ~30 min; pulsatile gonadotropin
#
# Estradiol:    MCR = 0.03470 min^-1, omega_n = 0.03470 rad/min, zeta = 0.7, H-type
#               (Table 1: MCR_min = 0.0347; Table 2: zeta = 0.7, type = H)
#               half-life ~20 min; cyclic HPG axis, female
#
# Progesterone: MCR = 0.00200 min^-1, omega_n = 0.00200 rad/min, zeta = 0.7, H-type
#               (Table 1: MCR_min = 0.0020; Table 2: zeta = 0.7, type = H)
#               half-life ~347 min; luteal phase maintenance
#
# LH is G-type (pulsatile, acute signalling); Estradiol and Progesterone are H-type.
# No additional external zero needed: G-type LH provides inherent phase lead.

HPG_elements <- list(
  list(type = "G", omegan = 0.02310, zeta = 0.5, label = "Pituitary (LH)"),
  list(type = "H", omegan = 0.03470, zeta = 0.7, label = "Gonad (Estradiol)"),
  list(type = "H", omegan = 0.00200, zeta = 0.7, label = "Luteal (Progesterone)")
)
HPG_zeros <- NULL

# ============================================================
# 9.  FREQUENCY GRID
# ============================================================
# Units: rad/min  (consistent with MCR values from Table 1)
# Range spans from T4 dynamics (~0.00007 rad/min) to
# epinephrine-class dynamics (~0.35 rad/min), covering all 36
# hormones in Table 1.

freq_range <- c(1e-5, 1e1)   # rad/min
n_pts      <- 3000
omega_vec  <- 10^seq(log10(freq_range[1]), log10(freq_range[2]), length.out = n_pts)

# ============================================================
# 10.  GENERATE ALL PANELS
# ============================================================
cat("=== Endocrine Cascade Bode Plot Generator ===\n")
cat("    Parameters from Boroujerdi (Tables 1 & 2)\n")
cat("    All frequencies in rad/min\n\n")

# HPA: Primary worked example from paper
cat("Processing HPA axis (primary worked example)...\n")
hpa <- bode_panel(HPA_elements, HPA_zeros, omega_vec,
                  title_text    = "HPA Axis — Worked Example",
                  subtitle_text = "ACTH (G-type, \u03c9n=0.0693) \u2192 Cortisol (H-type, \u03c9n=0.0077)",
                  x_label       = "Frequency (rad/min)")
cat(sprintf("  HPA  stable = %s | PM = %.1f\u00b0 | GM = %.1f dB | slope = %.0f dB/dec\n",
            hpa$metrics$stable, hpa$metrics$phase_margin,
            hpa$metrics$gain_margin, hpa$metrics$crossover_slope_dB_dec))

# HPT: Three H-type elements
cat("Processing HPT axis...\n")
hpt <- bode_panel(HPT_elements, HPT_zeros, omega_vec,
                  title_text    = "HPT Axis",
                  subtitle_text = "TSH (H, \u03c9n=0.0116) \u2192 T4 (H, \u03c9n=7e-5) \u2192 T3 (H, \u03c9n=4.8e-4) + stabilising zero",
                  x_label       = "Frequency (rad/min)")
cat(sprintf("  HPT  stable = %s | PM = %.1f\u00b0 | GM = %.1f dB | slope = %.0f dB/dec\n",
            hpt$metrics$stable, hpt$metrics$phase_margin,
            hpt$metrics$gain_margin, hpt$metrics$crossover_slope_dB_dec))

# HPG: Female reproductive axis
cat("Processing HPG axis (female)...\n")
hpg <- bode_panel(HPG_elements, HPG_zeros, omega_vec,
                  title_text    = "HPG Axis (Female)",
                  subtitle_text = "LH (G-type, \u03c9n=0.0231) \u2192 Estradiol (H, \u03c9n=0.0347) \u2192 Progesterone (H, \u03c9n=0.002)",
                  x_label       = "Frequency (rad/min)")
cat(sprintf("  HPG  stable = %s | PM = %.1f\u00b0 | GM = %.1f dB | slope = %.0f dB/dec\n",
            hpg$metrics$stable, hpg$metrics$phase_margin,
            hpg$metrics$gain_margin, hpg$metrics$crossover_slope_dB_dec))

# ============================================================
# 11.  COMBINED FIGURES
# ============================================================

# Figure 1: HPA worked example (full-width, primary figure)
fig1 <- hpa$panel +
  plot_annotation(
    title    = "HPA Axis Bode Plot — Open-loop (blue) & Closed-loop (red dashed)",
    subtitle = paste0("ACTH \u2192 Cortisol two-element cascade with long-loop negative feedback\n",
                      "Parameters from Boroujerdi Table 1 & 2 | Frequencies in rad/min"),
    theme    = theme(
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40")
    )
  )

# Figure 2: All three axes side by side
fig2 <- (hpa$panel | hpt$panel | hpg$panel) +
  plot_annotation(
    title    = "Endocrine Cascade Bode Plots — Open-loop (blue) & Closed-loop (red dashed)",
    subtitle = paste0("All \u03c9n values from MCR-derived parameters (Boroujerdi Tables 1 & 2) | Frequencies in rad/min\n",
                      "Vertical dotted line = gain-crossover frequency | Stability metrics in magnitude panel"),
    theme    = theme(
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40")
    )
  )

# ============================================================
# 12.  SAVE OUTPUTS
# ============================================================
out_dir <- "endocrine_bode_outputs"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cat("\nSaving figures...\n")

# Fig 1: HPA worked example
ggsave(file.path(out_dir, "Fig1_HPA_worked_example.png"),
       fig1, width = 9, height = 7, dpi = 300)
ggsave(file.path(out_dir, "Fig1_HPA_worked_example.pdf"),
       fig1, width = 9, height = 7)

# Fig 2: All three axes
ggsave(file.path(out_dir, "Fig2_three_endocrine_axes.png"),
       fig2, width = 20, height = 7, dpi = 300)
ggsave(file.path(out_dir, "Fig2_three_endocrine_axes.pdf"),
       fig2, width = 20, height = 7)

# Individual axis figures
for (nm in c("HPA", "HPT", "HPG")) {
  obj <- get(tolower(nm))
  ggsave(file.path(out_dir, paste0(nm, "_bode.png")),
         obj$panel, width = 9, height = 7, dpi = 300)
}

# ============================================================
# 13.  SUMMARY TABLE
# ============================================================
axes_names <- c("HPA (ACTH->Cortisol)", "HPT (TSH->T4->T3)", "HPG (LH->E2->Prog)")
axes_list  <- list(hpa$metrics, hpt$metrics, hpg$metrics)

summ_df <- do.call(rbind, lapply(seq_along(axes_names), function(i) {
  m <- axes_list[[i]]
  data.frame(
    Axis                  = axes_names[i],
    Stable                = m$stable,
    Phase_Margin_deg      = m$phase_margin,
    Gain_Margin_dB        = m$gain_margin,
    GC_Freq_rad_per_min   = signif(m$gc_freq, 3),
    Slope_dB_dec          = m$crossover_slope_dB_dec,
    Bandwidth_rad_per_min = signif(m$bandwidth, 3),
    stringsAsFactors      = FALSE
  )
}))

write.csv(summ_df, file.path(out_dir, "stability_summary.csv"), row.names = FALSE)

cat("\n=== STABILITY SUMMARY ===\n")
print(summ_df, row.names = FALSE)

cat(sprintf("\nAll outputs saved to: %s/\n", out_dir))
cat("  Fig1_HPA_worked_example.png/pdf    - HPA axis primary worked example\n")
cat("  Fig2_three_endocrine_axes.png/pdf  - HPA, HPT, HPG side-by-side panel\n")
cat("  HPA/HPT/HPG_bode.png               - individual axis plots\n")
cat("  stability_summary.csv              - gain margin, phase margin, slope\n")
cat("\nParameter source: Boroujerdi, Tables 1 & 2 (omega_n = MCR in rad/min)\n")
cat("Architecture rules (Table 2):\n")
cat("  G-type: acute signalling (ACTH, LH, GLP-1, Angiotensin II, etc.)\n")
cat("  H-type: maintenance/integrating (Cortisol, T4, T3, IGF-1, Testosterone, etc.)\n")
cat("Damping ratios (Table 2):\n")
cat("  0.5 = underdamped / oscillatory (G-type hormones, pulsatile axes)\n")
cat("  0.7 = lightly damped (general default)\n")
cat("  0.8-0.9 = well-damped / overdamped (very slow hormones: T4, T3, Vitamin D3)\n")
