#!/usr/bin/env Rscript
# =============================================================================
# endogenous_scaling_rerun.R
#
# Reruns ONLY the T-scaling experiment from endogenous_rsm_toeplitz_factor.R
# with higher resolution and more reps:
#   K_BY_scale = 1   (every k, not every 4th)
#   N_REP_scale = 500
#   K_RS_scale  = 500
#
# Saves to RSM_endogenous_results/scaling_v2/
# =============================================================================

setwd("/Users/boris/Desktop/RSShapley/Combination_Monte/Final/codess/JMP")

# ---- Load config + function definitions only (stop before main run block) ---
src_lines <- readLines("endogenous_rsm_toeplitz_factor.R")
cutoff    <- which(grepl("^# 9\\. Main run block", src_lines))[1] - 1L
cat(sprintf("Sourcing lines 1-%d (config + functions) from main script\n", cutoff))
eval(parse(text = paste(src_lines[seq_len(cutoff)], collapse = "\n")),
     envir = .GlobalEnv)

# ---- Override output directory ----------------------------------------------
OUTDIR_V2 <- file.path(getwd(), "RSM_endogenous_results", "scaling_v2")
dir.create(OUTDIR_V2, showWarnings = FALSE, recursive = TRUE)

# ---- Override T-scaling config with high-resolution settings ----------------
T_SCALING_CONFIG$N_REP_scale  <- 500L
T_SCALING_CONFIG$K_RS_scale   <- 500L
T_SCALING_CONFIG$K_BY_scale   <- 1L    # every single k value

cat("\n=== T-scaling rerun ===\n")
cat(sprintf("N_REP_scale : %d\n", T_SCALING_CONFIG$N_REP_scale))
cat(sprintf("K_RS_scale  : %d\n", T_SCALING_CONFIG$K_RS_scale))
cat(sprintf("K_BY_scale  : %d  (step size in k grid)\n", T_SCALING_CONFIG$K_BY_scale))
cat(sprintf("T_len_grid  : %s\n", paste(T_SCALING_CONFIG$T_len_grid, collapse = ", ")))
cat(sprintf("Output      : %s\n\n", OUTDIR_V2))

sink(file.path(OUTDIR_V2, "run_log.txt"), split = TRUE)
cat("Started:", format(Sys.time()), "\n\n")

scaling_results <- run_scaling_T_endo(T_SCALING_CONFIG, OUTDIR_V2)

cat("\nFinished:", format(Sys.time()), "\n")
sink()

cat("\nAll done. Results in:", OUTDIR_V2, "\n")
