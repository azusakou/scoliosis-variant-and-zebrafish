#!/usr/bin/env Rscript
# Mechanical analyses from the compact processed inputs using base R only.
# Usage: Rscript mechanics.R INPUT_DIR OUTPUT_DIR
# Original expression order, RNG kinds and seeds are retained.
# WT denotes the processed reference code; katnb1 controls are heterozygotes.
# Sampling uses animals, except the explicitly unclustered geometry segments.

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2L) stop("Usage: Rscript mechanics.R INPUT_DIR OUTPUT_DIR")

input_dir <- normalizePath(args[1], mustWork = TRUE)

output_dir <- args[2]

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read <- function(name) read.csv(file.path(input_dir, paste0(name, ".csv")), check.names = FALSE)

write <- function(x, name) write.csv(x, file.path(output_dir, paste0(name, ".csv")), row.names = FALSE, 
    na = "NA")

# Input validation and biological observation units.
s <- read("mechanics_animal_swe")

r <- read("mechanics_animal_redox")

g <- read("mechanics_geometry_segments")

c <- read("mechanics_cobb_paired")

stopifnot(nrow(s) == 209L, nrow(r) == 23L, nrow(g) == 140L, nrow(c) == 20L, all(is.finite(s$speed_m_s)), 
    all(s$speed_m_s > 0), max(abs(s$modulus_proxy_m2_s2 - s$speed_m_s^2)) < 1e-10, all(g$min_width_um > 
        0), all(g$max_width_um >= g$min_width_um), max(abs(g$width_ratio - g$max_width_um/g$min_width_um)) < 
        1e-09, all(r$roi_n == 3L), max(abs(r$ctf_mean - rowMeans(r[, c("ctf_roi1", "ctf_roi2", 
        "ctf_roi3")]))) < 1e-08, !anyDuplicated(s$animal_id), !anyDuplicated(r$animal_id), 
    !anyDuplicated(c$animal_id))

# Seeded untreated-to-treated redox calibration.
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

set.seed(260917)

B <- 10000L

gs <- c("WT", "Mut", "WT_NACET", "Mut_NACET")

s0 <- s[s$source_sheet == "Figure 5f", ]

rx <- setNames(lapply(gs, function(gr) r$ctf_mean[r$group == gr]), gs)

sx <- setNames(lapply(gs, function(gr) s0$modulus_proxy_m2_s2[s0$group == gr]), gs)

stopifnot(identical(as.integer(lengths(rx)), c(3L, 8L, 3L, 9L)), all(lengths(sx) == 15L))

bootmean <- function(x) replicate(B, mean(sample(x, length(x), replace = TRUE)))

rb <- lapply(rx, bootmean)

sb <- lapply(sx, bootmean)

rmn <- sapply(rx, mean)

smn <- sapply(sx, mean)

beta <- log(smn["Mut"]/smn["WT"])/log(rmn["Mut"]/rmn["WT"])

bb <- log(sb$Mut/sb$WT)/log(rb$Mut/rb$WT)

param <- data.frame(parameter = "redox_modulus_power_exponent", estimate = unname(beta), 
    lower = unname(quantile(bb, 0.025)), upper = unname(quantile(bb, 0.975)), bootstrap_n = B, 
    nonfinite_n = sum(!is.finite(bb)), redox_reversal_fraction = mean(rb$Mut <= rb$WT), scope = "ecological cross-cohort calibration; not individual mediation")

write(param, "mechanical_calibration_parameters")

pred <- do.call(rbind, lapply(gs, function(gr) {
    est <- smn["WT"] * (rmn[gr]/rmn["WT"])^beta
    pb <- sb$WT * (rb[[gr]]/rb$WT)^bb
    data.frame(group = gr, role = ifelse(grepl("NACET", gr), "held_out", "training"), redox_mean = unname(rmn[gr]), 
        redox_n = length(rx[[gr]]), observed = unname(smn[gr]), swe_n = length(sx[[gr]]), 
        observed_lower = unname(quantile(sb[[gr]], 0.025)), observed_upper = unname(quantile(sb[[gr]], 
            0.975)), predicted = unname(est), predicted_lower = unname(quantile(pb, 0.025, 
            na.rm = TRUE)), predicted_upper = unname(quantile(pb, 0.975, na.rm = TRUE)), 
        relative_error = unname(est/smn[gr] - 1))
}))

write(pred, "mechanical_heldout_predictions")

# Group-mean geometry calibration; segment clustering cannot be recovered.
gorder <- unique(g$group)

geo <- do.call(rbind, lapply(gorder, function(gr) {
    x <- g[g$group == gr, ]
    data.frame(group = gr, segments = nrow(x), published_animals = x$published_animal_n[1], 
        mean_ratio = mean(x$width_ratio), median_ratio = median(x$width_ratio), mean_max = mean(x$max_width_um), 
        mean_min = mean(x$min_width_um))
}))

gb <- log(smn["Mut"]/smn["WT"])/log(geo$mean_ratio[geo$group == "Mut"]/geo$mean_ratio[geo$group == 
    "WT"])

geo$predicted_modulus_proxy <- unname(smn["WT"] * (geo$mean_ratio/geo$mean_ratio[geo$group == 
    "WT"])^gb)

geo$observed_modulus_proxy <- unname(smn[geo$group])

geo$role <- ifelse(grepl("NACET", geo$group), "held_out", "training")

geo$fitted_exponent <- unname(gb)

geo$uncertainty_status <- "not estimable at animal level: original IDs absent"

write(geo, "mechanical_geometry_calibration")

# Within-source ratios, Welch comparisons and the specified BH family.
effects <- do.call(rbind, lapply(unique(s$source_sheet), function(sh) {
    q <- s[s$source_sheet == sh, ]
    groups <- unique(q$group)
    vals <- setNames(lapply(groups, function(gr) q$modulus_proxy_m2_s2[q$group == gr]), groups)
    do.call(rbind, lapply(setdiff(groups, "WT"), function(gr) {
        boot <- log(bootmean(vals[[gr]])/bootmean(vals$WT))
        tt <- t.test(log(vals[[gr]]), log(vals$WT))
        data.frame(source_sheet = sh, model = q$model[1], contrast = paste(gr, "vs WT"), 
            n_group = length(vals[[gr]]), n_control = length(vals$WT), log_modulus_ratio = log(mean(vals[[gr]])/mean(vals$WT)), 
            lower = unname(quantile(boot, 0.025)), upper = unname(quantile(boot, 0.975)), 
            p = tt$p.value)
    }))
}))

effects$q_BH <- p.adjust(effects$p, "BH")

write(effects, "mechanical_external_validation")

# Cobb-angle fits and fully refitted leave-one-animal-out comparisons.
models <- list(speed = cobb_degree ~ speed_m_s, log_modulus = cobb_degree ~ log(modulus_proxy_m2_s2))

fits <- list()

preds <- list()

for (nm in names(models)) {
    f <- lm(models[[nm]], data = c)
    cv <- vapply(seq_len(nrow(c)), function(i) unname(predict(lm(models[[nm]], data = c[-i, 
        ]), newdata = c[i, , drop = FALSE])), numeric(1))
    baseline <- vapply(seq_len(nrow(c)), function(i) mean(c$cobb_degree[-i]), numeric(1))
    ci <- confint(f)
    co <- coef(summary(f))
    fits[[nm]] <- data.frame(model = nm, n = nrow(c), intercept = unname(coef(f)[1]), slope = unname(coef(f)[2]), 
        slope_lower = ci[2, 1], slope_upper = ci[2, 2], p = co[2, 4], R2 = summary(f)$r.squared, 
        LOOCV_RMSE = sqrt(mean((cv - c$cobb_degree)^2)), LOOCV_MAE = mean(abs(cv - c$cobb_degree)), 
        baseline_LOOCV_RMSE = sqrt(mean((baseline - c$cobb_degree)^2)))
    preds[[nm]] <- data.frame(c, model = nm, fitted = unname(fitted(f)), heldout_prediction = cv)
}

write(do.call(rbind, fits), "mechanical_cobb_models")

write(do.call(rbind, preds), "mechanical_cobb_predictions")

# Four-group velocity model with genotype, treatment and interaction.
s0$mutant <- as.integer(s0$group %in% c("Mut", "Mut_NACET"))

s0$nacet <- as.integer(grepl("NACET", s0$group))

f <- lm(speed_m_s ~ mutant * nacet, data = s0)

ci <- confint(f)

interaction <- data.frame(term = names(coef(f)), estimate = unname(coef(f)), ci_low = ci[, 
    1], ci_high = ci[, 2], p = coef(summary(f))[, 4])

write(interaction, "stiffness_antioxidant_interaction")

set.seed(260917)

# Independent animal bootstrap for the velocity-difference interval.
va <- s0$speed_m_s[s0$group == "Mut"]

vb <- s0$speed_m_s[s0$group == "WT"]

vd <- replicate(B, mean(sample(va, length(va), replace = TRUE)) - mean(sample(vb, length(vb), 
    replace = TRUE)))

velocity_difference <- data.frame(source_sheet = "Figure 5f", contrast = "Mut vs untreated WT", 
    n_group = length(va), n_control = length(vb), mean_difference_m_s = mean(va) - mean(vb), 
    lower = unname(quantile(vd, 0.025)), upper = unname(quantile(vd, 0.975)), bootstrap_n = B, 
    seed = 260917, rng_design = "explicit independent seed reset for this contrast")

write(velocity_difference, "mechanical_velocity_difference")

gs <- c("WT", "Mut", "WT_NACET", "Mut_NACET")

# Comparison with the constant untreated-WT control-state baseline.
G <- setNames(vapply(gs, function(k) mean(s0$modulus_proxy_m2_s2[s0$group == k]), numeric(1)), 
    gs)

Q <- setNames(vapply(gs, function(k) mean(g$width_ratio[g$group == k]), numeric(1)), gs)

gamma <- log(G["Mut"]/G["WT"])/log(Q["Mut"]/Q["WT"])

geometry <- data.frame(group = gs, role = c("calibration", "calibration", "held_out", "held_out"), 
    segment_n = as.integer(table(factor(g$group, levels = gs))), published_animal_n = c(8, 
        9, 6, 7), geometry_ratio = as.numeric(Q), observed_proxy = as.numeric(G), geometry_prediction = as.numeric(G["WT"] * 
        (Q/Q["WT"])^gamma), constant_WT_prediction = rep(unname(G["WT"]), 4), exponent = unname(gamma))

geometry$geometry_error_percent <- 100 * (geometry$geometry_prediction/geometry$observed_proxy - 
    1)

geometry$constant_WT_error_percent <- 100 * (geometry$constant_WT_prediction/geometry$observed_proxy - 
    1)

write(geometry, "geometry_constant_baseline")

held <- geometry[geometry$role == "held_out", ]

write(data.frame(method = c("Geometry power law", "Untreated WT constant"), heldout_group_n = 2L, 
    RMSE = c(sqrt(mean((held$geometry_prediction - held$observed_proxy)^2)), sqrt(mean((held$constant_WT_prediction - 
        held$observed_proxy)^2))), mean_absolute_percent_error = c(mean(abs(held$geometry_error_percent)), 
        mean(abs(held$constant_WT_error_percent)))), "geometry_baseline_aggregate")

# Case-bootstrap uncertainty: hold out every copy of each animal identity.
cv_weighted <- function(x, y, w) {
    n <- sum(w)
    swx <- sum(w * x)
    swy <- sum(w * y)
    swxx <- sum(w * x * x)
    swxy <- sum(w * x * y)
    keep <- w > 0
    nt <- n - w[keep]
    xx <- swx - w[keep] * x[keep]
    yy <- swy - w[keep] * y[keep]
    xx2 <- swxx - w[keep] * x[keep]^2
    xy <- swxy - w[keep] * x[keep] * y[keep]
    denom <- xx2 - xx^2/nt
    if (any(nt <= 2) | any(denom <= 1e-12)) 
        return(c(RMSE = NA, baseline_RMSE = NA, improvement = NA, MAE = NA, distinct = sum(keep)))
    slope <- (xy - xx * yy/nt)/denom
    pred <- yy/nt + slope * (x[keep] - xx/nt)
    base <- yy/nt
    rm <- sqrt(sum(w[keep] * (pred - y[keep])^2)/n)
    bm <- sqrt(sum(w[keep] * (base - y[keep])^2)/n)
    c(RMSE = rm, baseline_RMSE = bm, improvement = bm - rm, MAE = sum(w[keep] * abs(pred - 
        y[keep]))/n, distinct = sum(keep))
}

x <- c$speed_m_s

y <- c$cobb_degree

n <- length(y)

p0 <- cv_weighted(x, y, rep(1, n))

B <- 20000L

set.seed(261028)

boots <- t(replicate(B, cv_weighted(x, y, tabulate(sample.int(n, n, replace = TRUE), nbins = n))))

qs <- apply(boots, 2, quantile, c(0.025, 0.5, 0.975), na.rm = TRUE)

write(data.frame(metric = names(p0), estimate = as.numeric(p0), lower95 = qs[1, ], bootstrap_median = qs[2, 
    ], upper95 = qs[3, ], resamples = B, nonfinite = sum(!is.finite(boots[, "RMSE"])), seed = 261028), 
    "cobb_cv_uncertainty")

# Regression diagnostics and HC3 coefficient-sensitivity calculations.
fc <- lm(cobb_degree ~ speed_m_s, data = c)

fs <- lm(speed_m_s ~ mutant * nacet, data = s0)

fits <- list(Cobb = fc, Velocity = fs)

diagrows <- list()

robust <- list()

resids <- list()

pdf(file.path(output_dir, "model_diagnostics.pdf"), width = 7.2, height = 4.6, family = "Helvetica", 
    useDingbats = FALSE)

par(mfrow = c(2, 3), mar = c(3.3, 3.6, 2.2, 0.6), mgp = c(2, 0.6, 0), cex = 0.7)

for (nm in names(fits)) {
    f <- fits[[nm]]
    nn <- nobs(f)
    e <- residuals(f)
    X <- model.matrix(f)
    inv <- solve(crossprod(X))
    h <- hatvalues(f)
    vc <- inv %*% crossprod(X, X * (e/(1 - h))^2) %*% inv
    se <- sqrt(diag(vc))
    est <- coef(f)
    df <- df.residual(f)
    tt <- est/se
    pp <- 2 * pt(abs(tt), df, lower.tail = FALSE)
    ci <- qt(0.975, df) * se
    robust[[nm]] <- data.frame(model = nm, term = names(est), estimate = as.numeric(est), 
        HC3_SE = se, HC3_t = tt, HC3_p = pp, HC3_lower95 = est - ci, HC3_upper95 = est + 
            ci)
    w <- e^2 - mean(e^2)
    aux <- lm.fit(X, w)
    bp <- nn * sum(aux$fitted.values^2)/sum(w^2)
    bpp <- pchisq(bp, ncol(X) - 1, lower.tail = FALSE)
    sw <- shapiro.test(e)
    cd <- cooks.distance(f)
    diagrows[[nm]] <- data.frame(model = nm, n = nn, parameter_n = ncol(X), residual_df = df, 
        Shapiro_W = unname(sw$statistic), Shapiro_p = sw$p.value, Breusch_Pagan_statistic = bp, 
        Breusch_Pagan_p = bpp, max_Cook_distance = max(cd), Cook_reference = 4/nn, above_Cook_reference = sum(cd > 
            4/nn), maximum_leverage = max(h))
    ids <- if (nm == "Cobb") 
        c$animal_id
    else s0$animal_id
    resids[[nm]] <- data.frame(model = nm, animal_id = ids, fitted = as.numeric(fitted(f)), 
        residual = as.numeric(e), standardized_residual = as.numeric(rstandard(f)), Cook_distance = as.numeric(cd), 
        leverage = as.numeric(h), index = seq_len(nn))
    plot(fitted(f), e, pch = 16, cex = 0.65, col = "#516C7C", xlab = "Fitted value", ylab = "Residual", 
        main = paste(nm, "residuals"))
    abline(h = 0, lty = 2, col = "grey60")
    qqnorm(rstandard(f), pch = 16, cex = 0.65, col = "#516C7C", main = paste("Normal Q-Q; p =", 
        signif(sw$p.value, 3)))
    qqline(rstandard(f), col = "#B57950")
    plot(cd, type = "h", col = "#516C7C", xlab = "Source animal index", ylab = "Cook's distance", 
        main = paste(nm, "influence"))
    abline(h = 4/nn, lty = 2, col = "#B57950")
}

dev.off()

write(do.call(rbind, diagrows), "model_diagnostics")

write(do.call(rbind, robust), "model_HC3_sensitivity")

write(do.call(rbind, resids), "model_residuals")

write(do.call(rbind, lapply(seq_len(n), function(i) {
    f <- lm(cobb_degree ~ speed_m_s, data = c[-i, ])
    data.frame(removed_animal = c$animal_id[i], n_training = n - 1, slope = unname(coef(f)[2]), 
        R2 = summary(f)$r.squared)
})), "cobb_delete_one_sensitivity")

