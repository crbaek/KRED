rm(list = ls())

############################################################
## KRED factor pipeline: transformation + basic factor analysis 
## (Section 3 of the paper)
##
## Workflow
## 1. Read raw KRED vintage file (with group-code and tcode rows)
## 2. Apply FRED-MD-style transformations and save transformed csv
## 3. Run factor-analysis, mR2, and diffusion-index routines
##
## Required external file:
##   - KRED-library.R  (must define at least ICselect() and getFactor())
############################################################

############################################################
## User settings
############################################################
#workdir <- "../KRED/github/2025"
setwd(workdir)

raw_file <- "kred-Dec2025.csv"
library_file <- "KRED-library.R"

## If NULL, transformed_file is built automatically from raw_file.
transformed_file <- NULL

output_tag <- "2025"
sample_start <- as.Date("2009-06-01")
sample_end   <- as.Date("2025-12-01")

n_factor_plot  <- 4
n_factor_total <- 20
n_factor_mr2   <- 16
n_factor_ic    <- 10
n_factor_gic   <- 5
n_top_table    <- 10

## Anchor variables used only to fix factor signs for interpretation.
anchor_f1 <- "T3YFFM"     # financial conditions
anchor_f2 <- "IPMANSICS"  # real activity
anchor_f3 <- "PERMIT"     # housing
anchor_f4 <- "CE16OV"     # labor / price pressure

## Start recursive factor re-estimation after this many observations.
recursive_start_n <- 24

############################################################
## Packages and library functions
############################################################
library(imputeTS)
library(ggplot2)
library(patchwork)

if (!file.exists(library_file)) {
  stop("Required library file not found: ", library_file)
}
source(library_file)

if (!exists("ICselect")) {
  stop("ICselect() was not found after sourcing ", library_file)
}
if (!exists("getFactor")) {
  stop("getFactor() was not found after sourcing ", library_file)
}

############################################################
## Helper functions
############################################################
build_transformed_name <- function(raw_name) {
  base <- sub("\\.[^.]*$", "", raw_name)
  paste0(base, "-transformed.csv")
}

parse_date_column <- function(x) {
  x <- trimws(as.character(x))
  out <- as.Date(x, format = "%Y.%m.%d")
  if (all(is.na(out))) {
    out <- as.Date(x)
  }
  out
}

apply_transformation <- function(value, tcode) {
  eps <- 1e-05
  value <- as.numeric(value)

  if (length(value) == 0) {
    return(numeric(0))
  }

  if (tcode == 1) {
    fx <- value
  } else if (tcode == 2) {
    fx <- diff(value)
  } else if (tcode == 3) {
    fx <- diff(diff(value))
  } else if (tcode == 4) {
    fx <- log(value + eps)
  } else if (tcode == 5) {
    fx <- diff(log(value + eps))
  } else if (tcode == 6) {
    fx <- diff(diff(log(value + eps)))
  } else if (tcode == 7) {
    tt <- length(value)
    growth <- value[-1] / value[-tt] - 1
    fx <- diff(growth)
  } else {
    stop("Unsupported transformation code: ", tcode)
  }

  as.numeric(fx)
}

fix_factor_signs <- function(factor_mat, data_mat, anchors) {
  if (nrow(factor_mat) < length(anchors)) {
    stop("factor_mat has fewer rows than the number of anchors.")
  }

  for (i in seq_along(anchors)) {
    nm <- anchors[i]
    if (nm %in% colnames(data_mat)) {
      this_cor <- suppressWarnings(cor(factor_mat[i, ], data_mat[, nm], use = "pairwise.complete.obs"))
      if (!is.na(this_cor) && this_cor < 0) {
        factor_mat[i, ] <- -factor_mat[i, ]
      }
    }
  }

  factor_mat
}

############################################################
## 1. Read raw file and extract metadata
############################################################
if (!file.exists(raw_file)) {
  stop("Raw data file not found: ", raw_file)
}

if (is.null(transformed_file)) {
  transformed_file <- build_transformed_name(raw_file)
}

raw_df <- read.csv(raw_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

if (nrow(raw_df) < 4) {
  stop("raw_file must contain at least 4 rows: header row, group-code row, tcode row, and data rows.")
}

raw_dates <- raw_df[-(1:3), 1]
ndata <- length(raw_dates)

raw_values <- raw_df[-(1:3), -1, drop = FALSE]
series_names <- colnames(raw_values)

for (j in seq_len(ncol(raw_values))) {
  raw_values[[j]] <- as.numeric(raw_values[[j]])
}

gcode_full <- as.numeric(raw_df[2, -1])
names(gcode_full) <- colnames(raw_df)[-1]

tcode_full <- as.numeric(raw_df[3, -1])
names(tcode_full) <- colnames(raw_df)[-1]

if (any(is.na(tcode_full))) {
  bad_names <- names(tcode_full)[is.na(tcode_full)]
  stop("Some transformation codes are missing or nonnumeric: ", paste(bad_names, collapse = ", "))
}

############################################################
## 2. Transform each series and save transformed csv
############################################################
cat("\nStarting transformation from raw file: ", raw_file, "\n", sep = "")

transformed_mat <- matrix(NA_real_, nrow = ndata, ncol = ncol(raw_values))
colnames(transformed_mat) <- series_names

for (j in seq_len(ncol(raw_values))) {
  data_j <- as.numeric(raw_values[[j]])

  first_non_na <- match(FALSE, is.na(data_j))
  if (is.na(first_non_na)) {
    next
  }

  data_trimmed <- data_j[first_non_na:length(data_j)]

  if (anyNA(data_trimmed)) {
    data_filled <- imputeTS::na_interpolation(data_trimmed)
  } else {
    data_filled <- data_trimmed
  }

  transformed_j <- apply_transformation(data_filled, tcode = tcode_full[j])
  transformed_mat[, j] <- c(rep(NA_real_, ndata - length(transformed_j)), transformed_j)

  cat("Transformed series ", j, "/", ncol(raw_values), ": ", series_names[j], "\n", sep = "")
}

transformed_df <- data.frame(date = raw_dates, transformed_mat, check.names = FALSE)
write.csv(transformed_df, transformed_file, row.names = FALSE)
cat("Saved transformed data to: ", transformed_file, "\n", sep = "")

############################################################
## 3. Read transformed data and select sample window
############################################################
dd <- transformed_df

for (j in 2:ncol(dd)) {
  dd[[j]] <- as.numeric(dd[[j]])
}

date_index <- parse_date_column(dd[[1]])
if (all(is.na(date_index))) {
  stop("Date column in transformed data could not be parsed as Date.")
}

sample_id <- which(date_index >= sample_start & date_index <= sample_end)
if (length(sample_id) == 0) {
  stop("No observations found in the requested sample window.")
}

dat <- dd[sample_id, , drop = FALSE]
date_used <- date_index[sample_id]

############################################################
## 4. Drop variables with no value at the first sample date
############################################################
na_first <- which(is.na(dat[1, ]))
if (length(na_first) > 0) {
  cat("\nDropped because first sample observation is NA:\n")
  print(colnames(dat)[na_first])
}

dat2 <- dat[, -na_first, drop = FALSE]
dat2 <- dat2[, -1, drop = FALSE]  # remove date column

############################################################
## 5. Match group codes to retained columns
############################################################
gcode <- gcode_full[colnames(dat2)]
if (any(is.na(gcode))) {
  stop("Group codes could not be matched for: ",
       paste(colnames(dat2)[is.na(gcode)], collapse = ", "))
}

############################################################
## 6. Impute missing values inside the selected sample
############################################################
dat3 <- dat2
for (j in seq_len(ncol(dat3))) {
  dat3[[j]] <- imputeTS::na_interpolation(as.numeric(dat2[[j]]))
}

save.image(file = paste0("KFRED-", output_tag, "-factor.Rdata"))

############################################################
## 7. Scale series and run factor analysis
############################################################
dat4 <- scale(dat3)
vv <- svd(cov(dat4))
v <- vv$d
lg <- length(unique(gcode)) - 1

cat("\nBy-group information criteria\n")
for (g in seq_len(lg)) {
  idg <- which(gcode == g)
  if (length(idg) > 0) {
    ic_g <- ICselect(Y = t(dat4[, idg, drop = FALSE]), maxr = n_factor_gic)
    cat("Group ", g, ": ", ic_g$order, "\n", sep = "")
  }
}

cat("\nOverall variance share\n")
print(v[1] / sum(v))
print(cumsum(v) / sum(v))

cat("\nOverall information criteria\n")
print(ICselect(Y = t(dat4), maxr = n_factor_ic))

out <- getFactor(Y = t(dat4), r = n_factor_total, isCorr = FALSE)
out <- fix_factor_signs(
  factor_mat = out,
  data_mat = dat4,
  anchors = c(anchor_f1, anchor_f2, anchor_f3, anchor_f4)
)

############################################################
## 8. Scree plot and cumulative variance plot
############################################################
scale_factor <- 12
cum_var <- 100 * cumsum(v) / sum(v)
cum_rescaled <- cum_var * scale_factor / 100

pdf(paste0("scree_plot_", output_tag, ".pdf"), width = 12, height = 8)
par(mar = c(5, 4, 4, 5) + 0.1)
plot(1:length(v), v,
     type = "h", col = "blue", pch = 16, cex = 2, lwd = 2,
     ylim = c(0, scale_factor), cex.lab = 1.5, cex.main = 1.5,
     xlab = "Number of PC", ylab = "Eigenvalue",
     main = "PCA Scree Plot with Cumulative Variance")
par(new = TRUE)
plot(1:length(v), cum_rescaled,
     type = "b", col = "red", pch = 16, cex = 1.2, lwd = 2,
     xlab = "", ylab = "", axes = FALSE, ylim = c(0, scale_factor))
axis(4, at = seq(0, scale_factor, length.out = 6), labels = seq(0, 100, 20))
mtext("Cumulative Variance (%)", side = 4, line = 3, cex = 1.5)
legend("bottomright",
       legend = c("eigenvalue", "cumulative variance"),
       col = c("blue", "red"),
       pch = c(16, 16),
       lty = c(1, 1),
       lwd = c(2, 2),
       bg = "white",
       box.col = "black",
       cex = 1.5)
dev.off()

############################################################
## 9. Plot the first four factors
############################################################
factor_plot_list <- vector("list", n_factor_plot)

for (i in seq_len(n_factor_plot)) {
  df_i <- data.frame(Time = date_used, Value = out[i, ])

  factor_plot_list[[i]] <- ggplot(df_i, aes(x = Time, y = Value)) +
    geom_line(color = "steelblue") +
    labs(
      title = bquote("Factor" ~ F[.(i)]),
      x = "Time",
      y = bquote(F[.(i)])
    ) +
    scale_x_date(date_breaks = "1 year", date_labels = "\n%Y") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.title.x = element_text(size = 14),
      axis.text = element_text(size = 12)
    )
}

pdf(paste0("factor_plots_", output_tag, ".pdf"), width = 12, height = 10)
print(factor_plot_list[[1]] / factor_plot_list[[2]] / factor_plot_list[[3]] / factor_plot_list[[4]])
dev.off()

############################################################
## 10. Plain diffusion index plots
############################################################
diffusion_plot_list <- vector("list", n_factor_plot)

for (i in seq_len(n_factor_plot)) {
  df_i <- data.frame(Time = date_used, Value = cumsum(out[i, ]))

  diffusion_plot_list[[i]] <- ggplot(df_i, aes(x = Time, y = Value)) +
    geom_line(color = "steelblue") +
    labs(
      title = bquote("Diffusion Index " ~ F[.(i)]),
      x = "Time",
      y = bquote(F[.(i)])
    ) +
    scale_x_date(date_breaks = "1 year", date_labels = "\n%Y") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.title.x = element_text(size = 14),
      axis.text = element_text(size = 12)
    )
}

pdf(paste0("diffusion_plots_", output_tag, ".pdf"), width = 12, height = 10)
print(diffusion_plot_list[[1]] / diffusion_plot_list[[2]] / diffusion_plot_list[[3]] / diffusion_plot_list[[4]])
dev.off()

############################################################
## 11. Calculate mR_i^2 using the first n_factor_mr2 factors
############################################################
n_series <- ncol(dat4)
S <- matrix(0, n_factor_mr2, n_series)
colnames(S) <- colnames(dat4)

for (k in seq_len(n_series)) {
  for (j in seq_len(n_factor_mr2)) {
    if (j == 1) {
      fit_j <- lm(dat4[, k] ~ 1 + out[1, ])
    } else {
      fit_j <- lm(dat4[, k] ~ ., data = as.data.frame(t(out[1:j, ])))
    }
    S[j, k] <- summary(fit_j)$r.squared
  }
}

Smean <- rowMeans(S)
print(diff(Smean))

s1 <- S[1, ]
s2 <- S[2, ]
s3 <- S[3, ]
s4 <- S[4, ]

group_labels <- as.factor(as.numeric(gcode))

df1 <- data.frame(index = seq_along(s1), value = s1, group = group_labels)
df2 <- data.frame(index = seq_along(s2), value = s2 - s1, group = group_labels)
df3 <- data.frame(index = seq_along(s3), value = s3 - s2, group = group_labels)
df4 <- data.frame(index = seq_along(s4), value = s4 - s3, group = group_labels)

y1 <- expression(italic(m) * R[i]^2 * (1))
y2 <- expression(italic(m) * R[i]^2 * (2) - italic(m) * R[i]^2 * (1))
y3 <- expression(italic(m) * R[i]^2 * (3) - italic(m) * R[i]^2 * (2))
y4 <- expression(italic(m) * R[i]^2 * (4) - italic(m) * R[i]^2 * (3))
y5 <- expression(R[i]^2 * (4))

make_mr2_plot <- function(df, title, ylab_expr) {
  ggplot(df, aes(x = index, y = value)) +
    geom_segment(aes(xend = index, yend = 0), color = "blue", linewidth = 1) +
    scale_x_continuous(
      breaks = seq(1, length(df$index), by = 5),
      labels = df$group[seq(1, length(df$index), by = 5)]
    ) +
    theme_minimal() +
    labs(title = title, x = "Group", y = ylab_expr) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 1, size = 12),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 16, face = "bold")
    )
}

p1 <- make_mr2_plot(df1, "First Factor", y1)
p2 <- make_mr2_plot(df2, "Second Factor", y2)
p3 <- make_mr2_plot(df3, "Third Factor", y3)
p4 <- make_mr2_plot(df4, "Fourth Factor", y4)

pdf(paste0("mR2_plots_", output_tag, ".pdf"), width = 12, height = 8)
print((p1 | p2) / (p3 | p4))
dev.off()

############################################################
## 12. Overall R^2 barplot and top-table export
############################################################
id_sorted <- order(s4, decreasing = TRUE)

pdf(paste0("barplot_output_", output_tag, ".pdf"), width = 8, height = 6)
par(mar = c(7, 5, 4, 2))
barplot(s4[id_sorted],
        names.arg = colnames(dat3)[id_sorted],
        las = 2,
        col = "blue",
        ylab = y5,
        main = "Sorted R²(4) Values",
        border = "blue",
        cex.names = 0.8,
        cex.lab = 1.2,
        cex.main = 1.3,
        ylim = c(0, 1),
        width = 0.2,
        space = 0.9)
dev.off()

t1 <- order(s1, decreasing = TRUE)[1:n_top_table]
t2 <- order(s2 - s1, decreasing = TRUE)[1:n_top_table]
t3 <- order(s3 - s2, decreasing = TRUE)[1:n_top_table]
t4 <- order(s4 - s3, decreasing = TRUE)[1:n_top_table]

mr2_table <- data.frame(
  factor1_name  = colnames(dat3)[t1],
  factor1_mR2   = s1[t1],
  factor1_group = gcode[t1],
  factor2_name  = colnames(dat3)[t2],
  factor2_mR2   = (s2 - s1)[t2],
  factor2_group = gcode[t2],
  factor3_name  = colnames(dat3)[t3],
  factor3_mR2   = (s3 - s2)[t3],
  factor3_group = gcode[t3],
  factor4_name  = colnames(dat3)[t4],
  factor4_mR2   = (s4 - s3)[t4],
  factor4_group = gcode[t4]
)

write.csv(mr2_table, paste0("mR2_table_", output_tag, ".csv"), row.names = FALSE)

############################################################
## 13. Recursively demeaned data for monitoring-style diffusion index
############################################################
xtilde <- matrix(NA_real_, nrow = nrow(dat3), ncol = ncol(dat3))
colnames(xtilde) <- colnames(dat3)

for (j in seq_len(ncol(dat3))) {
  xj <- as.numeric(dat3[[j]])
  xbar_cum <- cumsum(xj) / seq_along(xj)
  sj <- sd(xj, na.rm = TRUE)
  if (is.na(sj) || sj == 0) sj <- 1
  xtilde[, j] <- (xj - xbar_cum) / sj
}

recursive_factor <- matrix(NA_real_, nrow = n_factor_plot, ncol = nrow(xtilde))

for (tt in recursive_start_n:nrow(xtilde)) {
  xtmp <- xtilde[1:tt, , drop = FALSE]
  ftmp <- getFactor(Y = t(xtmp), r = n_factor_plot, isCorr = FALSE)
  ftmp <- fix_factor_signs(
    factor_mat = ftmp,
    data_mat = xtmp,
    anchors = c(anchor_f1, anchor_f2, anchor_f3, anchor_f4)
  )
  recursive_factor[, tt] <- ftmp[, ncol(ftmp)]
}

recursive_diffusion <- matrix(NA_real_, nrow = n_factor_plot, ncol = nrow(xtilde))
for (i in seq_len(n_factor_plot)) {
  recursive_diffusion[i, recursive_start_n:nrow(xtilde)] <-
    cumsum(recursive_factor[i, recursive_start_n:nrow(xtilde)])
}

cat("\nRecursive monitoring information criteria on recursively demeaned data\n")
print(ICselect(Y = t(xtilde), maxr = 20))

############################################################
## 14. Monitoring-style diffusion index plots
############################################################
diffusion_tilde_list <- vector("list", n_factor_plot)

for (i in seq_len(n_factor_plot)) {
  df_i <- data.frame(Time = date_used, Value = recursive_diffusion[i, ])

  diffusion_tilde_list[[i]] <- ggplot(df_i, aes(x = Time, y = Value)) +
    geom_line(color = "steelblue") +
    labs(
      title = bquote("Diffusion Index " ~ F[.(i)]),
      x = "Time",
      y = bquote(F[.(i)])
    ) +
    scale_x_date(date_breaks = "1 year", date_labels = "\n%Y") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.title.x = element_text(size = 14),
      axis.text = element_text(size = 12)
    )
}

pdf(paste0("diffusiontilde_plots_", output_tag, ".pdf"), width = 12, height = 10)
print(diffusion_tilde_list[[1]] / diffusion_tilde_list[[2]] / diffusion_tilde_list[[3]] / diffusion_tilde_list[[4]])
dev.off()

cat("\nDone: transformation and factor-analysis outputs saved for ", output_tag, "\n", sep = "")
