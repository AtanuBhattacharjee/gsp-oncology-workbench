# =============================================================================
# GSP in Oncology – Data workbench
# Good Statistical Practice in Oncology · HBCH & RC Muzaffarpur · 24 October 2026
# Upload a trial dataset and reproduce every analysis taught in the sessions:
# data checks, Table 1, Kaplan–Meier, follow-up, log-rank, Cox, proportional
# hazards, RMST, landmark differences, subgroups, sample size, diagnostic
# accuracy, ROC, ICC, Bland–Altman and kappa.
# Packages: shiny, bslib, DT, readxl, survival, ggplot2
#   install.packages(c("shiny", "bslib", "DT", "readxl", "survival", "ggplot2"))
# Run locally:   shiny::runApp("app.R")   (or open in RStudio and click Run App)
# Deploy:        copy app.R to a Shiny Server / Posit Connect / shinyapps.io folder.
# Data format:   one row per patient; a time column in months from randomisation,
#                an event indicator, a treatment arm, and any other columns
#                (baseline factors, test results, reader measurements).
# Settings:      time can be recorded in days, weeks, months or years and shown in
#                any of these; group colours come from a catalogue of palettes
#                or a custom colour picker (Data tab, Display settings).
#                The "Use the example dataset" button loads a simulated trial
#                with every column the tabs expect.
# =============================================================================

library(shiny)
library(bslib)
library(DT)
library(readxl)
library(survival)
library(ggplot2)

# ---- palette and theme -------------------------------------------------------
PL <- "#23527C"; AM <- "#E07A3F"; AMD <- "#A04F1E"; BL <- "#2E8B74"
GR <- "#8C96A3"; TX <- "#1F2933"; BG <- "#F7F9FB"; TP <- "#E6EEF6"
ARM_COLS <- c(AM, PL, BL, "#7A6FB0", "#B35C5C", "#8C96A3")
UNIT_DAYS <- c(Days = 1, Weeks = 7, Months = 30.4375, Years = 365.25)
PALETTES <- list(
  "Workbench (blue and orange)" = ARM_COLS,
  "Colour-blind safe (Okabe–Ito)" = c("#E69F00", "#0072B2", "#009E73", "#CC79A7", "#56B4E9", "#D55E00"),
  "Viridis" = c("#440154", "#2A788E", "#7AD151", "#414487", "#22A884", "#FDE725"),
  "Dark 2" = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02"),
  "Set 2 (soft)" = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#E5C494"),
  "Paired (strong)" = c("#1F78B4", "#E31A1C", "#33A02C", "#FF7F00", "#6A3D9A", "#B15928"),
  "Greyscale" = c("#1A1A1A", "#7A7A7A", "#B0B0B0", "#474747", "#999999", "#D0D0D0"))
swatches <- function(cols) div(style = "display:flex;gap:6px;margin:4px 0 10px 0", lapply(cols, function(c) span(style = sprintf("display:inline-block;width:26px;height:26px;border-radius:6px;background:%s;border:1px solid #D6DFE8", c), title = c)))

theme_gsp <- function(base = 15) {
  theme_minimal(base_size = base) +
    theme(plot.background = element_rect(fill = BG, colour = NA),
          panel.background = element_rect(fill = BG, colour = NA),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "#E3E8EE"),
          axis.text = element_text(colour = TX), axis.title = element_text(colour = TX),
          legend.position = "top", legend.title = element_blank(),
          plot.title = element_text(colour = PL, face = "bold"))
}

fmt <- function(x, d = 2) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = d))
fmtp <- function(p) ifelse(is.na(p), "NA", ifelse(p < 0.001, "< 0.001", formatC(p, format = "f", digits = 3)))

# ---- example dataset ---------------------------------------------------------
make_example <- function(n = 350, seed = 2410) {
  set.seed(seed)
  arm <- rep(c("Gefitinib", "Combination"), each = n / 2)
  x <- as.integer(arm == "Combination")
  ps2 <- rbinom(n, 1, 0.2)
  brain <- rbinom(n, 1, 0.18)
  mutation <- ifelse(runif(n) < 0.6, "Exon 19 deletion", "L858R")
  rate <- log(2) / 8 * exp(-0.69 * x + 0.45 * ps2 + 0.3 * brain)
  t_event <- rexp(n, rate)
  entry <- runif(n, 0, 30)
  cens <- pmin(42 - entry, rexp(n, 1 / 200))
  pfs_months <- pmax(0.1, round(pmin(t_event, cens), 1))
  pfs_event <- as.integer(t_event <= cens)
  # diagnostic study: tissue reference, plasma test and a continuous level
  tissue <- rbinom(n, 1, 0.30)
  plasma <- ifelse(tissue == 1, rbinom(n, 1, 0.75), rbinom(n, 1, 0.02))
  level <- round(ifelse(tissue == 1, rlnorm(n, log(12), 0.7), rlnorm(n, log(5), 0.7)), 1)
  # agreement study: two readers
  true_mm <- runif(n, 15, 85)
  reader_a <- round(true_mm + rnorm(n, 0, 2.2), 1)
  reader_b <- round(true_mm + 3 + rnorm(n, 0, 2.2), 1)
  cats <- c("CR", "PR", "SD", "PD")
  resp_a <- sample(cats, n, replace = TRUE, prob = c(0.1, 0.35, 0.35, 0.2))
  resp_b <- ifelse(runif(n) < 0.83, resp_a, sample(cats, n, replace = TRUE))
  data.frame(id = sprintf("P%03d", 1:n), arm, pfs_months, pfs_event,
             age = round(rnorm(n, 56, 10)), sex = sample(c("Female", "Male"), n, TRUE, c(0.48, 0.52)),
             ps = ifelse(ps2 == 1, "PS 2", "PS 0-1"), brain_mets = ifelse(brain == 1, "Yes", "No"),
             mutation, tissue_egfr = ifelse(tissue == 1, "Positive", "Negative"),
             plasma_egfr = ifelse(plasma == 1, "Positive", "Negative"), ctdna_level = level,
             reader_a_mm = reader_a, reader_b_mm = reader_b, recist_reader_a = resp_a,
             recist_reader_b = resp_b, stringsAsFactors = FALSE)
}

# ---- survival helpers ----------------------------------------------------------
km_steps <- function(fit) {
  s <- summary(fit, censored = TRUE)
  strata <- if (is.null(s$strata)) rep("All patients", length(s$time)) else sub(".*=", "", as.character(s$strata))
  d <- data.frame(time = s$time, surv = s$surv, lower = s$lower, upper = s$upper,
                  n.censor = s$n.censor, strata = strata)
  lev <- unique(strata)
  out <- lapply(lev, function(l) {
    di <- d[d$strata == l, ]
    rbind(data.frame(time = 0, surv = 1, lower = 1, upper = 1, n.censor = 0, strata = l), di)
  })
  do.call(rbind, out)
}

plot_km <- function(fit, xmax = NULL, show_ci = TRUE, landmark = NA, ylab = "Survival probability", pal = ARM_COLS, xlab = "Months from randomisation") {
  d <- km_steps(fit)
  if (is.null(xmax) || is.na(xmax)) xmax <- max(d$time)
  lev <- unique(d$strata); cols <- setNames(rep(pal, length.out = length(lev)), lev)
  d$strata <- factor(d$strata, levels = lev)
  brks <- pretty(c(0, xmax), 6); brks <- brks[brks <= xmax]
  risk <- summary(fit, times = brks, extend = TRUE)
  rs <- if (is.null(risk$strata)) rep(lev[1], length(risk$time)) else sub(".*=", "", as.character(risk$strata))
  rt <- data.frame(time = risk$time, n = risk$n.risk, strata = factor(rs, levels = lev))
  rowy <- -0.24 - 0.11 * (seq_along(lev) - 1)
  rt$y <- rowy[as.integer(rt$strata)]
  padl <- 0.06 + 0.021 * max(nchar(lev))
  lab <- data.frame(x = -0.065 * xmax, y = rowy, strata = factor(lev, levels = lev), label = lev)
  cens <- d[d$n.censor > 0, ]
  dci <- transform(d, lower = ifelse(is.na(lower), surv, lower), upper = ifelse(is.na(upper), surv, upper))
  g <- ggplot(d, aes(time, surv, colour = strata))
  if (show_ci) g <- g + geom_ribbon(data = dci, aes(ymin = lower, ymax = upper, fill = strata), alpha = 0.12, colour = NA)
  g <- g + geom_step(linewidth = 1.2) +
    geom_point(data = cens, shape = 3, size = 2.2, stroke = 1) +
    geom_hline(yintercept = 0.5, linetype = "dotted", colour = GR) +
    geom_hline(yintercept = -0.1, colour = "#DCE3EA") +
    geom_text(data = rt, aes(x = time, y = y, label = n), size = 4.2, show.legend = FALSE) +
    geom_text(data = lab, aes(x = x, y = y, label = label), hjust = 1, size = 4.2, fontface = "bold", show.legend = FALSE) +
    annotate("text", x = 0, y = -0.14, label = "Number at risk", hjust = 0, size = 4, colour = TX) +
    scale_colour_manual(values = cols) + scale_fill_manual(values = cols, guide = "none") +
    scale_x_continuous(breaks = brks, limits = c(-padl * xmax, xmax)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(min(rowy) - 0.05, 1)) +
    guides(colour = guide_legend(nrow = if (length(lev) > 2) 2 else 1)) +
    labs(x = xlab, y = ylab) + theme_gsp()
  if (!is.na(landmark)) g <- g + geom_vline(xintercept = landmark, linetype = "dashed", colour = BL)
  g
}

km_median_table <- function(fit, dg = 1) {
  tb <- summary(fit)$table
  if (is.null(dim(tb))) tb <- t(as.matrix(tb))
  nm <- if (nrow(tb) == 1 && is.null(rownames(tb))) "All patients" else sub(".*=", "", rownames(tb))
  data.frame(Group = nm, N = tb[, "records"], Events = tb[, "events"],
             `Median (95% CI)` = paste0(fmt(tb[, "median"], dg), " (", fmt(tb[, "0.95LCL"], dg), "–",
                                        ifelse(is.na(tb[, "0.95UCL"]), "NR", fmt(tb[, "0.95UCL"], dg)), ")"),
             check.names = FALSE, row.names = NULL)
}

landmark_table <- function(fit, t) {
  s <- summary(fit, times = t, extend = TRUE)
  g <- if (is.null(s$strata)) "All patients" else sub(".*=", "", as.character(s$strata))
  data.frame(Group = g, S = s$surv, SE = s$std.err, lower = s$lower, upper = s$upper)
}

rmst_fun <- function(time, status, tau) {
  o <- order(time); time <- time[o]; status <- status[o]
  ev <- sort(unique(time[status == 1 & time <= tau]))
  S <- 1; tt <- 0; ss <- 1; terms <- NULL
  for (ti in ev) {
    n <- sum(time >= ti); d <- sum(time == ti & status == 1)
    S <- S * (1 - d / n); tt <- c(tt, ti); ss <- c(ss, S); terms <- rbind(terms, c(n, d))
  }
  tt <- c(tt, tau)
  widths <- diff(tt); area <- sum(widths * ss)
  v <- 0
  if (!is.null(terms)) for (j in seq_len(nrow(terms))) {
    A <- sum(widths[(j + 1):length(widths)] * ss[(j + 1):length(ss)])
    n <- terms[j, 1]; d <- terms[j, 2]
    if (n > d) v <- v + A^2 * d / (n * (n - d))
  }
  c(rmst = area, se = sqrt(v))
}

# ---- sample size helpers -------------------------------------------------------
p_event <- function(median, accrual, follow) {
  l <- log(2) / median
  1 - (exp(-l * follow) - exp(-l * (accrual + follow))) / (l * accrual)
}

simon_design <- function(p0, p1, alpha, beta, nmax = 60) {
  best_opt <- NULL; best_mm <- NULL
  for (n in 5:nmax) {
    for (n1 in 1:(n - 1)) {
      n2 <- n - n1; x1 <- 0:n1; rr <- 0:n
      M0 <- outer(x1, rr, function(x, r) pbinom(r - x, n2, p0, lower.tail = FALSE))
      M1 <- outer(x1, rr, function(x, r) pbinom(r - x, n2, p1, lower.tail = FALSE))
      A0 <- apply(dbinom(x1, n1, p0) * M0, 2, function(v) rev(cumsum(rev(v))))
      A1 <- apply(dbinom(x1, n1, p1) * M1, 2, function(v) rev(cumsum(rev(v))))
      for (r1 in 0:(n1 - 1)) {
        a <- A0[r1 + 2, ]; ok <- which(a <= alpha & rr >= r1)
        if (!length(ok)) next
        j <- ok[1]
        if (A1[r1 + 2, j] >= 1 - beta) {
          pet <- pbinom(r1, n1, p0); en0 <- n1 + (1 - pet) * n2
          d <- c(r1 = r1, n1 = n1, r = rr[j], n = n, EN0 = en0, PET0 = pet, alpha = a[j], power = A1[r1 + 2, j])
          if (is.null(best_opt) || en0 < best_opt["EN0"] - 1e-9) best_opt <- d
          if (is.null(best_mm) || n < best_mm["n"] || (n == best_mm["n"] && en0 < best_mm["EN0"])) best_mm <- d
        }
      }
    }
    if (!is.null(best_mm) && n > best_mm["n"] + 15) break
  }
  list(optimal = best_opt, minimax = best_mm)
}

ahern <- function(p0, p1, alpha, beta, nmax = 200) {
  for (n in 5:nmax) for (r in 0:n) {
    if (pbinom(r - 1, n, p0, lower.tail = FALSE) <= alpha) {
      if (pbinom(r - 1, n, p1, lower.tail = FALSE) >= 1 - beta) return(c(n = n, r = r))
      break
    }
  }
  c(n = NA, r = NA)
}

# ---- diagnostic and agreement helpers --------------------------------------
wilson <- function(x, n, z = 1.96) {
  if (n == 0) return(c(NA, NA, NA))
  p <- x / n; den <- 1 + z^2 / n
  c(p, (p + z^2 / (2 * n) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / den,
    (p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / den)
}

roc_fun <- function(marker, truth) {
  ok <- !is.na(marker) & !is.na(truth); marker <- marker[ok]; truth <- truth[ok]
  thr <- sort(unique(marker), decreasing = TRUE)
  tpr <- sapply(thr, function(c) mean(marker[truth == 1] >= c))
  fpr <- sapply(thr, function(c) mean(marker[truth == 0] >= c))
  n1 <- sum(truth == 1); n0 <- sum(truth == 0)
  auc <- (sum(rank(marker)[truth == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  q1 <- auc / (2 - auc); q2 <- 2 * auc^2 / (1 + auc)
  se <- sqrt((auc * (1 - auc) + (n1 - 1) * (q1 - auc^2) + (n0 - 1) * (q2 - auc^2)) / (n1 * n0))
  j <- which.max(tpr - fpr)
  list(curve = data.frame(fpr = c(0, fpr, 1), tpr = c(0, tpr, 1)), auc = auc, se = se,
       cut = thr[j], sens = tpr[j], spec = 1 - fpr[j], n1 = n1, n0 = n0)
}

icc_fun <- function(M, conf = 0.95) {
  M <- as.matrix(M[complete.cases(M), ]); n <- nrow(M); k <- ncol(M); gm <- mean(M)
  MSR <- k * sum((rowMeans(M) - gm)^2) / (n - 1)
  MSC <- n * sum((colMeans(M) - gm)^2) / (k - 1)
  SSE <- sum((M - outer(rowMeans(M), rep(1, k)) - outer(rep(1, n), colMeans(M)) + gm)^2)
  MSE <- SSE / ((n - 1) * (k - 1)); a2 <- (1 - conf) / 2
  icc_c <- (MSR - MSE) / (MSR + (k - 1) * MSE)
  FL <- (MSR / MSE) / qf(1 - a2, n - 1, (n - 1) * (k - 1)); FU <- (MSR / MSE) * qf(1 - a2, (n - 1) * (k - 1), n - 1)
  ci_c <- c((FL - 1) / (FL + k - 1), (FU - 1) / (FU + k - 1))
  icc_a <- (MSR - MSE) / (MSR + (k - 1) * MSE + k * (MSC - MSE) / n)
  a <- k * icc_a / (n * (1 - icc_a)); b <- 1 + k * icc_a * (n - 1) / (n * (1 - icc_a))
  v <- (a * MSC + b * MSE)^2 / ((a * MSC)^2 / (k - 1) + (b * MSE)^2 / ((n - 1) * (k - 1)))
  Fs <- qf(1 - a2, n - 1, v); Fs2 <- qf(1 - a2, v, n - 1)
  ci_a <- c(n * (MSR - Fs * MSE) / (Fs * (k * MSC + (k * n - k - n) * MSE) + n * MSR),
            n * (Fs2 * MSR - MSE) / (k * MSC + (k * n - k - n) * MSE + n * Fs2 * MSR))
  list(n = n, k = k, MSR = MSR, MSC = MSC, MSE = MSE,
       table = data.frame(Form = c("ICC(2,1) absolute agreement", "ICC(3,1) consistency"),
                          ICC = c(icc_a, icc_c), Lower = c(ci_a[1], ci_c[1]), Upper = c(ci_a[2], ci_c[2])))
}

kappa_fun <- function(a, b) {
  ok <- !is.na(a) & !is.na(b); a <- a[ok]; b <- b[ok]
  lev <- sort(union(unique(a), unique(b)))
  tab <- table(factor(a, lev), factor(b, lev)); N <- sum(tab)
  po <- sum(diag(tab)) / N; pe <- sum(rowSums(tab) * colSums(tab)) / N^2
  k <- (po - pe) / (1 - pe); se <- sqrt(po * (1 - po) / (N * (1 - pe)^2))
  list(tab = tab, N = N, po = po, pe = pe, kappa = k, lower = k - 1.96 * se, upper = k + 1.96 * se)
}

read_any <- function(path, name) {
  ext <- tolower(tools::file_ext(name))
  if (ext %in% c("xlsx", "xls")) as.data.frame(read_excel(path))
  else if (ext == "tsv") read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  else read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

note <- function(...) div(class = "gsp-note", ...)
formula_card <- function(...) div(class = "gsp-formula", withMathJax(...))

# ---- UI -------------------------------------------------------------------------
gsp_theme <- bs_theme(version = 5, bg = BG, fg = TX, primary = PL, secondary = AM,
                      base_font = font_google("Source Sans 3", local = FALSE), heading_font = font_google("Source Sans 3", local = FALSE),
                      "navbar-bg" = "#FFFFFF")

css <- HTML(sprintf("
  .gsp-note{background:#E8F3EF;border-radius:12px;padding:12px 16px;margin:10px 0;color:%s}
  .gsp-formula{background:#FFFFFF;border:1px solid #D6DFE8;border-radius:12px;padding:12px 18px;margin:10px 0}
  .card-header{color:%s;font-weight:700;background:%s}
  h4,h5{color:%s;font-weight:700}
  .navbar-brand{color:%s !important;font-weight:700}
  .sentence{background:%s;border-radius:12px;padding:14px 18px;font-size:1.08rem;margin-bottom:10px}
  table.dataTable td, table.dataTable th{font-size:0.95rem}
", TX, PL, TP, PL, PL, TP))

ui <- page_navbar(
  title = "GSP in Oncology · data workbench", theme = gsp_theme, fillable = FALSE, header = tags$head(tags$style(css)),
  # 1 ---------------------------------------------------------------------------
  nav_panel("1 · Data",
    layout_sidebar(
      sidebar = sidebar(width = 330,
        h5("Load a dataset"),
        fileInput("file", "Upload CSV, TSV or Excel", accept = c(".csv", ".tsv", ".xlsx", ".xls")),
        actionButton("use_example", "Use the example dataset", class = "btn-primary"),
        downloadButton("dl_example", "Download the example CSV", class = "btn-outline-secondary mt-2"),
        hr(), h5("Map the variables"),
        selectInput("time", "Time from randomisation", NULL),
        selectInput("unit_in", "Time is recorded in", names(UNIT_DAYS), selected = "Months"),
        selectInput("status", "Event indicator", NULL),
        selectInput("event_value", "Value that means an event", NULL),
        selectInput("arm", "Treatment arm", NULL),
        selectInput("ref_arm", "Reference (control) arm", NULL),
        selectizeInput("t1vars", "Variables for Table 1", NULL, multiple = TRUE),
        hr(), h5("Display settings"),
        selectInput("unit_out", "Show results in", names(UNIT_DAYS), selected = "Months"),
        selectInput("palette", "Colour catalogue for groups", c(names(PALETTES), "Custom"), selected = names(PALETTES)[1]),
        uiOutput("pal_preview"),
        conditionalPanel("input.palette == 'Custom'",
          p("Click a square to pick any colour for groups 1–6 (reference group first).", style = "font-size:0.9rem"),
          div(style = "display:flex;gap:8px;flex-wrap:wrap",
              lapply(1:6, function(i) tags$input(type = "color", id = paste0("cc", i), value = ARM_COLS[i],
                style = "width:40px;height:34px;border:none;background:none;padding:0",
                oninput = sprintf("Shiny.setInputValue('ccol%d', this.value)", i)))))
      ),
      navset_card_underline(
        nav_panel("Data checks", tableOutput("checks"), uiOutput("check_notes")),
        nav_panel("Table 1", note("Describes the randomised arms: median (IQR) for numbers, n (%) with denominators, a row for missing values, no p-values."), tableOutput("table1")),
        nav_panel("Preview", DTOutput("preview"))
      )
    )
  ),
  # 2 ---------------------------------------------------------------------------
  nav_panel("2 · Survival",
    layout_sidebar(
      sidebar = sidebar(width = 300,
        numericInput("landmark", "Landmark time (months)", 12, min = 0),
        numericInput("xmax", "Maximum time shown (blank = all)", NA, min = 1),
        checkboxInput("show_ci", "Show 95% confidence bands", TRUE),
        downloadButton("dl_km", "Download the curve (PNG)", class = "btn-outline-secondary")
      ),
      layout_columns(col_widths = c(8, 4),
        card(card_header("Kaplan–Meier curves by arm"), plotOutput("km_plot", height = "600px")),
        card(card_header("Medians and landmark rates"), tableOutput("median_tab"), tableOutput("landmark_tab"),
             formula_card("$$\\hat S(t)=\\prod_{t_i\\le t}\\left(1-\\frac{d_i}{n_i}\\right)$$",
                          "$$\\widehat{\\mathrm{Var}}[\\hat S(t)]=\\hat S(t)^2\\sum_{t_i\\le t}\\frac{d_i}{n_i(n_i-d_i)}$$"))
      ),
      layout_columns(col_widths = c(6, 6),
        card(card_header("Median follow-up by reverse Kaplan–Meier"), plotOutput("fu_plot", height = "360px"), uiOutput("fu_text")),
        card(card_header("Kaplan–Meier by hand: first ten event times"), tableOutput("km_hand"))
      )
    )
  ),
  # 3 ---------------------------------------------------------------------------
  nav_panel("3 · Compare arms",
    layout_sidebar(
      sidebar = sidebar(width = 310,
        h5("Groups to compare"),
        selectInput("grp1", "Compare groups defined by", NULL),
        selectInput("grp2", "Split further by (optional)", NULL),
        selectInput("grp_ref", "Reference group", NULL),
        checkboxInput("cmp_ci", "Show 95% confidence bands", TRUE),
        downloadButton("dl_cmp_km", "Download the curve (PNG)", class = "btn-outline-secondary"),
        hr(),
        numericInput("tau", "RMST horizon τ (months)", 24, min = 1),
        numericInput("lm2", "Landmark for difference and NNT (months)", 12, min = 0),
        selectizeInput("adjust", "Adjust the Cox model for", NULL, multiple = TRUE),
        note("Every table in this tab follows the groups chosen above: treatment arm by default, or any baseline factor such as PS, sex or mutation type, alone or split by arm.")
      ),
      layout_columns(col_widths = c(8, 4),
        card(card_header(textOutput("cmp_title", inline = TRUE)), plotOutput("cmp_km", height = "600px")),
        card(card_header("Medians by group"), tableOutput("cmp_median"))
      ),
      layout_columns(col_widths = c(6, 6),
        card(card_header("Log-rank test"), tableOutput("lr_tab"), uiOutput("lr_text")),
        card(card_header("Cox model"), tableOutput("cox_tab"), uiOutput("cox_text"),
             formula_card("$$h(t\\mid x)=h_0(t)\\,e^{\\beta x},\\qquad \\mathrm{HR}=e^{\\beta},\\qquad 95\\%\\ \\mathrm{CI}=e^{\\beta\\pm1.96\\,\\mathrm{SE}}$$"))
      ),
      layout_columns(col_widths = c(6, 6),
        card(card_header("Proportional hazards check"), plotOutput("ll_plot", height = "340px"), tableOutput("ph_tab")),
        card(card_header("RMST and landmark difference"), tableOutput("rmst_tab"), tableOutput("lmdiff_tab"),
             formula_card("$$\\mathrm{RMST}(\\tau)=\\int_0^{\\tau}\\hat S(t)\\,dt \\qquad \\mathrm{NNT}(t)=\\frac{1}{\\hat S_1(t)-\\hat S_0(t)}$$"))
      )
    )
  ),
  # 4 ---------------------------------------------------------------------------
  nav_panel("4 · Subgroups",
    layout_sidebar(
      sidebar = sidebar(width = 300, selectizeInput("subvars", "Subgroup variables", NULL, multiple = TRUE),
                        note("Judge a subgroup by the interaction test, not by whether its own CI crosses 1.")),
      card(card_header("Forest plot of hazard ratios"), uiOutput("forest_ui"), tableOutput("forest_tab"))
    )
  ),
  # 5 ---------------------------------------------------------------------------
  nav_panel("5 · Sample size",
    navset_card_underline(
      nav_panel("Time-to-event",
        layout_columns(col_widths = c(4, 8),
          div(numericInput("ss_hr", "Target hazard ratio", 0.67, 0.1, 0.99, 0.01),
              numericInput("ss_alpha", "Two-sided α", 0.05, 0.001, 0.2, 0.005),
              numericInput("ss_power", "Power", 0.8, 0.5, 0.99, 0.05),
              numericInput("ss_k", "Allocation ratio (experimental : control = k : 1)", 1, 0.25, 4, 0.25),
              numericInput("ss_med", "Control median (months)", 10.2, 0.5),
              numericInput("ss_acc", "Accrual period (months)", 24, 1),
              numericInput("ss_fu", "Minimum follow-up (months)", 12, 0),
              numericInput("ss_drop", "Expected loss to follow-up", 0.10, 0, 0.5, 0.01)),
          div(uiOutput("ss_tte"), plotOutput("ss_curve", height = "340px"),
              formula_card("$$D=\\frac{(1+k)^2}{k}\\cdot\\frac{(z_{1-\\alpha/2}+z_{1-\\beta})^2}{(\\ln \\mathrm{HR})^2}$$")))),
      nav_panel("Non-inferiority",
        layout_columns(col_widths = c(4, 8),
          div(numericInput("ni_margin", "Non-inferiority margin (HR)", 1.3, 1.01, 3, 0.01),
              numericInput("ni_true", "Assumed true HR", 1, 0.5, 1.5, 0.01),
              numericInput("ni_alpha", "One-sided α", 0.025, 0.001, 0.1, 0.005),
              numericInput("ni_power", "Power", 0.8, 0.5, 0.99, 0.05)),
          div(uiOutput("ss_ni"), formula_card("$$D=\\frac{4\\,(z_{1-\\alpha}+z_{1-\\beta})^2}{(\\ln \\Delta-\\ln \\mathrm{HR})^2}$$")))),
      nav_panel("Two proportions",
        layout_columns(col_widths = c(4, 8),
          div(numericInput("bp0", "Control rate", 0.625, 0.01, 0.99, 0.005), numericInput("bp1", "Experimental rate", 0.75, 0.01, 0.99, 0.005),
              numericInput("b_alpha", "Two-sided α", 0.05, 0.001, 0.2, 0.005), numericInput("b_power", "Power", 0.8, 0.5, 0.99, 0.05)),
          div(uiOutput("ss_bin"), formula_card("$$n=\\frac{\\left[z_{1-\\alpha/2}\\sqrt{2\\bar p\\bar q}+z_{1-\\beta}\\sqrt{p_1q_1+p_2q_2}\\right]^2}{(p_1-p_2)^2}$$")))),
      nav_panel("Two means",
        layout_columns(col_widths = c(4, 8),
          div(numericInput("m_delta", "Difference worth detecting δ", 10), numericInput("m_sd", "Standard deviation σ", 20),
              numericInput("m_alpha", "Two-sided α", 0.05, 0.001, 0.2, 0.005), numericInput("m_power", "Power", 0.8, 0.5, 0.99, 0.05)),
          div(uiOutput("ss_mean"), formula_card("$$n=\\frac{2\\sigma^2(z_{1-\\alpha/2}+z_{1-\\beta})^2}{\\delta^2}$$")))),
      nav_panel("Single-arm phase II",
        layout_columns(col_widths = c(4, 8),
          div(numericInput("s_p0", "Uninteresting response rate p₀", 0.2, 0.01, 0.95, 0.01), numericInput("s_p1", "Promising rate p₁", 0.4, 0.02, 0.99, 0.01),
              numericInput("s_alpha", "One-sided α", 0.05, 0.001, 0.2, 0.005), numericInput("s_power", "Power", 0.8, 0.5, 0.99, 0.05),
              numericInput("s_nmax", "Largest total N searched", 60, 10, 90, 5),
              actionButton("s_go", "Find the designs", class = "btn-primary")),
          div(tableOutput("simon_tab"), note("Optimal minimises the expected N if the drug is inactive; minimax minimises the maximum N. Stop after stage 1 if responses ≤ r₁; declare promising if total responses > r."))))
    )
  ),
  # 6 ---------------------------------------------------------------------------
  nav_panel("6 · Diagnostics",
    layout_sidebar(
      sidebar = sidebar(width = 320,
        selectInput("ref_col", "Reference standard", NULL), selectInput("ref_pos", "Value meaning condition present", NULL),
        radioButtons("test_type", "Index test", c("Binary result" = "bin", "Continuous marker with threshold" = "cont")),
        conditionalPanel("input.test_type == 'bin'", selectInput("test_col", "Test result", NULL), selectInput("test_pos", "Value meaning positive", NULL)),
        conditionalPanel("input.test_type == 'cont'", selectInput("marker_col", "Marker", NULL), numericInput("thr", "Threshold (≥ = positive)", NA)),
        sliderInput("prev", "Prevalence for predictive values", 0.01, 0.8, 0.30, 0.01)
      ),
      layout_columns(col_widths = c(5, 7),
        card(card_header("2 × 2 table and accuracy"), tableOutput("tab22"), tableOutput("acc_tab"),
             formula_card("$$\\mathrm{PPV}=\\frac{\\mathrm{Se}\\,p}{\\mathrm{Se}\\,p+(1-\\mathrm{Sp})(1-p)}\\qquad \\mathrm{LR}^+=\\frac{\\mathrm{Se}}{1-\\mathrm{Sp}}$$")),
        card(card_header("Predictive values across prevalence"), plotOutput("prev_plot", height = "330px"),
             conditionalPanel("input.test_type == 'cont'", plotOutput("roc_plot", height = "360px"), uiOutput("roc_text")))
      )
    )
  ),
  # 7 ---------------------------------------------------------------------------
  nav_panel("7 · Agreement",
    layout_sidebar(
      sidebar = sidebar(width = 300,
        h5("Continuous measurements"), selectInput("ag_a", "Method or reader A", NULL), selectInput("ag_b", "Method or reader B", NULL),
        numericInput("clin_lim", "Clinically acceptable difference (optional)", NA),
        hr(), h5("Categories"), selectInput("kp_a", "Rater A", NULL), selectInput("kp_b", "Rater B", NULL)),
      layout_columns(col_widths = c(6, 6),
        card(card_header("Bland–Altman"), plotOutput("ba_plot", height = "360px"), tableOutput("ba_tab"),
             formula_card("$$\\bar d \\pm 1.96\\,s_d$$")),
        card(card_header("Scatter and ICC"), plotOutput("sc_plot", height = "360px"), tableOutput("icc_tab"),
             formula_card("$$\\mathrm{ICC}(2,1)=\\frac{MS_R-MS_E}{MS_R+(k-1)MS_E+\\frac{k}{n}(MS_C-MS_E)}$$"))
      ),
      card(card_header("Cohen's kappa"), layout_columns(col_widths = c(6, 6), tableOutput("kp_table"), div(tableOutput("kp_tab"),
           formula_card("$$\\kappa=\\frac{p_o-p_e}{1-p_e}$$"))))
    )
  ),
  # 8 ---------------------------------------------------------------------------
  nav_panel("8 · Report",
    card(card_header("Result sentences from your data"), uiOutput("sentences"),
         downloadButton("dl_results", "Download the results (CSV)", class = "btn-primary"),
         note("Check every number against the tables before it goes into a manuscript."))
  )
)

# ---- server -----------------------------------------------------------------------
server <- function(input, output, session) {
  rv <- reactiveValues(df = NULL)
  # display settings: time unit and colours
  unit <- reactive(if (is.null(input$unit_out)) "Months" else input$unit_out)
  u <- reactive(tolower(unit()))
  xlab_r <- reactive(paste(unit(), "from randomisation"))
  dg <- reactive(if (unit() == "Years") 2 else 1)
  ul <- function(v) if (isTRUE(v == 1)) sub("s$", "", u()) else u()
  pal <- reactive({
    if (is.null(input$palette) || input$palette != "Custom") return(PALETTES[[if (is.null(input$palette)) 1 else input$palette]])
    sapply(1:6, function(i) { v <- input[[paste0("ccol", i)]]; if (is.null(v) || !nzchar(v)) ARM_COLS[i] else v })
  })
  output$pal_preview <- renderUI(swatches(pal()))
  prev_unit <- reactiveVal("Months")
  observeEvent(input$unit_out, {
    f <- UNIT_DAYS[[prev_unit()]] / UNIT_DAYS[[input$unit_out]]; lab <- tolower(input$unit_out)
    sc <- function(v) if (is.null(v) || is.na(v)) NA else signif(v * f, 3)
    updateNumericInput(session, "landmark", label = paste0("Landmark time (", lab, ")"), value = sc(input$landmark))
    updateNumericInput(session, "xmax", label = paste0("Maximum time shown (", lab, "; blank = all)"), value = sc(input$xmax))
    updateNumericInput(session, "tau", label = paste0("RMST horizon τ (", lab, ")"), value = sc(input$tau))
    updateNumericInput(session, "lm2", label = paste0("Landmark for difference and NNT (", lab, ")"), value = sc(input$lm2))
    prev_unit(input$unit_out)
  }, ignoreInit = TRUE)
  observeEvent(input$use_example, rv$df <- make_example())
  observeEvent(input$file, {
    d <- tryCatch(read_any(input$file$datapath, input$file$name), error = function(e) NULL)
    if (is.null(d)) showNotification("The file could not be read. Save it as CSV or XLSX with one header row.", type = "error")
    else rv$df <- d
  })
  output$dl_example <- downloadHandler("gsp_example_dataset.csv", function(f) write.csv(make_example(), f, row.names = FALSE))

  # populate mappings
  observeEvent(rv$df, {
    d <- rv$df; nm <- names(d)
    num <- nm[sapply(d, is.numeric)]; lowcat <- nm[sapply(d, function(v) length(unique(na.omit(v))) <= 10)]
    pick <- function(cands, pool) { hit <- intersect(cands, pool); if (length(hit)) hit[1] else pool[1] }
    updateSelectInput(session, "time", choices = num, selected = pick(c("pfs_months", "time", "os_months"), num))
    updateSelectInput(session, "status", choices = lowcat, selected = pick(c("pfs_event", "status", "event"), lowcat))
    updateSelectInput(session, "arm", choices = lowcat, selected = pick(c("arm", "treatment", "group"), lowcat))
    updateSelectizeInput(session, "t1vars", choices = nm, selected = intersect(c("age", "sex", "ps", "brain_mets", "mutation"), nm))
    updateSelectizeInput(session, "adjust", choices = nm, selected = NULL)
    updateSelectInput(session, "grp1", choices = lowcat, selected = pick(c("arm", "treatment", "group"), lowcat))
    updateSelectInput(session, "grp2", choices = c("None", lowcat), selected = "None")
    updateSelectizeInput(session, "subvars", choices = lowcat, selected = intersect(c("sex", "ps", "brain_mets", "mutation"), nm))
    updateSelectInput(session, "ref_col", choices = lowcat, selected = pick(c("tissue_egfr", "reference"), lowcat))
    updateSelectInput(session, "test_col", choices = lowcat, selected = pick(c("plasma_egfr", "test"), lowcat))
    updateSelectInput(session, "marker_col", choices = num, selected = pick(c("ctdna_level", "marker"), num))
    updateSelectInput(session, "ag_a", choices = num, selected = pick(c("reader_a_mm"), num))
    updateSelectInput(session, "ag_b", choices = num, selected = pick(c("reader_b_mm"), num))
    updateSelectInput(session, "kp_a", choices = lowcat, selected = pick(c("recist_reader_a"), lowcat))
    updateSelectInput(session, "kp_b", choices = lowcat, selected = pick(c("recist_reader_b"), lowcat))
  })
  observeEvent(list(rv$df, input$status), {
    req(rv$df, input$status %in% names(rv$df)); v <- sort(unique(na.omit(rv$df[[input$status]])))
    updateSelectInput(session, "event_value", choices = v, selected = if (1 %in% v) 1 else v[length(v)])
  })
  observeEvent(list(rv$df, input$arm), {
    req(rv$df, input$arm %in% names(rv$df)); v <- sort(unique(na.omit(as.character(rv$df[[input$arm]]))))
    updateSelectInput(session, "ref_arm", choices = v, selected = if ("Gefitinib" %in% v) "Gefitinib" else v[1])
  })
  observeEvent(list(rv$df, input$ref_col), {
    req(rv$df, input$ref_col %in% names(rv$df)); v <- sort(unique(na.omit(as.character(rv$df[[input$ref_col]]))))
    updateSelectInput(session, "ref_pos", choices = v, selected = if ("Positive" %in% v) "Positive" else v[length(v)])
  })
  observeEvent(list(rv$df, input$test_col), {
    req(rv$df, input$test_col %in% names(rv$df)); v <- sort(unique(na.omit(as.character(rv$df[[input$test_col]]))))
    updateSelectInput(session, "test_pos", choices = v, selected = if ("Positive" %in% v) "Positive" else v[length(v)])
  })

  # analysis dataset
  sdat <- reactive({
    req(rv$df, input$time, input$status, input$arm, input$event_value, input$ref_arm)
    d <- rv$df; req(all(c(input$time, input$status, input$arm) %in% names(d)))
    fin <- if (is.null(input$unit_in)) "Months" else input$unit_in
    d$.time <- suppressWarnings(as.numeric(d[[input$time]])) * UNIT_DAYS[[fin]] / UNIT_DAYS[[unit()]]
    d$.status <- as.integer(as.character(d[[input$status]]) == as.character(input$event_value))
    a <- as.character(d[[input$arm]]); lv <- c(input$ref_arm, setdiff(sort(unique(na.omit(a))), input$ref_arm))
    d$.arm <- factor(a, levels = lv)
    d[!is.na(d$.time) & !is.na(d$.arm), ]
  })
  fit_arm <- reactive(survfit(Surv(.time, .status) ~ .arm, data = sdat(), conf.type = "log-log"))

  # 1 · checks and Table 1
  output$preview <- renderDT(datatable(req(rv$df), options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE))
  output$checks <- renderTable({
    d <- sdat(); raw <- rv$df
    by <- split(d, d$.arm)
    rbind(
      data.frame(Check = "Patients per arm", Result = paste(names(by), sapply(by, nrow), sep = ": ", collapse = " · ")),
      data.frame(Check = "Events per arm", Result = paste(names(by), sapply(by, function(x) sum(x$.status)), sep = ": ", collapse = " · ")),
      data.frame(Check = "Censored per arm", Result = paste(names(by), sapply(by, function(x) sum(1 - x$.status)), sep = ": ", collapse = " · ")),
      data.frame(Check = "Shortest and longest time", Result = paste(fmt(min(d$.time), dg()), "to", fmt(max(d$.time), dg()), u())),
      data.frame(Check = "Times of zero or below", Result = sum(d$.time <= 0)),
      data.frame(Check = "Values of the event indicator", Result = paste(sort(unique(raw[[input$status]])), collapse = ", ")),
      data.frame(Check = "Rows dropped (missing time or arm)", Result = nrow(raw) - nrow(d)))
  })
  output$check_notes <- renderUI({
    d <- sdat(); msgs <- NULL
    if (any(d$.time <= 0)) msgs <- c(msgs, "Some times are zero or negative: check the source dates.")
    if (nlevels(d$.arm) < 2) msgs <- c(msgs, "Only one arm found: comparisons need two.")
    if (length(msgs)) note(HTML(paste(msgs, collapse = "<br>"))) else note("No problems found in the six checks.")
  })
  output$table1 <- renderTable({
    d <- sdat(); req(length(input$t1vars) > 0); arms <- levels(d$.arm)
    out <- data.frame(Characteristic = character(), check.names = FALSE)
    rows <- list()
    for (v in input$t1vars) {
      x <- d[[v]]
      if (is.numeric(x) && length(unique(na.omit(x))) > 10) {
        r <- c(paste0(v, ", median (IQR)"), sapply(arms, function(a) { y <- x[d$.arm == a]
          paste0(fmt(median(y, na.rm = TRUE), 1), " (", fmt(quantile(y, .25, na.rm = TRUE), 1), "–", fmt(quantile(y, .75, na.rm = TRUE), 1), ")") }))
        rows[[length(rows) + 1]] <- r
      } else {
        rows[[length(rows) + 1]] <- c(paste0(v, ", n (%)"), rep("", length(arms)))
        for (l in sort(unique(na.omit(as.character(x))))) rows[[length(rows) + 1]] <- c(paste0("   ", l), sapply(arms, function(a) {
          y <- as.character(x[d$.arm == a]); n <- sum(!is.na(y)); k <- sum(y == l, na.rm = TRUE); paste0(k, " (", fmt(100 * k / n, 1), "%)") }))
      }
      miss <- sapply(arms, function(a) sum(is.na(x[d$.arm == a])))
      if (any(miss > 0)) rows[[length(rows) + 1]] <- c("   missing", miss)
    }
    t1 <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
    names(t1) <- c("Characteristic", paste0(arms, " (n = ", as.integer(table(d$.arm)), ")")); t1
  })

  # 2 · survival
  km_plot_obj <- reactive(plot_km(fit_arm(), input$xmax, input$show_ci, input$landmark, "Event-free probability", pal(), xlab_r()))
  output$km_plot <- renderPlot(km_plot_obj(), res = 100)
  output$dl_km <- downloadHandler("km_curve.png", function(f) ggsave(f, km_plot_obj(), width = 10, height = 6.5, dpi = 300, bg = BG))
  output$median_tab <- renderTable(km_median_table(fit_arm(), dg()), digits = 0)
  output$landmark_tab <- renderTable({
    lt <- landmark_table(fit_arm(), input$landmark)
    data.frame(Group = lt$Group, `Landmark` = paste(input$landmark, ul(input$landmark)),
               `Estimate (95% CI)` = paste0(fmt(100 * lt$S, 1), "% (", fmt(100 * lt$lower, 1), "–", fmt(100 * lt$upper, 1), ")"), check.names = FALSE)
  })
  fu_fit <- reactive(survfit(Surv(.time, 1 - .status) ~ 1, data = sdat()))
  output$fu_plot <- renderPlot({
    f <- fu_fit(); d <- data.frame(time = c(0, f$time), surv = c(1, f$surv)); m <- summary(f)$table["median"]
    ggplot(d, aes(time, surv)) + geom_step(colour = BL, linewidth = 1.2) +
      geom_segment(aes(x = 0, xend = m, y = 0.5, yend = 0.5), linetype = "dashed", colour = BL) +
      geom_segment(aes(x = m, xend = m, y = 0, yend = 0.5), linetype = "dashed", colour = BL) +
      geom_vline(xintercept = median(sdat()$.time), linetype = "dotted", colour = AM, linewidth = 1) +
      labs(x = xlab_r(), y = "Still under follow-up") + theme_gsp()
  }, res = 100)
  output$fu_text <- renderUI({
    m <- summary(fu_fit())$table
    note(HTML(paste0("Median follow-up by reverse Kaplan–Meier: <b>", fmt(m["median"], dg()), " ", u(), "</b> (95% CI ", fmt(m["0.95LCL"], dg()), "–", fmt(m["0.95UCL"], dg()),
                     "). The naive median of all observed times is ", fmt(median(sdat()$.time), dg()), " ", u(), " (dotted line) and understates follow-up.")))
  })
  output$km_hand <- renderTable({
    d <- sdat(); ev <- sort(unique(d$.time[d$.status == 1]))[1:min(10, length(unique(d$.time[d$.status == 1])))]
    S <- 1; rows <- NULL
    for (t in sort(unique(d$.time[d$.time <= max(ev)]))) {
      n <- sum(d$.time >= t); dd <- sum(d$.time == t & d$.status == 1); cc <- sum(d$.time == t & d$.status == 0)
      if (dd > 0) S <- S * (1 - dd / n)
      if (dd > 0) rows <- rbind(rows, data.frame(Time = fmt(t, dg() + 1), `At risk` = n, Events = dd, Censored = cc, `1 − d/n` = fmt(1 - dd / n, 3), `S(t)` = fmt(S, 3), check.names = FALSE))
    }
    rows
  })

  # 3 · compare (groups chosen by the user)
  grp_vars <- reactive({ req(input$grp1); v <- input$grp1; if (!is.null(input$grp2) && input$grp2 != "None" && input$grp2 != input$grp1) v <- c(v, input$grp2); v })
  grp_raw <- reactive({
    d <- sdat(); v <- grp_vars(); req(all(v %in% names(d)))
    g <- if (length(v) == 1) as.character(d[[v]]) else ifelse(is.na(d[[v[1]]]) | is.na(d[[v[2]]]), NA, paste(d[[v[1]]], d[[v[2]]], sep = " · "))
    d$.grp_chr <- g; d[!is.na(g), ]
  })
  observeEvent(grp_raw(), {
    lv <- sort(unique(grp_raw()$.grp_chr)); cur <- isolate(input$grp_ref)
    sel <- if (!is.null(cur) && cur %in% lv) cur else if (identical(grp_vars(), input$arm) && input$ref_arm %in% lv) input$ref_arm else lv[1]
    updateSelectInput(session, "grp_ref", choices = lv, selected = sel)
  })
  cdat <- reactive({
    d <- grp_raw(); lv <- sort(unique(d$.grp_chr)); ref <- if (!is.null(input$grp_ref) && input$grp_ref %in% lv) input$grp_ref else lv[1]
    d$.grp <- factor(d$.grp_chr, levels = c(ref, setdiff(lv, ref))); validate(need(nlevels(d$.grp) >= 2, "Choose a variable with at least two groups.")); d
  })
  fit_grp <- reactive(survfit(Surv(.time, .status) ~ .grp, data = cdat(), conf.type = "log-log"))
  output$cmp_title <- renderText(paste("Kaplan–Meier curves by", paste(grp_vars(), collapse = " and ")))
  cmp_plot_obj <- reactive(plot_km(fit_grp(), input$xmax, input$cmp_ci, NA, "Event-free probability", pal(), xlab_r()))
  output$cmp_km <- renderPlot(cmp_plot_obj(), res = 100)
  output$dl_cmp_km <- downloadHandler("km_by_group.png", function(f) ggsave(f, cmp_plot_obj(), width = 10, height = 6.5, dpi = 300, bg = BG))
  output$cmp_median <- renderTable(km_median_table(fit_grp(), dg()), digits = 0)
  lr <- reactive(survdiff(Surv(.time, .status) ~ .grp, data = cdat()))
  output$lr_tab <- renderTable({
    l <- lr(); data.frame(Group = sub(".*=", "", names(l$n)), N = as.integer(l$n), `Observed O` = as.integer(l$obs), `Expected E` = fmt(l$exp, 1), `O ÷ E` = fmt(l$obs / l$exp, 2), check.names = FALSE)
  })
  output$lr_text <- renderUI({ l <- lr(); df <- length(l$n) - 1; p <- pchisq(l$chisq, df, lower.tail = FALSE)
    note(HTML(paste0("χ² = ", fmt(l$chisq, 2), " on ", df, " df, P ", ifelse(p < 0.001, "< 0.001", paste("=", fmtp(p))),
      if (df > 1) ": a global test that the groups share one survival curve." else ".", " The log-rank test gives a P value only; the Cox model gives the size of the effect."))) })
  cox <- reactive({
    d <- cdat(); adj <- setdiff(input$adjust, c(input$time, input$status, grp_vars()))
    fml <- as.formula(paste("Surv(.time, .status) ~ .grp", if (length(adj)) paste("+", paste(sprintf("`%s`", adj), collapse = " + ")) else ""))
    coxph(fml, data = d)
  })
  output$cox_tab <- renderTable({
    s <- summary(cox()); co <- s$coefficients; ci <- s$conf.int
    data.frame(Term = sub("^\\.grp", "", rownames(co)), coef = fmt(co[, "coef"], 3), HR = fmt(co[, "exp(coef)"], 3), SE = fmt(co[, "se(coef)"], 3),
               z = fmt(co[, "z"], 2), P = fmtp(co[, "Pr(>|z|)"]), `95% CI` = paste0(fmt(ci[, "lower .95"], 3), " – ", fmt(ci[, "upper .95"], 3)), check.names = FALSE)
  })
  output$cox_text <- renderUI({ s <- summary(cox()); ci <- s$conf.int; rows <- grep("^\\.grp", rownames(ci)); ref <- levels(cdat()$.grp)[1]
    txt <- sapply(rows, function(i) paste0("<b>", sub("^\\.grp", "", rownames(ci)[i]), "</b> against <b>", ref, "</b>: HR ", fmt(ci[i, 1], 2), " (95% CI ", fmt(ci[i, 3], 2), "–", fmt(ci[i, 4], 2),
      "); at any time the event rate is ", fmt(100 * ci[i, 1], 0), "% of the reference rate."))
    note(HTML(paste(txt, collapse = "<br>"))) })
  output$ll_plot <- renderPlot({
    d <- km_steps(fit_grp()); d <- d[d$time > 0 & d$surv > 0 & d$surv < 1, ]
    lev <- unique(d$strata); d$strata <- factor(d$strata, levels = lev)
    ggplot(d, aes(log(time), log(-log(surv)), colour = strata)) + geom_step(linewidth = 1.1) +
      scale_colour_manual(values = setNames(rep(pal(), length.out = length(lev)), lev)) + guides(colour = guide_legend(nrow = if (length(lev) > 2) 2 else 1)) +
      labs(x = paste0("log(", u(), ")"), y = "log(−log S(t))") + theme_gsp()
  }, res = 100)
  output$ph_tab <- renderTable({ z <- cox.zph(cox()); tb <- z$table
    data.frame(Term = sub("^\\.grp", "Groups", rownames(tb)), `χ²` = fmt(tb[, "chisq"], 2), df = as.integer(tb[, "df"]), P = fmtp(tb[, "p"]), check.names = FALSE) })
  output$rmst_tab <- renderTable({
    d <- cdat(); gl <- levels(d$.grp); tau <- input$tau
    rr <- lapply(gl, function(a) rmst_fun(d$.time[d$.grp == a], d$.status[d$.grp == a], tau))
    out <- data.frame(Group = gl, RMST = sapply(rr, function(x) fmt(x["rmst"], 2)), SE = sapply(rr, function(x) fmt(x["se"], 2)), check.names = FALSE)
    for (j in seq_along(gl)[-1]) { df_ <- rr[[j]]["rmst"] - rr[[1]]["rmst"]; se <- sqrt(rr[[j]]["se"]^2 + rr[[1]]["se"]^2)
      out <- rbind(out, data.frame(Group = paste0("Difference: ", gl[j], " − ", gl[1]), RMST = paste0(fmt(df_, 2), " (", fmt(df_ - 1.96 * se, 2), " to ", fmt(df_ + 1.96 * se, 2), ")"), SE = fmt(se, 2), check.names = FALSE)) }
    names(out)[2] <- paste0("RMST (", u(), ")"); out
  })
  output$lmdiff_tab <- renderTable({
    lt <- landmark_table(fit_grp(), input$lm2); req(nrow(lt) >= 2)
    rows <- data.frame(Quantity = paste0("S(", input$lm2, ") ", lt$Group), Value = paste0(fmt(100 * lt$S, 1), "%"))
    for (j in 2:nrow(lt)) { d <- lt$S[j] - lt$S[1]; se <- sqrt(lt$SE[j]^2 + lt$SE[1]^2)
      rows <- rbind(rows, data.frame(Quantity = c(paste0("Difference: ", lt$Group[j], " − ", lt$Group[1]), paste0("NNT: ", lt$Group[j], " vs ", lt$Group[1])),
        Value = c(paste0(fmt(100 * d, 1), " points (", fmt(100 * (d - 1.96 * se), 1), " to ", fmt(100 * (d + 1.96 * se), 1), ")"), ifelse(d > 0, as.character(ceiling(1 / d)), "not defined (no benefit)")))) }
    rows
  })

  # arm-only models for the report tab
  lr_arm <- reactive(survdiff(Surv(.time, .status) ~ .arm, data = sdat()))
  cox_arm <- reactive(coxph(Surv(.time, .status) ~ .arm, data = sdat()))

  # 4 · subgroups
  forest_df <- reactive({
    d <- sdat(); req(length(input$subvars) > 0, nlevels(d$.arm) == 2)
    rows <- list(); ov <- summary(coxph(Surv(.time, .status) ~ .arm, data = d))$conf.int
    rows[[1]] <- data.frame(label = "Overall", n = nrow(d), hr = ov[1, 1], lo = ov[1, 3], hi = ov[1, 4], pint = NA, head = TRUE)
    for (v in input$subvars) {
      x <- as.character(d[[v]]); dd <- d[!is.na(x), ]; dd$.g <- factor(x[!is.na(x)])
      pint <- tryCatch({ m0 <- coxph(Surv(.time, .status) ~ .arm + .g, data = dd); m1 <- coxph(Surv(.time, .status) ~ .arm * .g, data = dd)
        anova(m0, m1)[2, "Pr(>|Chi|)"] }, error = function(e) NA)
      rows[[length(rows) + 1]] <- data.frame(label = v, n = NA, hr = NA, lo = NA, hi = NA, pint = pint, head = TRUE)
      for (l in levels(dd$.g)) { s <- dd[dd$.g == l, ]
        ci <- tryCatch(summary(coxph(Surv(.time, .status) ~ .arm, data = s))$conf.int, error = function(e) matrix(NA, 1, 4))
        rows[[length(rows) + 1]] <- data.frame(label = paste0("· ", l), n = nrow(s), hr = ci[1, 1], lo = ci[1, 3], hi = ci[1, 4], pint = NA, head = FALSE) }
    }
    f <- do.call(rbind, rows); f$y <- rev(seq_len(nrow(f))); f
  })
  output$forest_ui <- renderUI(plotOutput("forest", height = paste0(max(320, 34 * nrow(forest_df()) + 120), "px")))
  output$forest <- renderPlot({
    f <- forest_df(); p <- f[!is.na(f$hr), ]
    ggplot(p, aes(hr, y)) + geom_vline(xintercept = 1, linetype = "dashed", colour = GR) +
      geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, colour = PL, linewidth = 1.1) +
      geom_point(aes(shape = label == "Overall"), size = 3.6, colour = PL) +
      scale_shape_manual(values = c(15, 18), guide = "none") + scale_x_log10() +
      scale_y_continuous(breaks = f$y, labels = f$label) +
      labs(x = paste0("Hazard ratio, ", levels(sdat()$.arm)[2], " vs ", levels(sdat()$.arm)[1], " (log scale)"), y = NULL) + theme_gsp()
  }, res = 100)
  output$forest_tab <- renderTable({ f <- forest_df()
    data.frame(Subgroup = f$label, N = ifelse(is.na(f$n), "", f$n), `HR (95% CI)` = ifelse(is.na(f$hr), "", paste0(fmt(f$hr, 2), " (", fmt(f$lo, 2), "–", fmt(f$hi, 2), ")")),
               `Interaction P` = ifelse(is.na(f$pint), "", fmtp(f$pint)), check.names = FALSE) })

  # 5 · sample size
  output$ss_tte <- renderUI({
    z <- (qnorm(1 - input$ss_alpha / 2) + qnorm(input$ss_power))^2; k <- input$ss_k
    D <- ceiling((1 + k)^2 / k * z / log(input$ss_hr)^2)
    pc <- p_event(input$ss_med, input$ss_acc, input$ss_fu); pe <- p_event(input$ss_med / input$ss_hr, input$ss_acc, input$ss_fu)
    pavg <- (pc + k * pe) / (1 + k); N <- ceiling(D / pavg); Nd <- ceiling(N / (1 - input$ss_drop))
    note(HTML(paste0("<b>", D, " events</b> are needed. Probability of an event by the analysis: control ", fmt(pc, 2), ", experimental ", fmt(pe, 2),
      " → <b>", N, " patients</b>, or <b>", Nd, "</b> allowing for ", fmt(100 * input$ss_drop, 0), "% loss. Power comes from events, not patients.")))
  })
  output$ss_curve <- renderPlot({
    hr <- seq(0.5, 0.85, 0.01); k <- input$ss_k
    d <- do.call(rbind, lapply(c(0.8, 0.9), function(pw) data.frame(hr = hr, D = (1 + k)^2 / k * (qnorm(1 - input$ss_alpha / 2) + qnorm(pw))^2 / log(hr)^2, power = paste0(pw * 100, "% power"))))
    ggplot(d, aes(hr, D, colour = power)) + geom_line(linewidth = 1.3) + geom_vline(xintercept = input$ss_hr, linetype = "dashed", colour = GR) +
      scale_colour_manual(values = c(PL, AM)) + labs(x = "Target hazard ratio", y = "Events required") + theme_gsp()
  }, res = 100)
  output$ss_ni <- renderUI({ D <- ceiling(4 * (qnorm(1 - input$ni_alpha) + qnorm(input$ni_power))^2 / (log(input$ni_margin) - log(input$ni_true))^2)
    note(HTML(paste0("<b>", D, " events</b> to show non-inferiority with margin ", input$ni_margin, " if the true HR is ", input$ni_true, "."))) })
  output$ss_bin <- renderUI({ p1 <- input$bp0; p2 <- input$bp1; pb <- (p1 + p2) / 2
    n <- ceiling((qnorm(1 - input$b_alpha / 2) * sqrt(2 * pb * (1 - pb)) + qnorm(input$b_power) * sqrt(p1 * (1 - p1) + p2 * (1 - p2)))^2 / (p1 - p2)^2)
    note(HTML(paste0("<b>", n, " patients per arm</b> (", 2 * n, " in total)."))) })
  output$ss_mean <- renderUI({ n <- ceiling(2 * input$m_sd^2 * (qnorm(1 - input$m_alpha / 2) + qnorm(input$m_power))^2 / input$m_delta^2)
    note(HTML(paste0("<b>", n, " patients per arm</b>; standardised effect d = ", fmt(input$m_delta / input$m_sd, 2), "."))) })
  simon <- eventReactive(input$s_go, withProgress(message = "Searching designs", {
    list(s = simon_design(input$s_p0, input$s_p1, input$s_alpha, 1 - input$s_power, input$s_nmax), a = ahern(input$s_p0, input$s_p1, input$s_alpha, 1 - input$s_power)) }))
  output$simon_tab <- renderTable({ s <- simon(); f <- function(x, nm) if (is.null(x)) data.frame(Design = nm, `Stage 1: n₁ (stop if ≤ r₁)` = "none found", `Total N (promising if > r)` = "", `P(stop early | p₀)` = "", `Expected N | p₀` = "", check.names = FALSE) else
      data.frame(Design = nm, `Stage 1: n₁ (stop if ≤ r₁)` = paste0(x["n1"], " (", x["r1"], ")"), `Total N (promising if > r)` = paste0(x["n"], " (", x["r"], ")"), `P(stop early | p₀)` = fmt(x["PET0"], 2), `Expected N | p₀` = fmt(x["EN0"], 1), check.names = FALSE)
    rbind(f(s$s$optimal, "Simon optimal"), f(s$s$minimax, "Simon minimax"),
          data.frame(Design = "A'Hern one-stage", `Stage 1: n₁ (stop if ≤ r₁)` = "—", `Total N (promising if > r)` = paste0(s$a["n"], " (", s$a["r"] - 1, ")"), `P(stop early | p₀)` = "—", `Expected N | p₀` = s$a["n"], check.names = FALSE)) })

  # 6 · diagnostics
  diag_data <- reactive({
    d <- req(rv$df); req(input$ref_col %in% names(d), input$ref_pos)
    ref <- as.integer(as.character(d[[input$ref_col]]) == input$ref_pos)
    if (input$test_type == "bin") { req(input$test_col %in% names(d), input$test_pos); tst <- as.integer(as.character(d[[input$test_col]]) == input$test_pos); mk <- NULL }
    else { req(input$marker_col %in% names(d)); mk <- as.numeric(d[[input$marker_col]])
      thr <- if (is.na(input$thr)) roc_fun(mk, ref)$cut else input$thr; tst <- as.integer(mk >= thr) }
    ok <- !is.na(ref) & !is.na(tst); list(ref = ref[ok], tst = tst[ok], mk = if (is.null(mk)) NULL else mk[ok])
  })
  acc <- reactive({ x <- diag_data(); TP <- sum(x$tst == 1 & x$ref == 1); FP <- sum(x$tst == 1 & x$ref == 0); FN <- sum(x$tst == 0 & x$ref == 1); TN <- sum(x$tst == 0 & x$ref == 0)
    list(TP = TP, FP = FP, FN = FN, TN = TN, se = wilson(TP, TP + FN), sp = wilson(TN, TN + FP), ppv = wilson(TP, TP + FP), npv = wilson(TN, TN + FN), prev = (TP + FN) / (TP + FP + FN + TN)) })
  output$tab22 <- renderTable({ a <- acc(); data.frame(` ` = c("Test positive", "Test negative", "Total"), `Condition present` = c(a$TP, a$FN, a$TP + a$FN), `Condition absent` = c(a$FP, a$TN, a$FP + a$TN), check.names = FALSE) })
  output$acc_tab <- renderTable({ a <- acc(); r <- function(v) paste0(fmt(100 * v[1], 1), "% (", fmt(100 * v[2], 1), "–", fmt(100 * v[3], 1), ")")
    lrp <- a$se[1] / (1 - a$sp[1]); lrn <- (1 - a$se[1]) / a$sp[1]
    data.frame(Measure = c("Sensitivity", "Specificity", "PPV (sample prevalence)", "NPV (sample prevalence)", "Prevalence in sample", "LR+", "LR−"),
               `Estimate (95% CI)` = c(r(a$se), r(a$sp), r(a$ppv), r(a$npv), paste0(fmt(100 * a$prev, 1), "%"), fmt(lrp, 2), fmt(lrn, 2)), check.names = FALSE) })
  output$prev_plot <- renderPlot({ a <- acc(); se <- a$se[1]; sp <- a$sp[1]; p <- seq(0.01, 0.8, 0.005)
    d <- rbind(data.frame(p = p, v = se * p / (se * p + (1 - sp) * (1 - p)), m = "PPV"), data.frame(p = p, v = sp * (1 - p) / (sp * (1 - p) + (1 - se) * p), m = "NPV"))
    pp <- input$prev; pt <- data.frame(p = pp, v = c(se * pp / (se * pp + (1 - sp) * (1 - pp)), sp * (1 - pp) / (sp * (1 - pp) + (1 - se) * pp)), m = c("PPV", "NPV"))
    ggplot(d, aes(100 * p, 100 * v, colour = m)) + geom_line(linewidth = 1.3) + geom_point(data = pt, size = 4) +
      geom_label(data = pt, aes(label = paste0(m, " ", fmt(100 * v, 0), "%")), hjust = 0, nudge_x = 2, nudge_y = c(-6, 6), show.legend = FALSE, size = 4.3, fill = "#FFFFFF", label.size = 0) +
      scale_colour_manual(values = c(NPV = AM, PPV = PL)) + labs(x = "Prevalence (%)", y = "Predictive value (%)") + ylim(0, 100) + theme_gsp() }, res = 100)
  roc <- reactive({ x <- diag_data(); req(!is.null(x$mk)); roc_fun(x$mk, x$ref) })
  output$roc_plot <- renderPlot({ r <- roc(); ggplot(r$curve, aes(fpr, tpr)) + geom_abline(linetype = "dotted", colour = GR) +
      geom_step(colour = PL, linewidth = 1.3) + annotate("point", x = 1 - r$spec, y = r$sens, colour = AM, size = 5) +
      annotate("text", x = 0.62, y = 0.2, label = paste0("AUC ", fmt(r$auc, 2), " (", fmt(r$auc - 1.96 * r$se, 2), "–", fmt(min(1, r$auc + 1.96 * r$se), 2), ")"), colour = PL, size = 6, fontface = "bold") +
      coord_equal() + labs(x = "1 − specificity", y = "Sensitivity") + theme_gsp() }, res = 100)
  output$roc_text <- renderUI({ r <- roc(); note(HTML(paste0("Youden cut-off ", fmt(r$cut, 2), ": sensitivity ", fmt(100 * r$sens, 0), "%, specificity ", fmt(100 * r$spec, 0),
      "%. Leave the threshold blank to use it. CI by the Hanley–McNeil method."))) })

  # 7 · agreement
  ag <- reactive({ d <- req(rv$df); req(input$ag_a %in% names(d), input$ag_b %in% names(d)); a <- as.numeric(d[[input$ag_a]]); b <- as.numeric(d[[input$ag_b]])
    ok <- !is.na(a) & !is.na(b); data.frame(a = a[ok], b = b[ok]) })
  output$ba_plot <- renderPlot({ d <- ag(); df <- d$b - d$a; mn <- (d$a + d$b) / 2; bias <- mean(df); s <- sd(df)
    g <- ggplot(data.frame(mn, df), aes(mn, df)) + geom_point(colour = PL, alpha = 0.7, size = 2.4) +
      geom_hline(yintercept = bias, colour = PL, linewidth = 1.2) + geom_hline(yintercept = bias + c(-1.96, 1.96) * s, colour = AM, linetype = "dashed", linewidth = 1) +
      geom_hline(yintercept = 0, colour = GR, linetype = "dotted") + labs(x = "Mean of A and B", y = "B − A") + theme_gsp()
    if (!is.na(input$clin_lim)) g <- g + geom_hline(yintercept = c(-1, 1) * input$clin_lim, colour = BL, linetype = "dotdash")
    g }, res = 100)
  output$ba_tab <- renderTable({ d <- ag(); df <- d$b - d$a; n <- length(df); bias <- mean(df); s <- sd(df); se_b <- s / sqrt(n); se_l <- s * sqrt(3 / n); tq <- qt(0.975, n - 1)
    data.frame(Quantity = c("n pairs", "Bias (mean B − A)", "Lower limit of agreement", "Upper limit of agreement", "Pearson r (association only)"),
               `Estimate (95% CI)` = c(n, paste0(fmt(bias, 2), " (", fmt(bias - tq * se_b, 2), " to ", fmt(bias + tq * se_b, 2), ")"),
                                        paste0(fmt(bias - 1.96 * s, 2), " (", fmt(bias - 1.96 * s - tq * se_l, 2), " to ", fmt(bias - 1.96 * s + tq * se_l, 2), ")"),
                                        paste0(fmt(bias + 1.96 * s, 2), " (", fmt(bias + 1.96 * s - tq * se_l, 2), " to ", fmt(bias + 1.96 * s + tq * se_l, 2), ")"), fmt(cor(d$a, d$b), 3)), check.names = FALSE) })
  output$sc_plot <- renderPlot({ d <- ag(); lim <- range(c(d$a, d$b))
    ggplot(d, aes(a, b)) + geom_abline(linetype = "dashed", colour = GR) + geom_point(colour = PL, alpha = 0.7, size = 2.4) + coord_equal(xlim = lim, ylim = lim) +
      labs(x = input$ag_a, y = input$ag_b) + theme_gsp() }, res = 100)
  output$icc_tab <- renderTable({ i <- icc_fun(ag()); t <- i$table; data.frame(Form = t$Form, `ICC (95% CI)` = paste0(fmt(t$ICC, 3), " (", fmt(t$Lower, 3), "–", fmt(t$Upper, 3), ")"), check.names = FALSE) })
  kp <- reactive({ d <- req(rv$df); req(input$kp_a %in% names(d), input$kp_b %in% names(d)); kappa_fun(as.character(d[[input$kp_a]]), as.character(d[[input$kp_b]])) })
  output$kp_table <- renderTable({ k <- kp(); m <- as.data.frame.matrix(k$tab); cbind(`A \\ B` = rownames(m), m) }, rownames = FALSE)
  output$kp_tab <- renderTable({ k <- kp(); data.frame(Quantity = c("Observed agreement pₒ", "Chance agreement pₑ", "Cohen's κ (95% CI)"),
      Value = c(fmt(k$po, 3), fmt(k$pe, 3), paste0(fmt(k$kappa, 2), " (", fmt(k$lower, 2), "–", fmt(k$upper, 2), ")"))) })

  # 8 · report
  results <- reactive({
    fa <- fit_arm(); tb <- summary(fa)$table; arms <- levels(sdat()$.arm); hr <- summary(cox_arm())$conf.int[1, ]; l <- lr_arm()
    fu <- summary(fu_fit())$table; lt <- landmark_table(fa, input$landmark)
    list(tb = tb, arms = arms, hr = hr, p = pchisq(l$chisq, length(l$n) - 1, lower.tail = FALSE), fu = fu, lt = lt)
  })
  output$sentences <- renderUI({
    r <- results(); tb <- r$tb; a <- r$arms; md <- function(i) paste0(fmt(tb[i, "median"], dg()), " ", u(), " (95% CI ", fmt(tb[i, "0.95LCL"], dg()), "–", ifelse(is.na(tb[i, "0.95UCL"]), "NR", fmt(tb[i, "0.95UCL"], dg())), ")")
    s1 <- paste0("At a median follow-up of ", fmt(r$fu["median"], dg()), " ", u(), " (reverse Kaplan–Meier), the median was ", md(2), " with ", a[2], " and ", md(1), " with ", a[1],
                 "; HR ", fmt(r$hr[1], 2), " (95% CI ", fmt(r$hr[3], 2), "–", fmt(r$hr[4], 2), "), log-rank P ", ifelse(r$p < 0.001, "< 0.001", paste("=", fmtp(r$p))), ".")
    s2 <- paste0("At ", input$landmark, " ", ul(input$landmark), ", ", fmt(100 * r$lt$S[2], 1), "% of patients were event-free with ", a[2], " and ", fmt(100 * r$lt$S[1], 1), "% with ", a[1], ".")
    out <- list(div(class = "sentence", s1), div(class = "sentence", s2))
    ac <- tryCatch(acc(), error = function(e) NULL)
    if (!is.null(ac)) out <- c(out, list(div(class = "sentence", paste0("Sensitivity was ", fmt(100 * ac$se[1], 1), "% (95% CI ", fmt(100 * ac$se[2], 1), "–", fmt(100 * ac$se[3], 1),
      ") and specificity ", fmt(100 * ac$sp[1], 1), "% (", fmt(100 * ac$sp[2], 1), "–", fmt(100 * ac$sp[3], 1), "); at a prevalence of ", fmt(100 * ac$prev, 1), "%, PPV was ", fmt(100 * ac$ppv[1], 1), "% and NPV ", fmt(100 * ac$npv[1], 1), "%."))))
    g <- tryCatch(ag(), error = function(e) NULL)
    if (!is.null(g)) { df <- g$b - g$a; i <- icc_fun(g)$table
      out <- c(out, list(div(class = "sentence", paste0("The mean difference between ", input$ag_b, " and ", input$ag_a, " was ", fmt(mean(df), 2), " (95% limits of agreement ", fmt(mean(df) - 1.96 * sd(df), 2), " to ",
        fmt(mean(df) + 1.96 * sd(df), 2), "); ICC(2,1) ", fmt(i$ICC[1], 3), " (95% CI ", fmt(i$Lower[1], 3), "–", fmt(i$Upper[1], 3), ").")))) }
    tagList(out)
  })
  output$dl_results <- downloadHandler("gsp_results.csv", function(f) {
    r <- results(); tb <- r$tb
    out <- data.frame(item = c(paste("median", r$arms), "HR", "HR lower", "HR upper", "log-rank P", "median follow-up (reverse KM)", paste0("landmark ", input$landmark, " ", r$lt$Group)),
                      value = c(tb[, "median"], r$hr[1], r$hr[3], r$hr[4], r$p, r$fu["median"], r$lt$S))
    write.csv(out, f, row.names = FALSE) })
}

shinyApp(ui, server)
