setwd("C:\\Bayes\\Projeto500\\Reworked\\Resultados_para_salvar\\logs")

# ============================================================
# PIPELINE COMPLETO - EXTRAÇÃO, MÉTRICAS, TABELAS E FIGURAS
# ============================================================

rm(list = ls())

# ------------------------------------------------------------
# 0) PACOTES
# ------------------------------------------------------------
packs <- c("dplyr", "purrr", "stringr", "tidyr", "ggplot2", "readr")
to_install <- packs[!packs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)


library(dplyr)
library(purrr)
library(stringr)
library(tidyr)
library(ggplot2)
library(readr)

# ------------------------------------------------------------
# 1) CONFIGURAÇÕES
# ------------------------------------------------------------
base_dir  <- "C:/Bayes/Projeto500/Reworked/Resultados_para_salvar"
data_file <- file.path(base_dir, "ALLDFs.rda")

# ajuste para FALSE se quiser excluir arquivos com *_old
usar_prioris_old <- FALSE

# Como os dados foram gerados em Organization (1).R:
# beta = 1, eta = 1
beta_true <- 1
eta_true  <- 1

priors_main <- c("FLAT", "IJ", "Partial", "Weak", "Informative")

if (usar_prioris_old) {
  priors_all <- c("FLAT", "IJ", "Partial", "Weak_old", "Informative_old", "Weak", "Informative")
} else {
  priors_all <- priors_main
}

# Saídas
output_dir  <- file.path(base_dir, "DissertationResults")
fig_dir     <- file.path(output_dir, "figures")
table_dir   <- file.path(output_dir, "tables")
csv_dir     <- file.path(output_dir, "csv")
scenario_dir <- file.path(output_dir, "scenario_plots")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(scenario_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2) FUNÇÕES AUXILIARES
# ------------------------------------------------------------

filter_extreme_estimates <- function(df, estimator = c("mean", "median", "map"), mult = 100) {
  estimator <- match.arg(estimator)
  
  est_col <- paste0(estimator, "_est")
  
  df %>%
    mutate(
      true_value_ref = dplyr::case_when(
        parameter == "ypr"   ~ true_ypr,
        parameter == "mu"    ~ true_mu,
        parameter == "sigma" ~ true_sigma,
        parameter == "tpr1"  ~ calc_tpr(true_mu, true_sigma, 0.01),
        parameter == "tpr5"  ~ calc_tpr(true_mu, true_sigma, 0.05),
        parameter == "tpr10" ~ calc_tpr(true_mu, true_sigma, 0.10),
        parameter == "tpr50" ~ calc_tpr(true_mu, true_sigma, 0.50),
        TRUE ~ NA_real_
      )
    ) %>%
    filter(
      is.finite(.data[[est_col]]),
      is.finite(true_value_ref),
      .data[[est_col]] <= mult * true_value_ref
    ) %>%
    select(-true_value_ref)
}

safe_quantile <- function(x, probs = c(0.025, 0.5, 0.975)) {
  if (length(x) == 0 || all(is.na(x))) {
    return(rep(NA_real_, length(probs)))
  }
  as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE))
}

posterior_map <- function(samples) {
  samples <- samples[is.finite(samples)]
  if (length(samples) < 10) return(NA_real_)
  d <- density(samples, na.rm = TRUE)
  d$x[which.max(d$y)]
}

posterior_summary <- function(samples, true_value) {
  q <- safe_quantile(samples, probs = c(0.025, 0.5, 0.975))
  mean_est   <- mean(samples, na.rm = TRUE)
  median_est <- median(samples, na.rm = TRUE)
  map_est    <- posterior_map(samples)
  
  tibble(
    mean_est   = mean_est,
    median_est = median_est,
    map_est    = map_est,
    post_sd    = sd(samples, na.rm = TRUE),
    q025       = q[1],
    q975       = q[3],
    width      = q[3] - q[1],
    cover      = as.integer(!is.na(q[1]) && !is.na(q[3]) && true_value >= q[1] && true_value <= q[3]),
    bias_mean    = mean_est - true_value,
    bias_median  = median_est - true_value,
    bias_map     = map_est - true_value,
    abs_mean     = abs(mean_est - true_value),
    abs_median   = abs(median_est - true_value),
    abs_map      = abs(map_est - true_value),
    sq_mean      = (mean_est - true_value)^2,
    sq_median    = (median_est - true_value)^2,
    sq_map       = (map_est - true_value)^2,
    mad          = median(abs(samples - median_est))
  )
}

load_results_file <- function(file_path) {
  e <- new.env(parent = emptyenv())
  nm <- load(file_path, envir = e)
  
  if ("results" %in% nm) return(e$results)
  
  objs <- mget(nm, envir = e)
  is_candidate <- vapply(objs, function(x) is.list(x), logical(1))
  if (any(is_candidate)) return(objs[[which(is_candidate)[1]]])
  
  stop(sprintf("Não foi possível encontrar o objeto 'results' em %s", file_path))
}

get_true_values <- function(p_repar, beta = 1, eta = 1) {
  # ypr salvo no arquivo é t_pr no espaço original
  true_ypr <- qweibull(p_repar, shape = beta, scale = eta)
  
  # mu verdadeiro no Stan: mu = log(eta)
  true_mu <- log(eta)
  
  # sigma = 1 / beta
  true_sigma <- 1 / beta
  
  list(
    true_ypr = true_ypr,
    true_mu = true_mu,
    true_sigma = true_sigma
  )
}

get_sim_info <- function(List_Dataframes, prop_obs, n_fail_target, sim_id) {
  prp_name <- paste0("ExpectedProportion_", prop_obs * 100)
  exn_name <- paste0("ExpectedN_", n_fail_target)
  sim_name <- paste0("sim_", sim_id)
  
  sim_obj <- List_Dataframes[[prp_name]][[exn_name]][[sim_name]]
  if (is.null(sim_obj)) return(NULL)
  
  info_df <- sim_obj$informations
  dados   <- sim_obj$Dataframe
  
  pfailreal <- info_df$pfailreal[1]
  p_repar   <- pfailreal / 2
  
  truths <- get_true_values(p_repar = p_repar, beta = beta_true, eta = eta_true)
  
  n_fail <- sum(dados$Contagem[dados$Censura == "Failed"])
  n_cens <- sum(dados$Contagem[dados$Censura == "Censored"])
  n_obs  <- n_fail + n_cens
  tc     <- max(dados$Tempo)
  
  tibble(
    sim = sim_id,
    prop_obs = prop_obs,
    prop_obs_label = paste0(round(prop_obs * 100), "% observados"),
    n_fail_target = n_fail_target,
    pfailreal = pfailreal,
    p_repar = p_repar,
    tc = tc,
    n_obs = n_obs,
    n_fail = n_fail,
    n_cens = n_cens,
    seed_data = info_df$seed[1],
    tries = info_df$tries[1],
    true_ypr = truths$true_ypr,
    true_mu = truths$true_mu,
    true_sigma = truths$true_sigma
  )
}

calc_tpr <- function(mu, sigma, p_obs) {
  exp(mu + sigma * log(-log(1 - p_obs)))
}

extract_one_prior <- function(results_obj, prior_name, scenario_info, file_path) {
  if (!prior_name %in% names(results_obj)) return(NULL)
  
  obj <- results_obj[[prior_name]]
  
  needed <- c("ypr", "mu", "sigma")
  if (!all(needed %in% names(obj))) {
    return(NULL)
  }
  
  # mu foi salvo como exp(mu)=eta, então voltamos para o espaço log
  mu_draws_log <- log(obj$mu)
  
  # parâmetros já existentes
  ypr_sum <- posterior_summary(obj$ypr, scenario_info$true_ypr) %>%
    mutate(parameter = "ypr")
  
  mu_sum <- posterior_summary(mu_draws_log, scenario_info$true_mu) %>%
    mutate(parameter = "mu")
  
  sigma_sum <- posterior_summary(obj$sigma, scenario_info$true_sigma) %>%
    mutate(parameter = "sigma")
  
  # novos tpr de interesse
  tpr1_draws  <- calc_tpr(mu_draws_log, obj$sigma, 0.01)
  tpr5_draws  <- calc_tpr(mu_draws_log, obj$sigma, 0.05)
  tpr10_draws <- calc_tpr(mu_draws_log, obj$sigma, 0.10)
  tpr50_draws <- calc_tpr(mu_draws_log, obj$sigma, 0.50)
  
  # valores verdadeiros correspondentes
  true_tpr1  <- calc_tpr(scenario_info$true_mu, scenario_info$true_sigma, 0.01)
  true_tpr5  <- calc_tpr(scenario_info$true_mu, scenario_info$true_sigma, 0.05)
  true_tpr10 <- calc_tpr(scenario_info$true_mu, scenario_info$true_sigma, 0.10)
  true_tpr50 <- calc_tpr(scenario_info$true_mu, scenario_info$true_sigma, 0.50)
  
  tpr1_sum <- posterior_summary(tpr1_draws, true_tpr1) %>%
    mutate(parameter = "tpr1")
  
  tpr5_sum <- posterior_summary(tpr5_draws, true_tpr5) %>%
    mutate(parameter = "tpr5")
  
  tpr10_sum <- posterior_summary(tpr10_draws, true_tpr10) %>%
    mutate(parameter = "tpr10")
  
  tpr50_sum <- posterior_summary(tpr50_draws, true_tpr50) %>%
    mutate(parameter = "tpr50")
  
  bind_rows(
    ypr_sum, mu_sum, sigma_sum,
    tpr1_sum, tpr5_sum, tpr10_sum, tpr50_sum
  ) %>%
    mutate(
      prior = prior_name,
      prior_seed = if ("seed" %in% names(obj)) obj$seed else NA_integer_,
      file_path = file_path
    ) %>%
    bind_cols(scenario_info, .)
}

build_file_index <- function(base_dir) {
  all_files <- list.files(
    base_dir,
    pattern = "^Geral_results_sim[0-9]+\\.rda$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  if (length(all_files) == 0) {
    stop("Nenhum arquivo Geral_results_simX.rda encontrado.")
  }
  
  tibble(
    file_path = all_files,
    file_name = basename(all_files),
    sim = as.integer(str_match(basename(all_files), "Geral_results_sim([0-9]+)\\.rda$")[,2]),
    N_folder = as.integer(str_match(all_files, "N([0-9]+)")[,2]),
    prop_obs = as.numeric(str_match(all_files, "CENS([0-9.]+)")[,2]) / 100
  ) %>%
    arrange(prop_obs, N_folder, sim)
}

summarise_metrics <- function(df, estimator = c("mean", "median", "map")) {
  estimator <- match.arg(estimator)
  
  bias_col <- paste0("bias_", estimator)
  abs_col  <- paste0("abs_", estimator)
  sq_col   <- paste0("sq_", estimator)
  est_col  <- paste0(estimator, "_est")
  
  df %>%
    group_by(parameter, prior, prop_obs, prop_obs_label, n_fail_target) %>%
    summarise(
      n_sim = n(),
      true_value = mean(dplyr::case_when(
        parameter == "ypr" ~ true_ypr,
        parameter == "mu" ~ true_mu,
        parameter == "sigma" ~ true_sigma
      ), na.rm = TRUE),
      bias = mean(.data[[bias_col]], na.rm = TRUE),
      rmse = sqrt(mean(.data[[sq_col]], na.rm = TRUE)),
      mae = mean(.data[[abs_col]], na.rm = TRUE),
      coverage = mean(cover, na.rm = TRUE),
      width = mean(width, na.rm = TRUE),
      estimate_mean = mean(.data[[est_col]], na.rm = TRUE),
      estimate_sd = sd(.data[[est_col]], na.rm = TRUE),
      post_sd_mean = mean(post_sd, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(estimator = estimator)
}

save_plot <- function(plot_obj, filename, width = 10, height = 6) {
  ggsave(filename = filename, plot = plot_obj, width = width, height = height, dpi = 300)
}

# ------------------------------------------------------------
# 3) CARREGAR BASE DAS SIMULAÇÕES
# ------------------------------------------------------------
List_Dataframes <- readRDS(data_file)

# ------------------------------------------------------------
# 4) INDEXAR ARQUIVOS
# ------------------------------------------------------------
file_index <- build_file_index(base_dir)

write_csv(file_index, file.path(csv_dir, "file_index.csv"))

# ------------------------------------------------------------
# 5) EXTRAIR RESULTADOS POR SIMULAÇÃO
# ------------------------------------------------------------
all_results <- vector("list", nrow(file_index))
pb <- txtProgressBar(min = 0, max = nrow(file_index), style = 3)

for (i in seq_len(nrow(file_index))) {
  file_info <- file_index[i, ]
  
  scenario_info <- get_sim_info(
    List_Dataframes = List_Dataframes,
    prop_obs = file_info$prop_obs,
    n_fail_target = file_info$N_folder,
    sim_id = file_info$sim
  )
  
  if (is.null(scenario_info)) {
    all_results[[i]] <- NULL
    setTxtProgressBar(pb, i)
    next
  }
  
  results_obj <- tryCatch(
    load_results_file(file_info$file_path),
    error = function(e) NULL
  )
  
  if (is.null(results_obj)) {
    all_results[[i]] <- NULL
    setTxtProgressBar(pb, i)
    next
  }
  
  extracted <- map_dfr(priors_all, function(pr) {
    extract_one_prior(
      results_obj = results_obj,
      prior_name = pr,
      scenario_info = scenario_info,
      file_path = file_info$file_path
    )
  })
  
  all_results[[i]] <- extracted
  setTxtProgressBar(pb, i)
}
close(pb)

results_by_sim <- bind_rows(all_results)

# ------------------------------------------------------------
# 6) SALVAR BASE BRUTA RESUMIDA
# ------------------------------------------------------------
write_csv(results_by_sim, file.path(csv_dir, "results_by_simulation_summary2.csv"))

