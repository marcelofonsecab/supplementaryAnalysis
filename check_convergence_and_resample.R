# ============================================================
# REPARO DE CENÁRIOS EXTREMOS
# - identifica cenários problemáticos na FLAT
# - recria a base com novas seeds
# - tenta rodar FLAT primeiro
# - se FLAT rodar e a mediana ficar no intervalo, roda as demais prioris
# - salva base, resultados e logs
# - pula cenários já processados
# ============================================================

#rm(list = ls())

library(dplyr)
library(purrr)
library(readr)
library(tibble)
library(rstan)

# ============================================================
# 0) PRÉ-REQUISITOS
# ============================================================
# Este script supõe que já existam no ambiente ou carregados antes:
# - criar_base_simulada()
# - rodar_uma_simulacao()
# - analise_bayesiana_confiabilidade()
#
# Exemplo:
# source("funcoes_confiabilidade.R")
# source("Organization (1).R")
#
# Também supõe que exista o arquivo:
# - results_by_simulation_summary2.csv

# ============================================================
# 1) FUNÇÕES AUXILIARES
# ============================================================

calc_tpr <- function(mu, sigma, p_obs) {
  exp(mu + sigma * log(-log(1 - p_obs)))
}

# ------------------------------------------------------------
# resumo robusto do tpr50 usando apenas mu e sigma por nome
# ------------------------------------------------------------
resumo_tpr50_flat <- function(fit_flat, p_obs = 0.50) {
  cat("    -> entrando em resumo_tpr50_flat()\n")
  flush.console()
  
  post <- as.matrix(fit_flat, pars = c("mu", "sigma"))
  
  cat("    -> extração concluída; n_draws =", nrow(post), "\n")
  flush.console()
  
  if (!all(c("mu", "sigma") %in% colnames(post))) {
    stop("Os parâmetros 'mu' e 'sigma' não foram encontrados no objeto Stan.")
  }
  
  mu_draws_log <- post[, "mu"]
  sigma_draws  <- post[, "sigma"]
  
  tpr50_draws <- calc_tpr(mu_draws_log, sigma_draws, p_obs)
  tpr50_draws <- tpr50_draws[is.finite(tpr50_draws)]
  
  cat("    -> draws finitos =", length(tpr50_draws), "\n")
  flush.console()
  
  if (length(tpr50_draws) == 0) {
    return(tibble(
      mean_est   = NA_real_,
      median_est = NA_real_,
      sd_est     = NA_real_,
      q025       = NA_real_,
      q975       = NA_real_
    ))
  }
  
  tibble(
    mean_est   = mean(tpr50_draws),
    median_est = median(tpr50_draws),
    sd_est     = sd(tpr50_draws),
    q025       = unname(quantile(tpr50_draws, 0.025)),
    q975       = unname(quantile(tpr50_draws, 0.975))
  )
}

# ------------------------------------------------------------
# identifica cenários problemáticos
# ------------------------------------------------------------
identificar_cenarios_extremos <- function(results_by_sim,
                                          prior_name = "FLAT",
                                          parameter_name = "tpr50",
                                          lim_inf = -5,
                                          lim_sup = 5) {
  results_by_sim %>%
    filter(prior == prior_name, parameter == parameter_name) %>%
    filter(is.finite(median_est)) %>%
    filter(median_est <= lim_inf | median_est >= lim_sup) %>%
    distinct(prop_obs, n_fail_target, sim, .keep_all = TRUE) %>%
    mutate(
      prop_obs = as.numeric(prop_obs),
      n_fail_target = as.integer(n_fail_target),
      sim = as.integer(sim)
    ) %>%
    arrange(prop_obs, n_fail_target, sim)
}

# ------------------------------------------------------------
# gerador seguro de seeds
# ------------------------------------------------------------
gerar_novas_seeds <- function(prop_obs, n_fail_target, sim, tentativa,
                              n_priors = 5,
                              master_seed = 987654) {
  comp1 <- as.numeric(round(prop_obs * 1e6))
  comp2 <- as.numeric(n_fail_target)
  comp3 <- as.numeric(sim)
  comp4 <- as.numeric(tentativa)
  comp5 <- as.numeric(master_seed)
  
  seed_base <- (comp1 * 131 + comp2 * 137 + comp3 * 139 + comp4 * 149 + comp5) %% .Machine$integer.max
  seed_base <- as.integer(seed_base)
  
  if (is.na(seed_base) || seed_base <= 0) {
    seed_base <- 123456L
  }
  
  set.seed(seed_base)
  
  seed_dados <- sample.int(.Machine$integer.max - 1L, 1)
  seed_priors <- sample.int(.Machine$integer.max - 1L, n_priors)
  
  nomes_priors <- c("FLAT", "IJ", "PARTIAL", "WEAK", "INF")
  names(seed_priors) <- nomes_priors[seq_len(n_priors)]
  
  list(
    seed_base = seed_base,
    seed_dados = seed_dados,
    seed_priors = seed_priors
  )
}

# ------------------------------------------------------------
# checa reparos já existentes
# ------------------------------------------------------------
checar_reparos_existentes <- function(dir_saida = "Resultados_reparados") {
  pasta_resultados <- file.path(dir_saida, "resultados")
  
  if (!dir.exists(pasta_resultados)) {
    return(tibble(
      prop_obs = numeric(),
      n_fail_target = integer(),
      sim = integer()
    ))
  }
  
  arquivos <- list.files(
    pasta_resultados,
    pattern = "^Geral_results_sim[0-9]+_reparado_prop_.*_N_[0-9]+\\.rda$",
    full.names = FALSE
  )
  
  if (length(arquivos) == 0) {
    return(tibble(
      prop_obs = numeric(),
      n_fail_target = integer(),
      sim = integer()
    ))
  }
  
  tibble(arquivo = arquivos) %>%
    mutate(
      sim = as.integer(sub("^Geral_results_sim([0-9]+)_reparado_prop_.*$", "\\1", arquivo)),
      prop_obs = as.numeric(sub("^Geral_results_sim[0-9]+_reparado_prop_([^_]+)_N_.*$", "\\1", arquivo)),
      n_fail_target = as.integer(sub("^.*_N_([0-9]+)\\.rda$", "\\1", arquivo))
    ) %>%
    distinct(prop_obs, n_fail_target, sim) %>%
    mutate(
      prop_obs = as.numeric(prop_obs),
      n_fail_target = as.integer(n_fail_target),
      sim = as.integer(sim)
    )
}

# ============================================================
# 2) AJUSTE SOMENTE DA FLAT
# ============================================================

rodar_flat_apenas <- function(data_,
                              seed_flat,
                              stan_model_path = "testestan.stan",
                              chains = 1,
                              iter = 1000,
                              warmup = 500,
                              refresh = 50) {
  p_repar <- data_$informations$pfailreal / 2
  Dados__ <- data_$Dataframe
  tc <- max(Dados__$Tempo)
  
  cat("  -> Entrou em rodar_flat_apenas()\n")
  cat("  -> seed_flat:", seed_flat, "\n")
  cat("  -> p_repar:", p_repar, "\n")
  cat("  -> tc:", tc, "\n")
  cat("  -> n_linhas_agregadas:", nrow(Dados__), "\n")
  flush.console()
  
  set.seed(seed_flat)
  
  fit_flat <- analise_bayesiana_confiabilidade(
    times = Dados__,
    distribution = "weibull",
    p_r = p_repar,
    prior_tpr = list(type = "flat"),
    prior_beta_sigma = list(type = "flat"),
    t_c = tc,
    stan_model_path = stan_model_path,
    chains = chains,
    iter = iter,
    warmup = warmup,
    refresh = refresh
  )
  
  cat("  -> Saiu do stan() da FLAT\n")
  flush.console()
  
  fit_flat
}

# ============================================================
# 3) REFAZER UM CENÁRIO ATÉ ACEITAR
# ============================================================

refazer_cenario_ate_aceitar <- function(prop_obs,
                                        n_fail_target,
                                        sim,
                                        beta = 1,
                                        eta = 1,
                                        lim_inf = -5,
                                        lim_sup = 5,
                                        max_tentativas = 30,
                                        stan_model_path = "testestan.stan",
                                        master_seed = 987654,
                                        chains_flat = 1,
                                        iter_flat = 1000,
                                        warmup_flat = 500,
                                        refresh_flat = 50) {
  historico <- vector("list", max_tentativas)
  
  for (tentativa in seq_len(max_tentativas)) {
    cat("Tentativa:", tentativa, "\n")
    flush.console()
    
    novas_seeds <- gerar_novas_seeds(
      prop_obs = prop_obs,
      n_fail_target = n_fail_target,
      sim = sim,
      tentativa = tentativa,
      master_seed = master_seed
    )
    
    cat("  -> seed_dados:", novas_seeds$seed_dados, "\n")
    cat("  -> criando base simulada...\n")
    flush.console()
    
    base_simulada <- tryCatch(
      criar_base_simulada(
        N_esperado = n_fail_target,
        p_falha = prop_obs,
        beta = beta,
        eta = eta,
        seed_dados = novas_seeds$seed_dados,
        minimo_falhas = 3
      ),
      error = function(e) e
    )
    
    if (inherits(base_simulada, "error")) {
      cat("  -> erro na criação da base:", conditionMessage(base_simulada), "\n")
      flush.console()
      
      historico[[tentativa]] <- tibble(
        tentativa = tentativa,
        seed_base = novas_seeds$seed_base,
        seed_dados = novas_seeds$seed_dados,
        seed_flat = novas_seeds$seed_priors["FLAT"],
        aprovado = FALSE,
        median_tpr50_flat = NA_real_,
        erro = paste("erro base:", conditionMessage(base_simulada))
      )
      next
    }
    
    cat("  -> base criada\n")
    flush.console()
    
    base_simulada$informations_priorseeds <- novas_seeds$seed_priors
    
    cat("  -> rodando FLAT...\n")
    flush.console()
    
    fit_flat <- tryCatch(
      rodar_flat_apenas(
        data_ = base_simulada,
        seed_flat = novas_seeds$seed_priors["FLAT"],
        stan_model_path = stan_model_path,
        chains = chains_flat,
        iter = iter_flat,
        warmup = warmup_flat,
        refresh = refresh_flat
      ),
      error = function(e) e
    )
    
    cat("  -> terminou tryCatch da FLAT\n")
    flush.console()
    
    if (inherits(fit_flat, "error")) {
      cat("  -> erro na FLAT:", conditionMessage(fit_flat), "\n")
      flush.console()
      
      historico[[tentativa]] <- tibble(
        tentativa = tentativa,
        seed_base = novas_seeds$seed_base,
        seed_dados = novas_seeds$seed_dados,
        seed_flat = novas_seeds$seed_priors["FLAT"],
        aprovado = FALSE,
        median_tpr50_flat = NA_real_,
        erro = conditionMessage(fit_flat)
      )
      next
    }
    
    cat("  -> resumindo tpr50...\n")
    flush.console()
    
    resumo_flat <- tryCatch(
      resumo_tpr50_flat(fit_flat, p_obs = 0.50),
      error = function(e) e
    )
    
    if (inherits(resumo_flat, "error")) {
      cat("  -> erro no resumo da FLAT:", conditionMessage(resumo_flat), "\n")
      flush.console()
      
      historico[[tentativa]] <- tibble(
        tentativa = tentativa,
        seed_base = novas_seeds$seed_base,
        seed_dados = novas_seeds$seed_dados,
        seed_flat = novas_seeds$seed_priors["FLAT"],
        aprovado = FALSE,
        median_tpr50_flat = NA_real_,
        erro = paste("erro resumo:", conditionMessage(resumo_flat))
      )
      next
    }
    
    mediana_flat <- resumo_flat$median_est[[1]]
    
    aprovado <- is.finite(mediana_flat) &&
      mediana_flat > lim_inf &&
      mediana_flat < lim_sup
    
    cat("  -> mediana FLAT =", mediana_flat, "\n")
    cat("  -> aprovado =", aprovado, "\n")
    flush.console()
    
    historico[[tentativa]] <- tibble(
      tentativa = tentativa,
      seed_base = novas_seeds$seed_base,
      seed_dados = novas_seeds$seed_dados,
      seed_flat = novas_seeds$seed_priors["FLAT"],
      aprovado = aprovado,
      median_tpr50_flat = mediana_flat,
      erro = NA_character_
    )
    
    if (aprovado) {
      cat("  -> cenário aprovado; rodando todas as prioris...\n")
      flush.console()
      
      resultados_finais <- tryCatch(
        rodar_uma_simulacao(
          data_ = base_simulada,
          prior_seeds = novas_seeds$seed_priors,
          stan_model_path = stan_model_path
        ),
        error = function(e) e
      )
      
      if (inherits(resultados_finais, "error")) {
        cat("  -> erro ao rodar todas as prioris:", conditionMessage(resultados_finais), "\n")
        flush.console()
        
        return(list(
          aceito = FALSE,
          prop_obs = prop_obs,
          n_fail_target = n_fail_target,
          sim = sim,
          tentativa_aceita = tentativa,
          base_simulada = base_simulada,
          resultados = NULL,
          historico = bind_rows(historico[seq_len(tentativa)]),
          resumo_flat = resumo_flat,
          erro_final = conditionMessage(resultados_finais)
        ))
      }
      
      return(list(
        aceito = TRUE,
        prop_obs = prop_obs,
        n_fail_target = n_fail_target,
        sim = sim,
        tentativa_aceita = tentativa,
        base_simulada = base_simulada,
        resultados = resultados_finais,
        historico = bind_rows(historico[seq_len(tentativa)]),
        resumo_flat = resumo_flat,
        erro_final = NA_character_
      ))
    }
  }
  
  list(
    aceito = FALSE,
    prop_obs = prop_obs,
    n_fail_target = n_fail_target,
    sim = sim,
    tentativa_aceita = NA_integer_,
    base_simulada = NULL,
    resultados = NULL,
    historico = bind_rows(historico),
    resumo_flat = NULL,
    erro_final = "Nenhuma tentativa aprovada."
  )
}

# ============================================================
# 4) RODAR O REPARO PARA TODOS
# ============================================================

reparar_cenarios_extremos <- function(cenarios_ruins,
                                      beta = 1,
                                      eta = 1,
                                      lim_inf = -5,
                                      lim_sup = 5,
                                      max_tentativas = 30,
                                      stan_model_path = "testestan.stan",
                                      master_seed = 987654,
                                      dir_saida = "Resultados_reparados",
                                      pular_existentes = TRUE,
                                      chains_flat = 1,
                                      iter_flat = 1000,
                                      warmup_flat = 500,
                                      refresh_flat = 50) {
  dir.create(dir_saida, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dir_saida, "bases"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dir_saida, "resultados"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dir_saida, "logs"), recursive = TRUE, showWarnings = FALSE)
  
  cenarios_ruins <- cenarios_ruins %>%
    mutate(
      prop_obs = as.numeric(prop_obs),
      n_fail_target = as.integer(n_fail_target),
      sim = as.integer(sim)
    )
  
  if (pular_existentes) {
    reparos_existentes <- checar_reparos_existentes(dir_saida)
    
    cenarios_ruins <- cenarios_ruins %>%
      anti_join(
        reparos_existentes,
        by = c("prop_obs", "n_fail_target", "sim")
      )
    
    cat("Cenários pendentes após checagem:", nrow(cenarios_ruins), "\n")
    flush.console()
  }
  
  if (nrow(cenarios_ruins) == 0) {
    cat("Nenhum cenário pendente para reparar.\n")
    flush.console()
    return(invisible(list()))
  }
  
  saida <- vector("list", nrow(cenarios_ruins))
  
  for (i in seq_len(nrow(cenarios_ruins))) {
    linha <- cenarios_ruins[i, ]
    
    cat("\n----------------------------------------\n")
    cat("Reparando cenário:", i, "de", nrow(cenarios_ruins), "\n")
    cat("prop_obs =", linha$prop_obs,
        "| n_fail_target =", linha$n_fail_target,
        "| sim =", linha$sim, "\n")
    flush.console()
    
    obj <- refazer_cenario_ate_aceitar(
      prop_obs = linha$prop_obs,
      n_fail_target = linha$n_fail_target,
      sim = linha$sim,
      beta = beta,
      eta = eta,
      lim_inf = lim_inf,
      lim_sup = lim_sup,
      max_tentativas = max_tentativas,
      stan_model_path = stan_model_path,
      master_seed = master_seed,
      chains_flat = chains_flat,
      iter_flat = iter_flat,
      warmup_flat = warmup_flat,
      refresh_flat = refresh_flat
    )
    
    saida[[i]] <- obj
    
    hist_file <- file.path(
      dir_saida, "logs",
      paste0("log_reparo_prop_", linha$prop_obs,
             "_N_", linha$n_fail_target,
             "_sim_", linha$sim, ".csv")
    )
    
    write_csv(obj$historico, hist_file)
    
    if (isTRUE(obj$aceito) && !is.null(obj$base_simulada) && !is.null(obj$resultados)) {
      saveRDS(
        obj$base_simulada,
        file = file.path(
          dir_saida, "bases",
          paste0("base_reparada_prop_", linha$prop_obs,
                 "_N_", linha$n_fail_target,
                 "_sim_", linha$sim, ".rds")
        )
      )
      
      results <- obj$resultados
      
      save(
        results,
        file = file.path(
          dir_saida, "resultados",
          paste0("Geral_results_sim", linha$sim,
                 "_reparado_prop_", linha$prop_obs,
                 "_N_", linha$n_fail_target, ".rda")
        )
      )
      
      cat("  -> base e resultados salvos.\n")
      flush.console()
    } else {
      cat("  -> cenário não aceito ou sem resultados; nada salvo em resultados/.\n")
      flush.console()
    }
  }
  
  resumo_final <- bind_rows(lapply(saida, function(x) {
    tibble(
      prop_obs = x$prop_obs,
      n_fail_target = x$n_fail_target,
      sim = x$sim,
      aceito = x$aceito,
      tentativa_aceita = x$tentativa_aceita,
      erro_final = x$erro_final
    )
  }))
  
  write_csv(resumo_final, file.path(dir_saida, "resumo_reparos.csv"))
  
  invisible(saida)
}

# ============================================================
# 5) EXECUÇÃO
# ============================================================

results_by_sim <- read_csv("results_by_simulation_summary2.csv", show_col_types = FALSE)

cenarios_ruins <- identificar_cenarios_extremos(
  results_by_sim = results_by_sim,
  prior_name = "FLAT",
  parameter_name = "tpr50",
  lim_inf = -5,
  lim_sup = 5
)

reparos <- reparar_cenarios_extremos(
  cenarios_ruins = cenarios_ruins,
  beta = 1,
  eta = 1,
  lim_inf = -5,
  lim_sup = 5,
  max_tentativas = 30,
  stan_model_path = "testestan.stan",
  master_seed = 20260422,
  dir_saida = "Resultados_reparados",
  pular_existentes = TRUE,
  chains_flat = 1,
  iter_flat = 1000,
  warmup_flat = 500,
  refresh_flat = 50
)


