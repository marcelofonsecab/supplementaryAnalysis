# =========================================================
# SCRIPT DE SIMULACAO E AJUSTE BAYESIANO
# =========================================================
# Este script:
# 1. carrega as funções auxiliares;
# 2. gera bases simuladas com censura à direita;
# 3. utiliza seeds reprodutíveis para os dados e para as prioris;
# 4. ajusta o modelo Bayesiano sob cinco especificações a priori:
#    - FLAT
#    - IJ
#    - Partial
#    - Weak
#    - Informative
#
# Observação:
# Este script supõe que o arquivo "funcoes_confiabilidade.R"
# esteja disponível no diretório de trabalho e possa ser
# carregado via source().
# =========================================================


setwd("C:\\Bayes\\Projeto500\\Reworked\\Reworkdnewprior")

# =========================================================
# CARREGAR FUNCOES AUXILIARES
# =========================================================
source("funcoes_confiabilidade.R")


# =========================================================
# PARAMETROS GERAIS DO ESTUDO
# =========================================================
Er <- c(10, 25, 35, 50, 75, 100)
pfail <- c(0.01, 0.05, 0.10, 0.50)

beta_true <- 1
eta_true <- 1

n_sim <- 500


# =========================================================
# SEEDS UTILIZADOS NO ESTUDO
# =========================================================
# seedz:
#   array para geração das bases de dados
#   dimensões: simulacao x N esperado x proporcao de falha
#
# seedz_priors:
#   array para seeds das cinco prioris
#   dimensões: priori x simulacao x N esperado x proporcao de falha
#   ordem das prioris:
#   1 = FLAT
#   2 = IJ
#   3 = PARTIAL
#   4 = WEAK
#   5 = INFORMATIVE

set.seed(2026)
seedz <- sample(1:1e7, size = n_sim * length(Er) * length(pfail))
seedz <- array(
  seedz,
  dim = c(n_sim, length(Er), length(pfail))
)

set.seed(6202)
rnms <- c("FLAT", "IJ", "PARTIAL", "WEAK", "INF")
seedz_priors <- sample(1:1e9, size = 5 * n_sim * length(Er) * length(pfail))
seedz_priors <- array(
  seedz_priors,
  dim = c(5, n_sim, length(Er), length(pfail)),
  dimnames = list(rnms, NULL, NULL, NULL)
)


# =========================================================
# FUNCAO PARA GERAR DADOS CENSURADOS
# =========================================================
gerar_dados_weibull_censurados <- function(n,
                                           beta = 1,
                                           eta = 1,
                                           p_falha = 0.10,
                                           minimo_falhas = 3,
                                           seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  soma_falhas <- 0
  tentativas <- 0
  
  while (soma_falhas < minimo_falhas) {
    tentativas <- tentativas + 1
    
    tempos <- rweibull(n, shape = beta, scale = eta)
    
    tc <- exp(sevinf(
      x = p_falha,
      mu = log(eta),
      sigma = 1 / beta
    ))
    
    soma_falhas <- sum(tempos < tc)
  }
  
  status <- ifelse(tempos < tc, 1, 0)
  tempos[tempos > tc] <- tc
  
  list(
    tempos = tempos,
    status = status,
    tc = tc,
    n_falhas = sum(status == 1),
    proporcao_falhas_observada = mean(status),
    tentativas = tentativas
  )
}


# =========================================================
# FUNCAO PARA AGREGAR OS DADOS
# =========================================================
agregar_dados_simulados <- function(tempos, status) {
  tabela_falhas <- table(tempos[status == 1])
  tabela_cens   <- table(tempos[status == 0])
  
  df_falhas <- data.frame(
    Tempo = as.numeric(names(tabela_falhas)),
    Contagem = as.numeric(tabela_falhas),
    Censura = "Failed",
    stringsAsFactors = FALSE
  )
  
  df_cens <- data.frame(
    Tempo = as.numeric(names(tabela_cens)),
    Contagem = as.numeric(tabela_cens),
    Censura = "Censored",
    stringsAsFactors = FALSE
  )
  
  dados_agregados <- rbind(df_falhas, df_cens)
  dados_agregados <- dados_agregados[order(dados_agregados$Tempo), ]
  rownames(dados_agregados) <- NULL
  
  dados_agregados
}


# =========================================================
# FUNCAO PARA CRIAR UMA BASE SIMULADA
# =========================================================
criar_base_simulada <- function(N_esperado,
                                p_falha,
                                beta = 1,
                                eta = 1,
                                seed_dados = NULL,
                                minimo_falhas = 3) {
  n_planejado <- N_esperado / p_falha
  
  dados <- gerar_dados_weibull_censurados(
    n = n_planejado,
    beta = beta,
    eta = eta,
    p_falha = p_falha,
    minimo_falhas = minimo_falhas,
    seed = seed_dados
  )
  
  dados_agregados <- agregar_dados_simulados(
    tempos = dados$tempos,
    status = dados$status
  )
  
  list(
    informations = data.frame(
      seed = seed_dados,
      tries = dados$tentativas,
      pfailreal = dados$n_falhas / n_planejado
    ),
    Dataframe = dados_agregados
  )
}


# =========================================================
# FUNCAO AUXILIAR PARA EXTRAIR RESULTADOS DO STAN
# =========================================================
extrair_resultados_stan <- function(fit, seed_prior) {
  res_extr <- rstan::extract(fit, permuted = FALSE)
  
  list(
    ypr = c(exp(res_extr[, , 1])),
    sigma = c(res_extr[, , 4]),
    mu = c(exp(res_extr[, , 6])),
    seed = seed_prior
  )
}


# =========================================================
# FUNCAO PARA AJUSTAR AS CINCO PRIORIS
# =========================================================
rodar_uma_simulacao <- function(data_,
                                prior_seeds,
                                stan_model_path = "testestan.stan") {
  p_repar <- data_$informations$pfailreal / 2
  Dados__ <- data_$Dataframe
  tc <- max(Dados__$Tempo)
  
  # -------------------------------------------------------
  # FLAT
  # -------------------------------------------------------
  set.seed(prior_seeds[1])
  
  fit_flat <- analise_bayesiana_confiabilidade(
    times = Dados__,
    distribution = "weibull",
    p_r = p_repar,
    prior_tpr = list(type = "flat"),
    prior_beta_sigma = list(type = "flat"),
    t_c = tc,
    stan_model_path = stan_model_path,
    chains = 4,
    iter = 6000,
    warmup = 1000
  )
  
  results_flat <- extrair_resultados_stan(fit_flat, prior_seeds[1])
  
  # -------------------------------------------------------
  # IJ
  # -------------------------------------------------------
  set.seed(prior_seeds[2])
  
  fit_ij <- analise_bayesiana_confiabilidade(
    times = Dados__,
    distribution = "weibull",
    p_r = p_repar,
    prior_tpr = list(type = "cj"),
    prior_beta_sigma = list(type = "cj"),
    t_c = tc,
    stan_model_path = stan_model_path,
    chains = 4,
    iter = 6000,
    warmup = 1000
  )
  
  results_ij <- extrair_resultados_stan(fit_ij, prior_seeds[2])
  
  # -------------------------------------------------------
  # PARTIAL
  # -------------------------------------------------------
  set.seed(prior_seeds[3])
  
  fit_partial <- analise_bayesiana_confiabilidade(
    times = Dados__,
    distribution = "weibull",
    p_r = p_repar,
    prior_tpr = list(type = "cj"),
    prior_beta_sigma = list(type = "truncatednormal", range = c(0.3, 1.7)),
    t_c = tc,
    stan_model_path = stan_model_path,
    chains = 4,
    iter = 6000,
    warmup = 1000
  )
  
  results_partial <- extrair_resultados_stan(fit_partial, prior_seeds[3])
  
  # -------------------------------------------------------
  # WEAK
  # -------------------------------------------------------
  tprrange_weak <- c(
    qexp(p_repar * (1 - 0.40)),
    qexp(p_repar * (1 + 0.85))
  )
  
  set.seed(prior_seeds[4])
  
  fit_weak <- analise_bayesiana_confiabilidade(
    times = Dados__,
    distribution = "weibull",
    p_r = p_repar,
    prior_tpr = list(type = "lognormal", range = tprrange_weak),
    prior_beta_sigma = list(type = "lognormal", range = c(0.3, 1.7)),
    t_c = tc,
    stan_model_path = stan_model_path,
    chains = 4,
    iter = 6000,
    warmup = 1000
  )
  
  results_weak <- extrair_resultados_stan(fit_weak, prior_seeds[4])
  
  # -------------------------------------------------------
  # INFORMATIVE
  # -------------------------------------------------------
  tprrange_info <- c(
    qexp(p_repar * (1 - 0.25)),
    qexp(p_repar * (1 + 0.30))
  )
  
  set.seed(prior_seeds[5])
  
  fit_info <- analise_bayesiana_confiabilidade(
    times = Dados__,
    distribution = "weibull",
    p_r = p_repar,
    prior_tpr = list(type = "lognormal", range = tprrange_info),
    prior_beta_sigma = list(type = "lognormal", range = c(0.8, 1.2)),
    t_c = tc,
    stan_model_path = stan_model_path,
    chains = 4,
    iter = 6000,
    warmup = 1000
  )
  
  results_info <- extrair_resultados_stan(fit_info, prior_seeds[5])
  
  list(
    FLAT = results_flat,
    IJ = results_ij,
    Partial = results_partial,
    Weak = results_weak,
    Informative = results_info
  )
}


# =========================================================
# FUNCAO PARA RECONSTRUIR TODA A ESTRUTURA DAS BASES
# =========================================================
criar_lista_bases <- function(Er,
                              pfail,
                              n_sim,
                              beta = 1,
                              eta = 1,
                              seedz,
                              seedz_priors) {
  Data_List <- list()
  
  for (cs in seq_along(pfail)) {
    p_tmp <- pfail[cs]
    
    nome_prop <- paste0("ExpectedProportion_", p_tmp * 100)
    Data_List[[nome_prop]] <- list()
    
    for (Ni in seq_along(Er)) {
      N_exp <- Er[Ni]
      
      nome_N <- paste0("ExpectedN_", N_exp)
      Data_List[[nome_prop]][[nome_N]] <- list()
      
      for (m in 1:n_sim) {
        seed_data_choose <- seedz[m, Ni, cs]
        seeds_priors_choose <- seedz_priors[, m, Ni, cs]
        
        base_simulada <- criar_base_simulada(
          N_esperado = N_exp,
          p_falha = p_tmp,
          beta = beta,
          eta = eta,
          seed_dados = seed_data_choose,
          minimo_falhas = 3
        )
        
        base_simulada$informations_priorseeds <- seeds_priors_choose
        
        Data_List[[nome_prop]][[nome_N]][[paste0("sim_", m)]] <- base_simulada
      }
    }
  }
  
  Data_List
}

