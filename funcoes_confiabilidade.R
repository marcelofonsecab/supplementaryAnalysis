# =========================================================
# FUNCOES AUXILIARES PARA ANALISE BAYESIANA DE CONFIABILIDADE
# =========================================================
# Este arquivo contem as funcoes auxiliares utilizadas no
# ajuste Bayesiano dos modelos de confiabilidade.
#
# Para utilizacao, basta carregar com:
# source("funcoes_confiabilidade.R")
# =========================================================


# =========================================================
# PACOTES
# =========================================================
library(RSplida)
library(Hmisc)
library(rstan)
library(lsinf)

options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


# =========================================================
# CARREGAMENTO DA ROTINA COMPILADA DO PACOTE lsinf
# =========================================================
dllpath <- system.file("libs", "x64", "lsinf.dll", package = "lsinf")

if (file.exists(dllpath)) {
  dyn.load(dllpath)
}


# =========================================================
# FUNCAO PARA OBTER ELEMENTOS DA INFORMACAO DE FISHER
# =========================================================
outro <- function(z, censor.type, distribution) {
  if (!is.numeric(z)) {
    stop("z must be numeric")
  }
  
  if (!is.character(censor.type)) {
    stop("censor.type must be character")
  }
  
  if (!is.character(distribution)) {
    stop("distribution must be character string")
  }
  
  switch(
    tolower(distribution),
    weibull     = ,
    sev         = idist <- 1,
    frechet     = ,
    lev         = idist <- 2,
    lognormal   = ,
    normal      = idist <- 3,
    loglogistic = ,
    logistic    = idist <- 4,
    stop("Distribution must be sev, lev, normal, or logistic")
  )
  
  switch(
    censor.type,
    uncensored = icode <- 1,
    right      = icode <- 2,
    left       = icode <- 3,
    stop("censor.type must be uncensored, left, or right")
  )
  
  nrows <- length(z)
  
  zout <- .Fortran(
    "slsinf",
    as.integer(idist),
    as.integer(icode),
    as.double(z),
    as.double(z),
    f11 = double(nrows),
    f12 = double(nrows),
    f22 = double(nrows),
    as.integer(nrows),
    ifault = integer(1),
    irow = integer(1)
  )
  
  if (zout$ifault > 1) {
    warning("Evaluation error in Fortran routine.")
  }
  
  list(
    f11 = zout$f11,
    f12 = zout$f12,
    f22 = zout$f22
  )
}


# =========================================================
# FUNCAO AUXILIAR PARA VETORIZAR A FUNCAO 'outro'
# =========================================================
outro_vec <- function(v, censor.type, distribution) {
  res_lista <- lapply(v, outro,
                      censor.type = censor.type,
                      distribution = distribution)
  
  list(
    f11 = sapply(res_lista, `[[`, "f11"),
    f12 = sapply(res_lista, `[[`, "f12"),
    f22 = sapply(res_lista, `[[`, "f22")
  )
}


# =========================================================
# FUNCAO PARA OBTER PARAMETROS DE NORMAL TRUNCADA
# =========================================================
get_tnorm_params <- function(a, b, pL, pU) {
  if (!(a > 0 && b > 0 && a < b)) {
    stop("É necessário que 0 < a < b.")
  }
  
  if (!(pL >= 0 && pU <= 1 && pL < pU)) {
    stop("É necessário que 0 <= pL < pU <= 1.")
  }
  
  qL <- 1 - pL
  qU <- 1 - pU
  
  eq_gamma <- function(gamma) {
    eps <- 1e-12
    
    pL_star <- pL + qL * pnorm(-gamma)
    pU_star <- pU + qU * pnorm(-gamma)
    
    pL_star <- pmin(pmax(pL_star, eps), 1 - eps)
    pU_star <- pmin(pmax(pU_star, eps), 1 - eps)
    
    zL <- qnorm(pL_star)
    zU <- qnorm(pU_star)
    
    left  <- (gamma + zL) / a
    right <- (gamma + zU) / b
    
    left - right
  }
  
  grid <- seq(-15, 15, length.out = 500)
  vals <- sapply(grid, eq_gamma)
  ok <- is.finite(vals)
  
  grid <- grid[ok]
  vals <- vals[ok]
  
  idx <- which(diff(sign(vals)) != 0)
  
  if (length(idx) == 0) {
    stop("Não foi encontrada raiz para gamma.")
  }
  
  g1 <- grid[idx[1]]
  g2 <- grid[idx[1] + 1]
  
  gamma_sol <- uniroot(eq_gamma, interval = c(g1, g2))$root
  
  eps <- 1e-12
  pL_star <- pL + qL * pnorm(-gamma_sol)
  pL_star <- pmin(pmax(pL_star, eps), 1 - eps)
  
  zL <- qnorm(pL_star)
  sigma <- a / (gamma_sol + zL)
  mu <- gamma_sol * sigma
  
  list(
    mu = mu,
    sigma = sigma,
    gamma = gamma_sol
  )
}


# =========================================================
# FUNCAO PARA CONVERTER INTERVALOS EM PARAMETROS DE PRIORI
# =========================================================
get_params_from_range <- function(range, type = "lognormal") {
  lower <- range[1]
  upper <- range[2]
  z_995 <- qnorm(0.995)
  
  if (type == "lognormal") {
    log_mean <- (log(lower) + log(upper)) / 2
    log_sd   <- (log(upper) - log(lower)) / (2 * z_995)
    return(c(log_mean, log_sd))
  }
  
  if (type == "normal") {
    mean <- (lower + upper) / 2
    sd   <- (upper - lower) / (2 * z_995)
    return(c(mean, sd))
  }
  
  if (type == "truncatednormal") {
    parms <- get_tnorm_params(lower, upper, 1 - 0.995, 0.995)
    return(c(parms$mu, parms$sigma))
  }
  
  stop("Tipo inválido em get_params_from_range().")
}


# =========================================================
# FUNCAO INVERSA DA DISTRIBUICAO DE VALOR EXTREMO MINIMO
# =========================================================
# Para x em (0, 1), retorna:
# z = mu + sigma * log(-log(1 - x))
#
# No caso padronizado, basta usar mu = 0 e sigma = 1.
sevinf <- function(x, mu = 0, sigma = 1) {
  mu + sigma * log(-log(1 - x))
}


# =========================================================
# FUNCAO PARA PREPARAR OS DADOS PARA O STAN
# =========================================================
# A funcao aceita dois formatos:
# 1. Vetores 'times' e 'status'
# 2. Data frame agregado com colunas:
#    Tempo, Contagem e Censura
preparar_dados_confiabilidade <- function(times, status = NULL) {
  if (is.data.frame(times)) {
    if (!all(c("Tempo", "Contagem", "Censura") %in% names(times))) {
      stop("O data.frame deve conter as colunas 'Tempo', 'Contagem' e 'Censura'.")
    }
    
    data0 <- times
    times <- rep(data0$Tempo, data0$Contagem)
    status <- ifelse(
      rep(data0$Censura, data0$Contagem) == "Failed",
      1,
      0
    )
  } else {
    if (is.null(status)) {
      stop("Se 'times' for vetor, 'status' deve ser informado.")
    }
    
    if (length(times) != length(status)) {
      stop("'times' e 'status' devem ter o mesmo comprimento.")
    }
  }
  
  list(
    times = as.numeric(times),
    status = as.numeric(status)
  )
}


# =========================================================
# FUNCAO PRINCIPAL DE ANALISE BAYESIANA
# =========================================================
analise_bayesiana_confiabilidade <- function(
    times,
    status = NULL,
    distribution = "weibull",
    p_r = 0.10,
    prior_tpr = list(type = "flat"),
    prior_beta_sigma = list(type = "flat"),
    t_c = NULL,
    stan_model_path = "stanPRIORS.stan",
    ...
) {
  if (!file.exists(stan_model_path)) {
    stop("Arquivo do modelo Stan não encontrado.")
  }
  
  dados <- preparar_dados_confiabilidade(times, status)
  times <- dados$times
  status <- dados$status
  
  t_fail <- times[status == 1]
  t_cens <- times[status == 0]
  
  N_fail <- length(t_fail)
  N_cens <- length(t_cens)
  
  dist_map <- c(
    "weibull" = 1,
    "lognormal" = 2
  )
  
  dist_code <- dist_map[tolower(distribution)]
  
  if (is.na(dist_code)) {
    stop("Distribuição deve ser 'weibull' ou 'lognormal'.")
  }
  
  prior_type_map <- c(
    "flat" = 1,
    "lognormal" = 2,
    "truncatednormal" = 3,
    "cj" = 4
  )
  
  stan_data <- list(
    N_fail = N_fail,
    N_cens = N_cens,
    t_fail = t_fail,
    t_cens = t_cens,
    dist_code = dist_code,
    p_r = p_r
  )
  
  # -------------------------------------------------------
  # Priori para t_pr
  # -------------------------------------------------------
  prior_tpr$type <- tolower(prior_tpr$type)
  stan_data$prior_tpr_type <- prior_type_map[prior_tpr$type]
  
  if (is.na(stan_data$prior_tpr_type)) {
    stop("Tipo de priori para t_pr inválido.")
  }
  
  stan_data$prior_tpr_params <- c(0, 1)
  
  if (prior_tpr$type == "lognormal") {
    stan_data$prior_tpr_params <- get_params_from_range(
      prior_tpr$range,
      type = "lognormal"
    )
  } else if (prior_tpr$type == "truncatednormal") {
    stan_data$prior_tpr_params <- get_params_from_range(
      prior_tpr$range,
      type = "truncatednormal"
    )
  }
  
  # -------------------------------------------------------
  # Priori para beta/sigma
  # -------------------------------------------------------
  prior_beta_sigma$type <- tolower(prior_beta_sigma$type)
  stan_data$prior_beta_sigma_type <- prior_type_map[prior_beta_sigma$type]
  
  if (is.na(stan_data$prior_beta_sigma_type)) {
    stop("Tipo de priori para beta/sigma inválido.")
  }
  
  stan_data$prior_beta_sigma_params <- c(0, 1)
  
  if (prior_beta_sigma$type == "lognormal") {
    stan_data$prior_beta_sigma_params <- get_params_from_range(
      prior_beta_sigma$range,
      type = "lognormal"
    )
  } else if (prior_beta_sigma$type == "truncatednormal") {
    stan_data$prior_beta_sigma_params <- get_params_from_range(
      prior_beta_sigma$range,
      type = "truncatednormal"
    )
  }
  
  # -------------------------------------------------------
  # Quantidades adicionais para priori CJ/IJ
  # -------------------------------------------------------
  if (stan_data$prior_tpr_type == 4 || stan_data$prior_beta_sigma_type == 4) {
    if (is.null(t_c)) {
      stop("t_c deve ser fornecido para prioris do tipo CJ/IJ.")
    }
    
    stan_data$t_c <- t_c
    
    z_grid <- seq(-10, 10, length.out = 200)
    
    dist_lsinf <- switch(
      tolower(distribution),
      "weibull" = "sev",
      "lognormal" = "normal"
    )
    
    f_elements <- outro_vec(
      v = z_grid,
      censor.type = "right",
      distribution = dist_lsinf
    )
    
    stan_data$N_grid <- length(z_grid)
    stan_data$z_grid <- z_grid
    stan_data$f11_grid <- f_elements$f11
    stan_data$f12_grid <- f_elements$f12
    stan_data$f22_grid <- f_elements$f22
  } else {
    stan_data$t_c <- 0
    stan_data$N_grid <- 0
    stan_data$z_grid <- numeric(0)
    stan_data$f11_grid <- numeric(0)
    stan_data$f12_grid <- numeric(0)
    stan_data$f22_grid <- numeric(0)
  }
  
  fit <- stan(
    file = stan_model_path,
    data = stan_data,
    ...
  )
  
  return(fit)
}