// Salvar como: stanPRIORS.stan

functions {
  // Função de interpolação linear para usar as grades da Matriz de Informação de Fisher (FIM)
  real linear_interp(real x, vector x_grid, vector y_grid) {
    int n = size(x_grid);
    if (x <= x_grid[1]) return y_grid[1];
    if (x >= x_grid[n]) return y_grid[n];

    int i = 1;
    while (i < n && x > x_grid[i+1]) {
      i += 1;
    }

    real w = (x - x_grid[i]) / (x_grid[i+1] - x_grid[i]);
    return (1 - w) * y_grid[i] + w * y_grid[i+1];
  }
}

data {
  // --- Dados de Confiabilidade ---
  int<lower=0> N_fail;      // Número de falhas
  int<lower=0> N_cens;      // Número de observações censuradas à direita
  vector[N_fail] t_fail;    // Tempos de falha
  vector[N_cens] t_cens;    // Tempos de censura

  // --- Configuração do Modelo ---
  int<lower=1, upper=2> dist_code; // 1: Weibull, 2: Lognormal
  real<lower=0, upper=1> p_r;     // Quantil para reparametrização (ex: 0.10)
  real<lower=0> t_c;              // Tempo de censura para prioris CJ/IJ

  // --- Configuração das Prioris ---
  int<lower=1, upper=4> prior_tpr_type;   // 1:Flat, 2:Lognormal, 3:TruncNormal, 4:CJ
  int<lower=1, upper=4> prior_beta_sigma_type; // 1:Flat, 2:Lognormal, 3:TruncNormal, 4:CJ
  vector[2] prior_tpr_params;         // [log_mean, log_sd] ou [mean, sd]
  vector[2] prior_beta_sigma_params;  // [log_mean, log_sd] ou [mean, sd]

  // --- Dados para Prioris CJ/IJ (baseado na FIM) ---
  int<lower=0> N_grid;
  vector[N_grid] z_grid;
  vector[N_grid] f11_grid;
  vector[N_grid] f12_grid;
  vector[N_grid] f22_grid;
}

parameters {
  // Parâmetros no espaço não restrito para melhor performance do MCMC [cite: 180]
  real log_tpr;
  real log_sigma;
}

transformed parameters {
  // Parâmetros no espaço original (restrito)
  real<lower=0> t_pr = exp(log_tpr);
  real<lower=0> sigma = exp(log_sigma);
  real<lower=0> beta = 1.0 / sigma;
  real mu;

  // Reparametrização de (t_pr, sigma/beta) para (mu, sigma) [cite: 98, 99, 101]
  if (dist_code == 1) { // Weibull
    // mu = log(eta), eta = t_pr / [-log(1-p_r)]^(1/beta)
    mu = log_tpr - log(-log(1 - p_r)) / beta;
  } else { // Lognormal
    // mu = log(t_pr) - qnorm(p_r) * sigma
    mu = log_tpr - inv_Phi(p_r) * sigma;
  }
}

model {
  // --- DEFINIÇÃO DAS PRIORIS (Baseado na Tabela 2 do artigo) [cite: 351] ---
  
  // Priori para o quantil t_pr
  if (prior_tpr_type == 1) {
    // Priori Flat para log(t_pr): não adiciona nada ao target (equivale a target += 0)
  } else if (prior_tpr_type == 2) { // Lognormal
    // Priori Normal para log(t_pr)
    log_tpr ~ normal(prior_tpr_params[1], prior_tpr_params[2]);
  } else if (prior_tpr_type == 3) { // Normal Truncada para t_pr > 0
    // dltnorm na Tabela H1 do apêndice [cite: 746]
    target += normal_lpdf(t_pr | prior_tpr_params[1], prior_tpr_params[2]) - normal_lcdf(0 | prior_tpr_params[1], prior_tpr_params[2]);
  } else if (prior_tpr_type == 4) { // Conditional Jeffreys (CJ)
    real z_c = (log(t_c) - mu) / sigma;
    real f11 = linear_interp(z_c, z_grid, f11_grid);
    target += 0.5 * log(f11); // Proporcional a sqrt(f11)
  }

  // Priori para o parâmetro de forma (beta para Weibull, sigma para Lognormal)
  if (prior_beta_sigma_type == 1) {
    // Priori Flat para log(sigma): não adiciona nada ao target
  } else if (prior_beta_sigma_type == 2) { // Lognormal
    if (dist_code == 1) { // Weibull: Lognormal em beta=1/sigma
      log(beta) ~ normal(prior_beta_sigma_params[1], prior_beta_sigma_params[2]);
    } else { // Lognormal: Lognormal em sigma
      log_sigma ~ normal(prior_beta_sigma_params[1], prior_beta_sigma_params[2]);
    }
  } else if (prior_beta_sigma_type == 3) { // Normal Truncada
    if (dist_code == 1) { // Weibull: Truncada em beta
      target += normal_lpdf(beta | prior_beta_sigma_params[1], prior_beta_sigma_params[2]) - normal_lcdf(0 | prior_beta_sigma_params[1], prior_beta_sigma_params[2]);
    } else { // Lognormal: Truncada em sigma
      target += normal_lpdf(sigma | prior_beta_sigma_params[1], prior_beta_sigma_params[2]) - normal_lcdf(0 | prior_beta_sigma_params[1], prior_beta_sigma_params[2]);
    }
  } else if (prior_beta_sigma_type == 4) { // Conditional Jeffreys (CJ)
    real z_c = (log(t_c) - mu) / sigma;
    real z_pr = inv_Phi(p_r);
    real f11 = linear_interp(z_c, z_grid, f11_grid);
    real f12 = linear_interp(z_c, z_grid, f12_grid);
    real f22 = linear_interp(z_c, z_grid, f22_grid);
    // Formula da IJ prior para (log(t_pr), log(sigma)) em (12) [cite: 209]
    // A parte da CJ para log(sigma) é o segundo termo da raiz quadrada.
    target += 0.5 * log(f11 * z_pr^2 - 2 * f12 * z_pr + f22);
  }

  // --- VEROSSIMILHANÇA (LIKELIHOOD) ---
  if (dist_code == 1) { // Weibull(mu, sigma) = Weibull(shape=beta, scale=exp(mu))
    if (N_fail > 0) {
      target += weibull_lpdf(t_fail | beta, exp(mu));
    }
    if (N_cens > 0) {
      target += weibull_lccdf(t_cens | beta, exp(mu));
    }
  } else { // Lognormal(mu, sigma)
    if (N_fail > 0) {
      target += lognormal_lpdf(t_fail | mu, sigma);
    }
    if (N_cens > 0) {
      target += lognormal_lccdf(t_cens | mu, sigma);
    }
  }
}


