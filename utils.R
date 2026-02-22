library(data.table)

risk_lowest = function(r2,K,n)
{
  stopifnot(r2 >=0, r2 <=1, K > 0,  K < 1, n >= 0)
  if (n == 0) return(K)
  r = sqrt(r2)
  zk = qnorm(K, lower.tail=F)
  integrand_lowest = function(t)
  {
    arg = (zk-t*sqrt(1-r2/2)) / (r/sqrt(2))
    y = dnorm(t)*pnorm(arg, lower.tail=F)^n
    return(y)
  }
  
  risk = integrate(integrand_lowest,-Inf,Inf)$value
  return(risk)
}

risk_lowest <- Vectorize(risk_lowest)

# Make sure that for each variable combination we run the risk_lowest only once.
memory <- list()
risk_lowest_mem <- function(r2, K, transfers, live_births, memory) {
  if (length(memory) == 0) {
    value <- risk_lowest(r2, K, 0:transfers + live_births)
    memory[[rlang::hash(c(r2, K, transfers, live_births))]] <<- value
    return(value)
  }
  specific <- memory[[rlang::hash(c(r2, K, transfers, live_births))]]
  if (length(specific) == 0) {
    value <- risk_lowest(r2, K, 0:transfers + live_births)
    memory[[rlang::hash(c(r2, K, transfers, live_births))]] <<- value
    return(value)
  }
  return(specific)
}

estimate_risk <- function(data, index=NULL, p_lb, 
                          intention_to_screen=T, ...) {
  if (is.null(index)) {
    df <- data
  }
  else {
    df <- data[index,]
  }
  
  # Instead of working directly, find all combinations of sim transfers and live births
  # since live births is only 0-2, it reduces the sample size by a large margin.
  
  temp <- df[, .N, by = .(live_births, sim_transfers)]
  sim_transfers <- temp$sim_transfers
  live_births_all <- temp$live_births
  
  total_expected_risk <- 0
  total_expected_families <- 0
  
  for (i in 1:nrow(temp)) {
    n <- sim_transfers[i]
    live_births <- live_births_all[i]
    
    count <- temp[i, ]$N
    risk_vector <- risk_lowest_mem(r2 = r2, K = K, n, live_births, memory) 
    if (intention_to_screen) {
      # Here if n+live_birth=0, the risk reduction is K
      risk_per_family <- sum(risk_vector * dbinom(0:n, n, p_lb))
      
      total_expected_risk <- total_expected_risk + (count * risk_per_family)
      total_expected_families <- total_expected_families + count
    }
    else {
      # Here we condition on n+live_births>0
      # If both are 0, the risk is 0
      if (live_births > 0) {
        # We are already "ok", so the result is the same as before
        
        # P(N>0) is 1
        truncated_prob <- 1
        risk_per_family <- sum(risk_vector * dbinom(0:n, n, p_lb))
      } 
      else {
        # No known births
        # So we need to condition on >0
        truncated_prob <- pbinom(0, n, p_lb, lower.tail = F)
        if (truncated_prob == 0) {next}
        # Sum only on n>0
        
        risk_per_family <- sum(risk_vector[-1] * dbinom(1:n, n, p_lb))
      }
      total_expected_risk <- total_expected_risk + (count * risk_per_family)
      total_expected_families <- total_expected_families + (count * truncated_prob)
    }
  }
  risk <- total_expected_risk / total_expected_families
  (K-risk) / K
}