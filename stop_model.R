library(dplyr)

file_path <- "data.xlsx"
data <- read_xlsx(file_path)

estimate_probs <- function(df) {
  # Ignore cases without embryos
  df_valid <- df[!is.na(df$Last) & df$K > 0, ]
  if (nrow(df_valid) == 0) return(c(p_lb=NA, p_sat=NA, p_drop=NA))
  
  # P live birth
  total_lb <- sum(df_valid$X)
  total_transfer <- sum(df_valid$K)
  p_lb <- total_lb / total_transfer
  se_p <- sqrt(p_lb*(1-p_lb)/total_transfer)
  
  # Stop after live birth
  stop_after_lb <- sum((df_valid$K < df_valid$N) & 
                         (df_valid$Last == 1))
  no_stop_after_lb <- sum((df_valid$K == df_valid$N) & 
                            (df_valid$Last == 1))
  n_1 <- total_lb - no_stop_after_lb
  p1 <- if (n_1 > 0) stop_after_lb /
    n_1 else 0
  se_1 <- sqrt(p1*(1-p1)/n_1)
  
  # Stop after no live birth
  stop_after_no_lb <- sum((df_valid$K < df_valid$N) &
                            (df_valid$Last == 0))
  temp_2 <- total_transfer - total_lb
  no_stop_after_no_lb <- sum((df_valid$K == df_valid$N) & 
                               (df_valid$Last == 0))
  n_2 <- temp_2 - no_stop_after_no_lb
  p2 <- if (n_2 > 0) stop_after_no_lb / 
    n_2 else 0
  se_2 <- sqrt(p2*(1-p2)/n_2)
  
  return(c(p_lb = p_lb, 
           p1 = p1, 
           p2 = p2,
           se_lb = se_p,
           se1 = se_1,
           se2 = se_2))
}

simulate_data <- function(N, p_lb, p_sat, p_drop) {
  n <- length(N) 
  K <- numeric(n)
  X <- numeric(n)
  Last_Outcome <- rep(NA, n) 
  
  for (i in 1:n) {
    euploids <- N[i]
    if (euploids == 0) {
      K[i] <- 0 
      X[i] <- 0 
      next 
    }
    
    current_x <- 0
    current_k <- 0
    last_res <- NA 
    
    for (k in 1:euploids) {
      current_k <- k
      outcome <- rbinom(1, 1, p_lb)
      current_x <- current_x + outcome
      last_res <- outcome
      
      stop_now <- FALSE
      if (k < euploids) {
        if (outcome == 1) {
          if (rbinom(1, 1, p_sat) == 1) stop_now <- TRUE
        } else {
          if (rbinom(1, 1, p_drop) == 1) stop_now <- TRUE
        }
      }
      if (stop_now) break
    }
    K[i] <- current_k 
    X[i] <- current_x
    Last_Outcome[i] <- last_res
  }
  return(data.frame(N = N, K = K, X = X, Last_Res = Last_Outcome))
}

# The outcome in the last transfer (or NA)
last <- (data |> 
           rowwise() |> 
           mutate(temp=coalesce(LB_8, LB_7, LB_6, LB_5, 
                                LB_4, LB_3, LB_2, LB_1)) |> 
           ungroup() |> 
           select(temp))$temp

mle_est <- estimate_probs(data.frame(X=data$TOT_Lbscong,
                                     K=data$TOT_scong,
                                     N=data$eup,
                                     Last=last))
mle_est

sim_data <- simulate_data(data$eup, mle_est[1], 
                          mle_est[2], mle_est[3])
table(sim_data$K)
table(sim_data$X)

sim_data |>
  filter(K < N) |>
  summarise(sum(X)/sum(K))

sim_data |>
  filter(K == N) |>
  summarise(sum(X)/sum(K))

# And to make sure we recover the parameters under the simulated data:
estimate_probs(sim_data)