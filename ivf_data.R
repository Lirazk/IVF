library(dplyr)
source("utils.R")

# should probably refactor the loops so I wouldn't repeat the same thing
# each time with minimal change, but alas

file_path <- "data.xlsx"
data <- read_xlsx(file_path)
data_closed <- data |>
  filter(emb.residui == 0)

r2_range <- c(0.05, 0.1, 0.15)
K_range <- c(0.005, 0.01, 0.05, 0.1)
bootstrap_samples <- 1000

# For cases without untransferred embryos
prob_LB <- data |>
  filter(emb.residui == 0) |>
  summarise(sum(TOT_Lbscong) / sum(TOT_scong)) |>
  as.numeric()
# And for only > residual embryo?
prob_LB_open <- data |>
  filter(emb.residui > 0) |>
  summarise(sum(TOT_Lbscong) / sum(TOT_scong)) |>
  as.numeric()
prob_LB_all <- data |>
  summarise(sum(TOT_Lbscong) / sum(TOT_scong)) |>
  as.numeric()

data <- data |>
  mutate(age_group = factor(
    case_when(
      maternal.age < 35 ~ "<35",
      between(maternal.age, 35, 37) ~ "35-37",
      between(maternal.age, 38, 40) ~ "38-40",
      between(maternal.age, 41, 42) ~ "41-42",
      maternal.age > 42 ~ ">42"
    ),
    levels = c(
      "<35", "35-37", "38-40",
      "41-42", ">42"
    )
  ))

# Intention to screen
# Should be for closed, open, and all
prob_LBs <- c(prob_LB, prob_LB_open, prob_LB_all)
names <- c("Closed", "Open", "All")

data$transfers <- data$TOT_scong
data$sim_transfers <- data$eup-data$transfers
data$live_births <- data$TOT_Lbscong

result <- NULL
system.time({
  for (p_lb in prob_LBs) {
    for (r2 in r2_range) {
      for (K in K_range) {
        temp_data <- data.table(data)
        
        est <- estimate_risk(temp_data, p_lb = p_lb)
        bootstrap <- boot(temp_data,
                          estimate_risk,
                          R = bootstrap_samples,
                          p_lb = p_lb)
        CI <- try(boot.ci(bootstrap, type = "perc"), silent = T)
        if (is.null(CI)) {
          result <- rbind(result,
                          data.frame(
                            est = est,
                            p_lb = p_lb,
                            r2 = r2,
                            K = K,
                            lower = NA,
                            upper = NA
                          ))
        }
        else {
          result <- rbind(
            result,
            data.frame(
              est = est,
              p_lb = p_lb,
              r2 = r2,
              K = K,
              lower = CI$percent[4],
              upper = CI$percent[5]
            )
          )
        }
      }
    }
  }
})

knitr::kable(result)

# At least one live birth
result <- NULL
system.time({
  for (p_lb in prob_LBs) {
    for (r2 in r2_range) {
      for (K in K_range) {
        temp_data <- data.table(data)
        
        est <- estimate_risk(temp_data,
                             p_lb = p_lb,
                             intention_to_screen = F)
        bootstrap <- boot(
          temp_data,
          estimate_risk,
          R = bootstrap_samples,
          p_lb = p_lb,
          intention_to_screen = F
        )
        CI <- try(boot.ci(bootstrap, type = "perc"), silent = T)
        if (is.null(CI)) {
          result <- rbind(result,
                          data.frame(
                            est = est,
                            p_lb = p_lb,
                            r2 = r2,
                            K = K,
                            lower = NA,
                            upper = NA
                          ))
        }
        else {
          result <- rbind(
            result,
            data.frame(
              est = est,
              p_lb = p_lb,
              r2 = r2,
              K = K,
              lower = CI$percent[4],
              upper = CI$percent[5]
            )
          )
        }
      }
    }
  }
})

knitr::kable(result)

# All cycles
# Intention to screen
result <- NULL
system.time({
  for (p_lb in prob_LBs) {
    for (r2 in r2_range) {
      for (K in K_range) {
        total_attempts <- data |>
          group_by(ID) |>
          summarise(
            transfers = sum(transfers),
            live_births = sum(live_births),
            eup = sum(eup),
            sim_transfers = eup - transfers
          )
        
        temp_data <- data.table(total_attempts)
        
        est <- estimate_risk(temp_data, p_lb = p_lb)
        bootstrap <- boot(temp_data,
                          estimate_risk,
                          R = bootstrap_samples,
                          p_lb = p_lb)
        CI <- try(boot.ci(bootstrap, type = "perc"), silent = T)
        if (is.null(CI)) {
          result <- rbind(result,
                          data.frame(
                            est = est,
                            p_lb = p_lb,
                            r2 = r2,
                            K = K,
                            lower = NA,
                            upper = NA
                          ))
        }
        else {
          result <- rbind(
            result,
            data.frame(
              est = est,
              p_lb = p_lb,
              r2 = r2,
              K = K,
              lower = CI$percent[4],
              upper = CI$percent[5]
            )
          )
        }
      }
    }
  }}
)

knitr::kable(result)

# At least one live birth
result <- NULL
system.time({
  for (p_lb in prob_LBs) {
    for (r2 in r2_range) {
      for (K in K_range) {
        total_attempts <- data |>
          group_by(ID) |>
          summarise(
            transfers = sum(transfers),
            live_births = sum(live_births),
            eup = sum(eup),
            sim_transfers = eup - transfers
          )
        
        temp_data <- data.table(total_attempts)
        
        est <- estimate_risk(temp_data,
                             p_lb = p_lb,
                             intention_to_screen = F)
        bootstrap <- boot(
          temp_data,
          estimate_risk,
          R = bootstrap_samples,
          p_lb = p_lb,
          intention_to_screen = F
        )
        CI <- try(boot.ci(bootstrap, type = "perc"), silent = T)
        if (is.null(CI)) {
          result <- rbind(result,
                          data.frame(
                            est = est,
                            p_lb = p_lb,
                            r2 = r2,
                            K = K,
                            lower = NA,
                            upper = NA
                          ))
        }
        else {
          result <- rbind(
            result,
            data.frame(
              est = est,
              p_lb = p_lb,
              r2 = r2,
              K = K,
              lower = CI$percent[4],
              upper = CI$percent[5]
            )
          )
        }
      }
    }
  }
})

knitr::kable(result)

# By age group
# Intention to screen
result <- NULL
system.time({
  for (p_lb in prob_LBs) {
    for(age in levels(data$age_group)) {
      temp_data <- data |>
        filter(age_group == age)
      temp_data <- data.table(temp_data)
      
      est <- estimate_risk(temp_data, p_lb = p_lb)
      bootstrap <- boot(temp_data, estimate_risk, 
                        R=bootstrap_samples, p_lb = p_lb)
      CI <- try(boot.ci(bootstrap, type="perc"),
                silent = T)
      if (is.null(CI)) {
        result <- rbind(result,
                        data.frame(est = est,
                                   p_lb = p_lb,
                                   age_group=age,
                                   r2=r2,
                                   K=K,
                                   lower = NA,
                                   upper = NA))
      }
      else {
        result <- rbind(result,
                        data.frame(est = est,
                                   p_lb = p_lb,
                                   age_group=age,
                                   r2=r2,
                                   K=K,
                                   lower = CI$percent[4],
                                   upper = CI$percent[5]))
      }
    }
  }})

knitr::kable(result)

# At least one live birth
result <- NULL
system.time({
  for (p_lb in prob_LBs) {
    for(age in levels(data$age_group)) {
      temp_data <- data |>
        filter(age_group == age)
      temp_data <- data.table(temp_data)
      
      est <- estimate_risk(temp_data, p_lb = p_lb, 
                           intention_to_screen = F)
      bootstrap <- boot(temp_data, estimate_risk, 
                        R=bootstrap_samples, p_lb=p_lb,
                        intention_to_screen = F)
      CI <- try(boot.ci(bootstrap, type="perc"),
                silent = T)
      if (is.null(CI)) {
        result <- rbind(result,
                        data.frame(est = est,
                                   p_lb = p_lb,
                                   age_group=age,
                                   r2=r2,
                                   K=K,
                                   lower = NA,
                                   upper = NA))
      }
      else {
        result <- rbind(result,
                        data.frame(est = est,
                                   p_lb = p_lb,
                                   age_group=age,
                                   r2=r2,
                                   K=K,
                                   lower = CI$percent[4],
                                   upper = CI$percent[5]))
      }
    }
  }})

knitr::kable(result)