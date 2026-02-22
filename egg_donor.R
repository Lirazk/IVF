library(readxl)
library(boot)

source("utils.R")

file_path <- "egg_donor.xlsx"
data <- read_xlsx(paste0(path, "egg_donor.xlsx"))

r2_range <- c(0.05, 0.1, 0.15)
K_range <- c(0.005, 0.01, 0.05, 0.1)
bootstrap_sim <- 1000

df <- data |>
  select(id = ID_DONOR, 
         id_stim=`ID_DONOR STIMULATION`,
         id_recip=ID_Recipient,
         oocytes=`Number of oocytes retrieved`, 
         live_births = `Number of live-born`,
         transfers = `Number of transferred embryos`,
         blasto = `UsableBlastocyst`,
         age = `Edad`) |>
  # Not really needed as it is redefined after we combine
  mutate(sim_transfers = blasto - transfers)

p_lb <- sum(df$live_births) / sum(df$transfers)

# Combine all results per donor
df <- df |>
  group_by(id) |>
  select(id, oocytes, live_births, 
         transfers, blasto, age) |>
  summarise(oocytes=first(oocytes), 
            live_births=sum(live_births),
            transfers=sum(transfers), 
            blasto=sum(blasto),
            age=first(age)) |>
  mutate(sim_transfers = blasto - transfers)

df <- data.table(df)

# Intention to screen
result <- NULL
for (r2 in r2_range) {
  for (K in K_range) {
    est <- estimate_risk(df, p_lb = p_lb)
    bootstrap <- boot(df, estimate_risk, R = bootstrap_sim, p_lb = p_lb)
    CI <- boot.ci(bootstrap, type = "perc")
    result <- rbind(result,
                    data.frame(
                      est = est,
                      r2 = r2,
                      K = K,
                      lower = CI$percent[4],
                      upper = CI$percent[5]
                    ))
  }
}

knitr::kable(result)

# At least one live birth
for (r2 in r2_range) {
  for (K in K_range) {
    est <- estimate_risk(df, p_lb = p_lb, 
                         intention_to_screen = F)
    bootstrap <- boot(df, estimate_risk, R=bootstrap_sim,
                      p_lb = p_lb, intention_to_screen = F)
    CI <- boot.ci(bootstrap, type="perc")
    result <- rbind(result,
                    data.frame(est = est,
                               r2=r2,
                               K=K,
                               lower = CI$percent[4],
                               upper = CI$percent[5]))
  }
}

knitr::kable(result)