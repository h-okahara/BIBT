#
# Sourcing this R file (> source("database.R")) contains database used in main.R.
#
# Create an environment to store all datasets
  database <- new.env()
#
#########################  BEGIN import database  ##############################

###------------------------###
###    Import Real Data    ###
###------------------------###

## Import baseball data from retrosheet (2020-2025)
years <- 2020:2025
gamelogs.list <- lapply(years, function(y) {
  fname <- file.path(getwd(), paste0("baseball/gl", y, ".csv"))
  read.csv(fname, header = FALSE, stringsAsFactors = FALSE)
})
gamelog <- bind_rows(gamelogs.list) %>%
  mutate(
    V4 = if_else(V4 == "OAK", "ATH", V4),
    V7 = if_else(V7 == "OAK", "ATH", V7)
  )
results <- gamelog %>%
  mutate(
    winner = if_else(V10 > V11, V4, V7),
    loser  = if_else(V10 > V11, V7, V4)
  )
results <- results[, c("winner", "loser")]
database$name.mlb <- sort(unique(c(results$winner, results$loser)))
num.teams <- length(database$name.mlb)
mlb.matrix <- table(
  factor(results$winner, levels = database$name.mlb),
  factor(results$loser, levels = database$name.mlb)
)

# Convert to binomial format
database$mlb <- expand.grid(
  player1 = database$name.mlb,
  player2 = database$name.mlb,
  stringsAsFactors = FALSE
) %>%
  filter(player1 < player2)

database$mlb <- database$mlb %>%
  mutate(
    win1 = mapply(function(a, b) mlb.matrix[a, b], player1, player2),
    win2 = mapply(function(a, b) mlb.matrix[b, a], player1, player2),
    n_ij = win1 + win2,
    y_ij = win1
  )
database$freq.mlb <- nrow(gamelog)
database$mlb <- database$mlb %>%
  mutate(
    player1 = factor(player1, levels = database$name.mlb),
    player2 = factor(player2, levels = database$name.mlb)
  )
database$network.mlb <- plot.networks(compute.M(database$mlb), num.entities = num.teams, components = c("M"),
                                      draw.flag = FALSE, weight = "prop", layout = "fr", tie_mode = "thin")

##########################  END import database  ###############################
