#
# This script loads and formats the empirical datasets into a new environment (`database`) to be used in main.R.
# Sourcing this R file: > source("database.R")
#
database <- new.env()
#
#########################  BEGIN import database  ##############################

###---------------------------------------------------------------------###
###    Canary Dominance Data originated by Shoemaker (1939)             ###
###    DomArchive Data from https://github.com/DomArchive/DomArchive    ###
###---------------------------------------------------------------------###

local({
target_fileid <- "Shoemaker_1939"
raw_data <- dom.data[[target_fileid]]
metadata <- dom.metadata %>% filter(fileid == target_fileid)
if (metadata$matrix_edgelist == "Matrix") {
  mat <- raw_data$matrix
  database$dom1.long <- as.data.frame(as.table(mat)) %>%
    rename(winner = Var1, loser = Var2, count = Freq) %>%
    filter(count > 0) %>%
    uncount(count)
} else {
  database$dom1.long <- raw_data$edgelist
}
database$name.dom1 <- sort(unique(c(as.character(database$dom1.long$winner), as.character(database$dom1.long$loser))))
num.entities <- length(database$name.dom1)
dom.matrix <- table(
  factor(database$dom1.long$winner, levels = database$name.dom1),
  factor(database$dom1.long$loser, levels = database$name.dom1)
)

# Convert to binomial format
database$dom1 <- expand.grid(
  player1 = database$name.dom1,
  player2 = database$name.dom1,
  stringsAsFactors = FALSE
) %>%
  filter(player1 < player2) %>%
  mutate(
    win1 = mapply(function(a, b) dom.matrix[a, b], player1, player2),
    win2 = mapply(function(a, b) dom.matrix[b, a], player1, player2),
    n_ij = win1 + win2,
    y_ij = win1,
    player1 = factor(player1, levels = database$name.dom1),
    player2 = factor(player2, levels = database$name.dom1)
  )
database$network.dom1 <- plot.networks(compute.M(database$dom1), num.entities = num.entities, components = c("M"),
                                      draw.flag = FALSE, weight = "prop", layout = "fr", tie_mode = "thin")

## Covariates Information
database$covs.dom1 <- data.frame(
  entity.name = database$name.dom1,
  sex    = c("female", "female", "female", "female", "female", "male", "male", "male", "male", "male"),
  weight = c(19.1, 19.2, 18.3, 19.4, 16.1, 16.0, 20.4, 17.5, 16.4, 18.2)
)


## Build Covariate Matrix X_E for Canary Dominance Data
sex.idx <- ifelse(database$covs.dom1$sex == "male", 1, 0)
weights <- as.vector(scale(database$covs.dom1$weight))
m.mat <- matrix(0, nrow = num.entities, ncol = num.entities)
mate.idx <- list(c(1, 8), c(2, 7), c(3, 6), c(4, 10), c(5, 9))
for (mp in mate.idx) {
  i <- mp[1]
  j <- mp[2]
  m.mat[i, j] <- 1
  m.mat[j, i] <- 1
}

# Construct Covariate Matrix X_E
pairs <- t(combn(num.entities, 2))
num.pairs <- nrow(pairs)
X_E <- matrix(0, nrow = 2, ncol = num.pairs)
rownames(X_E) <- c("Sex Effect", "Mate Effect")
p1 <- pairs[, 1]
p2 <- pairs[, 2]
sex.diff <- sex.idx[p1] - sex.idx[p2]
mate.vec <- m.mat[cbind(p1, p2)]
X_E[1, ] <- sex.diff
X_E[2, ] <- mate.vec * sex.diff
database$X_E.dom1 <- X_E




###---------------------------------------------------------------------###
###    Guanaco Dominance Data originated by Correa et al. (2013)        ###
###    DomArchive Data from https://github.com/DomArchive/DomArchive    ###
###---------------------------------------------------------------------###

target_fileid <- "Correa_2013a" # Total Interactions
raw_data <- dom.data[[target_fileid]]
metadata <- dom.metadata %>% filter(fileid == target_fileid)
if (metadata$matrix_edgelist == "Matrix") {
  mat <- raw_data$matrix
  database$dom2.long <- as.data.frame(as.table(mat)) %>%
    rename(winner = Var1, loser = Var2, count = Freq) %>%
    filter(count > 0) %>%
    uncount(count)
} else {
  database$dom2.long <- raw_data$edgelist
}
database$name.dom2 <- sort(unique(c(as.character(database$dom2.long$winner), as.character(database$dom2.long$loser))))
num.entities <- length(database$name.dom2)
dom.matrix <- table(
  factor(database$dom2.long$winner, levels = database$name.dom2),
  factor(database$dom2.long$loser, levels = database$name.dom2)
)

# Convert to binomial format
database$dom2 <- expand.grid(
  player1 = database$name.dom2,
  player2 = database$name.dom2,
  stringsAsFactors = FALSE
) %>%
  filter(player1 < player2) %>%
  mutate(
    win1 = mapply(function(a, b) dom.matrix[a, b], player1, player2),
    win2 = mapply(function(a, b) dom.matrix[b, a], player1, player2),
    n_ij = win1 + win2,
    y_ij = win1,
    player1 = factor(player1, levels = database$name.dom2),
    player2 = factor(player2, levels = database$name.dom2)
  )
database$network.dom2 <- plot.networks(compute.M(database$dom2), num.entities = num.entities, components = c("M"),
                                      draw.flag = FALSE, weight = "prop", layout = "fr", tie_mode = "thin")
})
##########################  END import database  ###############################
