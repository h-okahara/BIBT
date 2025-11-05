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

# Import from BradleyTerry2 package
data(citations)
database$citations4 <- countsToBinomial(citations) 
database$citations4$n_ij <- database$citations4$win1 + database$citations4$win2
database$citations4$y_ij <- database$citations4$win1
database$name4 <- c("Biometrika", "Comm Statist", "JASA", "JRSS-B")
database$network.citations4 <- plot.networks(compute.M(database$citations4), num.entities = 4, components = c("M"), 
                                             draw.flag = FALSE, weight = "prop", layout = "fr", tie_mode = "thin")



# Cross-citations involving probability journals, giving total citations for the years 1987-88. 
# Rows correspond to citing journal, columns to cited journal.
citation.matrix <- matrix(c(
  468, 255,  33,  46,  72,   74,  # AnnPr
  333, 322,  47,  47,  72,   76,  # PrTh/ZW
  208, 155,  93,  76,  40,   41,  # StochPr
  121,  31,  37, 283,  26,   35,  # JAP
  101,  60,  23,  38, 344,   63,  # ThPrApp
   76,  81,  13,  14,  50, 1009   # AnnSt
), nrow = 6, byrow = TRUE)

# Name the rows and columns
rownames(citation.matrix) <- 
  colnames(citation.matrix) <-
  database$name6 <- c("AnnPr", "PrTh/ZW", "StochPr", "JAP", "ThPrApp", "AnnSt")

# Convert to binomial format
database$citations6 <- countsToBinomial(t(citation.matrix))
database$citations6$n_ij <- database$citations6$win1 + database$citations6$win2
database$citations6$y_ij <- database$citations6$win1
database$network.citations6 <- plot.networks(compute.M(database$citations6), num.entities = 6, components = c("M"), 
                                             draw.flag = FALSE, weight = "prop", layout = "fr", tie_mode = "thin")



# Cross-citations involving statistics journals, giving total citations for the years 1987-89. 
# Rows correspond to citing journal, columns to cited journal. 
citation.matrix <- matrix(c(
  1623,  42, 275,  47,  340, 179,  28,  57,  # AnnSt
   155, 770, 419,  37,  348, 163,  85,  66,  # Biometrics
   466, 141, 714,  33,  320, 284,  68,  81,  # Biometrika
  1025, 273, 730, 425,  813, 276,  94, 418,  # Comm Statist
   739, 264, 498,  68, 1072, 325, 104, 117,  # JASA
   182,  60, 221,  17,  142, 188,  43,  27,  # JRSS-B
    88, 134, 163,  19,  145, 104, 211,  62,  # JRSS-C
   112,  45, 147,  27,  181, 116,  41, 386   # Technometrics
), nrow = 8, byrow = TRUE)

# Name the rows and columns
rownames(citation.matrix) <- 
  colnames(citation.matrix) <- 
  database$name8 <- c(
    "AnnSt", 
    "Biometrics", 
    "Biometrika", 
    "Comm Statist", 
    "JASA", 
    "JRSS-B", 
    "JRSS-C", 
    "Technometrics"
  )

# Convert to binomial format
database$citations8 <- countsToBinomial(t(citation.matrix))
database$citations8$n_ij <- database$citations8$win1 + database$citations8$win2
database$citations8$y_ij <- database$citations8$win1
database$network.citations8 <- plot.networks(compute.M(database$citations8), num.entities = 8, components = c("M"), 
                                             draw.flag = FALSE, weight = "prop", layout = "fr", tie_mode = "thin")



# Cross-citations involving applied statistics journals, giving total citations for the years 1987-89.
# Rows correspond to citing journal, columns to cited journal.
citation.matrix <- matrix(c(
  578, 202,  69,  75, 129,  47,  51,  37,  33,  # ApplStat
  195, 696,  92, 134, 187,  83,  47,  76,  54,  # Biostats
  118, 177, 330, 118, 187,  76,  28,  33,  43,  # CommStat
   96, 171, 101, 251, 153,  89,  37,  38,  42,  # StatMed
   98, 165, 116, 108, 568, 101,  45,  50,  57,  # JCompGr
   28,  56,  44,  35,  62, 129,  18,  21,  17,  # StatSoft
   41,  57,  24,  38,  73,  37,  86,  43,  32,  # CompStat
   32,  49,  25,  39,  67,  33,  40, 116,  26,  # Technometrics
   25,  43,  30,  35,  64,  27,  37,  24,  81   # ApplMath
), nrow = 9, byrow = TRUE)

# Name the rows and columns
rownames(citation.matrix) <- 
  colnames(citation.matrix) <- 
  database$name9 <- c(
    "ApplStat", 
    "Biostats", 
    "CommStat", 
    "StatMed", 
    "JCompGr", 
    "StatSoft", 
    "CompStat", 
    "Technometrics", 
    "ApplMath"
  )

# Convert to binomial format using countsToBinomial
database$citations9 <- countsToBinomial(t(citation.matrix))
database$citations9$n_ij <- database$citations9$win1 + database$citations9$win2
database$citations9$y_ij <- database$citations9$win1
database$network.citations9 <- plot.networks(compute.M(database$citations9), num.entities = 9, components = c("M"), 
                                             draw.flag = FALSE, weight = "prop", layout = "fr", tie_mode = "thin")



# Real-world data (Sumo) Pairwise Comparison Data from 2005 to 2009
result.matrix <- matrix(c(
   0, 11, 14, 16, 21,  1, 16, 11, 19, 19,  # 朝青龍
  10,  0, 19, 20, 18,  4, 13, 14, 18, 15,  # 白鵬
   3,  2,  0,  9,  9,  0,  8,  7, 12, 14,  # 魁皇
   4,  4,  9,  0, 10,  2,  7, 13, 10, 10,  # 千代大海
   2,  8, 14, 15,  0,  3, 10, 10, 15, 13,  # 琴光喜
   0,  0,  1,  0,  1,  0,  4,  2,  3,  1,  # 高見盛
   4,  6, 12, 12, 14,  3,  0, 12,  9, 10,  # 曰馬富士
   4,  4, 12,  6, 10,  3,  9,  0, 12,  5,  # 稀勢ノ里
   2,  4,  8, 12,  9,  7,  6,  4,  0,  8,  # 雅山
   0,  2,  2,  6,  7,  7,  4,  8,  8,  0   # 旭天鵬 
), nrow = 10, byrow = TRUE)

# Name the rows and columns
rownames(result.matrix) <- 
  colnames(result.matrix) <- 
  database$name.sumo <- c(
    "朝青龍",
    "白鵬",
    "魁皇",
    "千代大海",
    "琴光喜",
    "高見盛",
    "曰馬富士",
    "稀勢ノ里",
    "雅山",
    "旭天鵬"
  )

# Convert to binomial format using countsToBinomial
database$sumo <- countsToBinomial(result.matrix)
database$sumo$n_ij <- database$sumo$win1 + database$sumo$win2
database$sumo$y_ij <- database$sumo$win1
database$network.sumo <- plot.networks(compute.M(database$sumo), num.entities = 10, components = c("M"), 
                                       draw.flag = FALSE, weight = "prop", layout = "fr", tie_mode = "thin")



# Paired comparisons consists in tasting experiments conducted in the National 
# Higher Institute for Education in the Food Industry in France in 2002.
result.matrix <- matrix(c(
   0, 39, 64, 40, 61, 76, 46, # 1
  61,  0, 65, 59, 55, 85, 60, # 2
  36, 35,  0, 31, 25, 41, 35, # 3
  60, 41, 69,  0, 41, 80, 28, # 4
  39, 45, 75, 59,  0, 71, 37, # 5
  24, 15, 59, 20, 29,  0, 18, # 6
  54, 40, 65, 72, 63, 82,  0  # 7
), nrow = 7, byrow = TRUE)

# Name the rows and columns
rownames(result.matrix) <- colnames(result.matrix) <- 
  database$name.cornflakes <- c("1","2","3","4","5","6","7")

# Convert to binomial format using countsToBinomial
database$cornflakes <- countsToBinomial(result.matrix)
database$cornflakes$n_ij <- database$cornflakes$win1 + database$cornflakes$win2
database$cornflakes$y_ij <- database$cornflakes$win1
database$network.cornflakes <- plot.networks(compute.M(database$cornflakes), num.entities = 7, components = c("M"), 
                                             draw.flag = FALSE, weight = "prop", layout = "fr", tie_mode = "thin")



# Sushi A (N = 10) preference data collected in Kamishima (2003)
sushiA.data <- read.csv("Sushi A.csv", header = FALSE, sep = " ")
sushiA.data <- sushiA.data[2:nrow(sushiA.data), 3:12]
sushiA.matrix <- matrix(0, nrow = 10, ncol = 10, dimnames = list(1:10, 1:10))

# Name the rows and columns
rownames(sushiA.matrix) <- colnames(sushiA.matrix) <-
  database$name.sushiA <- c("shrimp", "sea eel", "tuna", "squid", "sea urchin",
                            "salmon roe", "egg", "fatty tuna", "tuna roll", "cucumber  roll")

pairs.list <- lapply(seq_len(nrow(sushiA.data)), function(i) {
  row <- as.numeric(sushiA.data[i, ])
  t(combn(row, 2))
})
pairs.list <- as.data.frame(do.call(rbind, pairs.list))
colnames(pairs.list) <- c("winner", "loser")

pairs.counts <- table(pairs.list)
winner.idx <- as.numeric(rownames(pairs.counts)) + 1
loser.idx <- as.numeric(colnames(pairs.counts)) + 1
sushiA.matrix[winner.idx, loser.idx] <- pairs.counts

# Convert to binomial format using countsToBinomial
database$sushiA <- countsToBinomial(sushiA.matrix)
database$sushiA$n_ij <- database$sushiA$win1 + database$sushiA$win2
database$sushiA$y_ij <- database$sushiA$win1
database$network.sushiA <- plot.networks(compute.M(database$sushiA), num.entities = 10, components = c("M"), 
                                         draw.flag = FALSE, weight = "prop", layout = "fr", tie_mode = "thin")

##########################  END import database  ###############################
