#
# Sourcing this R file contains database used in main.R.
#
# Create an environment to store all datasets
  database <- new.env()
#
#########################  BEGIN import database  ##############################


###---------------------------------------###
###        Import real-world data         ###
###---------------------------------------###

####    BEGIN import citation data    ####
#
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

####    END import citation data    ####


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



####    BEGIN Sushi preference data collected in Kamishima (2003)    ####
#
# Sushi A (N = 10)
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

####    END Sushi preference data collected in Kamishima (2003)    ####



##########################  END import database  ###############################






#######################  BEGIN artificial database  ############################

## Generate artificial data for 5 entities
N <- 5
num.free <- choose(N-1,2)
database$artificial.name5 <- paste("Entity", 1:N)
num.pairs <- ncol(combn(N, 2))
database$freq.true5 <- rep(30, num.pairs)
database$s.true5 <- c(-2,-1,0,1,2)
#database$w.true5 <- mvrnorm(1, rep(0,num.free), diag(num.free))
database$w.true5 <- rep(0, num.free)
database$w.true5[1] <- 3
database$Phi.true5 <- compute.Phi.true(num.entities = N, weights = database$w.true5)
database$Phi.true5 <- round(database$Phi.true5, 2)
database$M.true5 <- compute.relations.true(num.entities = N,
                                           s = database$s.true5, 
                                           Phi = database$Phi.true5)

## Name the rows and columns
database$artificial5 <- generate.comparisons(num.entities = N, 
                                             freq.vec = database$freq.true5, 
                                             s = database$s.true5,
                                             Phi = database$Phi.true5)
rownames(database$artificial5) <- colnames(database$artificial5) <- database$artificial.name5

## Convert to binomial format
database$artificial5 <- countsToBinomial(database$artificial5)
database$artificial5$n_ij <- database$artificial5$win1 + database$artificial5$win2
database$artificial5$y_ij <- database$artificial5$win1



## Generate artificial data for 10 entities
N <- 10
num.free <- choose(N-1,2)
database$artificial.name10 <- paste("Entity", 1:N)
num.pairs <- ncol(combn(N, 2))
database$freq.true10 <- rep(50, num.pairs)
database$s.true10 <- seq(from = 0.5, to = 5, by = 0.5)
database$s.true10 <- database$s.true10 - mean(database$s.true10)
database$w.true10 <- rep(0, num.free)
database$w.true10[1] <- 7
database$w.true10[6] <- 5
database$w.true10[12] <-6
database$Phi.true10 <- compute.Phi.true(num.entities = N, weights = database$w.true10)
database$M.true10 <- compute.relations.true(num.entities = N,
                                            s = database$s.true10, 
                                            Phi = database$Phi.true10)

## Name the rows and columns
database$artificial10 <- generate.comparisons(num.entities = N, 
                                              freq.vec = database$freq.true10, 
                                              s = database$s.true10,
                                              Phi = database$Phi.true10)
rownames(database$artificial10) <- colnames(database$artificial10) <- database$artificial.name10

## Convert to binomial format
database$artificial10 <- countsToBinomial(database$artificial10)
database$artificial10$n_ij <- database$artificial10$win1 + database$artificial10$win2
database$artificial10$y_ij <- database$artificial10$win1



## Generate artificial data for 15 entities
N <- 15
num.free <- choose(N-1,2)
database$artificial.name15 <- paste("Entity", 1:N)
num.pairs <- ncol(combn(N, 2))
database$freq.true15 <- rep(100, num.pairs)
database$s.true15 <- seq(from = 0.5, to = 7.5, by = 0.5)
database$s.true15 <- database$s.true15 - mean(database$s.true15)
database$w.true15 <- rep(0, num.free)
database$w.true15[1] <- 4
database$w.true15[6] <- 3
database$w.true15[12] <-3
database$Phi.true15 <- compute.Phi.true(num.entities = N, weights = database$w.true15)
database$M.true15 <- compute.relations.true(num.entities = N,
                                            s = database$s.true15, 
                                            Phi = database$Phi.true15)

## Name the rows and columns
database$artificial15 <- generate.comparisons(num.entities = N, 
                                              freq.vec = database$freq.true15, 
                                              s = database$s.true15,
                                              Phi = database$Phi.true15)
rownames(database$artificial15) <- colnames(database$artificial15) <- database$artificial.name15

## Convert to binomial format
database$artificial15 <- countsToBinomial(database$artificial15)
database$artificial15$n_ij <- database$artificial15$win1 + database$artificial15$win2
database$artificial15$y_ij <- database$artificial15$win1



## Generate artificial data for 20 entities
N <- 20
num.free <- choose(N-1,2)
database$artificial.name20 <- paste("Entity", 1:N)
num.pairs <- ncol(combn(N, 2))
database$freq.true20 <- rep(100, num.pairs)
database$s.true20 <- seq(from = 0.5, to = 10, by = 0.5)
database$s.true20 <- database$s.true20 - mean(database$s.true20)
database$w.true20 <- rep(0, num.free)
database$w.true20[1] <- 4
database$w.true20[6] <- 3
database$w.true20[12] <-3
database$w.true20[20] <-2.5
database$Phi.true20 <- compute.Phi.true(num.entities = N, weights = database$w.true20)
database$M.true20 <- compute.relations.true(num.entities = N,
                                            s = database$s.true20, 
                                            Phi = database$Phi.true20)

## Name the rows and columns
database$artificial20 <- generate.comparisons(num.entities = N, 
                                              freq.vec = database$freq.true20, 
                                              s = database$s.true20,
                                              Phi = database$Phi.true20)
rownames(database$artificial20) <- colnames(database$artificial20) <- database$artificial.name20

## Convert to binomial format
database$artificial20 <- countsToBinomial(database$artificial20)
database$artificial20$n_ij <- database$artificial20$win1 + database$artificial20$win2
database$artificial20$y_ij <- database$artificial20$win1



## Generate artificial data for 30 entities
N <- 30
num.free <- choose(N-1,2)
database$artificial.name30 <- paste("Entity", 1:N)
num.pairs <- ncol(combn(N, 2))
database$freq.true30 <- rep(30, num.pairs)
database$s.true30 <- seq(from = 0.5, to = 15, by = 0.5)
database$s.true30 <- database$s.true30 - mean(database$s.true30)
database$w.true30 <- rep(0, num.free)
database$w.true30[1] <- 4
database$w.true30[6] <- 3
database$w.true30[12] <-3
database$w.true30[20] <-2.5
database$Phi.true30 <- compute.Phi.true(num.entities = N, weights = database$w.true30)
database$M.true30 <- compute.relations.true(num.entities = N,
                                            s = database$s.true30, 
                                            Phi = database$Phi.true30)

## Name the rows and columns
database$artificial30 <- generate.comparisons(num.entities = N, 
                                              freq.vec = database$freq.true30, 
                                              s = database$s.true30,
                                              Phi = database$Phi.true30)
rownames(database$artificial30) <- colnames(database$artificial30) <- database$artificial.name30

## Convert to binomial format
database$artificial30 <- countsToBinomial(database$artificial30)
database$artificial30$n_ij <- database$artificial30$win1 + database$artificial30$win2
database$artificial30$y_ij <- database$artificial30$win1



## Generate artificial data for 50 entities
N <- 50
num.free <- choose(N-1,2)
database$artificial.name50 <- paste("Entity", 1:N)
num.pairs <- ncol(combn(N, 2))
database$freq.true50 <- rep(50, num.pairs)
database$s.true50 <- seq(from = 0.5, to = 25, by = 0.5)
database$s.true50 <- database$s.true50 - mean(database$s.true50)
database$w.true50 <- rep(0, num.free)
database$w.true50[1] <- 4
database$w.true50[6] <- 3
database$w.true50[12] <-3
database$w.true50[20] <-2.5
database$w.true50[30] <-4
database$Phi.true50 <- compute.Phi.true(num.entities = N, weights = database$w.true50)
database$M.true50 <- compute.relations.true(num.entities = N,
                                            s = database$s.true50, 
                                            Phi = database$Phi.true50)

## Name the rows and columns
database$artificial50 <- generate.comparisons(num.entities = N, 
                                              freq.vec = database$freq.true50, 
                                              s = database$s.true50,
                                              Phi = database$Phi.true50)
rownames(database$artificial50) <- colnames(database$artificial50) <- database$artificial.name50

## Convert to binomial format
database$artificial50 <- countsToBinomial(database$artificial50)
database$artificial50$n_ij <- database$artificial50$win1 + database$artificial50$win2
database$artificial50$y_ij <- database$artificial50$win1

########################  END artificial database  #############################