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

# Import from BradleyTerry2 package
data(citations)
database$citations.4 <- countsToBinomial(citations) 
database$citations.4$n_ij <- database$citations.4$win1 + database$citations.4$win2
database$citations.4$y_ij <- database$citations.4$win1
database$name.4 <- c("Biometrika", "Comm Statist", "JASA", "JRSS-B")


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
  database$name.6 <- c(
    "AnnPr", 
    "PrTh/ZW", 
    "StochPr", 
    "JAP", 
    "ThPrApp", 
    "AnnSt"
  )

# Convert to binomial format
database$citations.6 <- countsToBinomial(t(citation.matrix))
database$citations.6$n_ij <- database$citations.6$win1 + database$citations.6$win2
database$citations.6$y_ij <- database$citations.6$win1




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
  database$name.8 <- c(
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
database$citations.8 <- countsToBinomial(t(citation.matrix))
database$citations.8$n_ij <- database$citations.8$win1 + database$citations.8$win2
database$citations.8$y_ij <- database$citations.8$win1




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
  database$name.9 <- c(
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
database$citations.9 <- countsToBinomial(t(citation.matrix))
database$citations.9$n_ij <- database$citations.9$win1 + database$citations.9$win2
database$citations.9$y_ij <- database$citations.9$win1



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

##########################  END import database  ###############################




#######################  BEGIN artificial database  ############################

###-----------------------###
###    Generate F.true    ###
###-----------------------###

generate.F.true <- function(N = NULL, K = NULL, decay = NULL, seed = 73) {
  ## Helper function: generate F.true
  generate.F.worths <- function(N = NULL, K = NULL, decay = NULL, seed = 73) {
    stopifnot(N >= 2, K >= 1, decay > 0, decay < 1)
    set.seed(seed)
    F.true <- matrix(NA_real_, nrow = K, ncol = N)
    
    for (k in 1:K){
      z <- rnorm(N)
      z <- scale(z, center = TRUE, scale = TRUE)
      F.true[k, ] <- sqrt(decay^(k-1)) * z
    }
    return(F.true)
  }
  F.true <- generate.F.worths(N = N, K = K, decay = decay, seed = seed)
  
  ## Sort the original column f_i into a descending index
  F.true[1, ] <- sort(F.true[1, ], decreasing = TRUE)
  F.true[2, ] <- sort(F.true[2, ], decreasing = FALSE)
  return(F.true)
}




###------------------------------###
###    Create artificial data    ###
###------------------------------###

Create.Artificial.Data <- function(num.freq = NULL, w0 = NULL, F0 = NULL, gamma0 = 0, seed = 73) {
  N <- ncol(F0) # number of entities
  K <- nrow(F0) # dimensionality of each entity
  
  ## Make skew-symmetric Gamma matrix
  J_skew.mat <- matrix(0, K, K)
  J_skew.mat[upper.tri(J_skew.mat)] <-  1
  J_skew.mat[lower.tri(J_skew.mat)] <- -1
  Gamma  <- gamma0 * J_skew.mat
  
  ## Initiate an N×N matrix storing results (diagonal is NA)
  result <- matrix(NA_integer_, nrow = N, ncol = N)
  rownames(result) <- colnames(result) <- 
    if (!is.null(colnames(F0))) colnames(F0) else paste0("Entity", 1:N)
  
  ## Simulate for each pair (i,j)
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      M_ij <- as.numeric(crossprod(w0, F0[, i] - F0[, j]) 
                         + crossprod(F0[, i], Gamma %*% F0[, j]))
      p_ij <- 1 / (1 + exp(-M_ij))
      win.freq <- rbinom(1, size = num.freq, prob = p_ij)
      result[i, j] <- win.freq
      result[j, i] <- num.freq - win.freq
    }
  }
  return(result)
}


# Generate artificial data for 10 entities
N <- 10
database$K0.true10 <- 2
database$artificial.name10 <- paste("Entity", 1:N)
database$gamma.true10 <- 1.5
database$w.true10 <- c(1, 0.8)
database$w.true10 <- database$w.true10 / norm(database$w.true10, type = "2")
database$F.true10 <- generate.F.true(N = N, K = database$K0.true10, decay = 0.7, seed = 73)
# worths <- crossprod(database$w.true10, database$F.true10)

# rowVars(database$F.true10)
# rowSds(database$F.true10)
# var10 <- rowVars(database$F.true10)
# var10 <- rowVars(database$w.true10 * database$F.true10)
# var10 / sum(var10)

# Name the rows and columns
database$artificial.10 <- Create.Artificial.Data(num.freq = 30, 
                                                 w0 = database$w.true10, 
                                                 F0 = database$F.true10,
                                                 gamma0 = database$gamma.true10)
rownames(database$artificial.10) <- 
  colnames(database$artificial.10) <- database$artificial.name10

# Convert to binomial format
database$artificial.10 <- countsToBinomial(database$artificial.10)
database$artificial.10$n_ij <- database$artificial.10$win1 + database$artificial.10$win2
database$artificial.10$y_ij <- database$artificial.10$win1



# Generate artificial data for 30 entities
N <- 15
database$K0.true15 <- 3
database$artificial.name15 <- paste("Entity", 1:N)
database$gamma.true15 <- 0.8
database$w.true15 <- c(1, 0.9, 0.8)
database$w.true15 <- database$w.true15 / norm(database$w.true15, type = "2")
database$F.true15 <- generate.F.true(N = N, K = database$K0.true15, decay = 0.3, seed = 73)

# rowVars(database$F.true15)
# var15 <- rowVars(database$F.true15)
# var15 <- rowVars(database$F.true15 * database$w.true15)
# var15 / sum(var15)

# Name the rows and columns
database$artificial.15 <- Create.Artificial.Data(num.freq = 30,
                                                 w0 = database$w.true15,
                                                 F0 = database$F.true15,
                                                 gamma0 = database$gamma.true15)
rownames(database$artificial.15) <- 
  colnames(database$artificial.15) <- database$artificial.name15

# Convert to binomial format
database$artificial.15 <- countsToBinomial(database$artificial.15)
database$artificial.15$n_ij <- database$artificial.15$win1 + database$artificial.15$win2
database$artificial.15$y_ij <- database$artificial.15$win1


########################### BEGIN Test #########################################

# J_skew.mat <- matrix(0, database$K0.true15, database$K0.true15)
# J_skew.mat[upper.tri(J_skew.mat)] <-  1
# J_skew.mat[lower.tri(J_skew.mat)] <- -1
# Gamma  <- database$gamma.true15 * J_skew.mat

# TRUE
# main.true <- crossprod(database$w.true15, database$F.true15[,1]-database$F.true15[,2])
# int.true <- crossprod(database$F.true15[,1], Gamma %*% database$F.true15[,2])
# main.true + int.true



# Generate artificial data for 30 entities
N <- 30
database$K0.true30 <- 4
database$artificial.name30 <- paste("Entity", 1:N)
database$w.true30 <- c(2.5, 2, 1.5, 1)
database$F.true30 <- generate.F.true(N = N, K = database$K0.true30, decay = 0.3, seed = 11)

# var30 <- rowVars(database$F.true30  * database$w.true30)
# var30 / sum(var30)

# Name the rows and columns
database$artificial.30 <- Create.Artificial.Data(num.freq = 10, 
                                                 w0 = database$w.true30,
                                                 F0 = database$F.true30)
rownames(database$artificial.30) <- 
  colnames(database$artificial.30) <- database$artificial.name30

# Convert to binomial format
database$artificial.30 <- countsToBinomial(database$artificial.30)
database$artificial.30$n_ij <- database$artificial.30$win1 + database$artificial.30$win2
database$artificial.30$y_ij <- database$artificial.30$win1



# Generate artificial data for 100 entities
N <- 100
database$K0.true100 <- 5
database$artificial.name100 <- paste("Entity", 1:N)
database$w.true100 <- c(3, 2.5, 2, 1.5, 1)
database$F.true100 <- generate.F.true(N = N, K = database$K0.true100, decay = 0.7, seed = 73)

# var100 <- rowVars(database$F.true100  * database$w.true100)
# var100 / sum(var100)

# Name the rows and columns
database$artificial.100 <- Create.Artificial.Data(num.freq = 30, 
                                                  w0 = database$w.true100,
                                                  F0 = database$F.true100)
rownames(database$artificial.100) <- 
  colnames(database$artificial.100) <- database$artificial.name100

# Convert to binomial format
database$artificial.100 <- countsToBinomial(database$artificial.100)
database$artificial.100$n_ij <- database$artificial.100$win1 + database$artificial.100$win2
database$artificial.100$y_ij <- database$artificial.100$win1


########################  END artificial database  #############################
