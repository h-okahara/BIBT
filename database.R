#
# Sourcing this R file contains database used in main.R.
#
#########################  BEGIN import database  ##############################

# Create an environment to store all datasets
database <- new.env()



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









###---------------------------------------###
###        Create artificial data         ###
###---------------------------------------###

# Generate artificial data for 4 entities
w.true <- c(0.9, 1.2, 1.1)
F.true <- matrix(c(
  0.24, -0.65, 0.08, 0.34,
  0.13, -0.37, 0.06, 0.18,
  0.05, -0.15, 0.02, 0.07 
), nrow=3, byrow = TRUE)

# Name the rows and columns
database$artificial.4 <- Create.Artificial.Data(num.freq = 800, w0 = w.true, F0 = F.true)
rownames(database$artificial.4) <- 
  colnames(database$artificial.4) <- database$name.4

# Convert to binomial format using countsToBinomial
database$artificial.4 <- countsToBinomial(database$artificial.4)
database$artificial.4$n_ij <- database$artificial.4$win1 + database$artificial.4$win2
database$artificial.4$y_ij <- database$artificial.4$win1



# Generate artificial data for 6 entities
w.true <- c(3.4, 3, 2.6)
F.true <- matrix(c(
  0.06, 0.02, -0.14, -0.01, 0.00, 0.04,
  0.08, 0.01, -0.16, -0.03, 0.00, 0.04,
  0.01, 0.00, -0.04, -0.00, 0.00, 0.02
), nrow = 3, byrow = TRUE)

# Name the rows and columns
database$artificial.6 <- Create.Artificial.Data(num.freq = 800, w0 = w.true, F0 = F.true)
rownames(database$artificial.6) <- 
  colnames(database$artificial.6) <- database$name.6

# Convert to binomial format
database$artificial.6 <- countsToBinomial(database$artificial.6)
database$artificial.6$n_ij <- database$artificial.6$win1 + database$artificial.6$win2
database$artificial.6$y_ij <- database$artificial.6$win1




##########################  END import database  ###############################