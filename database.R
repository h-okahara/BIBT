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
# plot.network(database$citations.6, weight = "prop", layout = "fr")


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
# plot.network(database$citations.8, weight = "prop", layout = "fr")



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
# g.9 <- plot.network(database$citations.9, weight = "prop", layout = "fr")


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
  #F.true[1, ] <- sort(F.true[1, ], decreasing = TRUE)
  #F.true[2, ] <- sort(F.true[2, ], decreasing = FALSE)
  return(F.true)
}




###-----------------------###
###    Generate M.true    ###
###-----------------------###

generate.M.true <- function(N = NULL, K = NULL, F.worths = NULL, w = NULL, gamma = NULL) {
  ## Preparation
  m <- floor(K/2)
  stopifnot(length(gamma) == m)
  pairs <- t(combn(N, 2))
  J.mat <- matrix(c(0, 1, -1, 0), nrow = 2, byrow = TRUE)
  
  # build block-diagonal Gamma (K x K)
  Gamma <- matrix(0, K, K)
  if (m > 0) {
    for (r in seq_len(m)) {
      idx <- (2*r - 1):(2*r)
      Gamma[idx, idx] <- gamma[r] * J.mat
    }
  }
  
  # transitive term: Phi %*% w
  Phi <- t(F.worths[, pairs[, 1], drop = FALSE] - F.worths[, pairs[, 2], drop = FALSE])
  trans.true <- as.vector(Phi %*% w)
  
  # intransitive term: for each pair (i,j), f_i^T Gamma f_j
  f_i <- F.worths[, pairs[, 1], drop = FALSE]
  f_j <- F.worths[, pairs[, 2], drop = FALSE]
  int.true <- as.vector(colSums(f_i * (Gamma %*% f_j)))
  
  # match-up function: M.true = trans.true + int.true
  M.true <- trans.true + int.true
  result <- cbind(trans.true, int.true, M.true)
  colnames(result) <- c("trans", "int", "M")
  return(result)
}





###--------------------------###
###    Plot latent score     ###
###--------------------------###

plot.scores <- function(F.worths = NULL, w = NULL, point.labels = NULL, add.arrows = FALSE,
                        point.size = 1, text.vjust = -0.8, seed.jitter = 73) {
  if (is.null(dim(F.worths)) || nrow(F.worths) != 2) {
    stop("F.worths must be a 2 × N matrix (K = 2).")
  }
  N <- ncol(F.worths)
  df <- data.frame(x = F.worths[1, ],
                   y = F.worths[2, ],
                   name = if (is.null(point.labels)) paste0("f", 1:N) else point.labels)

  if (!is.null(w)) {
    w <- w / norm(w, type = "2")
    df <- df %>%
      rowwise() %>%
      mutate(
        proj.scalar = c(x, y) %*% w,
        x.proj = proj.scalar * w[1],
        y.proj = proj.scalar * w[2]
      ) %>%
      ungroup()
  }
  
  p <- ggplot(df, aes(x, y)) +
    theme_bw() +
    coord_equal(xlim = c(-2, 2), ylim = c(-2, 2)) +
    xlab(expression(f[1])) +
    ylab(expression(f[2])) +
    geom_hline(yintercept = 0, colour = "grey75") +
    geom_vline(xintercept = 0, colour = "grey75") +
    geom_point(size = point.size, colour = "steelblue")
  
  if (add.arrows && !is.null(w)) {
    line.range <- 4
    line.data <- data.frame(
      x = -line.range * w[1],
      y = -line.range * w[2],
      xend = line.range * w[1],
      yend = line.range * w[2]
    )
    
    p <- p +
      geom_segment(aes(xend = x.proj, yend = y.proj),
                   arrow = arrow(length = unit(0.1, "cm")),
                   colour = "darkorange", alpha = 0.6) +
      geom_point(aes(x = x.proj, y = y.proj),
                 colour = "firebrick", size = 2) +
      geom_segment(data = line.data,
                   aes(x = x, y = y, xend = xend, yend = yend),
                   inherit.aes = FALSE, colour = "black", linetype = "dashed")
  }

  if (!is.null(point.labels)) {
    p <- p +
      geom_text(aes(label = name), vjust = text.vjust, size = 3)
  }
  return(p)
}




###--------------------------###
###    Plot match network    ###
###--------------------------###

plot.network <- function(bin_df, weight = c("diff", "prop"), edge.label = FALSE,  
                         layout = c("fr", "circle"), tie_mode = c("skip", "thin")) 
  {
  ## Preparation
  weight <- match.arg(weight)
  layout <- match.arg(layout)
  tie_mode <- match.arg(tie_mode)
  
  ## Setting nodes and edges
  nodes_df <- data.frame(name = sort(unique(c(bin_df$player1, bin_df$player2))))
  edges_df <- bin_df %>%
    mutate(
      is_tie = (win1 == win2),
      winner = if_else(win1 > win2, player1, player2),
      loser  = if_else(win1 > win2, player2, player1),
      w_win  = pmax(win1, win2),
      w_lose = pmin(win1, win2),
      metric = case_when(
        weight == "diff"  ~ w_win - w_lose,
        weight == "prop"  ~ w_win / w_lose
      ),
      label = paste(w_win, w_lose, sep = "-")
    ) 
  if (tie_mode == "skip") {
    edges_df <- edges_df %>% filter(!is_tie)
  }
  edges_df <- edges_df %>%
    select(
      from = winner,
      to   = loser,
      metric,
      label
    )
  
  ## Define graph object
  g <- graph_from_data_frame(vertices = nodes_df, d = edges_df, directed = TRUE)
  if (length(unique(E(g)$metric)) > 1) {
    E(g)$width <- rescale(E(g)$metric, to = c(0.5, 3)) # scaling width of all edges
  } else {
    E(g)$width <- 3
  }
  layout_coords <- switch(layout,
                          fr     = layout_with_fr(g),
                          circle = layout_in_circle(g))
  
  ## Detect cyclic structures and highlight them
  scc <- components(g, mode = "strong")
  memb <- scc$membership
  csize <- scc$csize
  eH <- as.integer(head_of(g, E(g)))  # 辺の終点
  eT <- as.integer(tail_of(g, E(g)))  # 辺の始点
  same_scc <- memb[eH] == memb[eT]
  scc_gt1  <- csize[memb[eH]] > 1
  loop_e   <- which_loop(g)           # 自己ループは常にサイクル
  on_cycle <- (same_scc & scc_gt1) | loop_e
  E(g)$color <- rgb(0.2, 0.5, 0.9, 1)
  E(g)[on_cycle]$color <- rgb(1, 0.1, 0.3, 1)

  ## Plot network graph
  plot(g,
       layout = layout_coords,
       vertex.size = 30,
       vertex.color = "grey95",
       vertex.frame.color = "grey40",
       vertex.label.color = "grey10",
       edge.width = E(g)$width,
       edge.color = E(g)$color,
       edge.arrow.size = 0.8,
       edge.curved = 0.1,
       edge.label = if (edge.label) E(g)$label else NA,
       edge.label.color = "grey20",
       main = sprintf("Pairwise Win Network (weight: %s)", weight)
  )
  return(g)
}


###------------------------------###
###    Create artificial data    ###
###------------------------------###

Create.Artificial.Data <- function(num.freq = NULL, w0 = NULL, F0 = NULL, 
                                   gamma0 = NULL, seed = 73) {
  N <- ncol(F0) # number of entities
  K <- nrow(F0) # dimensionality of each entity
  
  ## Generate intransitive matrix \Gamma
  m <- floor(K/2)
  J.mat <- matrix(c(0, 1, -1, 0), nrow = 2, byrow = TRUE)
  Gamma.mat <- matrix(0, K, K)
  if (m > 0) {
    for (r in seq_len(m)) {
      idx <- (2*r - 1):(2*r)
      Gamma.mat[idx, idx] <- gamma0[r] * J.mat
    }
  }
  
  ## Initiate an N×N matrix storing results (diagonal is NA)
  result <- matrix(NA_integer_, nrow = N, ncol = N)
  rownames(result) <- colnames(result) <- 
    if (!is.null(colnames(F0))) colnames(F0) else paste0("Entity", 1:N)
  
  ## Simulate for each pair (i,j)
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      M_ij <- as.numeric(crossprod(w0, F0[, i] - F0[, j]) 
                         + crossprod(F0[, i], Gamma.mat %*% F0[, j]))
      p_ij <- 1 / (1 + exp(-M_ij))
      win.freq <- rbinom(1, size = num.freq, prob = p_ij)
      result[i, j] <- win.freq
      result[j, i] <- num.freq - win.freq
    }
  }
  return(result)
}




## Generate artificial data for 7 entities
N <- 7
database$K.true7 <- 2
database$artificial.name7 <- paste("Entity", 1:N)
database$gamma.true7 <- 2
database$w.true7 <- c(1, 1)
database$F.true7 <- matrix(0, database$K.true7, N)
database$F.true7[,1] = c(-2,0)
database$F.true7[,2] = c(-1,0)
database$F.true7[,3] = c(-1/2, sqrt(3))
database$F.true7[,4] = c(-1/2, -sqrt(3))
database$F.true7[,5] = c(-1/2, 0)
database$F.true7[,6] = c(1/2, 0)
database$F.true7[,7] = c(2, 0)
database$M.true7 <- generate.M.true(N, K = database$K.true7, 
                                    F.worths = database$F.true7, 
                                    w = database$w.true7,
                                    gamma = database$gamma.true7)

## Name the rows and columns
database$artificial.7 <- Create.Artificial.Data(num.freq = 100,
                                                w0 = database$w.true7,
                                                F0 = database$F.true7,
                                                gamma0 = database$gamma.true7)
rownames(database$artificial.7) <- 
  colnames(database$artificial.7) <- database$artificial.name7

## Convert to binomial format
database$artificial.7 <- countsToBinomial(database$artificial.7)
database$artificial.7$n_ij <- database$artificial.7$win1 + database$artificial.7$win2
database$artificial.7$y_ij <- database$artificial.7$win1

## Plot scores and network
# plot.scores(database$F.true7, database$w.true7, point.labels = database$artificial.name7, add.arrows = TRUE)
# g.true.7 <- plot.network(database$artificial.7, weight = "prop", layout = "fr")




# Generate artificial data for 10 entities
N <- 10
database$K.true10 <- 2
database$artificial.name10 <- paste("Entity", 1:N)
database$gamma.true10 <- 2.1
database$w.true10 <- c(1, 1)
# database$w.true10 <- database$w.true10 / norm(database$w.true10, type = "2")
database$F.true10 <- generate.F.true(N = N, K = database$K.true10, decay = 0.6, seed = 73)
database$M.true10 <- generate.M.true(N, K = database$K.true10, 
                                     F.worths = database$F.true10, 
                                     w = database$w.true10, 
                                     gamma = database$gamma.true10)

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

# Plot database$F.true10
# plot.scores(database$F.true10, database$w.true10, point.labels = database$artificial.name10, add.arrows = TRUE)
# g.true.10 <- plot.network(database$artificial.10, weight = "prop", layout = "fr")



# Generate artificial data for 15 entities
N <- 15
database$K.true15 <- 4
database$artificial.name15 <- paste("Entity", 1:N)
database$gamma.true15 <- c(1, 0.5)
database$w.true15 <- rep(1, database$K.true15)
# database$w.true15 <- database$w.true15 / norm(database$w.true15, type = "2")
database$F.true15 <- generate.F.true(N = N, K = database$K.true15, decay = 1/3, seed = 73)
database$M.true15 <- generate.M.true(N, K = database$K.true15, 
                                     F.worths = database$F.true15, 
                                     w = database$w.true15, 
                                     gamma = database$gamma.true15)

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


########################  END artificial database  #############################
