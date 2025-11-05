/**
 * Bayesian Bradley-Terry (BBT) Model (Wainer,2023) 
 */
data {
  int<lower = 1> num_entities;    // players
  int<lower = 1> num_pairs;       // pairs
  array[num_pairs] int<lower=1, upper = num_entities> player1;    // player 1 for pairs n
  array[num_pairs] int<lower=1, upper = num_entities> player2;    // player 2 for pairs n
  array[num_pairs] int<lower = 0> win1;    // number of wins for player 1
  array[num_pairs] int<lower = 0> win2;    // number of wins for player 2
}

transformed data{
  array[num_pairs] int<lower = 0> n;
  for (i in 1:num_pairs) n[i] = win1[i] + win2[i];
}

parameters {
  real<lower = 0.001> sigma; // scale of ability variation
  vector[num_entities] s;    // intrinsic ability for player n
}

transformed parameters{
  vector[num_entities] w;
  vector[num_pairs] pi;
  vector[num_pairs] M;

  w = exp(s);
  for (i in 1:num_pairs){
    M[i] = s[player1[i]] - s[player2[i]];
    pi[i] =  w[player1[i]] / (w[player1[i]] + w[player2[i]]);
    if (is_nan(pi[i])) pi[i] = 0.0 ;    // get some nan errors sometimes
  }
}

model {
  sigma ~ lognormal(0, 0.5);
  s ~ normal(0, sigma);
  win1 ~ binomial(n, pi);
}

generated quantities {
  array[num_pairs] int win1_rep;
  array[num_pairs] real log_lik;
  
  win1_rep = binomial_rng(n, pi) ;
  for (i in 1:num_pairs) log_lik[i] = binomial_lpmf( win1[i] | n[i] , pi[i]);
}


