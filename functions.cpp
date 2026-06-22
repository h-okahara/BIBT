//
// This C++ script relies on Rcpp and RcppArmadillo to implement the core Gibbs sampling routines
// for the Covariate-Assisted BIBT (CA-BIBT) model, as well as subroutines 
// for calculating Stochastic Transitivity (ST) classes.
//
#include <RcppArmadillo.h>
#include <pg.h>   // Pólya-Gamma sampler from pg package
// [[Rcpp::depends(RcppArmadillo, pg)]]

// Helper function for vectorized rgamma
arma::vec rgamma_vec(int n, double shape, const arma::vec& rate) {
  arma::vec out(n);
  for (int i = 0; i < n; ++i) {
    out[i] = R::rgamma(shape, 1.0 / rate[i]);
  }
  return out;
}

// Helper function for getting M_ij from M_vec
inline double get_M(int i, int j, const arma::vec& M_vec, int N) {
  if (i == j) return 0.0;
  if (i > j) return -get_M(j, i, M_vec, N);
  int idx = (N - 1)*i - (i - 1)*i/2 + j - i - 1;
  return M_vec[idx];
}

// Helper function for computing V_S, V_M, V_W
inline void calc_ST_indicators(const arma::vec& M_vec, int N, double threshold, double& V_S, double& V_M, double& V_W) {
  V_S = -arma::datum::inf;
  V_M = -arma::datum::inf;
  V_W = -arma::datum::inf;
  double M_th = (threshold == 0.5) ? 0.0 : std::log(threshold / (1.0 - threshold));
  
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      if (i == j) continue;
      double M_ij = get_M(i, j, M_vec, N);
      if (M_ij < M_th) continue; // M_ij >= M_th
      
      for (int k = 0; k < N; ++k) {
        if (i == k || j == k) continue;
        double M_jk = get_M(j, k, M_vec, N);
        if (M_jk < M_th) continue; // M_jk >= M_th
        
        double M_ki = get_M(k, i, M_vec, N);
        double C_ijk = M_ij + M_jk + M_ki;
        
        double min_val = std::min(M_ij, M_jk);
        double max_val = std::max(M_ij, M_jk);
        double sum_val = M_ij + M_jk;
        
        double v_s = C_ijk - min_val;
        double v_m = C_ijk - max_val;
        double v_w = C_ijk - sum_val;
        
        if (v_s > V_S) V_S = v_s;
        if (v_m > V_M) V_M = v_m;
        if (v_w > V_W) V_W = v_w;
      }
    }
  }
}



// =====================================================//
//    Computing the ST class for the match-up matrix    //
// =====================================================//

// [[Rcpp::export]]
Rcpp::List calc_ST_class_cpp(const Rcpp::NumericMatrix& M_mat, double threshold) {
  int N = M_mat.nrow();
  double V_S = R_NegInf;
  double V_M = R_NegInf;
  double V_W = R_NegInf;
  double M_th = (threshold == 0.5) ? 0.0 : std::log(threshold / (1.0 - threshold));
  
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      if (i == j) continue;
      double M_ij = M_mat(i, j);
      if (M_ij < M_th) continue; // M_ij >= M_th
      
      for (int k = 0; k < N; ++k) {
        if (i == k || j == k) continue;
        double M_jk = M_mat(j, k);
        if (M_jk < M_th) continue; // M_jk >= M_th
        
        double M_ki = M_mat(k, i);
        double C_ijk = M_ij + M_jk + M_ki;
        
        double min_val = std::min(M_ij, M_jk);
        double max_val = std::max(M_ij, M_jk);
        double sum_val = M_ij + M_jk;
        
        double v_s = C_ijk - min_val;
        double v_m = C_ijk - max_val;
        double v_w = C_ijk - sum_val;
        
        if (v_s > V_S) V_S = v_s;
        if (v_m > V_M) V_M = v_m;
        if (v_w > V_W) V_W = v_w;
      }
    }
  }
  
  // Decide ST class deterministically
  int pi_S = 0, pi_M = 0, pi_W = 0, pi_I = 0;
  if (V_S <= 0) {
    pi_S = 1;
  } else if (V_M <= 0) {
    pi_M = 1;
  } else if (V_W <= 0) {
    pi_W = 1;
  } else {
    pi_I = 1;
  }
  
  return Rcpp::List::create(
    Rcpp::Named("V_S") = V_S,
    Rcpp::Named("V_M") = V_M,
    Rcpp::Named("V_W") = V_W,
    Rcpp::Named("pi_S") = pi_S,
    Rcpp::Named("pi_M") = pi_M,
    Rcpp::Named("pi_W") = pi_W,
    Rcpp::Named("pi_I") = pi_I
  );
}




// =======================================================================//
//    Covariate-Assisted Bayesian Intransitive Bradley-Terry (CA-BIBT)    //
// =======================================================================//

// [[Rcpp::export]]
Rcpp::List CA_BIBT_Gibbs_cpp(int mcmc, int burn, int thin, 
                             const arma::vec& n_ij, const arma::vec& kappa, 
                             const arma::sp_mat& G, const arma::sp_mat& C_ast, 
                             const arma::mat& D_g, const arma::mat& D_g_t, const arma::mat& B_g,
                             const arma::mat& D_c, const arma::mat& D_c_t, const arma::mat& B_c,
                             const arma::mat& X_E, int num_entities, int num_pairs, int num_triplets, 
                             int dim_u, int dim_z, int dim_cov,
                             arma::vec u, double sigma_u, arma::vec z, 
                             arma::vec lambda, double tau, arma::vec nu, double xi,
                             arma::vec beta, double sigma_beta,
                             double a = 0.5, double b = 0.5,
                             double threshold = 0.5) 
{
  // Initial values
  arma::sp_mat C      = C_ast.t();
  arma::vec grad_res_flow = D_g * u;
  arma::vec curl_res_flow = arma::zeros<arma::vec>(num_pairs);
  if (dim_z > 0) {
    curl_res_flow = D_c * z;
  }
  arma::vec cov_flow = arma::zeros<arma::vec>(num_pairs);
  if (dim_cov > 0) {
    cov_flow = X_E.t() * beta;
  }
  arma::vec M_vec     = grad_res_flow + curl_res_flow + cov_flow;
  arma::vec M_vec_abs = arma::abs(M_vec);
  arma::vec omega     = pg::rpg_hybrid(n_ij, M_vec_abs);
  
  // Define matrices for posterior samples
  int mcmc_row = (mcmc-burn)/thin;
  arma::mat u_pos(mcmc_row, dim_u);
  arma::mat s_pos(mcmc_row, num_entities);
  arma::vec sigma_u_pos(mcmc_row);
  arma::mat z_pos(mcmc_row, dim_z);
  arma::mat Phi_pos(mcmc_row, num_triplets);
  arma::mat lambda_pos(mcmc_row, dim_z);
  arma::vec tau_pos(mcmc_row);
  arma::mat nu_pos(mcmc_row, dim_z);
  arma::vec xi_pos(mcmc_row);
  
  arma::mat beta_pos(mcmc_row, dim_cov);
  arma::vec sigma_beta_pos(mcmc_row);
  
  arma::mat grad_res_pos(mcmc_row, num_pairs);
  arma::mat grad_cov_pos(mcmc_row, num_pairs);
  arma::mat grad_pos(mcmc_row, num_pairs);
  arma::mat curl_res_pos(mcmc_row, num_pairs);
  arma::mat curl_cov_pos(mcmc_row, num_pairs);
  arma::mat curl_pos(mcmc_row, num_pairs);
  arma::mat cov_pos(mcmc_row, num_pairs);
  arma::mat M_pos(mcmc_row, num_pairs);
  
  // Flow Contribution Ratios
  arma::vec R_gr_pos(mcmc_row);
  arma::vec R_cr_pos(mcmc_row);
  arma::vec R_gx_pos(mcmc_row);
  arma::vec R_cx_pos(mcmc_row);
  arma::vec R_g_pos(mcmc_row);
  arma::vec R_c_pos(mcmc_row);
  arma::vec R_x_pos(mcmc_row);
  arma::vec R_x_g_pos(mcmc_row);
  arma::vec R_x_c_pos(mcmc_row);
  
  arma::mat LV_pos(mcmc_row, num_triplets);
  arma::vec V_S_pos(mcmc_row);
  arma::vec V_M_pos(mcmc_row);
  arma::vec V_W_pos(mcmc_row);
  
  double count_S = 0.0;
  double count_M = 0.0;
  double count_W = 0.0;
  double count_I = 0.0;
  
  int sample_idx = 0;
  
  //=======================   BEGIN MCMC sampling   =============================
  for (int iter = 1; iter <= mcmc; ++iter) {
    
    // Updating omega: sample omega from Pólya-Gamma distribution
    arma::vec M_vec_abs_loop = arma::abs(M_vec);
    omega = pg::rpg_hybrid(n_ij, M_vec_abs_loop);
    
    // Updating beta: d x 1 covariate coefficient vector
    if (dim_cov > 0) {
      double sigma_beta_sq = sigma_beta * sigma_beta;
      double inv_sigma_beta_sq = 1.0 / sigma_beta_sq;
      
      arma::mat X_E_omega = X_E.each_row() % omega.t(); 
      arma::mat Prec_beta = X_E_omega * X_E.t(); 
      Prec_beta.diag() += inv_sigma_beta_sq;
      Prec_beta = arma::symmatu(Prec_beta);
      arma::mat U_beta = arma::chol(Prec_beta); 
      arma::vec B_beta = X_E * (kappa - omega % (grad_res_flow + curl_res_flow));
      arma::vec tmp_beta = arma::solve(arma::trimatl(U_beta.t()), B_beta);
      arma::vec mu_beta = arma::solve(arma::trimatu(U_beta), tmp_beta);
      arma::vec v_beta = Rcpp::rnorm(dim_cov);
      arma::vec z_beta = arma::solve(arma::trimatu(U_beta), v_beta);
      beta = mu_beta + z_beta;
      cov_flow = X_E.t() * beta;
      
      // Updating sigma_beta
      double a_sigma_beta = (1.0 + dim_cov) / 2.0;
      double b_sigma_beta_rate = (1.0 + arma::dot(beta, beta)) / 2.0;
      sigma_beta = std::sqrt(1.0 / R::rgamma(a_sigma_beta, 1.0 / b_sigma_beta_rate));
      sigma_beta_sq = sigma_beta * sigma_beta;
      inv_sigma_beta_sq = 1.0 / sigma_beta_sq;
    } else {
      cov_flow.zeros(); 
    }
    
    // Updating u: unconstrained score vector (q_u x 1)
    double sigma_u_sq = sigma_u * sigma_u;
    double inv_sigma_u_sq = 1.0 / sigma_u_sq;

    arma::mat D_g_omega = D_g;
    D_g_omega.each_col() %= omega;
    arma::mat Prec_u = D_g_t * D_g_omega;
    Prec_u.diag() += inv_sigma_u_sq;
    Prec_u = 0.5 * (Prec_u + Prec_u.t());
    arma::mat U_u = arma::chol(Prec_u);
    arma::vec B_u = D_g_t * (kappa - omega % (curl_res_flow + cov_flow));
    arma::vec tmp_u = arma::solve(arma::trimatl(U_u.t()), B_u);
    arma::vec mu_u = arma::solve(arma::trimatu(U_u), tmp_u);
    arma::vec v_u = Rcpp::rnorm(dim_u);
    arma::vec z_u = arma::solve(arma::trimatu(U_u), v_u);
    u = mu_u + z_u;
    grad_res_flow = D_g * u;
    
    // Updating sigma_u
    double a_sigma_u = (1.0 + dim_u) / 2.0;
    double b_sigma_u_rate = (1.0 + arma::dot(u, u)) / 2.0;
    sigma_u = std::sqrt(1.0 / R::rgamma(a_sigma_u, 1.0 / b_sigma_u_rate));
    sigma_u_sq = sigma_u * sigma_u;
    inv_sigma_u_sq = 1.0 / sigma_u_sq;
    
    // Updating z: unconstrained curl weights (q_z x 1)
    if (dim_z > 0) {
      arma::mat D_c_omega = D_c;
      D_c_omega.each_col() %= omega;
      arma::mat Prec_z = D_c_t * D_c_omega;
      arma::vec Prec_z_prior = 1.0 / arma::pow(tau * lambda, 2);
      Prec_z.diag() += Prec_z_prior;
      Prec_z = 0.5 * (Prec_z + Prec_z.t());
      arma::mat U_z = arma::chol(Prec_z);
      arma::vec B_z = D_c_t * (kappa - omega % (grad_res_flow + cov_flow));
      arma::vec tmp_z = arma::solve(arma::trimatl(U_z.t()), B_z);
      arma::vec mu_z = arma::solve(arma::trimatu(U_z), tmp_z);
      arma::vec v_z = Rcpp::rnorm(dim_z);
      arma::vec z_z = arma::solve(arma::trimatu(U_z), v_z);
      z = mu_z + z_z;
      curl_res_flow = D_c * z;
      
      // Updating lambda
      arma::vec b_lambda_rate = 1.0 / nu + arma::pow(z, 2) / (2 * tau * tau);
      lambda = arma::sqrt(1.0 / rgamma_vec(dim_z, b+0.5, b_lambda_rate));
      
      // Updating tau
      double a_tau = (dim_z + 1.0) / 2.0;
      double S = arma::sum(arma::pow(z / lambda, 2));
      double b_tau_rate = S / 2.0 + 1.0 / xi;
      tau = std::sqrt(1.0 / R::rgamma(a_tau, 1.0 / b_tau_rate));
      
      // Updating nu
      arma::vec b_nu_rate = 1.0 + 1.0 / arma::pow(lambda, 2);
      nu = 1.0 / rgamma_vec(dim_z, a+b, b_nu_rate);
      
      // Updating xi
      double b_xi_rate = 1.0 + 1.0 / (tau * tau);
      xi = 1.0 / R::rgamma(1.0, 1.0 / b_xi_rate);
    } else {
      curl_res_flow.zeros(); 
    }
    
    // Updating Match-up vector and Each Flows
    M_vec = grad_res_flow + curl_res_flow + cov_flow;
    arma::vec grad_flow = (G * G.t() * M_vec) / num_entities;
    arma::vec curl_flow = M_vec - grad_flow;
    arma::vec grad_cov_flow = grad_flow - grad_res_flow;
    arma::vec curl_cov_flow = curl_flow - curl_res_flow;
    
    // Flow Contribution Ratios
    double grad_res_norm = arma::dot(grad_res_flow, grad_res_flow);
    double curl_res_norm = arma::dot(curl_res_flow, curl_res_flow);
    double grad_cov_norm = arma::dot(grad_cov_flow, grad_cov_flow);
    double curl_cov_norm = arma::dot(curl_cov_flow, curl_cov_flow);
    double M_norm = arma::dot(M_vec, M_vec);
    
    double R_gr = (M_norm > 0) ? (grad_res_norm / M_norm) : 0.0;
    double R_cr = (M_norm > 0) ? (curl_res_norm / M_norm) : 0.0;
    double R_gx = (M_norm > 0) ? (grad_cov_norm / M_norm) : 0.0;
    double R_cx = (M_norm > 0) ? (curl_cov_norm / M_norm) : 0.0;
    double R_g = R_gr + R_gx;
    double R_c = R_cr + R_cx;
    double R_x = R_gx + R_cx;
    double R_x_g = (R_g > 0) ? (R_gx / R_g) : 0.0;
    double R_x_c = (R_c > 0) ? (R_cx / R_c) : 0.0;

    arma::vec LV = C * curl_flow; // Local Velocity
    
    if (iter > burn && (iter-burn) % thin == 0) { 
      u_pos.row(sample_idx)         = u.t();
      s_pos.row(sample_idx)         = (B_g * u).t();
      sigma_u_pos[sample_idx]       = sigma_u;
      z_pos.row(sample_idx)         = z.t();
      Phi_pos.row(sample_idx)       = (B_c * z).t();
      lambda_pos.row(sample_idx)    = lambda.t();
      tau_pos[sample_idx]           = tau;
      nu_pos.row(sample_idx)        = nu.t();
      xi_pos[sample_idx]            = xi;
      beta_pos.row(sample_idx)      = beta.t();       
      sigma_beta_pos[sample_idx]    = sigma_beta;     
      
      grad_res_pos.row(sample_idx)  = grad_res_flow.t();
      grad_cov_pos.row(sample_idx)  = grad_cov_flow.t();
      grad_pos.row(sample_idx)      = grad_flow.t();
      curl_res_pos.row(sample_idx)  = curl_res_flow.t();
      curl_cov_pos.row(sample_idx)  = curl_cov_flow.t();
      curl_pos.row(sample_idx)      = curl_flow.t();
      cov_pos.row(sample_idx)       = cov_flow.t();
      M_pos.row(sample_idx)         = M_vec.t();
      R_gr_pos[sample_idx]          = R_gr;
      R_cr_pos[sample_idx]          = R_cr;
      R_gx_pos[sample_idx]          = R_gx;
      R_cx_pos[sample_idx]          = R_cx;
      R_g_pos[sample_idx]           = R_g;
      R_c_pos[sample_idx]           = R_c;
      R_x_pos[sample_idx]           = R_x;
      R_x_g_pos[sample_idx]         = R_x_g;
      R_x_c_pos[sample_idx]         = R_x_c;
      LV_pos.row(sample_idx)        = LV.t();
      
      // Store posterior counts for V_S, V_M, V_W
      double v_s_val, v_m_val, v_w_val;
      calc_ST_indicators(M_vec, num_entities, threshold, v_s_val, v_m_val, v_w_val);
      
      V_S_pos[sample_idx] = v_s_val;
      V_M_pos[sample_idx] = v_m_val;
      V_W_pos[sample_idx] = v_w_val;
      
      if (v_s_val <= 0) count_S += 1;
      if (v_m_val <= 0) count_M += 1;
      if (v_w_val <= 0) count_W += 1;
      else              count_I += 1;
      
      sample_idx++;
    }
  }
  //=======================   END MCMC sampling   ==============================
  
  // Compute posterior probabilities
  double pi_S = count_S / mcmc_row;
  double pi_M = count_M / mcmc_row;
  double pi_W = count_W / mcmc_row;
  double pi_I = count_I / mcmc_row;
  
  return Rcpp::List::create(Rcpp::Named("s")          = s_pos,
                            Rcpp::Named("u")          = u_pos,
                            Rcpp::Named("sigma_u")    = sigma_u_pos,
                            Rcpp::Named("z")          = z_pos,
                            Rcpp::Named("Phi")        = Phi_pos,
                            Rcpp::Named("lambda")     = lambda_pos,
                            Rcpp::Named("tau")        = tau_pos,
                            Rcpp::Named("nu")         = nu_pos,
                            Rcpp::Named("xi")         = xi_pos,
                            Rcpp::Named("beta")       = beta_pos,
                            Rcpp::Named("sigma_beta") = sigma_beta_pos,
                            Rcpp::Named("grad_res")   = grad_res_pos,
                            Rcpp::Named("grad_cov")   = grad_cov_pos,
                            Rcpp::Named("grad")       = grad_pos,
                            Rcpp::Named("curl_res")   = curl_res_pos,
                            Rcpp::Named("curl_cov")   = curl_cov_pos,
                            Rcpp::Named("curl")       = curl_pos,
                            Rcpp::Named("cov")        = cov_pos,
                            Rcpp::Named("M")          = M_pos,
                            Rcpp::Named("R_gr")       = R_gr_pos,
                            Rcpp::Named("R_cr")       = R_cr_pos,
                            Rcpp::Named("R_gx")       = R_gx_pos,
                            Rcpp::Named("R_cx")       = R_cx_pos,
                            Rcpp::Named("R_g")        = R_g_pos,
                            Rcpp::Named("R_c")        = R_c_pos,
                            Rcpp::Named("R_x")        = R_x_pos,
                            Rcpp::Named("R_x_g")      = R_x_g_pos,
                            Rcpp::Named("R_x_c")      = R_x_c_pos,
                            Rcpp::Named("LV")         = LV_pos,
                            Rcpp::Named("V_S")        = V_S_pos,
                            Rcpp::Named("V_M")        = V_M_pos,
                            Rcpp::Named("V_W")        = V_W_pos,
                            Rcpp::Named("pi_S")       = pi_S,
                            Rcpp::Named("pi_M")       = pi_M,
                            Rcpp::Named("pi_W")       = pi_W,
                            Rcpp::Named("pi_I")       = pi_I);
}
