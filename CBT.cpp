#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <pg.h>   // Pólya-Gamma sampler from pg package

// [[Rcpp::depends(RcppArmadillo, pg, RcppEigen)]]


// Helper function to build sparse matrix
Eigen::SparseMatrix<double> arma_to_eigen_sparse(const arma::sp_mat& A) {
  Eigen::SparseMatrix<double> B(A.n_rows, A.n_cols);
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(A.n_nonzero);
  
  for (arma::sp_mat::const_iterator iter = A.begin(); iter != A.end(); ++iter) {
    triplets.push_back(Eigen::Triplet<double>(iter.row(), iter.col(), *iter));
  }
  
  B.setFromTriplets(triplets.begin(), triplets.end());
  return B;
}


// Helper function for vectorized rgamma
// R's rgamma(n, shape, rate) corresponds to Rcpp::rgamma(shape, scale = 1.0/rate)
arma::vec rgamma_vec(int n, double shape, const arma::vec& rate) {
  arma::vec out(n);
  for (int i = 0; i < n; ++i) {
    out[i] = R::rgamma(shape, 1.0 / rate[i]);
  }
  return out;
}




//========================================//
//    Cyclic Bradley-Terry (CBT) Model    //
//========================================//

// [[Rcpp::export]]
Rcpp::List CBT_Gibbs_cpp(int mcmc, int burn, int thin, 
                         const arma::vec& n_ij, const arma::vec& kappa, 
                         const arma::sp_mat& G, const arma::sp_mat& C_ast, 
                         const arma::mat& H, const arma::mat& D_ast, const arma::mat& D_ast_t, 
                         int num_entities, int num_pairs, int num_triplets, int num_free,
                         arma::vec s, double sigma, arma::vec Phi, 
                         arma::vec lambda, double tau, arma::vec nu, double xi) 
  {
  // Initial values
  arma::vec weights   = arma::zeros<arma::vec>(num_free);
  arma::vec grad_flow = G * s;
  arma::vec curl_flow = C_ast * Phi;
  arma::vec M_vec     = grad_flow + curl_flow;
  arma::vec M_vec_abs = arma::abs(M_vec);
  arma::vec omega     = pg::rpg_hybrid(n_ij, M_vec_abs);
  
  // Define matrices for posterior samples
  int mcmc_row = (mcmc-burn)/thin;
  arma::mat s_pos(mcmc_row, num_entities);
  arma::mat weights_pos(mcmc_row, num_free);
  arma::mat Phi_pos(mcmc_row, num_triplets);
  arma::mat lambda_pos(mcmc_row, num_free);
  arma::vec tau_pos(mcmc_row);
  arma::mat nu_pos(mcmc_row, num_free);
  arma::vec xi_pos(mcmc_row);
  
  arma::mat grad_pos(mcmc_row, num_pairs);
  arma::mat curl_pos(mcmc_row, num_pairs);
  arma::mat M_pos(mcmc_row, num_pairs);
  
  double sigma_sq = sigma * sigma;
  double inv_sigma_sq = 1.0 / sigma_sq;
  
  int sample_idx = 0;
  //=======================   BEGIN MCMC sampling   =============================
  for (int iter = 1; iter <= mcmc; ++iter) {
    // -----------------------  BEGIN Updating  ---------------------------------
    // Updating omega: sample omega from Pólya-Gamma distribution
    arma::vec M_vec_abs_loop = arma::abs(M_vec);
    omega = pg::rpg_hybrid(n_ij, M_vec_abs_loop);
    
    
    // Updating s: num_entities×1 score vector
    arma::sp_mat G_omega = G;
    for (arma::sp_mat::iterator iter = G_omega.begin(); iter != G_omega.end(); ++iter) {
      *iter *= std::sqrt(omega(iter.row()));
    }
    Eigen::SparseMatrix<double> G_omega_sparse = arma_to_eigen_sparse(G_omega);
    Eigen::SparseMatrix<double> Prec_s = G_omega_sparse.transpose() * G_omega_sparse;
    Eigen::VectorXd Prec_s_prior(num_entities);
    Prec_s_prior.fill(inv_sigma_sq);
    Prec_s.diagonal() += Prec_s_prior;  // Compute the precision matrix
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver; // Sparse Cholesky decomposition
    solver.compute(Prec_s);
    
    // Solve for mean mu_s
    arma::vec B_s = G.t() * (kappa - omega % (C_ast * Phi));
    Eigen::Map<Eigen::VectorXd> B_s_sparse(B_s.memptr(), B_s.n_elem);
    Eigen::VectorXd mu_s = solver.solve(B_s_sparse); // Solve for mean mu_s
    
    // Sample v_s ~ N(0, Prec_s^{-1})
    arma::vec v_s = Rcpp::rnorm(num_entities);
    Eigen::Map<Eigen::VectorXd> v_s_sparse(v_s.memptr(), v_s.n_elem);
    Eigen::VectorXd z_s = solver.matrixU().solve(v_s_sparse);
    s = arma::vec(mu_s.data(), mu_s.size()) + arma::vec(z_s.data(), z_s.size());
    s = s - arma::mean(s); // Identification
    
    
    // Updating w: num.free × 1 weight vector
    arma::mat D_ast_omega = D_ast;
    D_ast_omega.each_col() %= omega;
    arma::mat Prec_likelihood = D_ast_t * D_ast_omega;
    arma::mat Prec_w = Prec_likelihood;
    arma::vec Prec_w_prior = 1.0 / arma::pow(tau * lambda, 2);
    Prec_w.diag() += Prec_w_prior;
    Prec_w = arma::symmatu(Prec_w);
    arma::mat U_w = arma::chol(Prec_w);
    arma::vec B_w = D_ast_t * (kappa - omega % (G * s));
    arma::vec tmp_w = arma::solve(arma::trimatl(U_w.t()), B_w);
    arma::vec mu_w = arma::solve(arma::trimatu(U_w), tmp_w);
    arma::vec v_w = Rcpp::rnorm(num_free);
    arma::vec z_w = arma::solve(arma::trimatu(U_w), v_w);
    weights = mu_w + z_w;
    Phi = H * weights;
    
    
    // Updating lambda: num.free×1 vector
    arma::vec b_lambda_rate = 1.0 / nu + arma::pow(weights, 2) / (2 * tau * tau);
    lambda = arma::sqrt(1.0 / rgamma_vec(num_free, 1.0, b_lambda_rate));
    
      
    // Updating tau:
    double a_tau = (num_free + 1.0) / 2.0;
    double S = arma::sum(arma::pow(weights / lambda, 2));
    double b_tau_rate = S / 2.0 + 1.0 / xi;
    tau = std::sqrt(1.0 / R::rgamma(a_tau, 1.0 / b_tau_rate));
    
    
    // Updating nu: num.free×1 vector
    arma::vec b_nu_rate = 1.0 + 1.0 / arma::pow(lambda, 2);
    nu = 1.0 / rgamma_vec(num_free, 1.0, b_nu_rate);
    
    
    // Updating xi:
    double b_xi_rate = 1.0 + 1.0 / (tau * tau);
    xi = 1.0 / R::rgamma(1.0, 1.0 / b_xi_rate);
    
    
    // Updating other parameters
    grad_flow = G * s;
    curl_flow = C_ast * Phi;
    M_vec = grad_flow + curl_flow;
    
    // ------------------------  END Updating  ----------------------------------
    if (iter > burn && (iter-burn) % thin == 0) { // Store posterior samples
      s_pos.row(sample_idx)       = s.t();
      weights_pos.row(sample_idx) = weights.t();
      Phi_pos.row(sample_idx)     = Phi.t();
      lambda_pos.row(sample_idx)  = lambda.t();
      tau_pos[sample_idx]         = tau; // Store in vector
      nu_pos.row(sample_idx)      = nu.t();
      xi_pos[sample_idx]          = xi; // Store in vector
      grad_pos.row(sample_idx)    = grad_flow.t();
      curl_pos.row(sample_idx)    = curl_flow.t();
      M_pos.row(sample_idx)       = M_vec.t();
      
      sample_idx++;
    }
  }
  //=======================   END MCMC sampling   ==============================
  
  return Rcpp::List::create(Rcpp::Named("s")        = s_pos,
                            Rcpp::Named("weights")  = weights_pos,
                            Rcpp::Named("Phi")      = Phi_pos,
                            Rcpp::Named("lambda")   = lambda_pos,
                            Rcpp::Named("tau")      = tau_pos,
                            Rcpp::Named("nu")       = nu_pos,
                            Rcpp::Named("xi")       = xi_pos,
                            Rcpp::Named("grad")     = grad_pos,
                            Rcpp::Named("curl")     = curl_pos,
                            Rcpp::Named("M")        = M_pos);
}




//==========================================//
//    Bayesian Bradley-Terry (BBT) Model    //
//==========================================//

// [[Rcpp::export]]
Rcpp::List BBT_Gibbs_cpp(int mcmc, int burn, int thin, 
                         const arma::vec& n_ij, const arma::vec& kappa, 
                         const arma::sp_mat& G,
                         int num_entities, int num_pairs, 
                         arma::vec s, double sigma) 
{
  // Initial values
  arma::vec M_vec     = G * s;
  arma::vec M_vec_abs = arma::abs(M_vec);
  arma::vec omega     = pg::rpg_hybrid(n_ij, M_vec_abs);
  
  // Define matrices for posterior samples
  int mcmc_row = (mcmc-burn)/thin;
  arma::mat s_pos(mcmc_row, num_entities);
  arma::mat M_pos(mcmc_row, num_pairs);
  
  double sigma_sq = sigma * sigma;
  double inv_sigma_sq = 1.0 / sigma_sq;
  
  int sample_idx = 0;
  //=======================   BEGIN MCMC sampling   =============================
  for (int iter = 1; iter <= mcmc; ++iter) {
    // -----------------------  BEGIN Updating  ---------------------------------
    // Updating omega: sample omega from Pólya-Gamma distribution
    arma::vec M_vec_abs_loop = arma::abs(M_vec);
    omega = pg::rpg_hybrid(n_ij, M_vec_abs_loop);
    
    
    // Updating s: num_entities×1 score vector
    arma::sp_mat G_omega = G;
    for (arma::sp_mat::iterator iter = G_omega.begin(); iter != G_omega.end(); ++iter) {
      *iter *= std::sqrt(omega(iter.row()));
    }
    Eigen::SparseMatrix<double> G_omega_sparse = arma_to_eigen_sparse(G_omega);
    Eigen::SparseMatrix<double> Prec_s = G_omega_sparse.transpose() * G_omega_sparse;
    Eigen::VectorXd Prec_s_prior(num_entities);
    Prec_s_prior.fill(inv_sigma_sq);
    Prec_s.diagonal() += Prec_s_prior;  // Compute the precision matrix
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver; // Sparse Cholesky decomposition
    solver.compute(Prec_s);
    
    // Solve for mean mu_s
    arma::vec B_s = G.t() * kappa;
    Eigen::Map<Eigen::VectorXd> B_s_sparse(B_s.memptr(), B_s.n_elem);
    Eigen::VectorXd mu_s = solver.solve(B_s_sparse); // Solve for mean mu_s
    
    // Sample v_s ~ N(0, Prec_s^{-1})
    arma::vec v_s = Rcpp::rnorm(num_entities);
    Eigen::Map<Eigen::VectorXd> v_s_sparse(v_s.memptr(), v_s.n_elem);
    Eigen::VectorXd z_s = solver.matrixU().solve(v_s_sparse);
    s = arma::vec(mu_s.data(), mu_s.size()) + arma::vec(z_s.data(), z_s.size());
    s = s - arma::mean(s); // Identification
    

    // Updating other parameters
    M_vec = G * s;
    
    // ------------------------  END Updating  ----------------------------------
    if (iter > burn && (iter-burn) % thin == 0) { // Store posterior samples
      s_pos.row(sample_idx) = s.t();
      M_pos.row(sample_idx) = M_vec.t();
      
      sample_idx++;
    }
  }
  //=======================   END MCMC sampling   ==============================
  
  return Rcpp::List::create(Rcpp::Named("s")    = s_pos,
                            Rcpp::Named("grad") = M_pos,
                            Rcpp::Named("M")    = M_pos);
}