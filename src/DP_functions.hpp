#include "beta_rng.hpp"

Eigen::ArrayXXd initialize_mu(const int& L, const int& K, const double& mu_0,
                             const double& kappa_0, std::mt19937& rng){

  Eigen::ArrayXXd out(L,K);
  std::normal_distribution<double> rnorm(mu_0,sqrt(kappa_0));

  for(int l = 0; l < L; l ++){
    for( int k = 0; k < K; k ++)
      out(l,k) = rnorm(rng);
  }
  return(out);
}


Eigen::ArrayXXd initialize_tau(const int &L, const int &K, const double &sigma_0, const int &nu_0){

  Eigen::ArrayXXd out(L,K);
  for(int l = 0; l < L; l ++){
    for(int k = 0; k < K; k++)
      out(l,k) = (sigma_0 * nu_0) / R::rchisq(nu_0);
  }

  return(out);
}

//' returns stick breaking atoms (not weights) for a vector with constant beta parameters alpha, beta across all vector elements
//' @param n length of vector
//' @param alpha
Eigen::ArrayXd stick_break(const int& n, const double& alpha, const double& beta, std::mt19937& rng){

    Eigen::ArrayXd v(n);
    sftrabbit::beta_distribution<> rbeta(alpha,beta);
    for(int i = 0; i < (n-1); i++)
        v(i) = rbeta_log(alpha,beta,rng);
    v(n-1) = 0;
    return(v);
}

//' stick breaking atoms (not weights) for a vector with alpha fixed at 1 and constant beta across vector elements
//' @param n length of vector
//' @param beta parameter for stick breaking process
//' @param rng random number generator engine
Eigen::ArrayXd stick_break(const int n,const double beta, std::mt19937& rng){

    Eigen::ArrayXd v(n);
    for(int i = 0; i < (n-1); i++)
        v(i) = rbeta_log_const(1,beta,rng);
    v(n-1) = 0.0;
    return(v);
}

//' stick breaking atoms (not weights) for a vector with alpha,beta variable across vector elements
//' n length of vector
//' alpha vector of alpha parameters for posterior beta distribution
//' beta vector of beta parameters for posterior beta distribution
Eigen::ArrayXd stick_break(const int n, Eigen::ArrayXd& alpha,Eigen::ArrayXd& beta, std::mt19937& rng){

    Eigen::ArrayXd v(n);
    v = Eigen::ArrayXd::Zero(n);
    for(int i = 0; i < (n-1); i++){
        v(i) = rbeta_log(alpha(i),beta(i),rng);
    }
    v(n-1) = 0.0;

    return(v);
}

Eigen::ArrayXd stick_break_print(const int n, Eigen::ArrayXd& alpha,Eigen::ArrayXd& beta, std::mt19937& rng){

    Eigen::ArrayXd v(n);
    v = Eigen::ArrayXd::Zero(n);
    for(int i = 0; i < (n-1); i++){
        v(i) = rbeta_log_print(alpha(i),beta(i),rng);
		Rcpp::Rcout << "v(i): " << v(i) << std::endl;
    }
    v(n-1) = 0.0;

    return(v);
}

//' stick breaking atoms (not weights) for a vector with constant alpha (1) and beta across atoms - returns a matrix of size L x K
//' @param
Eigen::ArrayXXd stick_break(const int& rows, const int& cols, const double& beta, std::mt19937& rng){

    Eigen::ArrayXXd out(rows,cols);

    sftrabbit::beta_distribution<> rbeta(1,beta);
    for(int row_ix = 0; row_ix < rows; row_ix ++){
        for(int col_ix = 0; col_ix < cols; col_ix ++)
            out(row_ix,col_ix) = row_ix == (rows-1) ? 0.0 : rbeta_log(1,beta,rng);
    }

    return(out);
}

//' stick breaking
Eigen::ArrayXXd stick_break(Eigen::ArrayXXd& alpha, Eigen::ArrayXXd& beta, std::mt19937& rng){

    const int rows = alpha.rows();
    const int cols = alpha.cols();
    Eigen::ArrayXXd out(rows,cols);

    for(int row_ix = 0; row_ix < rows; row_ix ++){
        for(int col_ix = 0; col_ix < cols; col_ix ++){
            out(row_ix,col_ix) = row_ix == (rows-1) ? 0.0 : rbeta_log(alpha(row_ix,col_ix),beta(row_ix,col_ix), rng);
        }
    }

    return(out);
}


Eigen::ArrayXd stick_break_weights(const int n, Eigen::ArrayXd& v){

    Eigen::ArrayXd w(n);
    for(int i = 0; i < (n-1); i++)
        w(i) = i == 0 ? v(i) : v(i) * (Eigen::ArrayXd::Ones(i) - v.head(i)).prod();

    w(n-1) =  v(n-1) * (1 - v.head(n-1)).prod();
    return(w);
}

Eigen::ArrayXd stick_break_weights(Eigen::ArrayXd& v){

    const int n = v.rows();
    Eigen::ArrayXd w(n);
    for(int i = 0; i < (n-1); i++)
        w(i) = i == 0 ? v(i) : v(i) * (1 - v.head(i)).prod();

    w(n-1) =  v(n-1) * (1 - v.head(n-1)).prod();
    return(w);
}

Eigen::ArrayXXd stick_break_weights(Eigen::ArrayXXd& u){

    const int rows = u.rows();
    const int cols = u.cols();
    Eigen::ArrayXXd w(rows,cols);


    for(int col_ix = 0; col_ix < cols; col_ix ++){
        for(int row_ix = 0; row_ix < rows ; row_ix ++)
            w(row_ix,col_ix) =  row_ix == 0 ? u(row_ix,col_ix) : u(row_ix,col_ix) * (1 - u.block(0,col_ix,row_ix,1)).prod();
    }

    return(w);
}
