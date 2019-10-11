Eigen::ArrayXXd dnorm(const int &J, const Eigen::ArrayXd &r, const Eigen::ArrayXXi &n_j, Eigen::ArrayXd &pi,
                      Eigen::ArrayXXd &w, Eigen::ArrayXXd &mu, double &tau){


  const int L = mu.rows();
  const int K = mu.cols();
  Eigen::ArrayXXd q(J,K);
  Eigen::ArrayXd tmp;

  for(int j = 0; j < J; j++){
    for(int k = 0; k < K; k ++){
      tmp = Eigen::ArrayXd::Zero(n_j(j,1));
      for(int l = 0; l < L; l ++){
        tmp = tmp + w(l,k) * pow(2 * M_PI * tau,-.5) * exp(- 1.0 / ( 2.0 * tau) * pow(r.segment(n_j(j,0),n_j(j,1)) - mu(l,k),2 ));
      }
      q(j,k) = pi(k) + (tmp).prod();
    }
  }

  return(q);
}


Eigen::ArrayXXd dnorm(const int &J, const Eigen::ArrayXd &r, const Eigen::ArrayXXi &n_j, Eigen::ArrayXd &pi,
                      Eigen::ArrayXXd &w, Eigen::ArrayXXd &mu, Eigen::ArrayXXd &tau){


  const int L = mu.rows();
  const int K = mu.cols();
  Eigen::ArrayXXd q(J,K);
  Eigen::ArrayXd tmp;

  for(int j = 0; j < J; j++){
    for(int k = 0; k < K; k ++){
      tmp = Eigen::ArrayXd::Zero(n_j(j,1));
      for(int l = 0; l < L; l ++){
        tmp = tmp + w(l,k) * pow(2 * M_PI * tau(l,k),-.5) * exp(- 1.0 / ( 2.0 * tau(l,k)) * pow(r.segment(n_j(j,0),n_j(j,1)) - mu(l,k),2 ));
      }
      q(j,k) = pi(k) + tmp.prod();
    }
  }

  return(q);
}


Eigen::ArrayXXd dnorm(const Eigen::ArrayXd &r, const Eigen::ArrayXXi &n_j, Eigen::ArrayXXd &w, Eigen::ArrayXXd &mu, double &tau, Eigen::ArrayXi& zeta){

  const int L = mu.rows();
  const int J = zeta.size();
  Eigen::ArrayXXd b(r.rows(),L);


  for(int l = 0; l < L; l++){
    for (int j = 0; j < J; j ++){
      for(int i = 0; i < n_j(j,1) ; i ++)
        b(n_j(j,0) +i, l) = log(w(l,zeta(j))) + R::dnorm(r(n_j(j,0) +i),mu(l,zeta(j)), sqrt(tau), true );
    }
  }

  return(b);

}

Eigen::ArrayXXd dnorm(const Eigen::ArrayXd &r, const Eigen::ArrayXXi &n_j, Eigen::ArrayXXd &w, Eigen::ArrayXXd &mu, Eigen::ArrayXXd &tau, Eigen::ArrayXi& zeta){

  const int L = mu.rows();
  const int J = zeta.size();
  Eigen::ArrayXXd b(r.rows(),L);


  for(int l = 0; l < L; l++){
    for (int j = 0; j < J; j ++){
      for(int i = 0; i < n_j(j,1) ; i ++)
        b(n_j(j,0) +i, l) = log(w(l,zeta(j))) + R::dnorm(r(n_j(j,0) +i),mu(l,zeta(j)), sqrt(tau(l,zeta(j))), true );
    }
  }

  return(b);

}
