#ifndef AUX_FUN
#define AUX_FUN 
void print_progress(const int &iter_ix, const int &warm_up, const int &iter_max, const int &chain){

  if(iter_max > 20){
      if((iter_ix) % (int)round(.1 * iter_max) == 0 || iter_ix == 1 || iter_ix == (warm_up + 1) ){
          int progress = (int)round(iter_ix * 100 / iter_max);
          std::string str = (iter_ix) <= (warm_up) ? "\t [Warmup]" : "\t [Sampling]";
          Rcpp::Rcout << "[Chain " << chain << "] Beginning of iteration: " << (iter_ix) << " / " << iter_max << " (" << progress << "%)" << str  << std::endl;
      }
  }
  else{
          int progress = (int)round(iter_ix * 100 / iter_max);
          std::string str = (iter_ix) <= (warm_up) ? "\t [Warmup]" : "\t [Sampling]";
          Rcpp::Rcout << "[Chain " << chain << "] Beginning of iteration: " << (iter_ix) << " / " << iter_max << " (" << progress << "%)" << str  << std::endl;
  }

}

//' @param n length of vector
//' @param rng random number generator
Eigen::VectorXd rnorm_vector(const int& n, std::mt19937& rng){
    
    Eigen::VectorXd out(n);
    std::normal_distribution<double> rnorm_(0,1);
    for(int i = 0; i < n; i++)
        out(i) = rnorm_(rng);
    
    return(out);
}

//' @param n length of vector
//' @param rng random number generator
Eigen::ArrayXd rnorm(const int& n, std::mt19937& rng){
    
    Eigen::ArrayXd out(n);
    std::normal_distribution<double> rnorm_(0,1);
    for(int i = 0; i < n; i++)
        out(i) = rnorm_(rng);
    
    return(out);
}

//' @param L rows of matrix 
//' @param K cols of matrix
//' @param rng random number generator
Eigen::ArrayXd rnorm(const int& L, const int& K, std::mt19937& rng){
    
    Eigen::ArrayXXd out(L,K);
    std::normal_distribution<double> rnorm_(0,1);
    for(int l = 0; l < L; l ++){
      for(int k =0; k< K; k ++)
        out(l,k) = rnorm_(rng);
    }

    return(out);
}

#endif
