// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
#include <random>
#include "support/auxiliary_functions.hpp"
#include "support/DP_functions.hpp"
#include "support/NDP_functions.hpp"
#include<vector>
// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

//' Estimate the nonhomgogenous poisson process intensity function from grouped data
//' 
//' @template bda 
//'
//' @param r vector of distances associatd with different BEFs
//' @param n_j matrix of integers denoting the start and length of each observations associated BEF distances
//' @param d a 1D grid of positive real values over which the differing intensities are evaluated
//' @param L component truncation number
//' @param K intensity cluster truncation number
//' @param J number of rows in r matrix; number of groups
//' @param mu_0 normal base measure prior mean
//' @param kappa_0 normal base measure prior variance scale 
//' @param nu_0 inverse chi sqaure base measure prior degrees of freedom
//' @param sigma_0 inverse chi square base measure prior scale  
//' @param a_alpha hyperparameter for alpha gamma prior
//' @param b_alpha scale hyperparameter for alpha gamma prior
//' @param a_rho hyperparameter for rho gamma prior
//' @param b_rho scale hyperparameter for rho gamma prior
//' @param iter_max total number of iterations for which to run sampler
//' @param warm_up number of iterations for which to burn-in or "warm-up" sampler
//' @param thin number of iterations to thin by
//' @param seed integer with which to initialize random number generator
//' @param chain integer chain label
//' @param num_posterior_samples the total number of posterior samples after burn in 
//' @seealso the conjugate normal parameterization in the reference below
//'
// [[Rcpp::export]]
Rcpp::List nd_nhpp_fit(
        const Eigen::ArrayXd& r,
        const Eigen::MatrixXi& n_j,
        const Eigen::ArrayXd& d,
        const int& L,
        const int& K,
        const int& J,
        const double& mu_0,
        const double& kappa_0,
        const int& nu_0,
        const double& sigma_0,
        const double& a_alpha,
        const double& b_alpha,
        const double& a_rho,
        const double& b_rho,
        const int& iter_max,
        const int& warm_up,
        const int& thin,
        const int& seed,
        const int& chain,
        const int& num_posterior_samples
        )
{
    // set seed
    std::mt19937 rng;
    rng = std::mt19937(seed);

    //create sample containers
#include "support/create_sample_containers.hpp"

#include "support/during_sampling_containers.hpp"

    // create rng's
    std::gamma_distribution<double> rgam_alpha(a_alpha,b_alpha);
    std::gamma_distribution<double> rgam_rho(a_rho,b_rho);
    std::uniform_real_distribution<double> runif(0,1);
    std::normal_distribution<double> rnorm(0,1);

    //initialize concentration parameters
    alpha = rgam_alpha(rng);
    rho = rgam_rho(rng);
	// sample priors
#include "support/sample_concentration_priors.hpp"

    // initialize component weights
    u = stick_break(L,K,alpha,rng);
    v = stick_break(K,rho,rng);
    w = stick_break_weights(u);
    pi = stick_break_weights(v);

    mu = initialize_mu(L,K,mu_0,kappa_0,rng);
    tau = initialize_tau(L,K,sigma_0,nu_0);


    Rcpp::Rcout << "Beginning Sampling" << std::endl;
    Rcpp::Rcout << "----------------------------------------------------------------------" << std::endl;
    for(int iter_ix = 1; iter_ix <= iter_max; iter_ix ++){
      print_progress(iter_ix,warm_up,iter_max,chain);

      //calculate cluster probabilities
      q = dnorm(J,r,n_j,pi,w,mu,tau);

      for(int j = 0; j < J; j ++){
            probs = q.row(j);
            std::discrete_distribution<int> d(probs.data(),probs.data()+probs.size());
            iter_cluster_assignment(j) = d(rng);
        }

      // calculate w/in cluster component probabilities
        b = dnorm(r,n_j,w,mu,tau,iter_cluster_assignment);

        // assign distances to components within clusters
        component_count = Eigen::ArrayXXi::Zero(L,K);
        for(int j = 0 ; j < J; j ++){
            for(int i = 0 ; i < n_j(j,1) ; i++){
                for(int l = 0 ; l< L; l++)
                    prob(l) = b(n_j(j,0) + i,l);
                std::discrete_distribution<int> d(prob.data(),prob.data() + prob.size());
                iter_component_assignment(n_j(j,0) +i) = d(rng);
                component_count(iter_component_assignment(n_j(j,0)+i),iter_cluster_assignment(j)) += 1;
            }
        }

        // draw samples from intensity cluster parameters

        for(int k = 0; k < K; k++)
            cluster_count(k) = (iter_cluster_assignment == k).count();

        for(int k = 0; k < K; k++){
            v_posterior_beta_alpha(k) = 1 + cluster_count(k);
            v_posterior_beta_beta(k) = alpha + cluster_count.tail(K-k-1).sum();
        }

        v = stick_break(K, v_posterior_beta_alpha,v_posterior_beta_beta,rng);
        pi = stick_break_weights(v);



        for(int l = 0; l < L; l ++){
            for(int k = 0; k < K; k++){
                u_posterior_beta_alpha(l,k) = 1 + component_count(l,k);
                u_posterior_beta_beta(l,k) = rho + ((component_count.col(k)).tail(K-k-1)).sum();
            }
        }

        u = stick_break(u_posterior_beta_alpha,u_posterior_beta_beta,rng);
        w = stick_break_weights(u);

        // calculate mu hat (sq) sums

        ycount = Eigen::ArrayXXd::Zero(L,K);
        ycount_sq = Eigen::ArrayXXd::Zero(L,K);
        for(int l = 0 ; l < L; l ++){
          for(int k = 0; k < K; k++){
            for(int j = 0; j < J; j++){
              for(int i = 0; i < n_j(j,1); i ++){
                if(iter_cluster_assignment(j) == k && iter_component_assignment(n_j(j,0)+i) == l){
                  ycount(l,k) += r(n_j(j,0)+i);
                  ycount_sq(l,k) += pow(r(n_j(j,0)+i),2);
                }
              }
            }
          }
        }


        //sample mu via conjugacy

        for(int l = 0; l < L; l ++){
         for(int k = 0; k < K; k++){
           if(component_count(l,k) == 0){
             tau(l,k) = nu_0 * sigma_0 / R::rchisq(nu_0);
             mu(l,k) = rnorm(rng) * sqrt(tau(l,k) / kappa_0 ) + mu_0;
           }
           else{
             s_n  = nu_0 * sigma_0 + (ycount_sq(l,k) - (pow(ycount(l,k),2) / component_count(l,k) ) ) + (kappa_0 * component_count(l,k) / (kappa_0 + component_count(l,k))) * pow( (ycount(l,k) / component_count(l,k) - mu_0 ),2);
             tau(l,k) = s_n / R::rchisq(nu_0 + component_count(l,k));
             mu_n = ( (kappa_0  / tau(l,k)) * mu_0   + ycount(l,k) / tau(l,k)  ) / ( kappa_0 / tau(l,k) + component_count(l,k) / tau(l,k));
             s_n = 1.0 /( (kappa_0 /tau(l,k)) + (component_count(l,k) / tau(l,k)) ) ;
             mu(l,k) =  rnorm(rng) * sqrt(s_n) + mu_n;
           }
         }
       }

        // sample concentration parameters
        posterior_b_alpha = 1.0 / b_alpha - (log(1-v.array())).head(K-1).sum();
        posterior_b_rho =  1.0 / b_rho - log(1-u.block(0,0,L-1,K).array()).matrix().colwise().sum().sum();
        std::gamma_distribution<double> rgam_alpha(posterior_a_alpha, 1.0 / posterior_b_alpha);
        std::gamma_distribution<double> rgam_rho(posterior_a_rho, 1.0 / posterior_b_rho);
        alpha = rgam_alpha(rng);
        rho = rgam_rho(rng);


        // calculate adjacency matrix and store samples
        if((iter_ix > warm_up) && (iter_ix % thin ==0)){
          for(int j = 0; j < J; j++){
            for(int j_ = 0; j_ < j; j_++)
              cluster_matrix(j,j_) += (iter_cluster_assignment(j) == iter_cluster_assignment(j_)) ? 1: 0;
          }
          for(int k =0; k < K; k++){
             for(int d_ix = 0; d_ix < d_length; d_ix ++){
               for(int l = 0; l < L; l++){
                 intensities(sample_ix, k*d_length + d_ix) += w(l,k) * R::dnorm(d(d_ix),mu(l,k),sqrt(tau(l,k)),false);
			   }
				 global_intensity(sample_ix,d_ix) += pi(k) *  intensities(sample_ix,k*d_length+d_ix);
             }
          }

          cluster_assignment.row(sample_ix) = iter_cluster_assignment;
          cluster_component_assignment.row(sample_ix) = iter_component_assignment;
          pi_samps.row(sample_ix) = pi;
          Eigen::Map<Eigen::RowVectorXd> mu_samp(mu.data(),mu.size());
          Eigen::Map<Eigen::RowVectorXd> w_samp(w.data(),w.size());
          Eigen::Map<Eigen::RowVectorXd> tau_samp(tau.data(),tau.size());
          w_samps.row(sample_ix) = w_samp; // stored in column order
          mu_samps.row(sample_ix) = mu_samp;
          tau_samps.row(sample_ix) = tau_samp;
          alpha_samps(sample_ix,0) = alpha;
          rho_samps(sample_ix,0) = rho;
          sample_ix += 1;
        }
    }

    cluster_matrix = cluster_matrix / num_posterior_samples;

    return(Rcpp::List::create(Rcpp::Named("cluster_assignment") = cluster_assignment,
                              Rcpp::Named("component_assignment") =  cluster_component_assignment,
                              Rcpp::Named("cluster_pair_probability") =  cluster_matrix,
                              Rcpp::Named("pi_samples") = pi_samps,
                              Rcpp::Named("w_samples") =  w_samps,
                              Rcpp::Named("intensities") = intensities,
							  Rcpp::Named("global_intensity") = global_intensity,
                              Rcpp::Named("mu_samples") =  mu_samps,
                              Rcpp::Named("tau_samples") = tau_samps,
                              Rcpp::Named("alpha_samples") = alpha_samps,
                              Rcpp::Named("rho_samples") = rho_samps,
							  Rcpp::Named("alpha_prior") = alpha_prior,
							  Rcpp::Named("rho_prior") = rho_prior
    ));
}

//' Estimate the nonhomgogenous poisson process intensity function from grouped data with fixed concentration parameters
//'
//' @template bda
//'
//'
//' @param r vector of distances associatd with different BEFs
//' @param n_j matrix of integers denoting the start and length of each observations associated BEF distances
//' @param d a 1D grid of positive real values over which the differing intensities are evaluated
//' @param L component truncation number
//' @param K intensity cluster truncation number
//' @param J number of rows in r matrix; number of groups
//' @param mu_0 normal base measure prior mean
//' @param kappa_0 normal base measure prior variance scale 
//' @param nu_0 inverse chi sqaure base measure prior degrees of freedom
//' @param sigma_0 inverse chi square base measure prior scale  
//' @param alpha upper level concentration parameter (fixed)
//' @param rho lower level concentration parameter (fixed)
//' @param iter_max total number of iterations for which to run sampler
//' @param warm_up number of iterations for which to burn-in or "warm-up" sampler
//' @param thin number of iterations to thin by
//' @param seed integer with which to initialize random number generator
//' @param chain integer chain label
//' @param num_posterior_samples the total number of posterior samples after burn in 
//' @seealso the normal-inverse Chi-square conjugate parameterization in the reference
//'
// [[Rcpp::export]]
Rcpp::List nd_nhpp_fixed_fit(
        const Eigen::ArrayXd& r,
        const Eigen::MatrixXi& n_j,
        const Eigen::ArrayXd& d,
        const int& L,
        const int& K,
        const int& J,
        const double& mu_0,
        const double& kappa_0,
        const int& nu_0,
        const double& sigma_0,
		const double& alpha,
		const double& rho,
        const int& iter_max,
        const int& warm_up,
        const int& thin,
        const int& seed,
        const int& chain,
        const int& num_posterior_samples
        )
{
    // set seed
    std::mt19937 rng;
    rng = std::mt19937(seed);

    //create sample containers
    #include "support/create_sample_containers.hpp"

    #include "support/during_sampling_fixed.hpp"

    // create rng's
    std::uniform_real_distribution<double> runif(0,1);
    std::normal_distribution<double> rnorm(0,1);


    // initialize component weights
    u = stick_break(L,K,alpha,rng);
    v = stick_break(K,rho,rng);
    w = stick_break_weights(u);
    pi = stick_break_weights(v);

    mu = initialize_mu(L,K,mu_0,kappa_0,rng);
    tau = initialize_tau(L,K,sigma_0,nu_0);


    Rcpp::Rcout << "Beginning Sampling" << std::endl;
    Rcpp::Rcout << "----------------------------------------------------------------------" << std::endl;
    for(int iter_ix = 1; iter_ix <= iter_max; iter_ix ++){
      print_progress(iter_ix,warm_up,iter_max,chain);

      //calculate cluster probabilities
      q = dnorm(J,r,n_j,pi,w,mu,tau);

      for(int j = 0; j < J; j ++){
            probs = q.row(j);
            std::discrete_distribution<int> d(probs.data(),probs.data()+probs.size());
            iter_cluster_assignment(j) = d(rng);
        }

      // calculate w/in cluster component probabilities
        b = dnorm(r,n_j,w,mu,tau,iter_cluster_assignment);

        // assign distances to components within clusters
        component_count = Eigen::ArrayXXi::Zero(L,K);
        for(int j = 0 ; j < J; j ++){
            for(int i = 0 ; i < n_j(j,1) ; i++){
                for(int l = 0 ; l< L; l++)
                    prob(l) = b(n_j(j,0) + i,l);
                std::discrete_distribution<int> d(prob.data(),prob.data() + prob.size());
                iter_component_assignment(n_j(j,0) +i) = d(rng);
                component_count(iter_component_assignment(n_j(j,0)+i),iter_cluster_assignment(j)) += 1;
            }
        }

        // draw samples from intensity cluster parameters

        for(int k = 0; k < K; k++)
            cluster_count(k) = (iter_cluster_assignment == k).count();

        for(int k = 0; k < K; k++){
            v_posterior_beta_alpha(k) = 1 + cluster_count(k);
            v_posterior_beta_beta(k) = alpha + cluster_count.tail(K-k-1).sum();
        }

        v = stick_break(K, v_posterior_beta_alpha,v_posterior_beta_beta,rng);
        pi = stick_break_weights(v);



        for(int l = 0; l < L; l ++){
            for(int k = 0; k < K; k++){
                u_posterior_beta_alpha(l,k) = 1 + component_count(l,k);
                u_posterior_beta_beta(l,k) = rho + ((component_count.col(k)).tail(K-k-1)).sum();
            }
        }

        u = stick_break(u_posterior_beta_alpha,u_posterior_beta_beta,rng);
        w = stick_break_weights(u);

        // calculate mu hat (sq) sums

        ycount = Eigen::ArrayXXd::Zero(L,K);
        ycount_sq = Eigen::ArrayXXd::Zero(L,K);
        for(int l = 0 ; l < L; l ++){
          for(int k = 0; k < K; k++){
            for(int j = 0; j < J; j++){
              for(int i = 0; i < n_j(j,1); i ++){
                if(iter_cluster_assignment(j) == k && iter_component_assignment(n_j(j,0)+i) == l){
                  ycount(l,k) += r(n_j(j,0)+i);
                  ycount_sq(l,k) += pow(r(n_j(j,0)+i),2);
                }
              }
            }
          }
        }


        //sample mu via conjugacy

        for(int l = 0; l < L; l ++){
         for(int k = 0; k < K; k++){
           if(component_count(l,k) == 0){
             tau(l,k) = nu_0 * sigma_0 / R::rchisq(nu_0);
             mu(l,k) = rnorm(rng) * sqrt(tau(l,k) / kappa_0 ) + mu_0;
           }
           else{
             s_n  = nu_0 * sigma_0 + (ycount_sq(l,k) - (pow(ycount(l,k),2) / component_count(l,k) ) ) + (kappa_0 * component_count(l,k) / (kappa_0 + component_count(l,k))) * pow( (ycount(l,k) / component_count(l,k) - mu_0 ),2);
             tau(l,k) = s_n / R::rchisq(nu_0 + component_count(l,k));
             mu_n = ( (kappa_0  / tau(l,k)) * mu_0   + ycount(l,k) / tau(l,k)  ) / ( kappa_0 / tau(l,k) + component_count(l,k) / tau(l,k));
             s_n = 1.0 /( (kappa_0 /tau(l,k)) + (component_count(l,k) / tau(l,k)) ) ;
             mu(l,k) =  rnorm(rng) * sqrt(s_n) + mu_n;
           }
         }
       }


        // calculate adjacency matrix and store samples
        if((iter_ix > warm_up) && (iter_ix % thin ==0)){
          for(int j = 0; j < J; j++){
            for(int j_ = 0; j_ < j; j_++)
              cluster_matrix(j,j_) += (iter_cluster_assignment(j) == iter_cluster_assignment(j_)) ? 1: 0;
          }
          for(int k =0; k < K; k++){
             for(int d_ix = 0; d_ix < d_length; d_ix ++){
               for(int l = 0; l < L; l++){
                 intensities(sample_ix, k*d_length + d_ix) += w(l,k) * R::dnorm(d(d_ix),mu(l,k),sqrt(tau(l,k)),false);
			   }
				 global_intensity(sample_ix,d_ix) += pi(k) *  intensities(sample_ix,k*d_length+d_ix);
             }
          }

          cluster_assignment.row(sample_ix) = iter_cluster_assignment;
          cluster_component_assignment.row(sample_ix) = iter_component_assignment;
          pi_samps.row(sample_ix) = pi;
          Eigen::Map<Eigen::RowVectorXd> mu_samp(mu.data(),mu.size());
          Eigen::Map<Eigen::RowVectorXd> w_samp(w.data(),w.size());
          Eigen::Map<Eigen::RowVectorXd> tau_samp(tau.data(),tau.size());
          w_samps.row(sample_ix) = w_samp; // stored in column order
          mu_samps.row(sample_ix) = mu_samp;
          tau_samps.row(sample_ix) = tau_samp;
          alpha_samps(sample_ix,0) = alpha;
          rho_samps(sample_ix,0) = rho;
          sample_ix += 1;
        }
    }

    cluster_matrix = cluster_matrix / num_posterior_samples;

    return(Rcpp::List::create(Rcpp::Named("cluster_assignment") = cluster_assignment,
                              Rcpp::Named("component_assignment") =  cluster_component_assignment,
                              Rcpp::Named("cluster_pair_probability") =  cluster_matrix,
                              Rcpp::Named("pi_samples") = pi_samps,
                              Rcpp::Named("w_samples") =  w_samps,
                              Rcpp::Named("intensities") = intensities,
							  Rcpp::Named("global_intensity") = global_intensity,
                              Rcpp::Named("mu_samples") =  mu_samps,
                              Rcpp::Named("tau_samples") = tau_samps
    ));
}

// @param x double real number
double sigmoid(double& x){
    return((exp(x) / (1+exp(x))));
}

// initialize's matrix of mixture component means
// @param L number of mixture components
// @param K number of cluster components
// @param a_0 base measure hyperparameter
// @param b_0 base measure hyperparameter
// @param rng random number generator
Eigen::MatrixXd initialize_mu_beta(const int& L, const int& K, const double& a_0, const double& b_0, std::mt19937& rng){
    
    Eigen::MatrixXd out(L,K);
    sftrabbit::beta_distribution<> beta_dist(a_0,b_0);
    
    for(int l = 0; l < L; l ++){
        for(int k = 0; k < K; k++)
            out(l,k) = beta_dist(rng);
    }
    
    return(out);
    
}
// stick breaking atoms (not weights) for a vector with alpha,beta variable across vector elements
// @param n length of vector
// @param alpha vector of alpha parameters for posterior beta distribution
// @param beta vector of beta parameters for posterior beta distribution
Eigen::VectorXd stick_break(const int n, Eigen::VectorXd& alpha,Eigen::VectorXd& beta, std::mt19937& rng){
    
    Eigen::VectorXd v(n);
    v = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd w(n);
    w = Eigen::VectorXd::Zero(n);
    for(int i = 0; i < (n-1); i++){
        sftrabbit::beta_distribution<> beta_dist(alpha(i),beta(i));
        v(i) = beta_dist(rng);
    }
    v(n-1) = 1;
    
    return(v);
}

// stick breaking
Eigen::MatrixXd stick_break(Eigen::MatrixXd& alpha, Eigen::MatrixXd& beta, std::mt19937& rng){
    
    const int rows = alpha.rows();
    const int cols = alpha.cols();
    Eigen::MatrixXd out(rows,cols);
    
    for(int row_ix = 0; row_ix < rows; row_ix ++){
        for(int col_ix = 0; col_ix < cols; col_ix ++){
            sftrabbit::beta_distribution<> beta_dist(alpha(row_ix,col_ix),beta(row_ix,col_ix));
            out(row_ix,col_ix) = row_ix == (rows-1) ? 1.0 : beta_dist(rng);
        }
    }
    
    return(out);
}


Eigen::VectorXd stick_break_weights(const int n, Eigen::VectorXd& v){
    
    Eigen::VectorXd w(n);
    for(int i = 0; i < (n-1); i++)
        w(i) = i == 0 ? v(i) : v(i) * (Eigen::VectorXd::Ones(i) - v.head(i)).prod();
    
    w(n-1) =  v(n-1) * (Eigen::VectorXd::Ones(n-1) - v.head(n-1)).prod();
    return(w);
}

// stick break weights returning vector
Eigen::VectorXd stick_break_weights(Eigen::VectorXd& v){
    
    const int n = v.rows();
    Eigen::VectorXd w(n);
    for(int i = 0; i < (n-1); i++)
        w(i) = i == 0 ? v(i) : v(i) * (Eigen::VectorXd::Ones(i) - v.head(i)).prod();
    
    w(n-1) =  v(n-1) * (Eigen::VectorXd::Ones(n-1) - v.head(n-1)).prod();
    return(w);
}

// stick break weights returning matrix
Eigen::MatrixXd stick_break_weights(Eigen::MatrixXd& u){
    
    const int rows = u.rows();
    const int cols = u.cols();
    Eigen::MatrixXd w(rows,cols);
    
    
    for(int col_ix = 0; col_ix < cols; col_ix ++){
        for(int row_ix = 0; row_ix < rows ; row_ix ++)
            w(row_ix,col_ix) =  row_ix == 0 ? u(row_ix,col_ix) : u(row_ix,col_ix) * (Eigen::MatrixXd::Ones(row_ix,1) - u.block(0,col_ix,row_ix,1)).prod();
    }
    
    return(w);
}

// returns density from mixture of betas with global tau
Eigen::MatrixXd dbeta(const int& J, const Eigen::VectorXd& r,const Eigen::MatrixXi& n_j, Eigen::VectorXd& pi,Eigen::MatrixXd& w, Eigen::MatrixXd& mu,double& tau){
    
    const int rows = mu.rows();
    const int cols = mu.cols();
    Eigen::MatrixXd q(J,cols);     
    Eigen::VectorXd tmp;
    Eigen::MatrixXd blk;
    q = Eigen::MatrixXd::Zero(J,cols);
    
    for(int j = 0; j < J ; j++){
        blk = r.segment(n_j(j,0),n_j(j,1));
        for(int k = 0; k < cols; k ++){
            tmp = Eigen::RowVectorXd::Zero(n_j(j,1));
            for(int l = 0; l < rows; l ++){
                tmp = tmp.array() +  w(l,k) * (pow(blk.array(), mu(l,k)*tau-1) * pow( (Eigen::VectorXd::Ones(n_j(j,1)) - blk).array(), (1-mu(l,k))*tau - 1)).array() * exp(-lgamma(mu(l,k)*tau) - lgamma((1-mu(l,k))*tau) + lgamma(tau));
            }
            q(j,k) = pi(k) * tmp.prod();
        }
    }
    
    return(q);
}

// returns density from mixture of betas with cluster specific tau
Eigen::MatrixXd dbeta(const int& J, const Eigen::VectorXd& r,const Eigen::MatrixXi& n_j, Eigen::VectorXd& pi,Eigen::MatrixXd& w, Eigen::MatrixXd& mu,Eigen::MatrixXd& tau){
    
    const int rows = mu.rows();
    const int cols = mu.cols();
    Eigen::MatrixXd q(J,cols);     
    Eigen::VectorXd tmp;
    Eigen::MatrixXd blk;
    q = Eigen::MatrixXd::Zero(J,cols);
    
    for(int j = 0; j < J ; j++){
        blk = r.segment(n_j(j,0),n_j(j,1));
        for(int k = 0; k < cols; k ++){
            tmp = Eigen::RowVectorXd::Zero(n_j(j,1));
            for(int l = 0; l < rows; l ++){
                tmp = tmp.array() +  w(l,k) * (pow(blk.array(), mu(l,k)*tau(l,k)-1) * pow( (Eigen::VectorXd::Ones(n_j(j,1)) - blk).array(), (1-mu(l,k))*tau(l,k) - 1)).array() * exp(-lgamma(mu(l,k)*tau(l,k)) - lgamma((1-mu(l,k))*tau(l,k)) + lgamma(tau(l,k)));
            }
            q(j,k) = pi(k) * tmp.prod();
        }
    }
    
    return(q);
}

// Returns matrix of probabilities from mixture of betas with constant tau
Eigen::MatrixXd dbeta(const Eigen::VectorXd& r, const Eigen::MatrixXi& n_j, Eigen::MatrixXd& w, Eigen::MatrixXd& mu, double& tau, Eigen::VectorXi& zeta){
    
    const int L = mu.rows();
    const int J = zeta.rows();
    Eigen::MatrixXd b(r.rows(),L);
    
    for(int l = 0; l < L; l ++){
        for(int j = 0; j < J; j++){
            for(int i = 0; i < n_j(j,1) ; i++)
                b(n_j(j,0) + i,l) = w(l,zeta(j)) * R::dbeta(r(n_j(j,0) + i),mu(l,zeta(j))*tau,(1-mu(l,zeta(j)))*tau ,false);
        }
    }
    
    return(b);
}

// Returns matrix of probabilities from mixture of betas with cluster specific tau
Eigen::MatrixXd dbeta(const Eigen::VectorXd& r, const Eigen::MatrixXi& n_j, Eigen::MatrixXd& w, Eigen::MatrixXd& mu, Eigen::MatrixXd& tau, Eigen::VectorXi& zeta){
    
    const int L = mu.rows();
    const int J = zeta.rows();
    Eigen::MatrixXd b(r.rows(),L);
    
    for(int l = 0; l < L; l ++){
        for(int j = 0; j < J; j++){
            for(int i = 0; i < n_j(j,1) ; i++)
                b(n_j(j,0) + i,l) = w(l,zeta(j)) * R::dbeta(r(n_j(j,0) + i),mu(l,zeta(j))*tau(l,zeta(j)),(1-mu(l,zeta(j)))*tau(l,zeta(j)) ,false);
        }
    }
    
    return(b);
}

//' Estimate the nonhomgogenous poisson process intensity function from grouped data  
//' 
//' @param r vector of distances associatd with different BEFs 
//' @param n_j matrix of integers denoting the start and length of each school's associated BEF distances
//' @param d a 1D grid of positive real values over which the differing intensities are evaluated 
//' @param mu_sd scale for mu proposal dist'n
//' @param tau_sd scale for tau proposal dist'n
//' @param L component truncation number
//' @param K intensity cluster truncation number
//' @param J number of rows in r matrix; number of groups
//' @param a_0 hyperparameter for mu base measure
//' @param b_0 hyperparameter for mu base measure
//' @param a_alpha hyperparameter for alpha gamma prior
//' @param b_alpha hyperparameter for alpha gamma prior
//' @param a_rho hyperparameter for rho gamma prior
//' @param b_rho hyperparameter for rho gamma prior
//' @param iter_max total number of iterations for which to run sampler
//' @param warm_up number of iterations for which to burn-in or "warm-up" sampler
//' @param thin number of iterations to thin by
//' @param seed integer with which to initialize random number generator
//' @param chain integer chain label
//'
// [[Rcpp::export]]
Rcpp::List beta_nd_nhpp_fit(
        const Eigen::ArrayXd& r,
        const Eigen::MatrixXi& n_j,
        const Eigen::VectorXd& d,
        double& mu_sd,
        double& tau_sd,
        const int& L,
        const int& K,
        const int& J,
        const double& a_0,
        const double& b_0,
        const double& a_alpha,
        const double& b_alpha,
        const double& a_rho,
        const double& b_rho,
        const int& iter_max,
        const int& warm_up,
        const int& thin,
        const int& seed,
        const int& chain)
{
    // set seed
    std::mt19937 rng;
    rng = std::mt19937(seed);
    
    // create vector containers
    int num_posterior_samples = floor((double)((double)(iter_max - warm_up) / (double)thin));
    
    // post sampling
    Eigen::MatrixXd w_samps(num_posterior_samples,L*K);
    Eigen::MatrixXd mu_samps(num_posterior_samples,L*K);
    Eigen::MatrixXd pi_samps(num_posterior_samples,K);
    Eigen::MatrixXi cluster_assignment(num_posterior_samples,J);
    Eigen::MatrixXi cluster_component_assignment(num_posterior_samples,r.rows());
    Eigen::MatrixXd cluster_matrix(J,J);
    cluster_matrix = Eigen::MatrixXd::Zero(J,J);
    Eigen::MatrixXd alpha_samps(num_posterior_samples,1);
    Eigen::MatrixXd rho_samps(num_posterior_samples,1);
    Eigen::MatrixXd tau_samps(num_posterior_samples,1);
    Eigen::MatrixXd intensities(num_posterior_samples,K * d.size());
    Eigen::MatrixXd global_intensity(num_posterior_samples,d.size());
    Eigen::MatrixXd acceptance_rate(L,K);
    Eigen::MatrixXd batch_acceptance(L,K);
    acceptance_rate = Eigen::MatrixXd::Zero(L,K);
	global_intensity = Eigen::MatrixXd::Zero(num_posterior_samples,d.size());
    double tau_acceptance = 0;
    double tau_batch_acceptance = 0;
    int sample_ix = 0;
    
    // during sampling
    Eigen::MatrixXd u(L,K);
    Eigen::VectorXd v(K);
    Eigen::MatrixXd w(L,K);
    Eigen::VectorXd pi(K);
    Eigen::MatrixXd q(J,K);
    Eigen::MatrixXd b(r.rows(),L);
    Eigen::MatrixXd mu(L,K);
    Eigen::VectorXi iter_cluster_assignment(J);
    Eigen::VectorXi iter_component_assignment(r.rows()); 
    Eigen::VectorXi cluster_count(K);
    Eigen::MatrixXi component_count(L,K);
    Eigen::MatrixXd u_posterior_beta_alpha(L,K);
    Eigen::MatrixXd u_posterior_beta_beta(L,K);
    Eigen::VectorXd v_posterior_beta_alpha(K);
    Eigen::VectorXd v_posterior_beta_beta(K);
    Eigen::VectorXd probs(K);
    Eigen::VectorXd prob(L);
    const int d_length = d.size();
    double alpha;
    double rho;
    double posterior_a_alpha = a_alpha + (K-1);
    double posterior_b_alpha;
    double posterior_a_rho = a_rho + K*(L-1);
    double posterior_b_rho;
    double prob_temp;
    double prob_temp_prop;
    double ratio;
    double mu_temp;
    double mu_prop;
    double tau_prop;
    int progress = 0;
    
    // adaptation jigger
    Eigen::MatrixXd mu_sds(L,K);
    mu_sds = Eigen::MatrixXd::Ones(L,K) * mu_sd;
    Eigen::MatrixXd mu_lsi(L,K);
    mu_lsi = Eigen::MatrixXd::Ones(L,K) ;
    double tau_lsi = 1;
    double accept_temp;
    
    // create rng's 
    std::gamma_distribution<double> rgam_alpha(a_alpha,b_alpha);
    std::gamma_distribution<double> rgam_rho(a_rho,b_rho);
    std::uniform_real_distribution<double> runif(0,1);
    std::normal_distribution<double> rnorm(0,1);
    
    //initialize concentration parameters
    alpha = rgam_alpha(rng);
    rho = rgam_rho(rng);
    
    // initialize component weights
    u = stick_break(L,K,alpha,rng);
    v = stick_break(K,rho,rng);
    w = stick_break_weights(u);
    pi = stick_break_weights(v);
    
    // initialize variables from priors
    mu = initialize_mu_beta(L,K,a_0,b_0,rng);
    double tau = 1.5; // improper prior

    for(int iter_ix = 1; iter_ix <= iter_max ; iter_ix ++ ){
        if((iter_ix) % (int)round(.1 * iter_max) == 0 || iter_ix == 1 || iter_ix == (warm_up + 1) ){
            progress = (int)round( iter_ix * 100 / iter_max );
            std::string str = (iter_ix) <= (warm_up) ? "\t [Warmup]" : "\t [Sampling]";
            Rcpp::Rcout << "[Chain " << chain << "] Beginning of iteration: " << (iter_ix) << " / " << iter_max << " (" << progress << "%)" << str  << std::endl;
        }
        // calculate intensity cluster probabilities
        q = dbeta(J, r,n_j,pi,w,mu,tau);
        
        for(int j = 0; j < J; j ++){
            probs = q.row(j);
            std::discrete_distribution<int> d(probs.data(),probs.data()+probs.size());
            iter_cluster_assignment(j) = d(rng);
        }
        
        // calculate within cluster component probabilities
        
        b = dbeta(r,n_j,w,mu,tau,iter_cluster_assignment);
        
        for(int j = 0 ; j < J; j ++){
            for(int i = 0 ; i < n_j(j,1) ; i++){
                for(int l = 0 ; l< L; l++)
                    prob(l) = b(n_j(j,0) + i,l);
                std::discrete_distribution<int> d(prob.data(),prob.data() + prob.size());
                iter_component_assignment(n_j(j,0) +i) = d(rng);
            }
        }
        
        // draw samples from  intensity cluster  parameters 
        
        for(int k = 0; k < K; k++)
            cluster_count(k) = (iter_cluster_assignment.array() == k).count();
        
        for(int k = 0; k < K; k++){
            v_posterior_beta_alpha(k) = 1 + cluster_count(k);
            v_posterior_beta_beta(k) = alpha + cluster_count.tail(K-k-1).sum();
        }
        
        v = stick_break(K, v_posterior_beta_alpha,v_posterior_beta_beta,rng);
        pi = stick_break_weights(v);
        
        for(int l = 0; l < L; l++){
            for(int k = 0; k < K; k++)
                component_count(l,k) = (iter_component_assignment.col(k).array() == l).count();
        }
        
        for(int l = 0; l < L; l ++){
            for(int k = 0; k < K; k++){
                u_posterior_beta_alpha(l,k) = 1 + component_count(l,k);
                u_posterior_beta_beta(l,k) = rho + ((component_count.col(k)).tail(K-k-1)).sum();
            }
        }
        
        u = stick_break(u_posterior_beta_alpha,u_posterior_beta_beta,rng);
        w = stick_break_weights(u);
        
        
        // sample mus via MH step
        prob_temp = 0;
        prob_temp_prop = 0;
        for(int l = 0 ; l < L;  l ++){
            for(int k = 0; k < K; k++){
                mu_prop = rnorm(rng)*mu_sds(l,k) + mu(l,k);
                mu_prop = sigmoid(mu_prop);
                for(int j = 0; j < J; j ++){
                    for(int i = 0; i < n_j(j,1); i ++){
                        if(iter_cluster_assignment(j) == k && iter_component_assignment(n_j(j,0)+i) == l){
                            prob_temp_prop += R::dbeta(r(n_j(j,0)+i),mu_prop*tau,(1-mu_prop)*tau,true);
                            prob_temp += R::dbeta(r(n_j(j,0)+i),mu(l,k)*tau,(1-mu(l,k))*tau,true);
                        }
                    }
                }
                ratio = prob_temp_prop- prob_temp + log(mu_prop) + log(1-mu_prop) - log(mu(l,k)) - log(1-mu(l,k));
                if(prob_temp != 0){
                    mu(l,k) = runif(rng) <= exp(ratio) ?  mu_prop : mu(l,k); 
                    acceptance_rate(l,k) += mu(l,k) == mu_prop ? 1 : 0 ;
                    batch_acceptance(l,k) += mu(l,k) == mu_prop ? 1 : 0 ;
                }
                else
                    mu(l,k) = runif(rng); // sample from prior if no cluster members
                prob_temp = 0;
                prob_temp_prop = 0;
            }
        }

        //sort mu 
        //mu = sort_columns(mu);
        
        
        // sample concentration parameters
        posterior_b_alpha = 1.0 / b_alpha - (log(1-v.array())).head(K-1).sum() ;
        posterior_b_rho =  1.0 / b_rho - log(1-u.block(0,0,L-1,K).array()).matrix().colwise().sum().sum();
        std::gamma_distribution<double> rgam_alpha(posterior_a_alpha, 1.0 / posterior_b_alpha);
        std::gamma_distribution<double> rgam_rho(posterior_a_rho, 1.0 / posterior_b_rho);
        alpha = rgam_alpha(rng);
        rho = rgam_rho(rng);
        
        
        // sample tau
        prob_temp = 0;
        prob_temp_prop = 0;
        tau_prop = exp(log(tau) + rnorm(rng) * tau_sd);
        for(int j = 0; j < J; j++){
            for(int i = 0 ; i < n_j(j,1) ; i++){
                mu_temp = mu(iter_component_assignment(n_j(j,0)+i),iter_cluster_assignment(j));
                prob_temp_prop += R::dbeta(r(n_j(j,0)+i),mu_temp*tau_prop,(1-mu_temp)*tau_prop,true);
                prob_temp += R::dbeta(r(n_j(j,0)+i),mu_temp*tau,(1-mu_temp)*tau,true);
            }
        }

        if(runif(rng) <= exp(prob_temp_prop - prob_temp + log(tau_prop) - log(tau))){
            tau = tau_prop;
            tau_acceptance += 1;
            tau_batch_acceptance += 1;
        }
        
        
        // adapt RWM
        if( (iter_ix <= warm_up) && (iter_ix % 50 == 0) ){

            for(int l = 0; l < L; l ++){
                for(int k = 0; k < K; k ++){
                    accept_temp = batch_acceptance(l,k) / 50.0;
                    mu_lsi(l,k) = accept_temp >= .44 ? mu_lsi(l,k) + std::min(0.01,pow(iter_ix,-.5)) : mu_lsi(l,k) - std::min(0.01,pow(iter_ix,-.5));
                    mu_sds(l,k) =  pow(mu_sds(l,k),mu_lsi(l,k));
                    batch_acceptance(l,k) = 0;
                }
            }
            tau_sd = (tau_batch_acceptance / 50.0) >= .44 ? pow(tau_sd,tau_lsi + std::min(0.01,pow(iter_ix,-.5))) : pow(tau_sd,tau_lsi - std::min(0.01,pow(iter_ix,-.5)));
            tau_batch_acceptance = 0;
        }

        // store samples
        if((iter_ix > warm_up) && (iter_ix % thin == 0)){
            for(int j = 0; j < J; j++){
                for(int j_ = 0; j_ < j; j_++)
                    cluster_matrix(j,j_) +=  (iter_cluster_assignment(j) == iter_cluster_assignment(j_)) ? 1 : 0;
            }
            for(int k = 0; k < K ; k++){
                for(int d_ix = 0; d_ix < d_length; d_ix++){
                    for(int l = 0; l < L; l ++){
                        intensities(sample_ix, k*d_length + d_ix ) += w(l,k) * R::dbeta(d(d_ix), mu(l,k)*tau,(1-mu(l,k))*tau,false);
                    }
					 global_intensity(sample_ix,d_ix) += pi(k) *  intensities(sample_ix,k*d_length+d_ix);
                }
            }
            
            cluster_assignment.row(sample_ix) = iter_cluster_assignment;
            cluster_component_assignment.row(sample_ix) = iter_component_assignment;
            pi_samps.row(sample_ix) = pi;
            Eigen::Map<Eigen::RowVectorXd> mu_samp(mu.data(),mu.size());
            Eigen::Map<Eigen::RowVectorXd> w_samp(w.data(),w.size());
            w_samps.row(sample_ix) = w_samp; // stored in column order
            mu_samps.row(sample_ix) = mu_samp;
            tau_samps(sample_ix,0) = tau;
            alpha_samps(sample_ix,0) = alpha;
            rho_samps(sample_ix,0) = rho;
            sample_ix += 1;
        }
    }
    
    cluster_matrix = cluster_matrix / (num_posterior_samples);
    
    return(Rcpp::List::create(Rcpp::Named("cluster_assignment") = cluster_assignment, 
                              Rcpp::Named("component_assignment") =  cluster_component_assignment,
                              Rcpp::Named("cluster_pair_probability") =  cluster_matrix, 
                              Rcpp::Named("pi_samples") = pi_samps,
                              Rcpp::Named("w_samples") =  w_samps,
                              Rcpp::Named("intensities") = intensities,
                              Rcpp::Named("global_density") = global_intensity,
                              Rcpp::Named("mu_samples") =  mu_samps,
                              Rcpp::Named("mu_acceptance") = (acceptance_rate / iter_max),
                              Rcpp::Named("tau_samples") = tau_samps,
                              Rcpp::Named("tau_acceptance") = tau_acceptance / iter_max,
                              Rcpp::Named("alpha_samples") = alpha_samps,
                              Rcpp::Named("rho_samples") = rho_samps
    ));
}

//' Estimate the nonhomgogenous poisson process intensity function from grouped data using multiple taus 
//' 
//' @param r vector of distances associatd with different BEFs 
//' @param n_j matrix of integers denoting the start and length of each school's associated BEF distances
//' @param d a 1D grid of positive real values over which the differing intensities are evaluated 
//' @param mu_sd scale for mu proposal dist'n
//' @param tau_sd  not used
//' @param L component truncation number
//' @param K intensity cluster truncation number
//' @param J number of rows in r matrix; number of groups
//' @param a_0 hyperparameter for mu base measure
//' @param b_0 hyperparameter for mu base measure
//' @param a_alpha hyperparameter for alpha gamma prior
//' @param b_alpha hyperparameter for alpha gamma prior
//' @param a_rho hyperparameter for rho gamma prior
//' @param b_rho hyperparameter for rho gamma prior
//' @param iter_max total number of iterations for which to run sampler
//' @param warm_up number of iterations for which to burn-in or "warm-up" sampler
//' @param thin number of iterations to thin by
//' @param seed integer with which to initialize random number generator
//' @param chain integer chain label
//'
// [[Rcpp::export]]
Rcpp::List beta_nd_nhpp_fit_multiple_taus(
        const Eigen::VectorXd& r,
        const Eigen::MatrixXi& n_j,
        const Eigen::VectorXd& d,
        double& mu_sd,
        double& tau_sd,
        const int& L,
        const int& K,
        const int& J,
        const double& a_0,
        const double& b_0,
        const double& a_alpha,
        const double& b_alpha,
        const double& a_rho,
        const double& b_rho,
        const int& iter_max,
        const int& warm_up,
        const int& thin,
        const int& seed,
        const int& chain)
{
    // set seed
    std::mt19937 rng;
    rng = std::mt19937(seed);
    
    // create vector containers
    int num_posterior_samples = floor((double)((double)(iter_max - warm_up) / (double)thin));
    
    // post sampling
    Eigen::MatrixXd w_samps(num_posterior_samples,L*K);
    Eigen::MatrixXd mu_samps(num_posterior_samples,L*K);
    Eigen::MatrixXd pi_samps(num_posterior_samples,K);
    Eigen::MatrixXi cluster_assignment(num_posterior_samples,J);
    Eigen::MatrixXi cluster_component_assignment(num_posterior_samples,r.rows());
    Eigen::MatrixXd cluster_matrix(J,J);
    cluster_matrix = Eigen::MatrixXd::Zero(J,J);
    Eigen::MatrixXd alpha_samps(num_posterior_samples,1);
    Eigen::MatrixXd rho_samps(num_posterior_samples,1);
    Eigen::MatrixXd tau_samps(num_posterior_samples,L*K);
    Eigen::MatrixXd intensities(num_posterior_samples,K * d.size());
    Eigen::MatrixXd global_intensity(num_posterior_samples,d.size());
    Eigen::MatrixXd acceptance_rate(L,K);
    acceptance_rate = Eigen::MatrixXd::Zero(L,K);
	global_intensity = Eigen::MatrixXd::Zero(num_posterior_samples,d.size());
    int sample_ix = 0;
    
    // during sampling
    Eigen::MatrixXd u(L,K);
    Eigen::VectorXd v(K);
    Eigen::MatrixXd w(L,K);
    Eigen::VectorXd pi(K);
    Eigen::MatrixXd q(J,K);
    Eigen::MatrixXd b(r.rows(),L);
    Eigen::MatrixXd mu(L,K);
    Eigen::VectorXi iter_cluster_assignment(J);
    Eigen::VectorXi iter_component_assignment(r.rows()); 
    Eigen::VectorXi cluster_count(K);
    Eigen::MatrixXi component_count(L,K);
    Eigen::MatrixXd u_posterior_beta_alpha(L,K);
    Eigen::MatrixXd u_posterior_beta_beta(L,K);
    Eigen::VectorXd v_posterior_beta_alpha(K);
    Eigen::VectorXd v_posterior_beta_beta(K);
    Eigen::VectorXd probs(K);
    Eigen::VectorXd prob(L);
    const int d_length = d.size();
    double alpha;
    double rho;
    double posterior_a_alpha = a_alpha + (K-1);
    double posterior_b_alpha;
    double posterior_a_rho = a_rho + K*(L-1);
    double posterior_b_rho;
    double prob_temp;
    double prob_temp_prop;
    double ratio;
    double mu_prop;
    double tau_prop;
    int progress = 0;
    
    // adaptation jigger
    Eigen::MatrixXd batch_acceptance(L,K);
    batch_acceptance = Eigen::MatrixXd::Zero(L,K);
    Eigen::MatrixXd mu_sds(L,K);
    mu_sds = Eigen::MatrixXd::Ones(L,K) * mu_sd;
    Eigen::MatrixXd mu_lsi(L,K);
    mu_lsi = Eigen::MatrixXd::Ones(L,K) ;
    double accept_temp;

    // create rng's 
    std::gamma_distribution<double> rgam_alpha(a_alpha,b_alpha);
    std::gamma_distribution<double> rgam_rho(a_rho,b_rho);
    std::gamma_distribution<double> rgam_tau(1,1);
    std::uniform_real_distribution<double> runif(0,1);
    std::normal_distribution<double> rnorm(0,1);
    
    //initialize concentration parameters
    alpha = rgam_alpha(rng);
    rho = rgam_rho(rng);
    
    // initialize component weights
    u = stick_break(L,K,alpha,rng);
    v = stick_break(K,rho,rng);
    w = stick_break_weights(u);
    pi = stick_break_weights(v);
    
    // initialize variables from priors
    mu = initialize_mu_beta(L,K,a_0,b_0,rng);
    Eigen::MatrixXd tau(L,K); 
    tau = Eigen::MatrixXd::Ones(L,K) * 1.5;

    Eigen::VectorXd temp;
    
    for(int iter_ix = 1; iter_ix <= iter_max ; iter_ix ++ ){
        if((iter_ix) % (int)round(.1 * iter_max) == 0 || iter_ix == 1 || iter_ix == (warm_up + 1) ){
            progress = (int)round( iter_ix * 100 / iter_max );
            std::string str = (iter_ix) <= (warm_up) ? "\t [Warmup]" : "\t [Sampling]";
            Rcpp::Rcout << "[Chain " << chain << "] Beginning of iteration: " << (iter_ix) << " / " << iter_max << " (" << progress << "%)" << str  << std::endl;
        }
        // calculate intensity cluster probabilities
        q = dbeta(J, r,n_j,pi,w,mu,tau);
        
        for(int j = 0; j < J; j ++){
            probs = q.row(j);
            std::discrete_distribution<int> d(probs.data(),probs.data()+probs.size());
            iter_cluster_assignment(j) = d(rng);
        }
        
        // calculate within cluster component probabilities
        
        b = dbeta(r,n_j,w,mu,tau,iter_cluster_assignment);
        
        for(int j = 0 ; j < J; j ++){
            for(int i = 0 ; i < n_j(j,1) ; i++){
                for(int l = 0 ; l< L; l++)
                    prob(l) = b(n_j(j,0) + i,l);
                std::discrete_distribution<int> d(prob.data(),prob.data() + prob.size());
                iter_component_assignment(n_j(j,0) +i) = d(rng);
            }
        }
        
        // draw samples from  intensity cluster  parameters 
        
        for(int k = 0; k < K; k++)
            cluster_count(k) = (iter_cluster_assignment.array() == k).count();
        
        for(int k = 0; k < K; k++){
            v_posterior_beta_alpha(k) = 1 + cluster_count(k);
            v_posterior_beta_beta(k) = alpha + cluster_count.tail(K-k-1).sum();
        }
        
        v = stick_break(K, v_posterior_beta_alpha,v_posterior_beta_beta,rng);
        pi = stick_break_weights(v);
        
        for(int l = 0; l < L; l++){
            for(int k = 0; k < K; k++)
                component_count(l,k) = (iter_component_assignment.col(k).array() == l).count();
        }
        
        for(int l = 0; l < L; l ++){
            for(int k = 0; k < K; k++){
                u_posterior_beta_alpha(l,k) = 1 + component_count(l,k);
                u_posterior_beta_beta(l,k) = rho + ((component_count.col(k)).tail(K-k-1)).sum();
            }
        }
        
        u = stick_break(u_posterior_beta_alpha,u_posterior_beta_beta,rng);
        w = stick_break_weights(u);
        
        
        // sample mus/taus via MH step
        prob_temp = 0;
        prob_temp_prop = 0;
        for(int l = 0 ; l < L;  l ++){
            for(int k = 0; k < K; k++){
                mu_prop = rnorm(rng)*mu_sds(l,k) + mu(l,k);
                mu_prop = sigmoid(mu_prop);
                tau_prop = exp(rnorm(rng)*mu_sds(l,k) + log(tau(l,k)));
                for(int j = 0; j < J; j ++){
                    for(int i = 0; i < n_j(j,1); i ++){
                        if(iter_cluster_assignment(j) == k && iter_component_assignment(n_j(j,0)+i) == l){
                            prob_temp_prop += R::dbeta(r(n_j(j,0)+i),mu_prop*tau_prop,(1-mu_prop)*tau_prop,true);
                            prob_temp += R::dbeta(r(n_j(j,0)+i),mu(l,k)*tau(l,k),(1-mu(l,k))*tau(l,k),true);
                        }
                    }
                }
                ratio = prob_temp_prop - prob_temp + log(mu_prop) + log(1-mu_prop) - log(mu(l,k)) - log(1-mu(l,k)) + log(tau_prop) - log(tau(l,k))  ;
                if(prob_temp != 0){
                    if(runif(rng)<=exp(ratio)){
                        mu(l,k) = mu_prop;
                        tau(l,k) = tau_prop;
                        acceptance_rate(l,k) +=  1.0; 
                        batch_acceptance(l,k) += 1.0;
                    }
                }
                else{
                    mu(l,k) = runif(rng); // sample from prior if no cluster members
                    tau(l,k) = rgam_tau(rng);
                }
                prob_temp = 0;
                prob_temp_prop = 0;
            }
        }
        
        // sample concentration parameters
        posterior_b_alpha = 1.0 / b_alpha - (log(1-v.array())).head(K-1).sum() ;
        posterior_b_rho =  1.0 / b_rho - log(1-u.block(0,0,L-1,K).array()).matrix().colwise().sum().sum();
        std::gamma_distribution<double> rgam_alpha(posterior_a_alpha, 1.0 / posterior_b_alpha);
        std::gamma_distribution<double> rgam_rho(posterior_a_rho, 1.0 / posterior_b_rho);
        alpha = rgam_alpha(rng);
        rho = rgam_rho(rng);
        
        
        // adapt RWM
        if( (iter_ix <= warm_up) && (iter_ix % 50 == 0) ){
            for(int l = 0; l < L; l ++){
                for(int k = 0; k < K; k ++){
                    accept_temp = batch_acceptance(l,k) / 50.0;
                    mu_lsi(l,k) = accept_temp >= .44 ? mu_lsi(l,k) + std::min(0.01,pow(iter_ix,-.5)) : mu_lsi(l,k) - std::min(0.01,pow(iter_ix,-.5));
                    mu_sds(l,k) =  pow(mu_sds(l,k),mu_lsi(l,k));
                    batch_acceptance(l,k) = 0;
                }
            }
        }

        // store samples
        if((iter_ix > warm_up) && (iter_ix % thin == 0)){
            for(int j = 0; j < J; j++){
                for(int j_ = 0; j_ < j; j_++)
                    cluster_matrix(j,j_) +=  (iter_cluster_assignment(j) == iter_cluster_assignment(j_)) ? 1 : 0;
            }
            for(int k = 0; k < K ; k++){
                for(int d_ix = 0; d_ix < d_length; d_ix++){
                    for(int l = 0; l < L; l ++){
                        intensities(sample_ix, k*d_length + d_ix ) += w(l,k) * R::dbeta(d(d_ix), mu(l,k)*tau(k),(1-mu(l,k))*tau(k),false);
                    }
					 global_intensity(sample_ix,d_ix) += pi(k) *  intensities(sample_ix,k*d_length+d_ix);
                }
            }
            
            cluster_assignment.row(sample_ix) = iter_cluster_assignment;
            cluster_component_assignment.row(sample_ix) = iter_component_assignment;
            pi_samps.row(sample_ix) = pi;
            Eigen::Map<Eigen::RowVectorXd> mu_samp(mu.data(),mu.size());
            Eigen::Map<Eigen::RowVectorXd> w_samp(w.data(),w.size());
            Eigen::Map<Eigen::RowVectorXd> tau_samp(tau.data(),tau.size());
            w_samps.row(sample_ix) = w_samp; // stored in column order
            mu_samps.row(sample_ix) = mu_samp;
            tau_samps.row(sample_ix) = tau_samp;
            alpha_samps(sample_ix,0) = alpha;
            rho_samps(sample_ix,0) = rho;
            sample_ix += 1;
        }
    }
    
    cluster_matrix = cluster_matrix / (num_posterior_samples);

    
    return(Rcpp::List::create(Rcpp::Named("cluster_assignment") = cluster_assignment, 
                              Rcpp::Named("component_assignment") =  cluster_component_assignment,
                              Rcpp::Named("cluster_pair_probability") =  cluster_matrix, 
                              Rcpp::Named("pi_samples") = pi_samps,
                              Rcpp::Named("w_samples") =  w_samps,
                              Rcpp::Named("intensities") = intensities,
                              Rcpp::Named("global_density") = global_intensity,
                              Rcpp::Named("mu_samples") =  mu_samps,
                              Rcpp::Named("mu_acceptance") = (acceptance_rate / iter_max),
                              Rcpp::Named("tau_samples") = tau_samps,
                              Rcpp::Named("alpha_samples") = alpha_samps,
                              Rcpp::Named("rho_samples") = rho_samps
    ));
}

Eigen::VectorXd rnorm_draw(int &n, std::mt19937 &rng){

	Eigen::VectorXd out(n);
	std::normal_distribution<double> rnorm(0.0,1.0);

	for (int i = 0; i < n ; i++){
		out(i) = rnorm(rng);
	}
	return(out);
}

Eigen::VectorXd runif_draw(int &n, std::mt19937 &rng){

	Eigen::VectorXd out(n);
	std::uniform_real_distribution<double> runif(-2.0,2.0);
	for(int i =0; i < n ; i ++)
		out(i) = runif(rng);

	return(out);
}


class NHPP
{
	public:
		Eigen::ArrayXXd X;
		Eigen::ArrayXd n_j;
		Eigen::ArrayXd eta;
		Eigen::VectorXd beta_grad;
  	    NHPP(Eigen::ArrayXXd &input_X,
			 Eigen::ArrayXd &input_n_j){
			X = input_X;
			n_j = input_n_j;
	}

	double calculate_energy(Eigen::VectorXd &beta,
							Eigen::VectorXd &momenta){

		double out = 0;
		eta = X.matrix() * beta;
		out = (n_j * eta).sum()  - (exp(eta)).sum();
		out -=  beta.dot(beta) / (2.0 * 9.0) ;
		out -= momenta.dot(momenta) / 2.0 ;

		return(out);
	}

	double sample_u(Eigen::VectorXd &beta,Eigen::VectorXd &momenta,std::mt19937 &rng){

		double energy = calculate_energy(beta,momenta);
		std::uniform_real_distribution<double> runif(0.0,1.0);
		double out = energy + log(runif(rng));
		return(out);
	}

	Eigen::VectorXd calculate_gradient(Eigen::VectorXd &beta){

		beta_grad = (X.colwise() * n_j).colwise().sum() - (X.colwise() * exp( (X.matrix() * beta).array()) ).colwise().sum();
		beta_grad = beta_grad - beta / 9.0 ;

		return(beta_grad);

	}

	double FindReasonableEpsilon(int &p, Eigen::VectorXd &beta,std::mt19937 &rng){

		double epsilon = 1.0;
		int a;
		Eigen::VectorXd beta_temp;
		Eigen::VectorXd momenta_temp;
		Eigen::VectorXd momenta_static = rnorm_draw(p,rng);
		Eigen::VectorXd beta_grad;
		double ratio, energy_init, energy_prop;
		energy_init = calculate_energy(beta,momenta_static);
		beta_grad = calculate_gradient(beta);
		momenta_temp = momenta_static + epsilon * beta_grad / 2.0;
		beta_temp = beta + epsilon * momenta_temp;
		beta_grad = calculate_gradient(beta);
		momenta_temp  = momenta_temp + epsilon * beta_grad / 2.0;
		energy_prop = calculate_energy(beta_temp,momenta_temp);
		ratio  = energy_prop - energy_init;
		ratio = std::isinf(-ratio) ? - DBL_MAX: ratio;
		a = ratio > log(.5) ? 1 : - 1;
		int cntr = 0;
		while(a * ratio > -a * log(2)){
			if(cntr>100)
				break;
			epsilon = pow(2,a) * epsilon;
			beta_grad = calculate_gradient(beta);
			momenta_temp = momenta_static + epsilon * beta_grad / 2.0;
			beta_temp = beta + epsilon * momenta_temp;
			beta_grad = calculate_gradient(beta);
			momenta_temp  = momenta_temp + epsilon * beta_grad / 2.0;
			energy_prop = calculate_energy(beta_temp,momenta_temp);
			ratio = energy_prop - energy_init;
			cntr ++;
		}

		return(epsilon);
	}

};

#include "support/Tree.hpp"


//' Returns draw of beta regression coefficients from posterior distribution via No U-Turn Sampler
//'
//' @template nuts 
//' @param warm_up number of iterations for which to tune step-size 
//' @param iter_max total number of iterations for which to run the sampler
//' @param input_X design matrix
//' @param input_n_j outcome of counts
//' @param adapt_delta (0,1) scalar that denotes the average acceptance probability 
//' @param seed random number generator seed
//' @seealso the reference linked
//'
//[[Rcpp::export]]
Rcpp::List nhpp_gamma(
		const int &warm_up,
		const int &iter_max,
		Eigen::ArrayXXd &input_X,
		Eigen::ArrayXd &input_n_j,
		const double adapt_delta,
		const int &seed
		){

	// set seed
    std::mt19937 rng;
    rng = std::mt19937(seed);

	// tree tuning parameters
	double u;
	int s = 1;
	double n = 1.0;
	int j = 0;
	int p = input_X.cols();
	int chain = 1;
	// beta parameters
	Eigen::MatrixXd beta_draws(iter_max,p);
	beta_draws = Eigen::MatrixXd::Zero(iter_max, p);
	Eigen::ArrayXd epsilons(iter_max);
	Eigen::VectorXd acceptance(iter_max);
	acceptance = Eigen::VectorXd::Zero(iter_max);
	Eigen::VectorXi treedepth(iter_max);
	Eigen::VectorXd return_beta;
	Eigen::VectorXd current_beta = runif_draw(p,rng);
	// initialize model
	NHPP model(input_X,input_n_j);
	double epsilon = 0;
	epsilon = model.FindReasonableEpsilon(p,current_beta,rng);
	// epsilon tuning parameters
	double epsilon_bar = 1;
	double mu = log(10 * epsilon);
	double H_bar = 0;
	double kappa = 0.75;
	double gamma = 0.05;
	double t_0 = 10.0;
	// random number engine
	std::uniform_real_distribution<double> runif(0.0,1.0);

	for(int iter_ix = 1; iter_ix <= iter_max; iter_ix ++){
		print_progress(iter_ix,warm_up,iter_max,chain);
		Eigen::VectorXd beta_left = current_beta;
		Eigen::VectorXd beta_right = current_beta;
		Eigen::VectorXd momenta = rnorm_draw(p,rng);
		Eigen::VectorXd momenta_left = momenta;
		Eigen::VectorXd momenta_right = momenta;
		s = 1;
		n = 1.0;
		j = 0;
		u = model.sample_u(current_beta,momenta,rng);
		Tree tree(p,rng);
		while(s == 1){
			if(runif(rng)<=1){
				tree.buildTree(model,beta_left,momenta_left,u,-1,j,epsilon,current_beta,momenta,rng);
				beta_left = tree.get_bl();
				momenta_left = tree.get_ml();
			}
			else{
				tree.buildTree(model,beta_right,momenta_right,u,1,j,epsilon,current_beta,momenta,rng);
				beta_right = tree.get_br();
				momenta_right = tree.get_mr();
			}
			if(tree.get_s_prime() == 1){
				if(runif(rng)<= tree.get_n_prime() / n){
					acceptance(iter_ix-1) = 1;
					return_beta = tree.get_bn();
				}
			}
			n += tree.get_n_prime();
			s = tree.get_s_prime() * ((beta_right - beta_left).dot(momenta_left) >= 0.0 )  * ((beta_right - beta_left).dot(momenta_right) >= 0.0) ;
			j++;
			if(j>10)
				break;
		}
		treedepth(iter_ix-1) = j;
		if(iter_ix <= warm_up){
			epsilons(iter_ix-1) = epsilon;
			H_bar = ( 1 - 1.0 / (iter_ix + t_0)) * H_bar;
			H_bar +=  1.0 / (iter_ix + t_0) * (adapt_delta - tree.get_alpha_prime() / tree.get_n_alpha());
			epsilon = exp(mu - (sqrt(iter_ix) / gamma) * H_bar);
			epsilon_bar = exp(pow(iter_ix,-kappa) * log(epsilon) + (1.0 - pow(iter_ix,-kappa)) * log(epsilon_bar));
		}else{
			epsilon = epsilon_bar;
			epsilons(iter_ix-1) = epsilon;
		}
		// store sample
		if(acceptance(iter_ix-1) == 1){
			current_beta = return_beta;
			beta_draws.row(iter_ix - 1) = return_beta.transpose();
		}
	}

	return(Rcpp::List::create(
				Rcpp::Named("treedepth") = treedepth,
				Rcpp::Named("beta_samples") = beta_draws,
				Rcpp::Named("epsilons") = epsilons,
				Rcpp::Named("acceptance") = acceptance));
}

//' Computes Green and Lau loss function with unknown classification
//'
//' @param cluster_assignment iter_total x J cluster assignment matrix
//' @param pmat J x J pairwise probability of co-clustering matrix
//' @param tau penalty parameter 
// [[Rcpp::export]]
Eigen::ArrayXd green_loss_unknown(const Eigen::ArrayXXi &cluster_assignment,
							 const Eigen::ArrayXXd &pmat,
							 const double &tau){


	Eigen::ArrayXd loss(cluster_assignment.rows());
	loss = Eigen::ArrayXd::Zero(cluster_assignment.rows());

	for(int iter_ix = 0; iter_ix < cluster_assignment.rows() ; iter_ix ++){
		for(int j = 0; j < pmat.rows(); j ++){
			for(int j_ = 0; j_ < j ; j_++)
				loss(iter_ix) += cluster_assignment(iter_ix,j) == cluster_assignment(iter_ix,j_) ? (pmat(j,j_) - tau) : 0.0;
		}
	}

	return(loss);
}


//' Computes Green and Lau Loss function with known classification
//'
//' @param cluster_assignment iter_total x J cluster assignment matrix
//' @param pmat J x J pairwise probability of co-clustering matrix
//' @param true_cluster_assignment J x J true Adjacency Matrix
//' @param a mis-classification penalty parameter
//' @param b classification penalty parameter
// [[Rcpp::export]]
Eigen::ArrayXd green_loss_known(const Eigen::ArrayXXi &cluster_assignment,
								const Eigen::ArrayXXd &pmat,
								const Eigen::ArrayXXi &true_cluster_assignment,
								const double &a,
								const double &b
								 ){


	Eigen::ArrayXd loss(cluster_assignment.rows());
	loss = Eigen::ArrayXd::Zero(cluster_assignment.rows());

	for(int iter_ix = 0; iter_ix < cluster_assignment.rows() ; iter_ix ++){
		for(int j = 0; j < pmat.rows(); j ++){
			for(int j_ = 0; j_ < j ; j_++){
				loss(iter_ix) += ((cluster_assignment(iter_ix,j) != cluster_assignment(iter_ix,j_)) && (true_cluster_assignment(j,j_) == 1 )  ) ? a : 0.0;
				loss(iter_ix) += ((cluster_assignment(iter_ix,j) == cluster_assignment(iter_ix,j_)) && (true_cluster_assignment(j,j_) == 0 ) ) ? b : 0.0;
			}
		}
	}

	return(loss);
}


//' Computes Square loss with unknown classification
//'
//' @param cluster_assignment iter_total x J cluster assignment matrix
//' @param pmat J x J pairwise probability of co-clustering matrix
// [[Rcpp::export]]
Eigen::ArrayXd square_error(const Eigen::ArrayXXi &cluster_assignment,
							const Eigen::ArrayXXd &pmat){


	Eigen::ArrayXd loss(cluster_assignment.rows());
	loss = Eigen::ArrayXd::Zero(cluster_assignment.rows());

	for(int iter_ix = 0; iter_ix < cluster_assignment.rows() ; iter_ix ++){
		for(int j = 0; j < pmat.rows(); j ++){
			for(int j_ = 0; j_ < j ; j_++)
				loss(iter_ix) += pow(((cluster_assignment(iter_ix,j) == cluster_assignment(iter_ix,j_)) -  pmat(j,j_)),2) ;
		}
	}

	return(loss);
}


//' Estimate the Multivariate nonhomgogenous poisson process intensity function from grouped data via Nested Dirichlet Process
//' 
//' @template bda 
//'
//' @param r matrix of distances associatd with different BEFs
//' @param n_j matrix of integers denoting the start and length of each observations associated BEF distances
//' @param d a 1D grid of positive real values over which the differing intensities are evaluated
//' @param L component truncation number
//' @param K intensity cluster truncation number
//' @param J number of rows in r matrix; number of groups
//' @param mu_0 normal base measure prior mean
//' @param kappa_0 normal base measure prior variance scale 
//' @param nu_0 inverse chi sqaure base measure prior degrees of freedom
//' @param sigma_0 inverse chi square base measure prior scale  
//' @param a_alpha hyperparameter for alpha gamma prior
//' @param b_alpha scale hyperparameter for alpha gamma prior
//' @param a_rho hyperparameter for rho gamma prior
//' @param b_rho scale hyperparameter for rho gamma prior
//' @param iter_max total number of iterations for which to run sampler
//' @param warm_up number of iterations for which to burn-in or "warm-up" sampler
//' @param thin number of iterations to thin by
//' @param seed integer with which to initialize random number generator
//' @param chain integer chain label
//' @param num_posterior_samples the total number of posterior samples after burn in 
//' @seealso the conjugate normal parameterization in the reference below
//'
Rcpp::List nd_nhpp_multivariate_fit(
        const Eigen::ArrayXXd& r,
        const Eigen::MatrixXi& n_j,
        const Eigen::ArrayXd& d,
        const int& L,
        const int& K,
        const int& J,
        const double& mu_0,
        const double& kappa_0,
        const int& nu_0,
        const double& sigma_0,
        const double& a_alpha,
        const double& b_alpha,
        const double& a_rho,
        const double& b_rho,
        const int& iter_max,
        const int& warm_up,
        const int& thin,
        const int& seed,
        const int& chain,
        const int& num_posterior_samples
        )
{
	/*
    // set seed
    std::mt19937 rng;
    rng = std::mt19937(seed);

	const int dim_BEF = r.cols();

    //create sample containers
#include "support/create_mvsample_containers.hpp"

#include "support/during_sampling_containers.hpp"

    // create rng's
    std::gamma_distribution<double> rgam_alpha(a_alpha,b_alpha);
    std::gamma_distribution<double> rgam_rho(a_rho,b_rho);
    std::uniform_real_distribution<double> runif(0,1);
    std::normal_distribution<double> rnorm(0,1);

    //initialize concentration parameters
    alpha = rgam_alpha(rng);
    rho = rgam_rho(rng);
	// sample priors
#include "support/sample_concentration_priors.hpp"

    // initialize component weights
    u = stick_break(L,K,alpha,rng);
    v = stick_break(K,rho,rng);
    w = stick_break_weights(u);
    pi = stick_break_weights(v);

    mu = initialize_mu(L,K,mu_0,kappa_0,rng);
    tau = initialize_tau(L,K,sigma_0,nu_0);


    Rcpp::Rcout << "Beginning Sampling" << std::endl;
    Rcpp::Rcout << "----------------------------------------------------------------------" << std::endl;
    for(int iter_ix = 1; iter_ix <= iter_max; iter_ix ++){
      print_progress(iter_ix,warm_up,iter_max,chain);

      //calculate cluster probabilities
	  // sampler.calculate_q()
      q = dmvnorm(J,r,n_j,pi,w,mu,tau);

	  // sampler.assign_clusters();
      for(int j = 0; j < J; j ++){
            probs = q.row(j);
            std::discrete_distribution<int> d(probs.data(),probs.data()+probs.size());
            iter_cluster_assignment(j) = d(rng);
        }

	  // sampler.calculate_b();
      // calculate w/in cluster component probabilities
        b = dnorm(r,n_j,w,mu,tau,iter_cluster_assignment);

		// sampler.assign_components(rng);
        // assign distances to components within clusters
        component_count = Eigen::ArrayXXi::Zero(L,K);
        for(int j = 0 ; j < J; j ++){
            for(int i = 0 ; i < n_j(j,1) ; i++){
                for(int l = 0 ; l< L; l++)
                    prob(l) = b(n_j(j,0) + i,l);
                std::discrete_distribution<int> d(prob.data(),prob.data() + prob.size());
                iter_component_assignment(n_j(j,0) +i) = d(rng);
                component_count(iter_component_assignment(n_j(j,0)+i),iter_cluster_assignment(j)) += 1;
            }
        }

        // draw samples from intensity cluster parameters
		// sampler.count_clusters();

        for(int k = 0; k < K; k++)
            cluster_count(k) = (iter_cluster_assignment == k).count();

		// sampler.update_cluster_stick_posterior()
        for(int k = 0; k < K; k++){
            v_posterior_beta_alpha(k) = 1 + cluster_count(k);
            v_posterior_beta_beta(k) = alpha + cluster_count.tail(K-k-1).sum();
        }

        v = stick_break(K, v_posterior_beta_alpha,v_posterior_beta_beta,rng);
        pi = stick_break_weights(v);



		// sampler.update_componentr_stick_posterior()
        for(int l = 0; l < L; l ++){
            for(int k = 0; k < K; k++){
                u_posterior_beta_alpha(l,k) = 1 + component_count(l,k);
                u_posterior_beta_beta(l,k) = rho + ((component_count.col(k)).tail(K-k-1)).sum();
            }
        }

        u = stick_break(u_posterior_beta_alpha,u_posterior_beta_beta,rng);
        w = stick_break_weights(u);

        // calculate mu hat (sq) sums
		// sampler.update_mu();

        ycount = Eigen::ArrayXXd::Zero(L,K);
        ycount_sq = Eigen::ArrayXXd::Zero(L,K);
        for(int l = 0 ; l < L; l ++){
          for(int k = 0; k < K; k++){
            for(int j = 0; j < J; j++){
              for(int i = 0; i < n_j(j,1); i ++){
                if(iter_cluster_assignment(j) == k && iter_component_assignment(n_j(j,0)+i) == l){
                  ycount(l,k) += r(n_j(j,0)+i);
                  ycount_sq(l,k) += pow(r(n_j(j,0)+i),2);
                }
              }
            }
          }
        }


        //sample mu via conjugacy

        for(int l = 0; l < L; l ++){
         for(int k = 0; k < K; k++){
           if(component_count(l,k) == 0){
             tau(l,k) = nu_0 * sigma_0 / R::rchisq(nu_0);
             mu(l,k) = rnorm(rng) * sqrt(tau(l,k) / kappa_0 ) + mu_0;
           }
           else{
             s_n  = nu_0 * sigma_0 + (ycount_sq(l,k) - (pow(ycount(l,k),2) / component_count(l,k) ) ) + (kappa_0 * component_count(l,k) / (kappa_0 + component_count(l,k))) * pow( (ycount(l,k) / component_count(l,k) - mu_0 ),2);
             tau(l,k) = s_n / R::rchisq(nu_0 + component_count(l,k));
             mu_n = ( (kappa_0  / tau(l,k)) * mu_0   + ycount(l,k) / tau(l,k)  ) / ( kappa_0 / tau(l,k) + component_count(l,k) / tau(l,k));
             s_n = 1.0 /( (kappa_0 /tau(l,k)) + (component_count(l,k) / tau(l,k)) ) ;
             mu(l,k) =  rnorm(rng) * sqrt(s_n) + mu_n;
           }
         }
       }
		// sampler.update_concentration_parameters()

        // sample concentration parameters
        posterior_b_alpha = 1.0 / b_alpha - (log(1-v.array())).head(K-1).sum();
        posterior_b_rho =  1.0 / b_rho - log(1-u.block(0,0,L-1,K).array()).matrix().colwise().sum().sum();
        std::gamma_distribution<double> rgam_alpha(posterior_a_alpha, 1.0 / posterior_b_alpha);
        std::gamma_distribution<double> rgam_rho(posterior_a_rho, 1.0 / posterior_b_rho);
        alpha = rgam_alpha(rng);
        rho = rgam_rho(rng);

		// sampler.update_psi();)
		
		// sampler.calculate_adjacency_matrix()
		// sampler.store_samples(iter_ix);

        // calculate adjacency matrix and store samples
        if((iter_ix > warm_up) && (iter_ix % thin ==0)){
          for(int j = 0; j < J; j++){
            for(int j_ = 0; j_ < j; j_++)
              cluster_matrix(j,j_) += (iter_cluster_assignment(j) == iter_cluster_assignment(j_)) ? 1: 0;
          }
          for(int k =0; k < K; k++){
             for(int d_ix = 0; d_ix < d_length; d_ix ++){
               for(int l = 0; l < L; l++){
                 intensities(sample_ix, k*d_length + d_ix) += w(l,k) * R::dnorm(d(d_ix),mu(l,k),sqrt(tau(l,k)),false);
			   }
				 global_intensity(sample_ix,d_ix) += pi(k) *  intensities(sample_ix,k*d_length+d_ix);
             }
          }

          cluster_assignment.row(sample_ix) = iter_cluster_assignment;
          cluster_component_assignment.row(sample_ix) = iter_component_assignment;
          pi_samps.row(sample_ix) = pi;
          Eigen::Map<Eigen::RowVectorXd> mu_samp(mu.data(),mu.size());
          Eigen::Map<Eigen::RowVectorXd> w_samp(w.data(),w.size());
          Eigen::Map<Eigen::RowVectorXd> tau_samp(tau.data(),tau.size());
          w_samps.row(sample_ix) = w_samp; // stored in column order
          mu_samps.row(sample_ix) = mu_samp;
          tau_samps.row(sample_ix) = tau_samp;
          alpha_samps(sample_ix,0) = alpha;
          rho_samps(sample_ix,0) = rho;
          sample_ix += 1;
        }
    }

    cluster_matrix = cluster_matrix / num_posterior_samples;

    return(Rcpp::List::create(Rcpp::Named("cluster_assignment") = cluster_assignment,
                              Rcpp::Named("component_assignment") =  cluster_component_assignment,
                              Rcpp::Named("cluster_pair_probability") =  cluster_matrix,
                              Rcpp::Named("pi_samples") = pi_samps,
                              Rcpp::Named("w_samples") =  w_samps,
                              Rcpp::Named("intensities") = intensities,
							  Rcpp::Named("global_intensity") = global_intensity,
                              Rcpp::Named("mu_samples") =  mu_samps,
                              Rcpp::Named("tau_samples") = tau_samps,
                              Rcpp::Named("alpha_samples") = alpha_samps,
                              Rcpp::Named("rho_samples") = rho_samps,
							  Rcpp::Named("alpha_prior") = alpha_prior,
							  Rcpp::Named("rho_prior") = rho_prior
    ));
	*/
}
