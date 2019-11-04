//' Estimate the nonhomgogenous poisson process intensity function from grouped data with fixed concentration parameters
//'
//'
//'
//' @param r vector of distances associatd with different BEFs
//' @param n_j matrix of integers denoting the start and length of each observations associated BEF distances
//' @param d a 1D grid of positive real values over which the differing intensities are evaluated
//' @param L component truncation number
//' @param K intensity cluster truncation number
//' @param J number of rows in r matrix; number of groups
//' @param iter_max total number of iterations for which to run sampler
//' @param warm_up number of iterations for which to burn-in or "warm-up" sampler
//' @param thin number of iterations to thin by
//' @param seed integer with which to initialize random number generator
//' @param chain integer chain label
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
    #include "create_sample_containers.hpp"

    #include "during_sampling_fixed.hpp"

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
