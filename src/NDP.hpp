#ifndef NDP
#define NDP

//' Class to implement NDP... eventually
class NDP
{

  public:
    const Eigen::ArrayXd r;
    const Eigen::ArrayXXd n_j;
    const double mu_0;
    const double kappa_0;
    const double sigma_0;
    const int nu_0;
    const int L;
    const int K;
    const double a_alpha;
    const double b_alpha;
    const double a_rho;
    const double b_rho;
    Eigen::ArrayXXd mu;
    Eigen::ArrayXi iter_cluster_assignment;
    Eigen::ArrayXi iter_component_assignment;
    Eigen::ArrayXi cluster_count;
    Eigen::ArrayXXi component_count;
    Eigen::ArrayXXd u_posterior_beta_alpha;
    Eigen::ArrayXXd u_posterior_beta_beta;
    Eigen::ArrayXd v_posterior_beta_alpha;
    Eigen::ArrayXd v_posterior_beta_beta;
    Eigen::ArrayXd probs;
    Eigen::ArrayXd prob;
    Eigen::ArrayXd d;

    NDP(const Eigen::ArrayXd &input_r,
        const Eigen::ArrayXXi &input_n_j,
        const double &input_mu_0,
        const double &input_kappa_0,
        const int &input_nu_0,
        const double &input_sigma_0,
        const int &input_L,
        const int &input_K,
        const double &input_a_alpha,
        const double &input_b_alpha,
        const double &input_a_rho,
        const double &input_b_rho,
        const Eigen::ArrayXd &input_d){
      
      r = input_r;
      n_j = input_n_j;
      mu_0 = input_mu_0;
      kappa_0 = input_kappa_0;
      sigma_0 = input_sigma_0;
      nu_0 = input_nu_0;
      L = input_L;
      K = input_K;
      a_alpha = input_a_alpha;
      b_alpha = input_b_alpha;
      a_rho = input_a_rho;
      b_rho = input_b_rho;
      d = input_d;

    
    }

	void stick_break_pi();


}

#endif
