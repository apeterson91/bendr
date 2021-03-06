Eigen::ArrayXXd u(L,K);
Eigen::ArrayXd v(K);
Eigen::ArrayXXd w(L,K);
Eigen::ArrayXd pi(K);
Eigen::ArrayXXd q(J,K);
Eigen::ArrayXXd b(r.rows(),L);
Eigen::ArrayXXd mu(L,K);
Eigen::ArrayXXd tau(L,K);
Eigen::ArrayXi iter_cluster_assignment(J);
Eigen::ArrayXi iter_component_assignment(r.rows()); 
Eigen::ArrayXi cluster_count(K);
Eigen::ArrayXXi component_count(L,K);
Eigen::ArrayXXd u_posterior_beta_alpha(L,K);
Eigen::ArrayXXd u_posterior_beta_beta(L,K);
Eigen::ArrayXd v_posterior_beta_alpha(K);
Eigen::ArrayXd v_posterior_beta_beta(K);
Eigen::ArrayXd probs(K);
Eigen::ArrayXd prob(L);
Eigen::ArrayXXd ycount(L,K);
Eigen::ArrayXXd ycount_sq(L,K);
const int d_length = d.size();
double mu_n;
double s_n;
double alpha;
double rho;
double posterior_a_alpha = a_alpha + (K-1);
double posterior_b_alpha;
double posterior_a_rho = a_rho + K*(L-1);
double posterior_b_rho;
