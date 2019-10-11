
//create sample containers
Eigen::MatrixXd beta_samps(num_posterior_samples,X.cols());
Eigen::MatrixXd w_samps(num_posterior_samples,L*K);
Eigen::MatrixXd mu_samps(num_posterior_samples,L*K);
Eigen::MatrixXd pi_samps(num_posterior_samples,K);
Eigen::MatrixXi cluster_assignment(num_posterior_samples,J);
Eigen::MatrixXi cluster_component_assignment(num_posterior_samples,r.rows());
Eigen::MatrixXd cluster_matrix(J,J);
Eigen::MatrixXd alpha_samps(num_posterior_samples,1);
Eigen::MatrixXd rho_samps(num_posterior_samples,1);
Eigen::MatrixXd tau_samps(num_posterior_samples,L*K);
Eigen::MatrixXd intensities(num_posterior_samples,K * d.size());
Eigen::MatrixXd global_intensity(num_posterior_samples,d.size());
Eigen::ArrayXd alpha_prior(num_posterior_samples);
Eigen::ArrayXd rho_prior(num_posterior_samples);
cluster_matrix = Eigen::MatrixXd::Zero(J,J);
int sample_ix = 0;
intensities = Eigen::MatrixXd::Zero(num_posterior_samples,K * d.size());
global_intensity = Eigen::MatrixXd::Zero(num_posterior_samples,d.size());
