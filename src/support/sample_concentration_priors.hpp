Eigen::ArrayXd alpha_prior(num_posterior_samples);
Eigen::ArrayXd rho_prior(num_posterior_samples);

for(int i = 0; i < num_posterior_samples; i ++){
	alpha_prior(i) = rgam_alpha(rng);
	rho_prior(i) = rgam_rho(rng);
}
	

