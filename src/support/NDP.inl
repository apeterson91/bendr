void initialize_concentration(std::mt19937 &rng){

	std::gamma_distribution<double> rgamma(a_alpha,b_alpha);

	alpha = rgamma(rng);

	std::gamma_distribution<double> rgamma(a_rho,b_rho);

	rho = rgamma(rng);

}

void NDP::initialize_sigma(std::mt19937 &rng){

	for(int i = 0; i < dim_BEF; i ++){
		for(int l= 0; l < L; l ++){
			for(int k = 0; k < K; k++)
				sigma[i](l,k) = (sigma_0(i) * nu_0(i)) / R::rchisq(nu_0(i));
		}
	}

}


void NDP::initialize_mu(std::mt19937 &rng){

	std::normal_distribution<double> rnorm(0,1);

	for(int i = 0; i < dim_BEF; i ++){
		for(int l = 0; l < L; l++){
			for(int k = 0; k < K; k++)
				mu[0](l,k) = rnorm(rng) * sigma[i](l,k) + mu(i)
		}
	}

}

void cluster_stick_break(std::mt19937 &rng){

	u = Eigen::ArrayXd::Zero(K);
	sftrabbit::beta_distribution<double> rbeta(1.0,alpha);
	for(int k = 0 ; k < K-1; k ++)
		u(k) = beta_distribution(rng);

	u(K-1) = 1;

	for(int k = 0; k < K-1; k++ )
		pi(k) = k ==0 ? u(k) : u(k) * (1-u.head(k)).prod();

	pi(k-1) = u(k-1) * (1 - u.head(K-1)).prod();

}

void component_stick_break(std::mt19937 &rng){

	v = Eigen::ArrayXXd::Zero(L*K);
	sftrabbit::beta_distribution<double> rbeta(1.0,rho);

	for(int l = 0; l < L; l ++){
		for(int k = 0; k < K; k ++)
			v(l,k) = l == (L-1) ? 1.0 : rbeta(rng);
	}

	for(int l = 0; l < L; l ++){
		for(int k = 0; k < K; k ++)
			w(l,k) =  row_ix == 0 ? v(l,k) : v(l,k) * (1 -v.block(0,k,l,1)).prod(); 
	}

}


void NDP::initialize_pars(std::mt19937 &rng){

	initialize_concentration();

	initialize_sigma();

	initialize_mu();

	cluster_stick_break();

	component_stick_break();

}

void NDP::create_Sigma(const int j){

	Eigen::MatrixXd Sigma(n_j_arr(j),n_j_arr(j)); 
	int row_start = 0;
	int col_start = 0;
	// block diagonal independence
	for(int i = 0; i < dim_BEF; i++){
		Sigma.block(row_start,col_start,n_j(j,(i+1)*2),n_j(j,(i+1)*2)) = Eigen::Matrix::Identity(n_j(j,(i+1)*2),n_j(j,(i+1)*2));
		row_start = n_j(j,(i+1)*2);
		col_start = n_j(j,(i+1)*2);
	}
	// off diagonal dependence
	// to be implemented later
}

void NDP::adjust_Sigma(const int j, const int k, const int l){

	Eigen::MatrixXd Sigma(n_j_arr(j),n_j_arr(j)); 
	int row_start = 0;
	int col_start = 0;
	for(int i = 0; i < dim_BEF; i ++){
		Sigma.block(row_start,col_start,n_j(j,(i+1)*2),n_j(j,(i+1)*2)) = pow(sigma[i](l,k),2) * Eigen::Matrix::Identity(n_j(j,(i+1)*2),n_j(j,(i+1)*2));
		row_start = n_j(j,(i+1)*2);
		col_start = n_j(j,(i+1)*2);
	}




}

void NDP::calculate_q(){

	q = Eigen::ArrayXXd::Zero(J,K);
	for(int j = 0; j < J; j ++){
			create_r_vector(j);
		for(int k = 0; k < K; k ++){
			for(int l = 0; l < L; l++){
				create_mu_vector(j,k,l);
				adjust_Sigma(j,k,l);
				q(j,k) += w(l,k) * Sigma.inv().det() + exp(- .5 * ( r_vector - mu_vector).transpose() * Sigma.inv() * ( r_vector - mu_vector )  );
			}
			q(j,k) *= pi(k);
		}
	}

}


void NDP::assign_clusters(){



}

void NDP::calculate_b(){

	b = Eigen::ArrayXXd::Zero(r.rows(),L);
	
	for(int l = 0; l < L; l ++){
		for(int j = 0; j < J; j ++){
			for(int i = 0; i < n_j(j,2*dim_BEF) ; i ++){
				create_mu_vector(j,k,l);
				adjust_Sigma(l,zeta(j),k);
				b(n_j(j,0) +i, l) = w(l,zeta(j)) * Sigma.inv().det() * exp( - .5 * (r_vector - mu_vector).transpose() * Sigma.inv() * (r_vector - mu_vector) );
			}
		}
	}
	
}

