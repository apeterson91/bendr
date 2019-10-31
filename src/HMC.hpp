Eigen::VectorXd rnorm_draw(int &n, std::mt19937 &rng){

	Eigen::VectorXd out(n);
	std::normal_distribution<double> rnorm(0.0,1.0);

	for (int i = 0; i < n ; i++){
		out(i) = rnorm(rng);
	}
	return(out);
}

Eigen::VectorXd runif_draww(int &n std::mt19937 &rng){

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
		ratio = isinf(-ratio) ? - DBL_MAX: ratio;
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

#include "Tree.hpp"


//' Returns draw of beta from posterior distribution via No U-Turn Sampler
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
