Eigen::VectorXd rnorm_draw(int n, std::mt19937 &rng){

	Eigen::VectorXd out(n);
	std::normal_distribution<double> rnorm(0.0,1.0);

	for (int i = 0; i < n ; i++){
		out(i) = rnorm(rng);
	}
	return(out);
}

class NHPP
{
	public:
		Eigen::MatrixXd X;
		Eigen::ArrayXd n_j;
		Eigen::ArrayXd eta;
		Eigen::VectorXd beta_grad;
  	    NHPP(Eigen::VectorXd &input_X,
			 Eigen::ArrayXd &input_n_j){
			X = input_X;
			n_j = input_n_j;
	}

	double calculate_energy(Eigen::VectorXd &beta,
							Eigen::VectorXd &momenta){

		double out = 0;
		eta = X * beta;
		out = (n_j * eta  + exp(eta)).sum();
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

		beta_grad = (n_j * X.array()).rowwise().sum() + (X.array() * exp((X*beta).array())).rowwise().sum();

		return(beta_grad);

	}

	double FindReasonableEpsilon(Eigen::VectorXd &beta){

		double epsilon = 1.0;


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
		Eigen::VectorXd &input_X,
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
	int n = 1;
	int v = 1;
	int j = 1;
	int accept = 0;
	// beta parameters
	Eigen::MatrixXd beta_draws(iter_max - warm_up,input_X.cols());
	Eigen::VectorXd acceptance(iter_max - warm_up);
	acceptance = Eigen::VectorXd::Zero(iter_max - warm_up);
	Eigen::VectorXd return_beta;
	Eigen::VectorXd current_beta = rnorm_draw(input_X.cols(),rng);
	// initialize model
	NHPP model(input_X,input_n_j);
	double epsilon = 0;
	epsilon = model.FindReasonableEpsilon(current_beta);
	// epsilon tuning parameters
	double epsilon_bar = 1;
	double mu = log(10 *epsilon);
	double H_bar = 0;
	double kappa = 0.75;
	double gamma = 0.05;
	double t_0 = 10;
	double sample_ix = 0;
	// random number engine
	std::uniform_real_distribution<double> runif(0.0,1.0);

	for(int iter_ix = 1; iter_ix <= iter_max; iter_ix ++){
		Eigen::VectorXd beta_left = current_beta;
		Eigen::VectorXd beta_right = current_beta;
		Eigen::VectorXd momenta = rnorm_draw(current_beta.size(),rng);
		Eigen::VectorXd momenta_left = momenta;
		Eigen::VectorXd momenta_right = momenta;

		accept = 0;
		n = 1;
		s = 1;
		u = model.sample_u(current_beta,momenta,rng);
		Tree tree;
		while(s == 1){
			if(runif(rng)<=1){
				tree.buildTree(model,beta_left,momenta_left,u,v,j,epsilon,current_beta,momenta,rng);
				beta_left = tree.get_bl();
				momenta_left = tree.get_ml();
			}else{
				tree.buildTree(model,beta_right,momenta_right,u,v,j,epsilon,current_beta,momenta,rng);
				beta_right = tree.get_br();
				momenta_right = tree.get_mr();
			}
			if(tree.get_s_prime() == 1){
				if(runif(rng)<= tree.get_n_prime() / n){
					accept = 1;
					return_beta = tree.get_bn();
				}
			}
			n += tree.get_n_prime();
			s = tree.get_s_prime() * ((beta_right - beta_left).dot(momenta_left) >=0 )  * ((beta_right - beta_left).dot(momenta_right) >= 0) ;
			j += 1;
		}
		if(iter_ix <= warm_up){
			H_bar = ( 1 - 1.0 / (iter_ix + t_0)) * H_bar;
			H_bar +=  1.0 / (iter_ix + t_0) * (adapt_delta - tree.get_alpha_prime() / tree.get_n_alpha());
			epsilon = exp(mu - sqrt(iter_ix)/gamma * H_bar);
			epsilon_bar = exp(pow(iter_ix,-kappa) * log(epsilon)) + (1 - pow(iter_ix,-kappa))*log(epsilon_bar);
		}else{
			epsilon = epsilon_bar;
			// store sample
			if(accept == 1){
				current_beta = return_beta;
				beta_draws.row(sample_ix) = return_beta.transpose();
				sample_ix += 1;
			}
		}
	}

	return(Rcpp::List::create(Rcpp::Named("beta_samples") = beta_draws,
							  Rcpp::Named("acceptance") = acceptance));
}
