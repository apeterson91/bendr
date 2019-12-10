
class NDP
{
	private:
		const Eigen::ArrayXXd r;
		Eigen::ArrayXd q;
		Eigen::ArrayXXd zeta;
		Eigen::ArrayXXd b;
		Eigen::ArrayXd pi;
		Eigen::ArrayXXd w;
		Eigen::VectorXd mu_vector;
		Eigen::VectorXd r_vector;
		std::vector mu;
		std::vector sigma;
		Eigen::ArrayXXd v;
		Eigen::ArrayXd u;
		double alpha;
		double rho;
		double psi;
		const Eigen::MatrixXi n_j;
		const Eigen::ArrayXi n_j_arr;
		const Eigen::ArrayXd d;
		const int L;
		const int K;
		const int J;
		const int dim_BEF;
		const Eigen::VectorXd mu_0;
		const Eigen::ArrayXd nu_0;
		const Eigen::ArrayXd sigma_0;
		const double a_alpha;
		const double b_alpha;
		const double a_rho;
		const double b_rho;
		const int iter_max;
		const int warm_up;
		const int thin;
		const int num_posterior_samples;

	public:

		NDP(const Eigen::ArrayXXd &input_r,
			const Eigen::MatrixXi &input_n_j,
			const Eigen::ArrayXi &input_n_j_arr,
			const Eigen::ArrayXd &input_d,
			const int &input_L,
			const int &input_K,
			const int &input_J,
			const int &input_dim_BEF,
			const Eigen::VectorXd &input_mu_0,
			const Eigen::ArrayXd &input_nu_0,
			const Eigen::ArrayXd &input_sigma_0,
			const double &input_a_alpha,
			const double &input_b_alpha,
			const double &input_a_rho,
			const double &input_b_rho,
			const int &input_iter_max,
			const int &input_warm_up,
			const int &input_thin, 
			const int &input_num_posterior_samples):
			r(input_r),
			n_j(input_n_j),
			n_j_arr(input_n_j_arr),
			d(input_d),
			L(input_L),
			K(input_K),
			J(input_J),
			dim_BEF(input_dim_BEF),
			mu_0(input_mu_0),
			nu_0(input_nu_0),
			sigma_0(input_sigma_0),
			a_alpha(input_a_alpha),
			b_alpha(input_b_alpha),
			a_rho(input_a_rho),
			b_rho(input_b_rho),
			iter_max(input_iter_max),
			warm_up(input_warm_up),
			thin(input_thin),
			num_posterior_samples(input_num_posterior_samples){

				for(int i = 0; i < dim_BEF; i++){
					mu.push_back(Eigen::ArrayXXd::Zero(L,K));
					sigma.push_back(Eigen::ArrayXXd::Zero(L,K));
				}

			}

		void initialize_pars(std::mt19937 &rng);

		void calculate_q();

		void assign_clusters();

		void calculate_b();

		void assign_components(std::mt19937 &rng);

		void count_clusters();

		void update_cluster_stick_posterior();

		void update_component_stick_posterior();

		void update_mu();

		void update_psi();

		void calculate_adjacency_matrix();

		void store_samples(const int iter_ix);

		//auxiliary functions
		
		void cluster_stick_break(std::mt19937 &rng);

		void component_stick_break(std::mt199937 &rng);

		void initialize_mu(std::mt199937 &rng);

		void initialize_sigma();

		void initialize_concentration(std::mt19937 &rng);

		void covariance_function(double &r_1, double &r_2);

};

#include "NDP.inl"
