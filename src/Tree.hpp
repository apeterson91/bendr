#include<random>

class Tree
{
    private :
		Eigen::VectorXd beta_left;
		Eigen::VectorXd beta_right;
		Eigen::VectorXd beta_new;
		Eigen::VectorXd momenta_right;
		Eigen::VectorXd momenta_left;
		Eigen::VectorXd momenta_new;
		Eigen::VectorXd beta_grad;
        double n_prime;
		double energy_init;
		double energy_new;
        int s_prime;
        double alpha_prime;
        double n_alpha;

    public:
		Tree(int p,std::mt19937 &rng){
			momenta_new = rnorm_draw(p,rng);
			momenta_left = momenta_new;
			momenta_right = momenta_new;
		}
        void buildTree(
				NHPP &model,
				Eigen::VectorXd &new_beta,
				Eigen::VectorXd &new_momenta,
                double& u, int v, int j,
                double &epsilon,
				Eigen::VectorXd &beta_naught,
				Eigen::VectorXd &momenta_naught,
                std::mt19937 &rng);

        void leapfrog(
				NHPP &model,
				Eigen::VectorXd &beta,
				Eigen::VectorXd &momenta,
				double epsilon);

        const int get_s_prime() const{
            return(s_prime);
        }

        const double get_n_prime() const{
            return(n_prime);
        }

        const double get_alpha_prime() const{
            return(alpha_prime);
        }

        const double get_n_alpha() const{
            return(n_alpha);
        }

		Eigen::VectorXd get_bl() const{
            return(beta_left);
        }

		Eigen::VectorXd get_ml() const{
			return(momenta_left);
		}

        Eigen::VectorXd get_br() const{
            return(beta_right);
        }

		Eigen::VectorXd get_mr() const{
			return(momenta_right);
		}

		Eigen::VectorXd get_bn() const{
            return(beta_new);
        }

};


#include "Tree.inl"
