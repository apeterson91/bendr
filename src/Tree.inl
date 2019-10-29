
void Tree::buildTree(
		NHPP &model,
		Eigen::VectorXd &new_beta,
		Eigen::VectorXd &new_momenta,
		double &u, int v, int j, 
        double &epsilon, 
		Eigen::VectorXd &beta_naught,
		Eigen::VectorXd &momenta_naught,
		std::mt19937 &rng){

    if( j == 0 ){
        energy_init = model.calculate_energy(beta_naught,momenta_naught);
        leapfrog(model,beta_new,momenta_new,v*epsilon);
        energy_new = model.calculate_energy(beta_new,momenta_new);
        n_prime = u <= energy_new ? 1: 0;
        s_prime = u < (1000 + energy_new) ? 1:0;
		beta_left = beta_new;
		beta_right = beta_new;
        alpha_prime = std::min(1.0,exp(energy_new - energy_init));
        n_alpha = 1.0;
    }else{
		std::uniform_real_distribution<double> die(0.0,1.0);
        Tree subtree;
        subtree.buildTree(model,new_beta,new_momenta,u,v,j-1,epsilon,beta_naught,momenta_naught,rng);
        s_prime = subtree.get_s_prime();
        n_prime = subtree.get_n_prime();
        n_alpha = subtree.get_n_alpha();
		beta_new = subtree.get_bn();
		beta_left = subtree.get_bl();
		momenta_left = subtree.get_ml();
		beta_right = subtree.get_br();
		momenta_right = subtree.get_mr();
        alpha_prime = subtree.get_alpha_prime();
        if(subtree.get_s_prime() == 1){
            Tree subsubtree;
            if( v == -1 ){
                subsubtree.buildTree(model, beta_left, momenta_left, u, v, j-1, epsilon, beta_naught, momenta_naught, rng);
				beta_left = subsubtree.get_bl();
				momenta_left = subsubtree.get_ml();
            }else{
                subsubtree.buildTree(model, beta_right, momenta_right, u, v, j-1, epsilon, beta_naught, momenta_naught, rng);
				beta_right = subsubtree.get_br();
				momenta_right = subsubtree.get_ml();
            }
            double p = (subsubtree.get_n_prime() == 0.0 && subtree.get_n_prime() == 0.0) ? 0.0 : subsubtree.get_n_prime() / (subtree.get_n_prime() + subsubtree.get_n_prime());
            std::uniform_real_distribution<double> die(0.0,1.0);
            if(die(rng) <= p)
				beta_new = subsubtree.get_bn();
            alpha_prime = subsubtree.get_alpha_prime() + subtree.get_alpha_prime();
            n_alpha = subtree.get_n_alpha() +  subsubtree.get_n_alpha();
            bool UTI_one = (beta_right - beta_left).dot(momenta_left) >= 0;
            bool UTI_two = (beta_right- beta_left).dot(momenta_right) >= 0; 
            s_prime = (UTI_one && UTI_two ) ? subsubtree.get_s_prime() : 0 ;
            n_prime = subtree.get_n_prime() + subsubtree.get_n_prime();
        }
    }
}

void Tree::leapfrog(NHPP &model,
		Eigen::VectorXd &beta,
		Eigen::VectorXd &momenta,
		double epsilon){

    beta_grad = model.calculate_gradient(beta);

	momenta_new = momenta_new + epsilon * beta_grad / 2.0;

	beta_new = beta + epsilon * momenta_new;
	
	beta_grad = model.calculate_gradient(beta);

	momenta_new  = momenta_new + epsilon * beta_grad / 2.0;

}
