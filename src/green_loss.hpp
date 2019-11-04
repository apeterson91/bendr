//' Computes Green and Lau loss function with unknown classification
//'
//' @param cluster_assignment iter_total x J cluster assignment matrix
//' @param pmat J x J pairwise probability of co-clustering matrix
//' @param tau penalty parameter 
// [[Rcpp::export]]
Eigen::ArrayXd green_loss_unknown(const Eigen::ArrayXXi &cluster_assignment,
							 const Eigen::ArrayXXd &pmat,
							 const double &tau){


	Eigen::ArrayXd loss(cluster_assignment.rows());
	loss = Eigen::ArrayXd::Zero(cluster_assignment.rows());

	for(int iter_ix = 0; iter_ix < cluster_assignment.rows() ; iter_ix ++){
		for(int j = 0; j < pmat.rows(); j ++){
			for(int j_ = 0; j_ < j ; j_++)
				loss(iter_ix) += cluster_assignment(iter_ix,j) == cluster_assignment(iter_ix,j_) ? (pmat(j,j_) - tau) : 0.0;
		}
	}

	return(loss);
}


//' Computes Green and Lau Loss function with known classification
//'
//' @param cluster_assignment iter_total x J cluster assignment matrix
//' @param pmat J x J pairwise probability of co-clustering matrix
//' @param true_cluster_assignment J x J true Adjacency Matrix
//' @param a mis-classification penalty parameter
//' @param b classification penalty parameter
// [[Rcpp::export]]
Eigen::ArrayXd green_loss_known(const Eigen::ArrayXXi &cluster_assignment,
								const Eigen::ArrayXXd &pmat,
								const Eigen::ArrayXXi &true_cluster_assignment,
								const double &a,
								const double &b
								 ){


	Eigen::ArrayXd loss(cluster_assignment.rows());
	loss = Eigen::ArrayXd::Zero(cluster_assignment.rows());

	for(int iter_ix = 0; iter_ix < cluster_assignment.rows() ; iter_ix ++){
		for(int j = 0; j < pmat.rows(); j ++){
			for(int j_ = 0; j_ < j ; j_++){
				loss(iter_ix) += ((cluster_assignment(iter_ix,j) != cluster_assignment(iter_ix,j_)) && (true_cluster_assignment(j,j_) == 1 )  ) ? a : 0.0;
				loss(iter_ix) += ((cluster_assignment(iter_ix,j) == cluster_assignment(iter_ix,j_)) && (true_cluster_assignment(j,j_) == 0 ) ) ? b : 0.0;
			}
		}
	}

	return(loss);
}


//' Computes Square loss with unknown classification
//'
//' @param cluster_assignment iter_total x J cluster assignment matrix
//' @param pmat J x J pairwise probability of co-clustering matrix
// [[Rcpp::export]]
Eigen::ArrayXd square_error(const Eigen::ArrayXXi &cluster_assignment,
							const Eigen::ArrayXXd &pmat){


	Eigen::ArrayXd loss(cluster_assignment.rows());
	loss = Eigen::ArrayXd::Zero(cluster_assignment.rows());

	for(int iter_ix = 0; iter_ix < cluster_assignment.rows() ; iter_ix ++){
		for(int j = 0; j < pmat.rows(); j ++){
			for(int j_ = 0; j_ < j ; j_++)
				loss(iter_ix) += pow(((cluster_assignment(iter_ix,j) == cluster_assignment(iter_ix,j_)) -  pmat(j,j_)),2) ;
		}
	}

	return(loss);
}
