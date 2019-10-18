
//' computes green loss function
//' @param cluster_assignment iter_total x J cluster assignment matrix
//' @param pmat J x J pairwise probability of co-clustering matrix
//' @param chain
// [[Rcpp::export]]
Eigen::ArrayXd green_loss_engine(const Eigen::ArrayXXi &cluster_assignment,
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
