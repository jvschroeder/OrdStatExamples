//#include<vector>
{
	const T zero = v1[0]*0;
	const T one = zero+1;
	std::vector< std::vector<T> > res(n1+1, std::vector<T>(n2+1,one));
	std::vector< std::vector<T> > M = res;
	std::vector<T> M0(n1+1,one);
	std::vector<T> partial_result(n2+1,zero);
	double i=0;
	for(int m1=0; m1<=n1; m1++) {
	  //Rcpp::Rcout << i/((n1+1)*(n2 + 1)) << std::endl;
		for(int m2=0; m2<=n2; m2++) {
		  i++;
		  if(m1+m2<=max_idx) {
  			for(int k1=0; k1 <= m1; k1++) {
				#ifdef ALGORITHM_PARALLEL
				tbb::parallel_for(0, m2+1, [&] (int k2)
				#else
				for(int k2=0; k2 <= m2; k2++)
				#endif
				{
  					if(k1 < m1 || k2 < m2){
						#ifdef ALGORITHM_PARALLEL
							partial_result[k2] = M[k1][k2] * res[k1][k2];
						#else
							res[m1][m2] -= M[k1][k2] * res[k1][k2];
						#endif
  					}
  					//add column to M - Takes M^{(m1,m2)} and returns M^{(m1,m2+1)}
					if(m2 < n2 && k1+k2 < max_idx) M[k1][k2] = M[k1][k2] * (m2+1)/(m2+1-k2) * (1-v2[k1+k2]);
  				}
				#ifdef ALGORITHM_PARALLEL
				);
				typedef typename std::vector<T>::const_iterator vec_it;
				typedef tbb::blocked_range< vec_it > range_type;
				vec_it end = partial_result.cbegin();
				end += m2;
				res[m1][m2] -= tbb::parallel_reduce(
				  range_type(partial_result.cbegin(),end),
				  zero,
				  [](range_type const& r, T value)->T {
				    return std::accumulate(r.begin(),r.end(),value);
				  },
				  std::plus<T>());
				#endif
  				//add column to M - Takes M^{(m1,m2)} and returns M^{(m1,m2+1)}
  				if(m2 < n2 && k1 < m1 && k1+m2+1 < max_idx) {
  					M[k1][m2+1] = M[k1+1][m2] * (k1+1)/(m1-k1) * (1-v1[k1+m2+1]);
  				}
  			}
		  } else {
		    res[m1][m2] = -1;
		  }
		}
		//add row to M0 going from M^{(m1,0)} to M^{(m1+1,0)}
		if(m1 < n1) {
			for(int k1=0; k1 <= m1; k1++) {
				M0[k1] = M0[k1] * (m1+1)/(m1+1-k1) * (1-v1[k1]);
				M[k1][0] = M0[k1];
			}
		}
	}
	return res;
}