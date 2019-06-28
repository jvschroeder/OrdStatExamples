#ifdef FAST_POW
	#define pow PairArithmetic::fast_pow<T>
#else
	#define pow PairArithmetic::pow<T>
#endif
//#include<vector>
{
	const T zero = v1[0]*0;
	std::vector< std::vector<T> > res(n1+1, std::vector<T>(n2+1,zero));
	for(int m1=0; m1<=n1; m1++) {
		for(int m2=0; m2<=n2; m2++) {
			if(m1+m2<=max_idx) {
				T v1_ = zero, v2_ = zero;
				if(m1+m2>0){
					v1_ = v1[m1+m2-1];
					v2_ = v2[m1+m2-1];
				}
				res[m1][m2] = pow(v1_,m1) * pow(v2_,m2);
				#ifdef ALGORITHM_PARALLEL
				std::vector<T> partial_result(std::max(0,m1+m2-1),zero);
				tbb::parallel_for(0, m1+m2-1, [&] (int j)
				#else
				for(int j = 0; j <= m1+m2-2; j++)
				#endif
				{
					int k1 = j > m2 ? j-m2 : 0;
					int k2 = j > m2 ? m2 : j;
					T d1 = v1_-v1[k1+k2];
					const T d2 = v2_-v2[k1+k2];
					T m = pow(d1,m1-k1) * pow(d2,m2-k2);
					T a = binomial<T>(m1,k1)*binomial<T>(m2,k2);
					if(d1 == 0) d1 = 1;
					for(; k1 <= m1 && k2 >= 0; k1++, k2--) {
						/*res[m1][m2] -= binomial<T>(m1,k1)*binomial<T>(m2,k2) *
											PairArithmetic::pow<T>(d1,m1-k1) *
											PairArithmetic::pow<T>(d2,m2-k2) *
											res[k1][k2];*/
						#ifdef ALGORITHM_PARALLEL
						partial_result[j] += a * m * res[k1][k2];
						#else
						res[m1][m2] -= a * m * res[k1][k2];
						#endif
						m *= d2 / d1;
						a /= k1+1;
						a *= m1-k1;
						a /= m2-k2+1;
						a *= k2;
					}
				}
				#ifdef ALGORITHM_PARALLEL
				);
				typedef typename std::vector<T>::iterator vec_it;
				typedef tbb::blocked_range< vec_it > range_type;
				res[m1][m2] -= tbb::parallel_reduce(
				  range_type(partial_result.begin(),partial_result.end()),
				  zero,
				  [](range_type const& r, T value)->T {
				    return std::accumulate(r.begin(),r.end(),value);
				  },
				  std::plus<T>());
				#endif
			}
		}
	}
	return res;
}

#undef pow