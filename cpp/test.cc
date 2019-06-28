#include<Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include<RcppParallel.h>
// [[Rcpp::plugins(cpp11)]]

#include "gmpxx.h"

#include<cmath>

#include<map>

#undef ERROR
#undef WARN
#define CATCH_CONFIG_RUNNER
#include<catch.hpp>

#include "misc/helpers.hpp"
#define NO_K
#include "PairArithmetic/PairArithmetic.hpp"
#include "misc/OpCounter.hpp"

#include "algorithms/bolshev_lower.hpp"
#include "algorithms/steck_lower.hpp"

template<typename T>
std::vector< std::vector<T> > bolshev2_lower(const T* v1, const T* v2, const int n1, const int n2, const int max_idx)
#include "algorithms/bolshev2_lower.hpp"

template<typename T>
std::vector< std::vector<T> > noe2_lower(const T* v1, const T* v2, const int n1, const int n2, const int max_idx)
#include "algorithms/noe2_lower.hpp"

template<typename T>
std::vector< std::vector<T> > steck2_lower_slow(const T* v1, const T* v2, const int n1, const int n2, const int max_idx)
#include "algorithms/steck2_lower.hpp"

#define FAST_POW
template<typename T>
std::vector< std::vector<T> > steck2_lower(const T* v1, const T* v2, const int n1, const int n2, const int max_idx)
#include "algorithms/steck2_lower.hpp"

#define ALGORITHM_PARALLEL
template<typename T>
std::vector< std::vector<T> > noe2_lower_p(const T* v1, const T* v2, const int n1, const int n2, const int max_idx)
#include "algorithms/noe2_lower.hpp"

template<typename T>
std::vector< std::vector<T> > noe2_lower_approx(const T* v1_, const T* v2_, const T* v1, const T* v2, const int n1, const int n2, const int max_idx)
#include "algorithms/noe2_lower_approx.hpp"

template<typename T>
std::vector< std::vector<T> > bolshev2_lower_p(const T* v1, const T* v2, const int n1, const int n2, const int max_idx)
#include "algorithms/bolshev2_lower.hpp"

template<typename T>
std::vector< std::vector<T> > steck2_lower_p(const T* v1, const T* v2, const int n1, const int n2, const int max_idx)
#include "algorithms/steck2_lower.hpp"

template<typename T>
std::vector< std::vector<T> > steck2_lower_b(const T* v1, const T* v2, const int n1, const int n2, const int max_idx)
#include "algorithms/steck2_lower_b.hpp"

template<typename T>
T val_for_linear_threshold(T alpha, const int i, const int n) {
	if(i == 0) return 1;
	T ret = alpha;
	ret = ret / n;
	ret *= PairArithmetic::fast_pow<T>(ret * (i+1), i-1);
	return ret;
}

template<typename T>
std::vector<T> linear_thresholds(T alpha, const int n) {
	std::vector<mpq_class> v(n);
	for(int i = 1; i <= n; i++) {
		v[i-1] = alpha;
		v[i-1] = v[i-1]  * i;
		v[i-1] = v[i-1]  / n;
	}
	return v;
}
#include<limits>
TEST_CASE( "Exact values of the recursions", "[exact]" ) {
	std::vector<int> ns(20);
	std::iota(ns.begin(), ns.end(), 1);
	std::vector<mpq_class> alphas;
	alphas.push_back(mpq_class(1,2));
	alphas.push_back(mpq_class(1,4));
	alphas.push_back(mpq_class(1,10));
	alphas.push_back(mpq_class(1,20));
	alphas.push_back(mpq_class(1,33));
	for(auto alpha : alphas) {
		for(auto n : ns) {
			std::vector<mpq_class> v = linear_thresholds(alpha,n);
			std::vector<PairArithmetic::DoublePair> v1(n);
			for(int i=0; i < n; i++) v1[i] = v[i].get_d();
			std::vector<mpq_class> bolshev_res = bolshev_lower(v.data(),n);
			std::vector<mpq_class> steck_res = steck_lower(v.data(),n);
			std::vector< std::vector<mpq_class> > bolshev2_res = bolshev2_lower(v.data(),v.data(),n,n,n);
			std::vector< std::vector<mpq_class> > noe2_res = noe2_lower(v.data(),v.data(),n,n,n);
			std::vector< std::vector<PairArithmetic::DoublePair> > noe2f_res = noe2_lower(v1.data(),v1.data(),n,n,n);
			for(int i = 0; i <= n; i++) {
				mpq_class exact = val_for_linear_threshold<mpq_class>(alpha,i,n);
				REQUIRE( (exact == bolshev_res[i]) );
				REQUIRE( (exact == steck_res[i]) );
				REQUIRE( (exact == bolshev2_res[i][0]) );
				REQUIRE( (exact == bolshev2_res[0][i]) );
				REQUIRE( (exact == noe2_res[i][0]) );
				REQUIRE( (exact == noe2_res[0][i]) );
				for(int j = 0; j <= n; j++) {
					if(i+j>n) {
						REQUIRE( bolshev2_res[i][j].get_d() == -1 );
						REQUIRE( noe2_res[i][j].get_d() == -1 );
					}
					else {
						mpq_class exact1 = val_for_linear_threshold<mpq_class>(alpha,i+j,n);
						REQUIRE( ( exact1 == bolshev2_res[i][j] ) );
						REQUIRE( ( exact1 == noe2_res[i][j] ) );
						REQUIRE( std::abs(exact1.get_d() - (double)noe2f_res[i][j]) / exact1.get_d() <= std::numeric_limits<double>::round_error() );
					}
				}
			}
		}
	}
}

Catch::Session session;

// [[Rcpp::export]]
int runTests() {
	return session.run();
}

// [[Rcpp::export]]
Rcpp::NumericVector steck_double(Rcpp::NumericVector v) {
	return toNumericVector(steck_lower<double>(fromNumericVector<double>(v).data(),v.length()));
}

// [[Rcpp::export]]
Rcpp::NumericVector bolshev_double(Rcpp::NumericVector v) {
	return toNumericVector(bolshev_lower<double>(fromNumericVector<double>(v).data(),v.length()));
}

// [[Rcpp::export]]
Rcpp::NumericVector bolshev_exact(Rcpp::NumericVector v) {
	return toNumericVector(bolshev_lower<mpq_class>(fromNumericVector<mpq_class>(v).data(),v.length()));
}

// [[Rcpp::export]]
Rcpp::NumericMatrix bolshev2_exact(Rcpp::NumericVector v1,Rcpp::NumericVector v2,const int n1,const int n2) {
	return toNumericMatrix<mpq_class>(bolshev2_lower<mpq_class>(fromNumericVector<mpq_class>(v1).data(),fromNumericVector<mpq_class>(v2).data(),n1,n2,std::min(v1.length(),v2.length())));
}

// [[Rcpp::export]]
Rcpp::NumericMatrix steck2_exact(Rcpp::NumericVector v1,Rcpp::NumericVector v2,const int n1,const int n2) {
	return toNumericMatrix<mpq_class>(steck2_lower<mpq_class>(fromNumericVector<mpq_class>(v1).data(),fromNumericVector<mpq_class>(v2).data(),n1,n2,std::min(v1.length(),v2.length())));
}

// [[Rcpp::export]]
Rcpp::NumericMatrix steck2_exact_p(Rcpp::NumericVector v1,Rcpp::NumericVector v2,const int n1,const int n2) {
	return toNumericMatrix<mpq_class>(steck2_lower_p<mpq_class>(fromNumericVector<mpq_class>(v1).data(),fromNumericVector<mpq_class>(v2).data(),n1,n2,std::min(v1.length(),v2.length())));
}

// [[Rcpp::export]]
Rcpp::NumericMatrix noe2_faithful_p(Rcpp::NumericVector v1,Rcpp::NumericVector v2,const int n1,const int n2) {
	return toNumericMatrix<PairArithmetic::DoublePair>(noe2_lower_p<PairArithmetic::DoublePair>(fromNumericVector<PairArithmetic::DoublePair>(v1).data(),fromNumericVector<PairArithmetic::DoublePair>(v2).data(),n1,n2,std::min(v1.length(),v2.length())));
}

// [[Rcpp::export]]
Rcpp::NumericMatrix noe2_p(Rcpp::NumericVector v1,Rcpp::NumericVector v2,const int n1,const int n2) {
  return toNumericMatrix<double>(noe2_lower_p<double>(fromNumericVector<double>(v1).data(),fromNumericVector<double>(v2).data(),n1,n2,std::min(v1.length(),v2.length())));
}

// [[Rcpp::export]]
Rcpp::NumericMatrix bolshev2_exact_p(Rcpp::NumericVector v1,Rcpp::NumericVector v2,const int n1,const int n2) {
	return toNumericMatrix<mpq_class>(bolshev2_lower_p<mpq_class>(fromNumericVector<mpq_class>(v1).data(),fromNumericVector<mpq_class>(v2).data(),n1,n2,std::min(v1.length(),v2.length())));
}


// [[Rcpp::export]]
Rcpp::NumericMatrix steck2_double_p(Rcpp::NumericVector v1,Rcpp::NumericVector v2,const int n1,const int n2) {
	return toNumericMatrix<double>(steck2_lower_p<double>(fromNumericVector<double>(v1).data(),fromNumericVector<double>(v2).data(),n1,n2,std::min(v1.length(),v2.length())));
}

// [[Rcpp::export]]
Rcpp::NumericMatrix bolshev2_double_p(Rcpp::NumericVector v1,Rcpp::NumericVector v2,const int n1,const int n2) {
	return toNumericMatrix<double>(bolshev2_lower_p<double>(fromNumericVector<double>(v1).data(),fromNumericVector<double>(v2).data(),n1,n2,std::min(v1.length(),v2.length())));
}

// [[Rcpp::export]]
Rcpp::NumericMatrix bolshev2_double(Rcpp::NumericVector v1,Rcpp::NumericVector v2,const int n1,const int n2) {
	return toNumericMatrix<double>(bolshev2_lower<double>(fromNumericVector<double>(v1).data(),fromNumericVector<double>(v2).data(),n1,n2,std::min(v1.length(),v2.length())));
}

// [[Rcpp::export]]
long int bolshev2_opcount(const int n1,const int n2) {
  OpCounter::reset();
	std::vector<OpCounter> v(n1+n2);
	std::vector< std::vector<OpCounter> > res = bolshev2_lower<OpCounter>(v.data(),v.data(),n1,n2,n1+n2);
	return OpCounter::global_k;
}

// [[Rcpp::export]]
long int noe2_opcount(const int n1,const int n2) {
  OpCounter::reset();
	std::vector<OpCounter> v(n1+n2);
	std::vector< std::vector<OpCounter> > res = noe2_lower<OpCounter>(v.data(),v.data(),n1,n2,n1+n2);
	return OpCounter::global_k;
}

// [[Rcpp::export]]
long int steck2_opcount(const int n1,const int n2) {
  OpCounter::reset();
	std::vector<OpCounter> v(n1+n2);
	std::vector< std::vector<OpCounter> > res = steck2_lower<OpCounter>(v.data(),v.data(),n1,n2,n1+n2);
	return OpCounter::global_k;
}

// [[Rcpp::export]]
long int steck2_slow_opcount(const int n1,const int n2) {
  OpCounter::reset();
	std::vector<OpCounter> v(n1+n2);
	std::vector< std::vector<OpCounter> > res = steck2_lower_slow<OpCounter>(v.data(),v.data(),n1,n2,n1+n2);
	return OpCounter::global_k;
}

// Code for examples

template<typename T>
std::vector< std::vector<T> > fd_fm_(const std::vector< std::vector<T> >& order_stat,const Rcpp::NumericVector& v1,const Rcpp::NumericVector& v2,const int m0,const int m,const int k_offset=0) {
	std::vector< std::vector<T> > res(m0+1,std::vector<T>(m+1,-1));
	// *Rf_choose(m0,j)
	T c1 = 1;
	for(int j = 0; j <= m0; j++) {
		// *Rf_choose(m-m0,k-j)
		T c2 = 1;
		for(int k = j; k <= m; k++) {
			if(j>=k-m+m0) {
				if(k>=j+k_offset) {
					T val = c1 * c2 * order_stat[m0-j][m-k-(m0-j)]; 
					if(k>0) val *= PairArithmetic::pow(T(1-v1[m-k]),j) * PairArithmetic::pow(T(1-v2[m-k]),k-j);
					res[j][k] = val;
				}
				c2 *= m-m0-(k-j);
				c2 /= k-j+1;
			} else res[j][k] = 0;
		}
		c1 *= m0-j;
		c1 /= j+1;
	}
	return res;
}

template<typename T>
std::vector< std::vector<T> > fd_fm_partial_(const std::vector< std::vector<T> >& order_stat,const Rcpp::NumericVector& v1,const Rcpp::NumericVector& v2,const int m0,const int m,const int k_offset=0) {
	std::vector< std::vector<T> > res(m0+1,std::vector<T>(m+1,-1));
	// *Rf_choose(m0,j)
	T c1 = 1;
	//T mx = -1;
	for(int j = 0; j <= m0; j++) {
		// *Rf_choose(m-m0,k-j)
		T c2 = 1;
		for(int k = j; k <= m; k++) {
			if(j>=k-m+m0) {
				if(k>=j+k_offset) {
					T val = c1 * c2;
					bool flag = true;
					if(m0-j<order_stat.size() && m-k-(m0-j)<order_stat[0].size()) val *= order_stat[m0-j][m-k-(m0-j)];
					else flag = false;
					if(k>0) val *= PairArithmetic::pow(T(1-v1[m-k]),j) * PairArithmetic::pow(T(1-v2[m-k]),k-j);
					if(flag) res[j][k] = val;
					else {
						res[j][k] = -1;
            //if((double) val > (double) mx) mx = val;
					}
				}
				c2 *= m-m0-(k-j);
				c2 /= k-j+1;
			} else res[j][k] = 0;
		}
		c1 *= m0-j;
		c1 /= j+1;
	}
	//Rcpp::Rcout << "test " << (double) mx << std::endl;
	return res;
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector fd_fm(Rcpp::NumericVector v1,Rcpp::NumericVector v2,const int m0,const int m) {
	if(v1.length() < m || v2.length() < m) ::Rf_error("v1 or v2 too short!");
	std::vector< std::vector<PairArithmetic::DoublePair> > order_stat = noe2_lower_p(std::vector<PairArithmetic::DoublePair>(v1.begin(),v1.end()).data(), std::vector<PairArithmetic::DoublePair>(v2.begin(),v2.end()).data(),m0,m-m0, std::min(v1.length(),v2.length()));
	return toNumericMatrix(fd_fm_(order_stat,v1,v2,m0,m));
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector fd_fm_double(Rcpp::NumericVector v1,Rcpp::NumericVector v2,const int m0,const int m) {
	if(v1.length() < m || v2.length() < m) ::Rf_error("v1 or v2 too short!");
	std::vector< std::vector<double> > order_stat = noe2_lower_p(std::vector<double>(v1.begin(),v1.end()).data(), std::vector<double>(v2.begin(),v2.end()).data(),m0,m-m0, std::min(v1.length(),v2.length()));
	return toNumericMatrix(fd_fm_(order_stat,v1,v2,m0,m));
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector fd_fm_b(Rcpp::NumericVector v1,Rcpp::NumericVector v2,const int m0,const int m) {
	if(v1.length() < m || v2.length() < m) ::Rf_error("v1 or v2 too short!");
	std::vector< std::vector<mpf_class> > order_stat = steck2_lower_b(std::vector<mpf_class>(v1.begin(),v1.end()).data(), std::vector<mpf_class>(v2.begin(),v2.end()).data(),m0,m-m0, std::min(v1.length(),v2.length()));
	return toNumericMatrix(fd_fm_(order_stat,v1,v2,m0,m));
}

// [[Rcpp:export(rng = false)]]
int mpf_default_prec() {
	return mpf_get_default_prec();
}

template<typename T>
std::vector< std::vector<T> > toMatrix(const Rcpp::NumericMatrix& m) {
  std::vector< std::vector<T> > ret;
  for(int i=0; i < m.nrow(); i++){
    Rcpp::NumericMatrix::ConstRow tmp = m(i,Rcpp::_);
    ret.push_back(std::vector<T>(tmp.begin(),tmp.end()));
  } 
  return ret;
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix fd_fm_order_stat(Rcpp::NumericVector v1,Rcpp::NumericVector v2,const int m0,const int m, Rcpp::NumericMatrix order_stat) {
	if(v1.length() < m || v2.length() < m) ::Rf_error("v1 or v2 too short!");
	return toNumericMatrix(fd_fm_(toMatrix<PairArithmetic::DoublePair>(order_stat),v1,v2,m0,m));
}

// [[Rcpp::export(rng = false)]]
double calc_pwr(Rcpp::NumericMatrix mat, const int m, const int m0) {
	std::vector< std::vector<double> > fd_fm = toMatrix<double>(mat);
	double res = 0;
	for(int j = 0; j <= m0; j++) {
		for(int k = j; k <= (m-m0+j); k++) {
			double val = fd_fm[j][k];
			val *= (k-j);
			val /= (double)(m-m0);
			res += val;
		}
	}
	return res;
}


template<typename T>
Rcpp::NumericVector fd_fm_pwr_(Rcpp::NumericVector v1,Rcpp::NumericVector v2, std::vector< std::vector<T> > order_stat, const int m, int m0_max=-1) {
	if(m0_max<0) m0_max=m;
	if(v1.length() < m || v2.length() < m) ::Rf_error("v1 or v2 too short!");
	Rcpp::NumericVector res(m,0.0);
	tbb::parallel_for(0, m0_max, [&] (int m0) {
	//for(int m0 = 0; m0 < m0_max; m0++) {
		std::vector< std::vector<T> > tmp = fd_fm_partial_(order_stat,v1,v2,m0,m);
		for(int j = 0; j <= m0; j++) {
			for(int k = j; k <= (m-m0+j); k++) {
				T val = tmp[j][k];
				val *= (k-j);
				val /= (m-m0);
				res[m0] += toDouble(val);
			}
		}
	//}
	});
	return res;
}

// [[Rcpp::export(rng = false)]]
Rcpp::List fdp_dist(Rcpp::NumericMatrix order_stat, Rcpp::NumericVector v1,Rcpp::NumericVector v2,const int m, const int m0) {
	Rcpp::NumericMatrix dist = fd_fm_order_stat(v1,v2,m0,m,order_stat);
	std::map<double, double> vals;
	for(int k=1; k<=m+1; k++) {
		for(int j=1; j<=std::min(m0+1,k); j++) {
			double v = j-1;
			v /= std::max(k-1,1);
			auto res = vals.find(v);
			double p = dist(j-1,k-1);
			if (res != vals.end()) p+=res->second;
			vals[v] = p;
		}
	}
	Rcpp::NumericVector supp(vals.size(),0.0);
	Rcpp::NumericVector p(vals.size(),0.0);
	
	int i = 0;
	for (auto it = vals.begin(); it != vals.end(); ++it) {
		supp[i] = it->first;
		p[i] = it->second;
		++i;
	}
	
	return Rcpp::List::create(Rcpp::Named("supp") = supp , Rcpp::Named("p") = p);
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector fd_fm_pwr(Rcpp::NumericVector v1,Rcpp::NumericVector v2,const int m) {
	std::vector< std::vector<PairArithmetic::DoublePair> > order_stat = noe2_lower_p(std::vector<PairArithmetic::DoublePair>(v1.begin(),v1.end()).data(), std::vector<PairArithmetic::DoublePair>(v2.begin(),v2.end()).data(),m,m, std::min(v1.length(),v2.length()));
	return fd_fm_pwr_(v1,v2,order_stat,m);
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector fd_fm_pwr_double(Rcpp::NumericVector v1,Rcpp::NumericVector v2,const int m) {
	std::vector< std::vector<double> > order_stat = noe2_lower_p(std::vector<double>(v1.begin(),v1.end()).data(), std::vector<double>(v2.begin(),v2.end()).data(),m,m, std::min(v1.length(),v2.length()));
	return fd_fm_pwr_(v1,v2,order_stat,m);
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector fd_fm_pwr_lim(Rcpp::NumericVector v1,Rcpp::NumericVector v2,const int m,const int m0) {
	std::vector< std::vector<PairArithmetic::DoublePair> > order_stat = noe2_lower_p(std::vector<PairArithmetic::DoublePair>(v1.begin(),v1.end()).data(), std::vector<PairArithmetic::DoublePair>(v2.begin(),v2.end()).data(),m0,m, std::min(v1.length(),v2.length()));
	return fd_fm_pwr_(v1,v2,order_stat,m,m0);
}

template<typename T,typename T1>
std::vector< std::vector<T> > toMatrix(const std::vector< std::vector<T1> >& m) {
	std::vector< std::vector<T> > ret;
	for(auto &v : m) ret.push_back(std::vector<T>(v.begin(),v.end()));
	return ret;
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector fd_fm_pwr_bolshev(Rcpp::NumericVector v1,Rcpp::NumericVector v2,const int m) {
	std::vector< std::vector<double> > order_stat = bolshev2_lower(std::vector<double>(v1.begin(),v1.end()).data(), std::vector<double>(v2.begin(),v2.end()).data(),m,m, std::min(v1.length(),v2.length()));
	return fd_fm_pwr_(v1,v2,toMatrix<PairArithmetic::DoublePair,double>(order_stat),m);
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector fd_fm_pwr_bolshev_exact(Rcpp::NumericVector v1,Rcpp::NumericVector v2,const int m) {
	std::vector< std::vector<mpq_class> > order_stat = bolshev2_lower(std::vector<mpq_class>(v1.begin(),v1.end()).data(), std::vector<mpq_class>(v2.begin(),v2.end()).data(),m,m, std::min(v1.length(),v2.length()));
	return fd_fm_pwr_(v1,v2,order_stat,m);
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector fd_fm_pwr_order_stat(Rcpp::NumericVector v1,Rcpp::NumericVector v2,const int m, Rcpp::NumericMatrix order_stat) {
	return fd_fm_pwr_(v1,v2,toMatrix<PairArithmetic::DoublePair>(order_stat),m);
}

// [[Rcpp::export(rng = false)]]
double fd_fm_lambda_pwr(Rcpp::NumericVector v1,Rcpp::NumericVector v2,const int m,const double pr, const double lambda) {
	if(v1.length() < m || v2.length() < m) ::Rf_error("v1 or v2 too short!");
	PairArithmetic::DoublePair res = 0;
	std::vector< std::vector<PairArithmetic::DoublePair> > order_stat = noe2_lower_p(std::vector<PairArithmetic::DoublePair>(v1.begin(),v1.end()).data(), std::vector<PairArithmetic::DoublePair>(v2.begin(),v2.end()).data(),m,std::ceil((1-lambda)*m), std::min(v1.length(),v2.length()));
	for(int m0 = 0; m0 <= m; m0++) {
		const int k_offset = std::ceil(lambda*(m-m0));
		std::vector< std::vector<PairArithmetic::DoublePair> > tmp = fd_fm_(order_stat,v1,v2,m0,m,k_offset);
		for(int j = 0; j <= m0; j++) {
			for(int k = j+k_offset; k <= (m-m0+j); k++) {
				res += PairArithmetic::fast_pow(PairArithmetic::DoublePair(pr),m0) * PairArithmetic::fast_pow(PairArithmetic::DoublePair(1-pr),m-m0) * PairArithmetic::choose<PairArithmetic::DoublePair>(m,m0) * tmp[j][k];
			}
		}
	}
	return (double)res;
}

// [[Rcpp::export(rng = false)]]
double fd_fm_lambda_pwr_double(Rcpp::NumericVector v1,Rcpp::NumericVector v2,const int m,const double pr, const double lambda) {
	if(v1.length() < m || v2.length() < m) ::Rf_error("v1 or v2 too short!");
	PairArithmetic::DoublePair res = 0;
	std::vector< std::vector<double> > order_stat = noe2_lower_p(std::vector<double>(v1.begin(),v1.end()).data(), std::vector<double>(v2.begin(),v2.end()).data(),m,std::ceil((1-lambda)*m), std::min(v1.length(),v2.length()));
	for(int m0 = 0; m0 <= m; m0++) {
		const int k_offset = std::ceil(lambda*(m-m0));
		std::vector< std::vector<double> > tmp = fd_fm_(order_stat,v1,v2,m0,m,k_offset);
		for(int j = 0; j <= m0; j++) {
			for(int k = j+k_offset; k <= (m-m0+j); k++) {
				res += PairArithmetic::fast_pow(pr,m0) * PairArithmetic::fast_pow(1-pr,m-m0) * PairArithmetic::choose<double>(m,m0) * tmp[j][k];
			}
		}
	}
	return (double)res;
}

// [[Rcpp::export(rng = false)]]
double fd_fm_lambda_pwr_steck(Rcpp::NumericVector v1,Rcpp::NumericVector v2,const int m,const double pr, const double lambda) {
	if(v1.length() < m || v2.length() < m) ::Rf_error("v1 or v2 too short!");
	PairArithmetic::DoublePair res = 0;
	std::vector< std::vector<double> > order_stat = steck2_lower_p(std::vector<double>(v1.begin(),v1.end()).data(), std::vector<double>(v2.begin(),v2.end()).data(),m,std::ceil((1-lambda)*m), std::min(v1.length(),v2.length()));
	for(int m0 = 0; m0 <= m; m0++) {
		const int k_offset = std::ceil(lambda*(m-m0));
		std::vector< std::vector<double> > tmp = fd_fm_(order_stat,v1,v2,m0,m,k_offset);
		for(int j = 0; j <= m0; j++) {
			for(int k = j+k_offset; k <= (m-m0+j); k++) {
				res += PairArithmetic::fast_pow(pr,m0) * PairArithmetic::fast_pow(1-pr,m-m0) * PairArithmetic::choose<double>(m,m0) * tmp[j][k];
			}
		}
	}
	return (double)res;
}
