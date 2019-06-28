#include<Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include<RcppParallel.h>
// [[Rcpp::plugins(cpp11)]]

#define NO_OP

#include "misc/helpers.hpp"
#include "PairArithmetic/PairArithmetic.hpp"

#define ALGORITHM_PARALLEL
template<typename T>
std::vector< std::vector<T> > noe2_lower_k(const T* v1, const T* v2, const int n1, const int n2, const int max_idx)
#include "algorithms/noe2_lower.hpp"

#define LOGN_SUMMATION
template<typename T>
std::vector< std::vector<T> > noe2_lower_k_large(const T* v1, const T* v2, const int n1, const int n2, const int max_idx)
#include "algorithms/noe2_lower.hpp"

// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix noe_faithful_k(Rcpp::NumericVector v1,Rcpp::NumericVector v2,const int n1,const int n2) {
	std::vector< std::vector<PairArithmetic::DoublePair> > res =
		noe2_lower_k<PairArithmetic::DoublePair>(std::vector<PairArithmetic::DoublePair>(v1.begin(),v1.end()).data(), std::vector<PairArithmetic::DoublePair>(v2.begin(),v2.end()).data(), n1, n2, std::min(v1.length(),v2.length()) );
	Rcpp::NumericMatrix m(n1+1, n2+1);
	for(int i=0;i < n1+1; i++) {
	  for(int j=0; j < n2+1; j++) {
		double d = (double)res[i][j];
	    if(d == -1) m(i,j) = NA_REAL;
		else m(i,j) = res[i][j].getK();
	  }
	}
	return m;
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix noe_faithful_k_large(Rcpp::NumericVector v1,Rcpp::NumericVector v2,const int n1,const int n2) {
	std::vector< std::vector<PairArithmetic::DoublePair> > res =
		noe2_lower_k_large<PairArithmetic::DoublePair>(std::vector<PairArithmetic::DoublePair>(v1.begin(),v1.end()).data(), std::vector<PairArithmetic::DoublePair>(v2.begin(),v2.end()).data(), n1, n2, std::min(v1.length(),v2.length()) );
	Rcpp::NumericMatrix m(n1+1, n2+1);
	for(int i=0;i < n1+1; i++) {
	  for(int j=0; j < n2+1; j++) {
		double d = (double)res[i][j];
	    if(d == -1) m(i,j) = NA_REAL;
		else m(i,j) = res[i][j].getK();
	  }
	}
	return m;
}