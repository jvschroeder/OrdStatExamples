#ifndef HELPERS_HPP
#define HELPERS_HPP

template<typename T>
T binomial(int n, int k) {
	T ret = 1;
	for(int i=1; i <= k; i++) {
		ret *= n+1-i;
		ret /= i;
	}
	return ret;
}

template<typename T>
std::vector<T> rev_vec(const T* v_, const int l) {
	std::vector<T> v(v_,v_+l);
	std::transform(v.begin(),v.end(),v.begin(),[](const T& val) -> T {return 1-val;});
	std::reverse(v.begin(),v.end());
	return v;
}

template<typename T>
std::vector<T> rev_vec1(const T* v_, const int l) {
	std::vector<T> v(v_,v_+l);
	std::transform(v.begin(),v.end(),v.begin(),[](const T& val) -> T {return 1-val;});
	return v;
}

template<typename T>
double toDouble(const T& v) {
	return (double) v;
}

#ifdef Rcpp_hpp
template<typename T>
Rcpp::NumericMatrix toNumericMatrix(const std::vector< std::vector<T> >& m) {
	int n1 = m.size(), n2=m[0].size();
	Rcpp::NumericMatrix ret(n1,n2);
	for(int i=0;i < n1; i++) {
		for(int j=0; j < n2; j++) {
			double d = toDouble(m[i][j]);
			if(d == -1) d = NA_REAL;
			ret(i,j) = d;
		}
	}
	return ret;
}
template<typename T>
std::vector<T> fromNumericVector(const Rcpp::NumericVector& v) {
	return std::vector<T>(v.begin(),v.end());
}
template<typename T>
Rcpp::NumericVector toNumericVector(const std::vector<T>& v) {
	return Rcpp::NumericVector(v.begin(),v.end());
}
#ifdef __GMP_PLUSPLUS__
template<>
Rcpp::NumericVector toNumericVector<mpq_class>(const std::vector<mpq_class>& v) {
	Rcpp::NumericVector ret(v.size());
	std::transform(v.begin(),v.end(),ret.begin(),[](const mpq_class& val) -> double { return val.get_d(); });
	return ret;
}
template<>
Rcpp::NumericMatrix toNumericMatrix<mpq_class>(const std::vector< std::vector<mpq_class> >& m) {
	int n1 = m.size(), n2=m[0].size();
	Rcpp::NumericMatrix ret(n1,n2);
	for(int i=0;i < n1; i++) {
		for(int j=0; j < n2; j++) {
			double d = m[i][j].get_d();
			if(d == -1) d = NA_REAL;
			ret(i,j) = d;
		}
	}
	return ret;
}
template<>
double toDouble(const mpq_class& q) {
	return q.get_d();
}
template<>
double toDouble(const mpf_class& f) {
  return f.get_d();
}
#endif
#endif
#endif