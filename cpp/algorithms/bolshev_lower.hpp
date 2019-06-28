#include<vector>

template<typename T>
std::vector<T> bolshev_lower(const T* v, const int l)
{
	//This line may look strange but it allows us to keep the precision if the input is of type T=mpf_class
	T one = v[0]*0+1;
	std::vector<T> Psi(l+1,one);
	std::vector<T> summands_vec(l,1-v[0]);
	T* summands = summands_vec.data();
	for(int k = 2; k <= l; k++) {
		T* Fk = &Psi[k-1];
		const T* curr_v = &v[0];
		T* curr_summand = &summands[0];
		Rcpp::checkUserInterrupt();
		for(int j = k; j > 1; --j) {
			*Fk -= *curr_summand;
			*curr_summand *= ((1-*curr_v) * k) / j;
			++curr_summand;
			++curr_v;
		}
		summands[k-1] = k * *Fk * (1-v[k-1]);
	}
	T* ret = &Psi[l];
	for(int j = 0; j < l; j++) *ret -= summands[j];
	return Psi;
}