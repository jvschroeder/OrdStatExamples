#include<vector>

template<typename T>
std::vector<T> steck_lower(const T* v, const int l)
{
	//This line may look strange but it allows us to keep the precision if the input is of type T=mpf_class
	T one = v[0]*0+1;
	std::vector<T> summands(l+1,one);
	std::vector<T> Psi(l+1,one);
	T* summands_ptr = summands.data();;
	for(int i = 1; i <= l; ++i) {
		T* Fi = &Psi[i];
		*Fi = PairArithmetic::fast_pow<T>(v[i-1],i);
		for(int j=0; j <= i-1; ++j) {
			*Fi -= summands_ptr[j] * PairArithmetic::fast_pow<T>(v[i-1]-v[j],i-j);
			summands_ptr[j] *= i+1;
			summands_ptr[j]	/= i+1-j;
		}
		summands_ptr[i] = *Fi;
		summands_ptr[i] *= (i+1);
	}
	return Psi;
}