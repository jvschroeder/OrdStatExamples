#define CATCH_CONFIG_MAIN
#include<catch.hpp>

#include<limits>
#include "./PairArithmetic.hpp"

TEMPLATE_TEST_CASE( "No modification by arithmetic operations", "[no-mod]", float, double ) {
	typedef TestType T;
	const T zeroT = 0;
	const T oneT = 1;
	const T eps2 = std::numeric_limits<T>::epsilon()/((T)2);
	PairArithmetic::Pair<T> zero(0);
	PairArithmetic::Pair<T> oneEps2(1);
	oneEps2 += eps2;
	
	double val = (T) (zero + oneEps2);
	REQUIRE( oneEps2.getFl1() == oneT );
	REQUIRE( oneEps2.getFl2() == eps2 );
	REQUIRE( zero.getFl1() == zeroT );
	REQUIRE( zero.getFl2() == zeroT );
	REQUIRE( val == oneT);
	
	val = (T) (oneEps2 * oneEps2 / oneEps2 * oneEps2);
	REQUIRE( oneEps2.getFl1() == oneT );
	REQUIRE( oneEps2.getFl2() == eps2 );
	REQUIRE( val == oneT + std::numeric_limits<T>::epsilon());
	
	val = (T) (oneEps2 / oneEps2 * oneEps2 / oneEps2);
	REQUIRE( oneEps2.getFl1() == oneT );
	REQUIRE( oneEps2.getFl2() == eps2 );
	REQUIRE( val == oneT);
	
	val = (T) (oneEps2 - oneEps2 + oneEps2);
	REQUIRE( oneEps2.getFl1() == oneT );
	REQUIRE( oneEps2.getFl2() == eps2 );
	REQUIRE( val == oneT);
	
	PairArithmetic::Pair<T> val1 = oneEps2;
	val1 += oneEps2;
	val1 -= oneEps2;
	val1 *= oneEps2;
	val1 /= oneEps2;
	REQUIRE( oneEps2.getFl1() == oneT );
	REQUIRE( oneEps2.getFl2() == eps2 );
	REQUIRE( (T)val1 == oneT);
}

TEMPLATE_TEST_CASE( "Calculations involving 1+eps/2 are faithful", "[faithful]", float, double ) {
	typedef TestType T;
	const T eps = std::numeric_limits<T>::epsilon();
	const T eps2 = eps/((T)2);
	const PairArithmetic::Pair<T> one(1);
	
	REQUIRE( ((T)1) + eps2 == ((T)1)); //make sure we are actually testing something
	
	PairArithmetic::Pair<T> val = one;
	val += eps2;
	REQUIRE((T)(val+val) == (((T)2)+eps));
	REQUIRE((T)(val*val) == (((T)1)+eps));
	REQUIRE((val*val).getFl1() == ((T)1));
	REQUIRE((val*val).getFl2() == (eps + eps2 * eps2));
	REQUIRE((T)(one/val)<1);
	REQUIRE((T)((one/val) * val) == ((T)1));
	
	val -= eps;
	REQUIRE(val.getFl2() != ((T)0));
	REQUIRE((T)val != ((T)1));
	REQUIRE((T)(val+val) == (((T)2)-eps));
	REQUIRE((T)(val*val) == (((T)1)-eps));
	REQUIRE((T)(one/val)>((T)1));
	REQUIRE((T)((one/val) * val) == ((T)1));
	REQUIRE((T)(val+eps2-((T)1))==((T)0));
	
	REQUIRE((T) (((val/eps2)/val)*eps2) == ((T)1));
	
	val = one - eps2;
	REQUIRE((T) (((val/eps2)/val)*eps2) == ((T)1));
}