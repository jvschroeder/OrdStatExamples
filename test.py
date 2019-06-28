import crossprob
import math
res = []
n = 2
def difference(a,b):
	if a == b:
		return -17
	else:
		return math.log10(abs(a-b) / b)

while True:
	v = [0.05 * i / n for i in range(1,n+1)]
	probs = (crossprob.ecdf_lower_noncrossing_probability(n,v),
		crossprob.ecdf_lower_noncrossing_probability_n2(n,v),
		crossprob.ecdf_lower_noncrossing_probability_n2logn(n,v))
	print((n,probs))
	v = 0.05 / n
	v = v*math.pow((n+1)*v,(n-1))
	if v == 0:
		break
	res.append((n,tuple(difference(probs[i],v) for i in range(0,3))))
	n += 1

import json
json.dump(res,open('./test.json','w'))