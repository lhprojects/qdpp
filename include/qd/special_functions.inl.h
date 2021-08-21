

inline qd_real factorial(int n) {
	if(n<0) {
		//qd_real::error("factorial: n<0");
		return qd_real::_nan;
	}
	qd_real f = 1.0;
	for(int i = 1; i <= n; ++i) {
		f *= i;
	}
	return f;
}

inline qd_real factorial(int n, int m) {
	if(n<0 || m>n) {
		//qd_real::error("factorial: n<0 || m>n");
		return qd_real::_nan;
	}
	qd_real f = 1.0;
	for(int i = m; i <= n; ++i) {
		f *= i;
	}
	return f;
}

inline qd_real combinal(int n, int k) {
	if(n<0 || k<0 || k>n) {
		//qd_real::error("combinal: n<0 || k<0 || k>n");
		return qd_real::_nan;
	}
	if(2*k > n) { k = n - k; }
	if(k==0) return qd_real(1.0);
	return factorial(n, n-k+1)/factorial(k);
}

inline qd_real tgamma(qd_real a) {
	qd_real p = 1.0;
	for(qd_real i = a-1; i != 0; i -= 1) {
		p *= i;
	}
	return p;
}

inline qd_real tgamma_lower(qd_real a, qd_real x) {
	if(a - 1 < 0) {
		//qd_real::error("gamma(int,qd_real): a-1 < 0");
		return qd_real::_nan;
	}
	qd_real sum = 0.0;
	qd_real k = 1.0;
	qd_real p = 1.0;
	for(int n = 0; n < 100000000; ++n) {
		qd_real oldsum = sum;
		k /= (a+n);
		sum += k*p;
		p *= x;
		if(oldsum==sum) {
			break;
		}
	}
	return exp(-x+a*log(x))*sum;
	//return exp(-x)*pow(x,a)*sum;
}

inline qd_real tgamma(qd_real a, qd_real x) {
	if(x==0) return tgamma(a);
	if(x < a+1) {
		return tgamma(a)-tgamma_lower(a, x);
	} else {
		static qd_real FPMIN = std::numeric_limits<qd_real>::min()
			/std::numeric_limits<qd_real>::epsilon();
		qd_real b = x+1.0-a;
		qd_real c = 1.0/FPMIN;
		qd_real d = 1.0/b;
		qd_real h = d;
		qd_real an;
		int count = 0;
		for(int i=1; count<2; ++i) {
			an = -i*(i-a);
			b += 2;
			d = an*d+b;
			if(fabs(d) < FPMIN) d = FPMIN;
			c=b+an/c;
			if(fabs(c) < FPMIN) c = FPMIN;
			d = 1.0/d;
			qd_real del = d*c;
			h *= del;
			if(fabs(del-1) <= 3*std::numeric_limits<qd_real>::epsilon()) ++count;
		}
		return exp(-x+a*log(x))*h;
	}
}



