
#include <stdio.h>
#include <stdint.h>
#include <qd/qd_real.h>
#include <iostream>
#include <cmath>
int main() {
	dd_real dd = 1;
	qd_real qd = 1;
	qd_real sqd =sin(qd);
	dd_real sdd =sin(dd);
	std::cout << sqd.to_string(qd_real::_ndigits) << std::endl;
	std::cout << sdd.to_string(dd_real::_ndigits) << std::endl;
	printf("0.84147098480789650665250232163029899962256306079837106567275170999\n");


}

