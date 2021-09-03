
#include <stdio.h>
#include <stdint.h>
#include <qd/dd.h>
#include <qd/qd_real.h>
#include <qd/dd_math.h>
#include <qd/special_functions.h>
#include <iostream>
#include <cmath>
#include <sstream>


int nerr = 0;
#define QdAssert(x) do {\
	bool v = bool(x);\
	if (!v) {\
		++nerr;\
		char const *expr = #x;\
		printf("%2d: failed: %s\n", __LINE__, expr);\
	}\
} while(false)

bool qd_close(double a, double b)
{
	if (std::isnan(a) || std::isnan(b)) {
		return std::isnan(a) == std::isnan(b);
	} else if (std::isinf(a) || std::isinf(b)) {
		return a == b;
	} else {
		if (std::abs(a - b) < 1E-8) {
			return true;
		} else if (std::abs(a - b) < a * 1E-8) {

		}
	}
	return false;
}



void test_cmp()
{
	QdAssert(abs(dd_real(0.)) == dd_real(0.));
	QdAssert(abs(dd_real(-0.)) == dd_real(0.));
	QdAssert(abs(dd_real(1.)) == dd_real(1.));
	QdAssert(isnan(abs(dd_real::_nan)));
	QdAssert(abs(dd_real::_inf) == dd_real::_inf);
	QdAssert(abs(-dd_real::_inf) == dd_real::_inf);

	QdAssert(dd_real(0.).is_zero());
	QdAssert(dd_real(-0.).is_zero());
	QdAssert(!dd_real(1.).is_zero());
	QdAssert(!dd_real::_nan.is_zero());
	QdAssert(!(dd_real::_inf).is_zero());
	QdAssert(!((-dd_real::_inf).is_zero()));

	QdAssert(!dd_real(0.).is_one());
	QdAssert(dd_real(1.).is_one());
	QdAssert(!dd_real::_nan.is_one());
	QdAssert(!(dd_real::_inf).is_one());
	QdAssert(!((-dd_real::_inf).is_one()));

	QdAssert(!dd_real(0.).is_positive());
	QdAssert(!dd_real(-0.).is_positive());
	QdAssert(dd_real(1.).is_positive());
	QdAssert(!dd_real::_nan.is_positive());
	QdAssert((dd_real::_inf).is_positive());
	QdAssert(!((-dd_real::_inf).is_positive()));

	QdAssert(!dd_real(0.).is_negative());
	QdAssert(!dd_real(-0.).is_negative());
	QdAssert(!dd_real(1.).is_negative());
	QdAssert(!dd_real::_nan.is_negative());
	QdAssert(!(dd_real::_inf).is_negative());
	QdAssert(((-dd_real::_inf).is_negative()));

	QdAssert(dd_real(0.).isfinite());
	QdAssert(dd_real(1.).isfinite());
	QdAssert(!dd_real::_nan.isfinite());
	QdAssert(!(dd_real::_inf).isfinite());
	QdAssert(!(-dd_real::_inf).isfinite());

	QdAssert(!dd_real(0.).isinf());
	QdAssert(!dd_real(1.).isinf());
	QdAssert(!dd_real::_nan.isinf());
	QdAssert((dd_real::_inf).isinf());
	QdAssert((-dd_real::_inf).isinf());

	QdAssert(!dd_real(0.).isnan());
	QdAssert(!dd_real(1.).isnan());
	QdAssert(dd_real::_nan.isnan());
	QdAssert(!(dd_real::_inf).isnan());
	QdAssert(!(-dd_real::_inf).isnan());

	QdAssert(dd_real(0.) == dd_real(0.));
	QdAssert(!(dd_real(0.) == dd_real(1.)));
	QdAssert(!(dd_real(0.) == dd_real::_nan));
	QdAssert(!(dd_real::_nan == dd_real::_nan));

	QdAssert(!(dd_real(0.) != dd_real(0.)));
	QdAssert((dd_real(0.) != dd_real(1.)));
	QdAssert((dd_real(0.) != dd_real::_nan));
	QdAssert((dd_real::_nan != dd_real::_nan));

	QdAssert(!(dd_real(0.) < 0.));
	QdAssert(dd_real(0.) < 1.);
	QdAssert(!(dd_real(1.) < dd_real(0.)));
	QdAssert(!(dd_real(0.) < dd_real::_nan));
	QdAssert(!(dd_real::_nan < dd_real::_nan));

	QdAssert((dd_real(0.) <= 0.));
	QdAssert(dd_real(0.) <= 1.);
	QdAssert(!(dd_real(1.) <= dd_real(0.)));
	QdAssert(!(dd_real(0.) <= dd_real::_nan));
	QdAssert(!(dd_real::_nan <= dd_real::_nan));

	QdAssert(!(dd_real(0.) > dd_real(0.)));
	QdAssert(!(dd_real(0.) > dd_real(1.)));
	QdAssert((dd_real(1.) > dd_real(0.)));
	QdAssert(!(dd_real(0.) > dd_real::_nan));
	QdAssert(!(dd_real::_nan > dd_real::_nan));

	QdAssert((dd_real(0.) >= dd_real(0.)));
	QdAssert(!(dd_real(0.) >= dd_real(1.)));
	QdAssert((dd_real(1.) >= dd_real(0.)));
	QdAssert(!(dd_real(0.) >= dd_real::_nan));
	QdAssert(!(dd_real::_nan >= dd_real::_nan));

}

void test_1() {
	printf("test 1\n");
char const *sin_table[] = {"0",
	"0.433883739117558120475768332848358754609990727787459876444547303532",
	"0.781831482468029808708444526674057750232334518708687528980634958045",
	"0.974927912181823607018131682993931217232785800619997437648079575088",
	"0.974927912181823607018131682993931217232785800619997437648079575088",
	"0.781831482468029808708444526674057750232334518708687528980634958045",
	"0.433883739117558120475768332848358754609990727787459876444547303532",
	"0","-0.433883739117558120475768332848358754609990727787459876444547303532",
	"-0.781831482468029808708444526674057750232334518708687528980634958045",
	"-0.974927912181823607018131682993931217232785800619997437648079575088",
	"-0.974927912181823607018131682993931217232785800619997437648079575088",
	"-0.781831482468029808708444526674057750232334518708687528980634958045",
	"-0.433883739117558120475768332848358754609990727787459876444547303532",
	"0","0.433883739117558120475768332848358754609990727787459876444547303532",
	"0.781831482468029808708444526674057750232334518708687528980634958045",
	"0.974927912181823607018131682993931217232785800619997437648079575088",
	"0.974927912181823607018131682993931217232785800619997437648079575088",
	"0.781831482468029808708444526674057750232334518708687528980634958045",
	"0.433883739117558120475768332848358754609990727787459876444547303532",
	"0","-0.433883739117558120475768332848358754609990727787459876444547303532",
	"-0.781831482468029808708444526674057750232334518708687528980634958045",
	"-0.974927912181823607018131682993931217232785800619997437648079575088",
	"-0.974927912181823607018131682993931217232785800619997437648079575088",
	"-0.781831482468029808708444526674057750232334518708687528980634958045",
	"-0.433883739117558120475768332848358754609990727787459876444547303532",
	"0","0.433883739117558120475768332848358754609990727787459876444547303532",
	"0.781831482468029808708444526674057750232334518708687528980634958045",
	"0.974927912181823607018131682993931217232785800619997437648079575088",
	"0.974927912181823607018131682993931217232785800619997437648079575088",
	"0.781831482468029808708444526674057750232334518708687528980634958045",
	"0.433883739117558120475768332848358754609990727787459876444547303532",
	"0","-0.433883739117558120475768332848358754609990727787459876444547303532",
	"-0.781831482468029808708444526674057750232334518708687528980634958045",
	"-0.974927912181823607018131682993931217232785800619997437648079575088",
	"-0.974927912181823607018131682993931217232785800619997437648079575088",
	"-0.781831482468029808708444526674057750232334518708687528980634958045",
	"-0.433883739117558120475768332848358754609990727787459876444547303532",
	"0","0.433883739117558120475768332848358754609990727787459876444547303532",
	"0.781831482468029808708444526674057750232334518708687528980634958045",
	"0.974927912181823607018131682993931217232785800619997437648079575088",
	"0.974927912181823607018131682993931217232785800619997437648079575088",
	"0.781831482468029808708444526674057750232334518708687528980634958045",
	"0.433883739117558120475768332848358754609990727787459876444547303532",
	"0","-0.433883739117558120475768332848358754609990727787459876444547303532",
	"-0.781831482468029808708444526674057750232334518708687528980634958045",
	"-0.974927912181823607018131682993931217232785800619997437648079575088",
	"-0.974927912181823607018131682993931217232785800619997437648079575088",
	"-0.781831482468029808708444526674057750232334518708687528980634958045",
	"-0.433883739117558120475768332848358754609990727787459876444547303532",
	"0","0.433883739117558120475768332848358754609990727787459876444547303532",
	"0.781831482468029808708444526674057750232334518708687528980634958045",
	"0.974927912181823607018131682993931217232785800619997437648079575088",
	"0.974927912181823607018131682993931217232785800619997437648079575088",
	"0.781831482468029808708444526674057750232334518708687528980634958045",
	"0.433883739117558120475768332848358754609990727787459876444547303532",
	"0","-0.433883739117558120475768332848358754609990727787459876444547303532",
	"-0.781831482468029808708444526674057750232334518708687528980634958045",
	"-0.974927912181823607018131682993931217232785800619997437648079575088",
	"-0.974927912181823607018131682993931217232785800619997437648079575088",
	"-0.781831482468029808708444526674057750232334518708687528980634958045",
	"-0.433883739117558120475768332848358754609990727787459876444547303532",
	"0","0.433883739117558120475768332848358754609990727787459876444547303532",
	"0.781831482468029808708444526674057750232334518708687528980634958045",
	"0.974927912181823607018131682993931217232785800619997437648079575088",
	"0.974927912181823607018131682993931217232785800619997437648079575088",
	"0.781831482468029808708444526674057750232334518708687528980634958045"};
	for(int i = 0; i <= 75; ++i) {

		qd_real x = i*qd_real::_pi/7;
		qd_real sin_x =sin(x);
		qd_real std;
		qd_real::read(sin_table[i], std);
		if(abs(std-sin_x) >= qd_real::_eps*4) {
			printf("%s %s %s\n", x.to_string().c_str(), std.to_string().c_str(), sin_x.to_string().c_str());
			printf("Check Faild!");
			exit(1);
		}
	}
}

void test_2() {
	printf("test 2\n");
	qd_real x;
	qd_real v;
	int Gamma_table[] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880};
	for(int i = 1; i <= sizeof(Gamma_table)/sizeof(int); ++i) {
		if(tgamma(i*qd_real(1.0))!=Gamma_table[i-1]) {
			printf("Check Faild!");
			exit(1);
		}
	}
}

void test_3() {
	printf("test 3\n");
	char const * Gamma_table[] = {
		"0.909795989568950135405699302486770680162877203230780433524338",
		"5.98949026466225505808752040803790697773894158793597118736856",
		"119.998300207521318901114214250911357835649491117793144263387",
		"5039.99968652758046807612556879628380805666635172881446378363",
		"362879.999937959493974986995000580802993417116770655450487250",
		"3.99167999999871679570681830940518713882588664527623911142446E7",
		"6.22702079999999726475627927537195797507482382486371277330191E9",
		"1.30767436799999999940406923045860835853602503711402057287992E12",
		"3.55687428095999999999867990123387178372056559459418776719846E14",
		"1.21645100408831999999999970373756323449594486565301305938580E17",
		"5.10909421717094399999999999932809698876963067879323032610300E19",
		"2.58520167388849766399999999999984629478244866489625464412093E^22",
		"1.55112100433309859839999999999999996458298896379207603586115E^25",
		"1.08888694504183521607679999999999999999178885840099319961921E28",
		"8.84176199373970195454361599999999999999998086229530652606821E30",
		"8.22283865417792281772556287999999999999999999551904605326863E33",
		"8.68331761881188649551819440127999999999999999999894658279108E36",
		"1.03331479663861449296666513375231999999999999999999997514710E^40",
		"1.37637530912263450463159795815809023999999999999999999999412E^43",
		"2.03978820811974433586402817399028973568000000000000000000000E^46"
	};
	qd_real v;
	for(int i = 1; i <= 10; ++i) {
		qd_real::read(Gamma_table[i-1],v);
		qd_real v1 = tgamma(2*i, 0.5);
		qd_real err = fabs(v1/v-1)/std::numeric_limits<qd_real>::epsilon();
		if(err>10000) {
			printf("Check Faild!");
			exit(1);
		}
	}
}

void test_4() {
	printf("test 4\n");
	char const * Gamma_table[] = {
		"0.0404276819945128025798162905388905454930975101641305125818332",
		"1.59015549178417023480610742786302812272850206645580016155211",
		"73.9152785799675740500974512019382140098799477500852871704517",
		"4367.80676268716319072335202982173453507425499813266057934126",
		"351330.234564724397916809735961257831636242527983899685159554", 
		"3.96991300207267880046984535218553027688694257387531974660362E7", 
		"6.22267440188170373530997263408146620569345769319774283247154E9", 
		"1.30758412769094458459778908328688676357756576218131055114597E12", 
		"3.55685501572571064015622738482602118353601227595085545210693E14", 
		"1.21645058415291143993479997683647258622218074572476594879302E17", 
		"5.10909412415204761466403100158688597544330994695119310981074E19", 
		"2.58520167180158508136063467792648702585574197048172768615716E22", 
		"1.55112100428578986100407638687527228843169769489594082524589E25", 
		"1.08888694504075351741653350572295794842353851548646790501732E28", 
		"8.84176199373945283634191491016293081671563703438108496025333E30", 
		"8.22283865417791704507526605937279804330921461432052280757564E33", 
		"8.68331761881188636104152042906643053275051642615938342104468E36", 
		"1.03331479663861449265194629808619143286840615416937955167338E40", 
		"1.37637530912263450462420260993527325131683323321434440209888E43", 
		"2.03978820811974433586385377161390621771368322922953468779619E46"};
	qd_real v;
	for(int i = 1; i <= 10; ++i) {
		qd_real::read(Gamma_table[i-1],v);
		qd_real v1 = tgamma(2*i, 5);
		qd_real err = fabs(v1/v-1)/std::numeric_limits<qd_real>::epsilon();
		if(err>10000) {
			printf("Check Faild!");
			exit(1);
		}
	}
}

void test_5() {
	printf("test 5\n");
	char const * Gamma_table[] = {
		"9.83662422461598069338844836428776413123944652127434084557883E-21",
		"2.56149552308696065091401230091686432002624693079286296144694E-17",
		"6.68130751283770659214039617782994301211338717721988601544016E-14",
		"1.74585432489678264318681060991625633504011780614627642628321E-10",
		"4.57086717662191011585723263637856974949855748292477845515219E-7",
		"0.00119924185071822172186093260400160730493904264525207527360010",
		"3.15366474682060753807545679094109136968349009180743320036135",
		"8314.17023877902288915595841842281335260753940914571631926568", 
		"2.19798136515028398224937874969566721185730366739427274642909E7", 
		"5.82845139894218535390090062826551152667833677541004050001044E10", 
		"1.55076838562628180105943394887534397381557298298187802227562E14", 
		"4.14159232994245602733179600565339066158663229573636237269006E17", 
		"1.11071391141013231951552139754531643905374734690878758827372E21", 
		"2.99274831961811428018063768990583550523600745438579055814308E24", 
		"8.10638258198600976818780091581044900936208832996158375213490E27", 
		"2.20888708222262352485656349795568744428822634265169965392557E31", 
		"6.05988292778165041639831917301346526131209623454135552260383E34", 
		"1.67540421407209350136685778318390972215703398028270642310843E38", 
		"4.67346778587781406905432600069484245512781114315976510513678E41", 
		"1.31709877119268661901602570405284590691920347189519906591034E45"};
	qd_real v;
	for(int i = 1; i <= 10; ++i) {
		qd_real::read(Gamma_table[i-1],v);
		qd_real v1 = tgamma(2*i, 50);
		qd_real err = fabs(v1/v-1)/std::numeric_limits<qd_real>::epsilon();
		if(err>10000) {
			printf("Check Faild!");
			exit(1);
		}
	}
}

void test_log()
{
	{
		dd_real a1 = 2.;
		dd_real a = log(a1);
		dd_real ref = dd_real::read("0.693147180559945309417232121458176568075500134360255");
		QdAssert(fabs(to_double(a / ref - 1.)) / dd_real::_eps < 4);
	}

	{
		dd_real a1 = 1. + dd_real(ldexp(0.9990234375, -1));
		dd_real a = log(a1);
		dd_real ref = dd_real::read("0.40513953428142396425401455828287536809888916586806444");
		QdAssert(fabs(to_double(a / ref - 1.)) / dd_real::_eps < 4);
	}
	{
		dd_real a1 = 1. + dd_real(ldexp(0.9990234375, -10));
		dd_real a = log(a1);
		dd_real ref = dd_real::read("0.00097513322869915911147780634863041169586551321001853317");
		QdAssert(fabs(to_double(a / ref - 1.)) / dd_real::_eps < 4);
	}
	{
		dd_real a1 = 1. + dd_real(ldexp(0.9990234375, -20));
		dd_real a = log(a1);
		dd_real ref = dd_real::read("0.00000095274253997231664801869224953591596729340018178035198");
		QdAssert(fabs(to_double(a / ref - 1.)) / dd_real::_eps < 4);
	}
}

void test_asinh()
{
	{
		dd_real a1 = 2;
		dd_real a = asinh(a1);
		dd_real ref = dd_real::read("1.4436354751788103424932767402731052694055530031569805");
		QdAssert(fabs(to_double(a / ref - 1.)) / dd_real::_eps < 4);
	}

	{
		dd_real a1 = -2;
		dd_real a = asinh(a1);
		dd_real ref = dd_real::read("-1.4436354751788103424932767402731052694055530031569805");
		QdAssert(fabs(to_double(a / ref - 1.)) / dd_real::_eps < 4);
	}
	{
		dd_real a1 = 1. / 1024;
		dd_real a = asinh(a1);
		dd_real ref = dd_real::read("0.00097656234477963751076391095890851381048161859275417364");
		QdAssert(fabs(to_double(a / ref - 1.)) / dd_real::_eps < 4);
	}
	{
		dd_real a1 = -1. / 1024;
		dd_real a = asinh(a1);
		dd_real ref = dd_real::read("-0.00097656234477963751076391095890851381048161859275417364");
		QdAssert(fabs(to_double(a / ref - 1.)) / dd_real::_eps < 4);
	}
}

#if 0
#include "qd/unused/DecimalFloat.h"
void foo()
{
	QdAssert(to_string(DecimalFloat(0)) == "0");
	QdAssert(to_string(DecimalFloat(1)) == "1");
	QdAssert(to_string(DecimalFloat(12)) == "12");
	QdAssert(to_string(DecimalFloat(123)) == "123");
	QdAssert(to_string(DecimalFloat(-1)) == "-1");
	QdAssert(to_string(DecimalFloat(-12)) == "-12");
	QdAssert(to_string(DecimalFloat(-123)) == "-123");
	QdAssert(to_string(DecimalFloat(INT_MAX)) == std::to_string(INT_MAX));
	QdAssert(to_string(DecimalFloat(INT_MIN)) == std::to_string(INT_MIN));

	{
		DecimalFloat a(1);
		a.mul10();
		QdAssert(to_string(a) == "10");
	}
	{
		DecimalFloat a(1);
		a.mul2();
		QdAssert(to_string(a) == "2");
	}
	{
		DecimalFloat a(5);
		a.mul2();
		QdAssert(to_string(a) == "10");
	}
	{
		DecimalFloat a(6);
		a.mul2();
		QdAssert(to_string(a) == "12");
	}
	{
		DecimalFloat a(55);
		a.mul2();
		QdAssert(to_string(a) == "110");
	}
	{
		DecimalFloat a(1);
		a.div2();
		QdAssert(to_string(a) == "0.5");
	}
	{
		DecimalFloat a(1);
		a.div2(); a.div2();
		QdAssert(to_string(a) == "0.25");
	}
	{
		DecimalFloat a(2);
		a.div2();
		QdAssert(to_string(a) == "1");
	}
	{
		DecimalFloat a(3);
		a.div2();
		QdAssert(to_string(a) == "1.5");
	}
	{
		DecimalFloat a(10);
		a.div2();
		QdAssert(to_string(a) == "5");
	}
	{
		DecimalFloat a(1100);
		a.div2();
		QdAssert(to_string(a) == "550");
	}



}
#endif



void test_constexpr()
{
	QD_CONSTEXPR dd_real ceil1 = ceil(dd_real(1.));
	QD_CONSTEXPR dd_real sqrt1 = sqrt(dd_real(1.));

    QD_CONSTEXPR dd_real sin1 = sin(dd_real(1.));
    QD_CONSTEXPR dd_real cos1 = cos(dd_real(1.));
    QD_CONSTEXPR dd_real tan1 = tan(dd_real(1.));
    QD_CONSTEXPR dd_real asin1 = asin(dd_real(1.));
    QD_CONSTEXPR dd_real acos1 = acos(dd_real(1.));
    QD_CONSTEXPR dd_real atan1 = atan(dd_real(1.));

    QD_CONSTEXPR dd_real exp2 = exp(dd_real(2.));
    QD_CONSTEXPR dd_real log2 = log(dd_real(2.));
    QD_CONSTEXPR dd_real sinh2 = sinh(dd_real(2.));
	QD_CONSTEXPR dd_real sinh0d01 = sinh(dd_real(0.01));
	QD_CONSTEXPR dd_real cosh2 = cosh(dd_real(2.));
    QD_CONSTEXPR dd_real tanh2 = tanh(dd_real(2.));
	QD_CONSTEXPR dd_real asinh2 = asinh(dd_real(2.));
	QD_CONSTEXPR dd_real acosh2 = acosh(dd_real(2.));
    QD_CONSTEXPR dd_real atanh2 = atanh(dd_real(2.));
}

void test_constexpr_qd()
{
	QD_CONSTEXPR qd_real sqrt1 = sqrt(qd_real(1.));

	QD_CONSTEXPR qd_real sin1 = sin(qd_real(1.));
	QD_CONSTEXPR qd_real cos1 = cos(qd_real(1.));
	QD_CONSTEXPR qd_real tan1 = tan(qd_real(1.));
	QD_CONSTEXPR qd_real asin1 = asin(qd_real(1.));
	QD_CONSTEXPR qd_real acos1 = acos(qd_real(1.));
	QD_CONSTEXPR qd_real atan1 = atan(qd_real(1.));

	QD_CONSTEXPR qd_real exp2 = exp(qd_real(2.));
	QD_CONSTEXPR qd_real log2 = log(qd_real(2.));
	QD_CONSTEXPR qd_real sinh2 = sinh(qd_real(2.));
	QD_CONSTEXPR qd_real sinh0d01 = sinh(qd_real(0.01));
	QD_CONSTEXPR qd_real cosh2 = cosh(qd_real(2.));
	QD_CONSTEXPR qd_real tanh2 = tanh(qd_real(2.));
	QD_CONSTEXPR qd_real asinh2 = asinh(qd_real(2.));
	QD_CONSTEXPR qd_real acosh2 = acosh(qd_real(2.));
	QD_CONSTEXPR qd_real atanh2 = atanh(qd_real(2.));
}

void test_drand()
{
	std::mt19937 gen;

	{
		double acc = 0;
		for (int i = 0; i < 1000 * 1000; ++i)
			acc += drand(gen);
		acc /= double(1000) * 1000;
		QdAssert(acc < 0.505 && acc > 0.495);
	}
	{
		double acc = 0;
		for (int i = 0; i < 1000 * 1000; ++i)
			acc += drand_fine(gen);
		acc /= double(1000) * 1000;
		QdAssert(acc < 0.505 && acc > 0.495);
	}

}

void test_ddrand()
{
	std::mt19937 gen;

	{
		dd_real acc = 0;
		for (int i = 0; i < 1000 * 1000; ++i)
			acc += ddrand(gen);
		acc /= double(1000) * 1000;
		QdAssert(acc < 0.505 && acc > 0.495);
	}

}

void test_qdrand()
{
	std::mt19937 gen;
	{
		qd_real acc = 0;
		for (int i = 0; i < 1000 * 1000; ++i)
			acc += qdrand(gen);
		acc /= double(1000) * 1000;
		QdAssert(acc < 0.505 && acc > 0.495);
	}

}

void test_read()
{
	QdAssert(to_int("0") == 0);
	QdAssert(to_int("1") == 1);
	QdAssert(to_int("2") == 2);
	QdAssert(to_int("11") == 11);
	QdAssert(to_int("-1") == -1);
	QdAssert(to_int("-11") == -11);
	QdAssert(to_int("2147483647") == 2147483647);
	QdAssert(to_int("-2147483648") == -2147483648);

	constexpr dd_real a0 = dd_real::read("0");
	constexpr dd_real a1 = dd_real::read("-0");
	constexpr dd_real a2 = dd_real::read("1");
	constexpr dd_real a3 = dd_real::read("-1");
	constexpr dd_real a4 = dd_real::read("1.");
	constexpr dd_real a5 = dd_real::read("-1.");
	constexpr dd_real a6 = dd_real::read("10");
	constexpr dd_real a7 = dd_real::read("-10");

	{
		constexpr dd_real dd_1 = dd_real::read("1");
		constexpr dd_real dd_10 = dd_real::read("1E1");
		constexpr qd_real qd_1 = qd_real::read("1");
		constexpr qd_real qd_10 = qd_real::read("1E1");
		constexpr dd_real dd_m1 = dd_real::read("-1");
		constexpr dd_real dd_m10 = dd_real::read("-1E1");
		constexpr qd_real qd_m1 = qd_real::read("-1");
		constexpr qd_real qd_m10 = qd_real::read("-1E1");

		constexpr dd_real dd_m1Em1 = dd_real::read("-1E-1");
		constexpr dd_real dd_m10Em1 = dd_real::read("-10E-1");
		constexpr qd_real qd_m1Em1 = qd_real::read("-1E-1");
		constexpr qd_real qd_m10Em1 = qd_real::read("-10E-1");
		qd_real qd_m10Em1_ = qd_real::read("-10E-1");
	}

	dd_real a;
	a = dd_real::read("0");
	a = dd_real::read("1");
	a = dd_real::read("2");
	a = dd_real::read("3");

	using namespace qd_literals;
	const auto vv = 1ull;

	QdAssert(1_dd == 1.);
	QdAssert(1._dd == 1.);

	QdAssert(1_qd == 1.);
	QdAssert(1._dd == 1.);
}

void test_dd_read()
{
	{
        using namespace qd_literals;
        std::stringstream ss("1.1E333");
        double a = 1;
        ss >> a;
        QdAssert(a == 0. && ss.eof() && ss.fail());
	}
	{
		using namespace qd_literals;
		std::stringstream ss("1.1E-333");
		double a = 1;
		ss >> a;
		QdAssert(a == 0. && ss.eof() && ss.fail());
	}
	{
		using namespace qd_literals;
		std::stringstream ss("a");
		std::string v;
		ss >> v;
		double a = 1;
		ss >> a;
		QdAssert(a == 1. && ss.eof() && ss.fail());
	}
	{
		using namespace qd_literals;
		std::stringstream ss("a");
		double a = 1;
		ss >> a;
		QdAssert(a == 0. && ss.fail());
	}
	{
		using namespace qd_literals;
		std::stringstream ss(" a");
		double a = 1;
		ss >> a;
		QdAssert(a == 0. && ss.fail());
	}
	{
		using namespace qd_literals;
		std::stringstream ss(" ");
		double a = 1;
		ss >> a;
		// check sentry failed
		QdAssert(a == 1. && ss.eof() && ss.fail());
	}
	{
		using namespace qd_literals;
		std::stringstream ss(".E");
		double a;
		ss >> a;
		QdAssert(!ss.eof() && ss.fail());
		ss.clear();
		QdAssert(ss.peek() == 'E');
	}
	{
		using namespace qd_literals;
		std::stringstream ss("1.1E1111111111111a");
		double a;
		ss >> a;
		QdAssert(!ss.eof() && ss.fail());
		ss.clear();
		QdAssert(ss.peek() == 'a');
	}
	{
		using namespace qd_literals;
		std::stringstream ss("1.1a");
		dd_real a;
		ss >> a;
		QdAssert(!a.isnan() && (ss.peek() == 'a') && abs(a / 1.1 - 1) < 0.0001);
	}
	{
		using namespace qd_literals;
		std::stringstream ss("11E-1");
		dd_real a;
		ss >> a;
		QdAssert(!a.isnan() && !ss.fail() && ss.eof() && abs(a / 1.1 - 1) < 0.0001);
	}
	{
		using namespace qd_literals;
		std::stringstream ss("0.11E+1");
		dd_real a;
		ss >> a;
		QdAssert(!a.isnan() && !ss.fail() && ss.eof() && abs(a / 1.1 - 1) < 0.0001);
	}
	{
		using namespace qd_literals;
		std::stringstream ss("0.11E+1111111111111111111111111111111111");
		dd_real a;
		ss >> a;
		QdAssert(a.isnan() && ss.fail() && ss.eof());
	}
	{
		using namespace qd_literals;
		std::stringstream ss("0.11E+1111111111111111111111111111111111a");
		dd_real a;
		ss >> a;
		QdAssert(a.isnan() && ss.fail() && !ss.eof());
		ss.clear();
		QdAssert(ss.peek() == 'a');
	}
	{
		using namespace qd_literals;
		std::stringstream ss("-0.11E1");
		dd_real a;
		ss >> a;
		bool ssf = ss.fail();
		QdAssert(!a.isnan() && !ss.fail() && ss.eof() && abs(a / -1.1 - 1) < 0.0001);
	}
	{
		using namespace qd_literals;
		std::stringstream ss(" 1 1");
		dd_real a;
		ss >> a;
		bool ssf = ss.fail();
		QdAssert(!a.isnan() && !ss.fail() && (ss.peek() == ' ') && abs(a - 1) < 0.0001);
	}
	{
		using namespace qd_literals;
		std::stringstream ss(".E");
		dd_real a;
		ss >> a;
		QdAssert(!ss.eof() && ss.fail());
		ss.clear();
		QdAssert(ss.peek() == 'E');
	}
	{
		std::stringstream ss("1 2 3");
		dd_real a, b, c;
		ss >> a >> b >> c;
		QdAssert(a == 1 && b == 2 && c == 3 && ss.eof());

	}
	{
		using namespace qd_literals;
		std::stringstream ss(" ");
		dd_real a = 1;
		ss >> a;
		// check sentry failed
		QdAssert(a == 1. && ss.eof() && ss.fail());
	}
}

void test_qd_read()
{

	{
		using namespace qd_literals;
		std::stringstream ss("1.1a");
		qd_real a;
		ss >> a;
		QdAssert(!a.isnan() && (ss.peek() == 'a') && abs(a / 1.1 - 1) < 0.0001);
	}
	{
		using namespace qd_literals;
		std::stringstream ss("11E-1");
		qd_real a;
		ss >> a;
		QdAssert(!a.isnan() && !ss.fail() && ss.eof() && abs(a / 1.1 - 1) < 0.0001);
	}
	{
		using namespace qd_literals;
		std::stringstream ss("0.11E+1");
		qd_real a;
		ss >> a;
		QdAssert(!a.isnan() && !ss.fail() && ss.eof() && abs(a / 1.1 - 1) < 0.0001);
	}
	{
		using namespace qd_literals;
		std::stringstream ss("0.11E+1111111111111111111111111111111111");
		qd_real a;
		ss >> a;
		QdAssert(a.isnan() && ss.fail() && ss.eof());
	}
	{
		using namespace qd_literals;
		std::stringstream ss("-0.11E1");
		qd_real a;
		ss >> a;
		bool ssf = ss.fail();
		QdAssert(!a.isnan() && !ss.fail() && ss.eof() && abs(a / -1.1 - 1) < 0.0001);
	}
	{
		using namespace qd_literals;
		std::stringstream ss(" 1 1");
		qd_real a;
		ss >> a;
		bool ssf = ss.fail();
		QdAssert(!a.isnan() && !ss.fail() && (ss.peek() == ' ') && abs(a - 1) < 0.0001);
	}
}

void test_nan()
{
    using namespace qd_literals;
	1. / dd_real::_inf;

    QdAssert((dd_real::_nan + 0.).isnan());
    QdAssert((dd_real::_nan + 0_dd).isnan());
    QdAssert((dd_real::_nan - 0.).isnan());
    QdAssert((dd_real::_nan - 0_dd).isnan());
    QdAssert((dd_real::_nan * 0.).isnan());
    QdAssert((dd_real::_nan * 0_dd).isnan());
    QdAssert((dd_real::_nan / 0.).isnan());
    QdAssert((dd_real::_nan / 0_dd).isnan());
    QdAssert((0. / dd_real::_nan).isnan());
    QdAssert((0_dd / dd_real::_nan).isnan());

	QdAssert((qd_real::_nan + 0.).isnan());
	QdAssert((qd_real::_nan + 0_dd).isnan());
	QdAssert((qd_real::_nan + 0_qd).isnan());
	QdAssert((qd_real::_nan - 0.).isnan());
	QdAssert((qd_real::_nan - 0_dd).isnan());
	QdAssert((qd_real::_nan - 0_qd).isnan());
	QdAssert((qd_real::_nan * 0.).isnan());
	QdAssert((qd_real::_nan * 0_dd).isnan());
	QdAssert((qd_real::_nan * 0_qd).isnan());
	QdAssert((qd_real::_nan / 0.).isnan());
	QdAssert((qd_real::_nan / 0_dd).isnan());
	QdAssert((qd_real::_nan / 0_qd).isnan());
	QdAssert((0. / qd_real::_nan).isnan());
	QdAssert((0_qd / qd_real::_nan).isnan());
	QdAssert((0_dd / qd_real::_nan).isnan());
}

//#include "qd/unused/DecimalFloat.h"

#include <type_traits>
int main() {

#ifdef QD_HAS_MPFR
	test_accuracy_set();
#endif

	test_dd_read();
	test_qd_read();
	test_read();
	test_constexpr();
	test_constexpr_qd();

	test_nan();

	test_drand();
	test_ddrand();
	test_qdrand();

	test_log();
	test_asinh();
	test_cmp();
	test_1();
	test_2();
	test_3();
	test_4();
	test_5();
	return (int)bool(nerr);
	return (int)bool(nerr);
}


