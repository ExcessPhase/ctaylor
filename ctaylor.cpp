#include "ctaylor.h"
/// the max number of terms of a polynomial in n variables and a maximum degree of m is (n+m)!/n!/m!
/// for n = 2 and m=4 this would be 6!/2!/4!=720/2/24=720/48=15
int main(int argc, char**argv)
{	using namespace taylor;
	using namespace boost::mp11;

		/// maximum order of derivatives
		/// careful -- the number of coefficients explodes
		/// and so does the compile time
		/// for MAX=1 it is cheaper to use a dual number class dedicated to first oder derivatives
	static constexpr std::size_t MAX = 2;
	if (argc != 6)
	{	std::cerr << "missing arguments -- must be 4 floating point numbers and one integer!" << std::endl;
		return 1;
	}
		/// create an independent variable for x0 (this is what the unused boolean is for)
	const auto s0 = ctaylor<makeIndependent<0>, MAX>(std::atof(argv[1]), false);
		/// create an independent variable for x1 (this is what the unused boolean is for)
	const auto s1 = ctaylor<makeIndependent<1>, MAX>(std::atof(argv[2]), false);
	const auto s2 = ctaylor<makeIndependent<2>, MAX>(std::atof(argv[3]), false);
	const auto s3 = ctaylor<makeIndependent<3>, MAX>(std::atof(argv[4]), false);
		/// some calculation
	auto s4 = -s0 + s1 - s2 + s1*s2 - s0*s1 + s2*s3;
	if (std::atoi(argv[5]))
	{
		std::cerr << "s4_0=" << s4 << "\n";
		s4 = -s0 + s1 - s2 + s1*s2;
		std::cerr << "s4_1=" << s4 << "=" << -s0 + s1 - s2 + s1*s2 << "\n";
		s4 -= s0*s1;
		std::cerr << "s4_2=" << s4 << "=" << -s0 + s1 - s2 + s1*s2 - s0*s1 << "\n";
		s4 += s2*s3;
		std::cerr << "s4_3=" << s4 << "=" << -s0 + s1 - s2 + s1*s2 - s0*s1 + s2*s3 << "\n";
	}
	const auto s41 = s4.convert2Independent(mp_size_t<4>());
	const auto s5 = fmod(s4*s4, 1.0 - s4*s4);
	const auto s52 = fmod(s41*s41, 1.0 - s41*s41);
	const auto s51 = s52.chainRule(s4, mp_size_t<4>());
		/// print the entire polynomial
	std::cout << "s4=" << s4 << "\n";
	std::cout << "s41=" << s41 << "\n";
	std::cout << "s5=" << s5 << "\n";
	std::cout << "s52=" << s52 << "\n";
	std::cout << "s51=" << s51 << "\n";
	typedef mp_list<
		pair<
			mp_size_t<0>,	// wrt x0
			mp_size_t<MAX>	// to the order of MAX
		>
	> LIST_OF_PAIRS;
	std::cout << "s5.getDer(LIST_OF_PAIRS())=" << s5.getDer(LIST_OF_PAIRS()) << "\n";
	typedef mp_list<
			/// elements must be sorted by variable-enum
		pair<
			mp_size_t<0>,	// wrt x0
			mp_size_t<1>	// to the order of 1
		>,
		pair<
			mp_size_t<1>,	// wrt x1
			mp_size_t<1>	// to the order of 1
		>
	> LIST_OF_PAIRS_2;
	std::cout << "s5.getDer(LIST_OF_PAIRS_2())=" << s5.getDer(LIST_OF_PAIRS_2()) << "\n";
	auto s6 = s5;
	s6 += s0;
	std::cout << s6 << "=" << s5 + s0 << "\n";
}
