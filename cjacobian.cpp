#include "cjacobian.h"


int main(int argc, char**argv)
{	using namespace jacobian;
	using namespace boost::mp11;

	if (argc != 6)
	{	std::cerr << argv[0] << ": Need 4 floating point arguments and one integer!" << std::endl;
		std::cerr << argv[0] << ": with testing -= and +=: " << argv[0] << " 1.2 1.3 1.4 1.5 1" << std::endl;
		std::cerr << argv[0] << ": without testing -= and +=: " << argv[0] << " 1.2 1.3 1.4 1.5 0" << std::endl;
		return 1;
	}
	else
	{	std::cout << "the 0th derivative (value) is printed last for cjacobian!" << std::endl;
		std::cout << "derivatives are printed sorted by enum!" << std::endl;
	}
		/// creation of 4 independent variables
		/// this is what the boolean argument is foe
		/// without we would simply create a value with derivative equal 0.0
	const cjacobian<mp_list<mp_size_t<0> > > s0(std::atof(argv[1]), true);
	const cjacobian<mp_list<mp_size_t<1> > > s1(std::atof(argv[2]), true);
	const cjacobian<mp_list<mp_size_t<2> > > s2(std::atof(argv[3]), true);
	const cjacobian<mp_list<mp_size_t<3> > > s3(std::atof(argv[4]), true);
		/// to create some intermediate value with derivatives vs the 4 independent variables
	auto s4 = -s0 + s1 - s2 + s1*s2 - s0*s1 + s2*s3;
		/// testing of the -= and ++ operator
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
		/// demonstration of chain=rule optimization
		/// create a copy of s4 with only a single derivative with enum=4
	const auto s41 = s4.convert2Independent(mp_size_t<4>());
		/// original expensive calculation (without chain-rule optimization)
	const auto s5 = fmod(s4*s4, 1.0 - s4*s4);
		/// do the same calculation like for s5 AND undo chain-rule optimization
		/// until chainRule() is being called, only a single derivative is being carried and calculated
		/// the last argument must be the same as passed above to convert2Independent()
		/// and it must be different than any other used ENUM -- see s0..s3 above which are using 0..3
	const auto s51 = fmod(s41*s41, 1.0 - s41*s41).chainRule(s4, mp_size_t<4>());;
		/// s4 is an intermediate result
	std::cerr << "s4=" << s4 << "\n";
		/// this is the value being printed by maxima.txt
	std::cerr << "compare to output of maxima.txt: s5=" << s5 << "\n";
		/// s41 is an intermediate result
	std::cerr << "s41=" << s41 << "\n";
		/// this is the value being printed by maxima.txt
	std::cerr << "compare to output of maxima.txt: s51=" << s51 << "\n";
}
