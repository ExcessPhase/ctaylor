#include "cjacobian.h"


int main(int argc, char**argv)
{	using namespace jacobian;
	using namespace boost::mp11;

	if (argc != 6)
	{	std::cerr << argv[0] << ": Need 4 floating point arguments and one integer!" << std::endl;
		return 1;
	}
	const cjacobian<mp_list<mp_size_t<0> > > s0(std::atof(argv[1]), true);
	const cjacobian<mp_list<mp_size_t<1> > > s1(std::atof(argv[2]), true);
	const cjacobian<mp_list<mp_size_t<2> > > s2(std::atof(argv[3]), true);
	const cjacobian<mp_list<mp_size_t<3> > > s3(std::atof(argv[4]), true);
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
		/// create a copy of s4 with only a single derivative with enum=4
	const auto s41 = s4.convert2Independent(mp_size_t<4>());
	const auto s5 = fmod(s4*s4, 1.0 - s4*s4);
		/// do the same calculation like for s5 AND undo chain-rule optimization
		/// until chainRule() is being called, only a single derivative is being carried and calculated
		/// the last argument must be the same as passed above to convert2Independent()
		/// and it must be different than any other used ENUM -- see s0..s3 above which are using 0..3
	const auto s51 = fmod(s41*s41, 1.0 - s41*s41).chainRule(s4, mp_size_t<4>());;
	std::cerr << "s4=" << s4 << "\n";
	std::cerr << "s5=" << s5 << "\n";
	std::cerr << "s41=" << s41 << "\n";
	std::cerr << "s51=" << s51 << "\n";
}
