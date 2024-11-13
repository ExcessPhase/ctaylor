#include "cjacobian.h"


int main(int argc, char**argv)
{	using namespace jacobian;
	using namespace boost::mp11;

	if (argc != 5)
	{	std::cerr << argv[0] << ": Need 4 floating point arguments!" << std::endl;
		return 1;
	}
	const cjacobian<mp_list<mp_size_t<0> > > s0(std::atof(argv[1]), true);
	const cjacobian<mp_list<mp_size_t<1> > > s1(std::atof(argv[2]), true);
	const cjacobian<mp_list<mp_size_t<2> > > s2(std::atof(argv[3]), true);
	const cjacobian<mp_list<mp_size_t<3> > > s3(std::atof(argv[4]), true);
	const auto s4 = if_(
		s0 > s1,
		[&](void)
		{	return s0 / s1;
		},
		[&](void)
		{	return s0 / s2;
		}
	);
	std::cerr << "s4=" << s4 << "\n";
	std::cerr << "s0/s2=" << s0/s2 << "\n";
	std::cerr << "s4.getDer(mp_size_t<0>())=" << s4.getDer(mp_size_t<0>()) << std::endl;
	std::cerr << "s4.getDer(mp_size_t<1>())=" << s4.getDer(mp_size_t<1>()) << std::endl;
	std::cerr << "s4.getDer(mp_size_t<2>())=" << s4.getDer(mp_size_t<2>()) << std::endl;
	std::cerr << "value(s4)=" << value(s4) << std::endl;
	std::cerr << "s3=" << s3 << std::endl;
	auto s5 = s4 + 0.0*s3;
	std::cerr << "s5=" << s5 << "\n";
	s5 += s3;
	std::cerr << "s5=" << s5 << "\n";
}
