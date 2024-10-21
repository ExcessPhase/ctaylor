#include "cjacobian.h"


int main(int argc, char**argv)
{	using namespace jacobian;
	using namespace boost::mp11;

	if (argc != 4)
	{	std::cerr << argv[0] << ": Need 3 floating point arguments!" << std::endl;
		return 1;
	}
	const cjacobian<mp_list<mp_size_t<0> > > s0(std::atof(argv[1]), true);
	const cjacobian<mp_list<mp_size_t<1> > > s1(std::atof(argv[2]), true);
	const cjacobian<mp_list<mp_size_t<2> > > s2(std::atof(argv[3]), true);
	const auto s3 = if_(
		s0 > s1,
		[&](void)
		{	return s0 / s1;
		},
		[&](void)
		{	return s0 / s2;
		}
	);
	std::cerr << "s3=" << s3 << "\n";
	std::cerr << "s0/s2=" << s0/s2 << "\n";
	std::cerr << "s3.getDer(mp_size_t<0>())=" << s3.getDer(mp_size_t<0>()) << std::endl;
	std::cerr << "s3.getDer(mp_size_t<1>())=" << s3.getDer(mp_size_t<1>()) << std::endl;
	std::cerr << "s3.getDer(mp_size_t<2>())=" << s3.getDer(mp_size_t<2>()) << std::endl;
	std::cerr << "value(s3)=" << value(s3) << std::endl;
}
