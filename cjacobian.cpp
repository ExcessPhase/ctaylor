#include "cjacobian.h"


int main(int, char**argv)
{	using namespace jacobian;
	using namespace boost::mp11;

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
	std::cerr << s3 << "\n";
	std::cerr << s0/s2 << "\n";
}
