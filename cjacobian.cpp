#include "cjacobian.h"


int main(int, char**)
{	using namespace jacobian;
	using namespace boost::mp11;

	const cjacobian<mp_list<mp_size_t<0> > > s0(1.2, true);
	const cjacobian<mp_list<mp_size_t<1> > > s1(1.3, true);

	std::cout << s0 * s1 << "\n";
}
