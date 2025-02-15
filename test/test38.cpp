#include "test.h"
#include <boost/test/unit_test.hpp>


template<typename X>
static auto f(const X&_rX)
{	return sqr(sqr(sqr(_rX) - 1.0));
}
static constexpr auto X0 = 1.23;

BOOST_AUTO_TEST_CASE(taylor_38)
{	using namespace taylor;
	using namespace boost::mp11;
	const auto sMap = read("data38.txt");
	constexpr std::size_t MAX = 10;
	auto const sX = ctaylor<makeIndependent<0>, MAX>(X0, true);
	const auto s5 = f(sX);
	__EQUAL_TAYLOR__();
}
BOOST_AUTO_TEST_CASE(jacobian_38)
{	using namespace jacobian;
	using namespace boost::mp11;
	const auto sMap = read("data38.txt", true);
	auto const sX = cjacobian<mp_list<mp_size_t<0> > >(X0, true);
	const auto s5 = f(sX);
	__EQUAL_JACOBIAN__();
}
