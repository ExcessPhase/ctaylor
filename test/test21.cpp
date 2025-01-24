#include "test.h"
#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_CASE(taylor_21)
{	using namespace taylor;
	using namespace boost::mp11;
	const auto sMap = read("data21.txt");
	const auto s4 = getS4Taylor(false);
	const auto s5 = cbrt(s4);
	__EQUAL_TAYLOR__();
}
BOOST_AUTO_TEST_CASE(jacobian_21)
{	using namespace taylor;
	using namespace boost::mp11;
	const auto sMap = read("data21.txt", true);
	const auto s4 = getS4Jacobian();
	const auto s5 = cbrt(s4);
	__EQUAL_JACOBIAN__();
}