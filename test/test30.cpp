#include "test.h"
#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_CASE(taylor_30)
{	using namespace taylor;
	using namespace boost::mp11;
	const auto sMap = read("data30.txt");
	const auto s4 = getS4Taylor(false);
	const auto s5 = log1p(s4);
	__EQUAL_TAYLOR__();
}
BOOST_AUTO_TEST_CASE(jacobian_30)
{	using namespace taylor;
	using namespace boost::mp11;
	const auto sMap = read("data30.txt", true);
	const auto s4 = getS4Jacobian();
	const auto s5 = log1p(s4);
	__EQUAL_JACOBIAN__();
}
