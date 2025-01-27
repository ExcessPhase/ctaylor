#include "test.h"
#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_CASE(taylor_25)
{	using namespace taylor;
	using namespace boost::mp11;
	const auto sMap = read("data25.txt");
	const auto s4 = getS4Taylor(false);
	const auto s5 = abs(s4);
	__EQUAL_TAYLOR__();
}
BOOST_AUTO_TEST_CASE(jacobian_25)
{	using namespace taylor;
	using namespace boost::mp11;
	const auto sMap = read("data25.txt", true);
	const auto s4 = getS4Jacobian();
	const auto s5 = abs(s4);
	__EQUAL_JACOBIAN__();
}
