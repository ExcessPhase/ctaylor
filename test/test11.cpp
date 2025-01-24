#include "test.h"
#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_CASE(taylor_11)
{	using namespace taylor;
	using namespace boost::mp11;
	const auto sMap = read("data11.txt");
	const auto s4 = getS4Taylor(false);
	const auto s5 = atanh(1.0/s4);
	__EQUAL_TAYLOR__();
}
BOOST_AUTO_TEST_CASE(jacobian_11)
{	using namespace taylor;
	using namespace boost::mp11;
	const auto sMap = read("data11.txt", true);
	const auto s4 = getS4Jacobian();
	const auto s5 = atanh(1.0/s4);
	__EQUAL_JACOBIAN__();
}
