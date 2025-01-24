#include "test.h"
#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_CASE(taylor_5)
{	using namespace taylor;
	using namespace boost::mp11;
	const auto sMap = read("data5.txt");
	const auto s4 = getS4Taylor();
	const auto s5 = erf(s4);
	__EQUAL_TAYLOR__();
}
BOOST_AUTO_TEST_CASE(jacobian_5)
{	using namespace taylor;
	using namespace boost::mp11;
	const auto sMap = read("data5.txt", true);
	const auto s4 = getS4Jacobian();
	const auto s5 = erf(s4);
	__EQUAL_JACOBIAN__();
}
