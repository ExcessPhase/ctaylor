#include "test.h"
#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_CASE(taylor_1)
{	using namespace taylor;
	using namespace boost::mp11;
	const auto sMap = read("data1.txt");
	const auto s4 = getS4Taylor();
	const auto s5 = fmod(s4*s4, 1.0 - s4*s4);
	__EQUAL_TAYLOR__();
}
BOOST_AUTO_TEST_CASE(taylor_chain_1)
{	using namespace taylor;
	using namespace boost::mp11;
	const auto sMap = read("data1.txt");
	const auto s4 = getS4Taylor();
	const auto s41 = s4.convert2Independent(mp_size_t<4>());
	const auto s51 = fmod(s41*s41, 1.0 - s41*s41);
	const auto s5 = s51.chainRule(s4, mp_size_t<4>());
	__EQUAL_TAYLOR__();
}
BOOST_AUTO_TEST_CASE(jacobian_1)
{	using namespace jacobian;
	using namespace boost::mp11;
	const auto sMap = read("data1.txt", true);
	const auto s4 = getS4Jacobian();
	const auto s5 = fmod(s4*s4, 1.0 - s4*s4);
	__EQUAL_JACOBIAN__();
}
BOOST_AUTO_TEST_CASE(jacobian_chain_1)
{	using namespace jacobian;
	using namespace boost::mp11;
	const auto sMap = read("data1.txt", true);
	const auto s4 = getS4Jacobian();
	const auto s41 = s4.convert2Independent(mp_size_t<4>());
	const auto s51 = fmod(s41*s41, 1.0 - s41*s41);
	const auto s5 = s51.chainRule(s4, mp_size_t<4>());
	__EQUAL_JACOBIAN__();
}
