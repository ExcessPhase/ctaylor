#include "../ctaylor.h"
#include "../cjacobian.h"
#include <regex>
#include <map>
#include <fstream>
#include <numeric>


#define BOOST_TEST_MODULE mytests
#include <boost/test/included/unit_test.hpp>

typedef std::map<std::size_t, std::size_t> ind2Expr;
static std::size_t order(const ind2Expr&_r)
{	return std::accumulate(
		_r.cbegin(),
		_r.cend(),
		std::size_t(),
		[&](const std::size_t _i, const ind2Expr::value_type&_r)
		{	return _i + _r.second;
		}
	);
}
struct compare
{	bool operator()(const ind2Expr&_r0, const ind2Expr&_r1) const
	{	const auto i0 = order(_r0);
		const auto i1 = order(_r1);
		if (i0 < i1)
			return true;
		else
		if (i1 < i0)
			return false;
		else
			return std::lexicographical_compare(
				_r0.crbegin(),
				_r0.crend(),
				_r1.crbegin(),
				_r1.crend()
			);
	}
};
typedef std::map<ind2Expr, double, compare> taylorMap;
static taylorMap read(const char *const _p, const bool _bOnlyFirstOrder = false)
{	taylorMap sMap;
	std::ifstream sFile(_p);
	if (!sFile)
		throw std::runtime_error(std::string("Cannot open file \"") + _p + "\"!");
	std::string sLine;
	std::regex sRegEx(R"(^([+-]?(?:\d+(?:\.\d*)?)(?:[eE][+-]?\d+)?)((?:\*X\d+(?:\^\d+)?)*)\s*$)");
	std::regex sVarRegEx(R"((?:\*X(\d+)(?:\^(\d+))?)\s*)");

	//const auto sPath = std::filesystem::current_path();
	//std::cerr << sPath << "\n";
	while (std::getline(sFile, sLine))
	{	std::istringstream sSS(sLine);
		std::smatch sMatch;
		if (std::regex_match(sLine, sMatch, sRegEx))
		{ // Extract the numerical coefficient
			const double d = std::stod(sMatch[1]);
			ind2Expr sI2E;
			const std::string sVars = sMatch[2]; // match.position(5) gives position after the coefficient
			std::smatch sVarsMatch;
			auto pStart = sVars.cbegin();
			while (pStart != sVars.cend())
				if (std::regex_search(pStart, sVars.cend(), sVarsMatch, sVarRegEx))
				{	//const std::string &sVar = sVarsMatch[1];
					const auto iID = std::atoi(std::string(sVarsMatch[1]).c_str());
					const std::string sExp = sVarsMatch[2];
					const auto iExp = sExp.empty() ? 1 : std::atoi(sExp.c_str());
					sI2E.emplace(iID, iExp);
					pStart = sVarsMatch.suffix().first;
					//sVars = sVars.substr(sVarsMatch.str(0).size());
				}
				else
					break;
			if (!_bOnlyFirstOrder || order(sI2E) < 2)
				sMap.emplace(sI2E, d);
		}
	}
	return sMap;
}
static auto sIsIdentical(const taylorMap::value_type&_r, const double _d)
{	static constexpr auto epsilon = 1e-12;
	return std::abs(_r.second - _d) < std::abs(epsilon*_r.second);
}
template<typename THREE = boost::mp11::mp_size_t<3> >
static auto getS4Taylor(const bool _b = false, const THREE& = THREE())
{	using namespace taylor;
	//using namespace boost::mp11;
	constexpr const std::size_t MAX = THREE::value;
	const auto s0 = ctaylor<makeIndependent<0>, MAX>(1.2, false);
		/// create an independent variable for x1 (this is what the unused boolean is for)
	const auto s1 = ctaylor<makeIndependent<1>, MAX>(1.3, false);
	const auto s2 = ctaylor<makeIndependent<2>, MAX>(1.4, false);
	const auto s3 = ctaylor<makeIndependent<3>, MAX>(1.5, false);
		/// some calculation
	if (_b)
	{	decltype(-s0 + s1 - s2 + s1*s2 - s0*s1 + s2*s3) s4 = 0.0;
		s4 -= s0;
		s4 += s1;
		s4 -= s2;
		s4 += s1*s2;
		s4 -= s0*s1;
		s4 += s2*s3;
		return s4;
	}
	else
		return -s0 + s1 - s2 + s1*s2 - s0*s1 + s2*s3;
}
static auto getS4Jacobian(const bool _b = false)
{	using namespace jacobian;
	using namespace boost::mp11;
	const auto s0 = cjacobian<mp_list<mp_size_t<0> > >(1.2, false);
		/// create an independent variable for x1 (this is what the unused boolean is for)
	const auto s1 = cjacobian<mp_list<mp_size_t<1> > >(1.3, false);
	const auto s2 = cjacobian<mp_list<mp_size_t<2> > >(1.4, false);
	const auto s3 = cjacobian<mp_list<mp_size_t<3> > >(1.5, false);
		/// some calculation
	if (_b)
	{	decltype(-s0 + s1 - s2 + s1*s2 - s0*s1 + s2*s3) s4 = 0.0;
		s4 -= s0;
		s4 += s1;
		s4 -= s2;
		s4 += s1*s2;
		s4 -= s0*s1;
		s4 += s2*s3;
		return s4;
	}
	else
		return -s0 + s1 - s2 + s1*s2 - s0*s1 + s2*s3;
}
#define __EQUAL_TAYLOR__()\
do\
{	BOOST_CHECK(s5.m_s.size() == sMap.size());\
	BOOST_CHECK(\
		std::equal(\
			sMap.cbegin(),\
			sMap.cend(),\
			s5.m_s.cbegin(),\
			sIsIdentical\
		)\
	);\
} while (false)
#define __EQUAL_JACOBIAN__()\
do\
{	BOOST_CHECK(s5.m_s.size() == sMap.size());\
	BOOST_CHECK(\
		std::equal(\
			std::next(sMap.cbegin()),\
			sMap.cend(),\
			s5.m_s.cbegin(),\
			sIsIdentical\
		)\
	);\
	BOOST_CHECK(\
		sIsIdentical(*sMap.cbegin(), value(s5))\
	);\
} while (false)
BOOST_AUTO_TEST_CASE(taylor_0)
{	using namespace taylor;
	using namespace boost::mp11;
	const auto sMap = read("data0.txt");
	const auto s4 = getS4Taylor();
	const auto s5 = exp(-1.0/(s4*s4));
	__EQUAL_TAYLOR__();
}
BOOST_AUTO_TEST_CASE(taylor_1_0)
{	using namespace taylor;
	using namespace boost::mp11;
	const auto sMap = read("data0.txt");
	const auto s4 = getS4Taylor(true);
	const auto s5 = exp(-1.0/(s4*s4));
	__EQUAL_TAYLOR__();
}
BOOST_AUTO_TEST_CASE(taylor_chain_0)
{	using namespace taylor;
	using namespace boost::mp11;
	const auto sMap = read("data0.txt");
	const auto s4 = getS4Taylor();
	const auto s41 = s4.convert2Independent(mp_size_t<4>());
	const auto s51 = exp(-1.0/(s41*s41));
	const auto s5 = s51.chainRule(s4, mp_size_t<4>());
	__EQUAL_TAYLOR__();
}
BOOST_AUTO_TEST_CASE(jacobian_0)
{	using namespace jacobian;
	using namespace boost::mp11;
	const auto sMap = read("data0.txt", true);
	const auto s4 = getS4Jacobian();
	const auto s5 = exp(-1.0/(s4*s4));
	__EQUAL_JACOBIAN__();
}
BOOST_AUTO_TEST_CASE(jacobian_1_0)
{	using namespace jacobian;
	using namespace boost::mp11;
	const auto sMap = read("data0.txt", true);
	const auto s4 = getS4Jacobian(true);
	const auto s5 = exp(-1.0/(s4*s4));
	__EQUAL_JACOBIAN__();
}
BOOST_AUTO_TEST_CASE(jacobian_chain_0)
{	using namespace jacobian;
	using namespace boost::mp11;
	const auto sMap = read("data0.txt", true);
	const auto s4 = getS4Jacobian();
	const auto s41 = s4.convert2Independent(mp_size_t<4>());
	const auto s51 = exp(-1.0/(s41*s41));
	const auto s5 = s51.chainRule(s4, mp_size_t<4>());
	__EQUAL_JACOBIAN__();
}
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
BOOST_AUTO_TEST_CASE(taylor_2)
{	using namespace taylor;
	using namespace boost::mp11;
	const auto sMap = read("data2.txt");
	const auto s4 = getS4Taylor();
	const auto s5 = hypot(s4*s4, 1.0 - s4*s4);
	__EQUAL_TAYLOR__();
}
BOOST_AUTO_TEST_CASE(taylor_chain_2)
{	using namespace taylor;
	using namespace boost::mp11;
	const auto sMap = read("data2.txt");
	const auto s4 = getS4Taylor();
	const auto s41 = s4.convert2Independent(mp_size_t<4>());
	const auto s51 = hypot(s41*s41, 1.0 - s41*s41);
	const auto s5 = s51.chainRule(s4, mp_size_t<4>());
	__EQUAL_TAYLOR__();
}
BOOST_AUTO_TEST_CASE(jacobian_2)
{	using namespace jacobian;
	using namespace boost::mp11;
	const auto sMap = read("data2.txt", true);
	const auto s4 = getS4Jacobian();
	const auto s5 = hypot(s4*s4, 1.0 - s4*s4);
	__EQUAL_JACOBIAN__();
}
BOOST_AUTO_TEST_CASE(jacobian_chain_2)
{	using namespace jacobian;
	using namespace boost::mp11;
	const auto sMap = read("data2.txt", true);
	const auto s4 = getS4Jacobian();
	const auto s41 = s4.convert2Independent(mp_size_t<4>());
	const auto s51 = hypot(s41*s41, 1.0 - s41*s41);
	const auto s5 = s51.chainRule(s4, mp_size_t<4>());
	__EQUAL_JACOBIAN__();
}
BOOST_AUTO_TEST_CASE(taylor_3)
{	using namespace taylor;
	using namespace boost::mp11;
	const auto sMap = read("data3.txt");
	const auto s4 = getS4Taylor();
	const auto s5 = atan(1.0/s4 - s4*s4);
	__EQUAL_TAYLOR__();
}
BOOST_AUTO_TEST_CASE(taylor_chain_3)
{	using namespace taylor;
	using namespace boost::mp11;
	const auto sMap = read("data3.txt");
	const auto s4 = getS4Taylor();
	const auto s41 = s4.convert2Independent(mp_size_t<4>());
	const auto s51 = atan(1.0/s41 - s41*s41);
	const auto s5 = s51.chainRule(s4, mp_size_t<4>());
	__EQUAL_TAYLOR__();
}
BOOST_AUTO_TEST_CASE(jacobian_3)
{	using namespace jacobian;
	using namespace boost::mp11;
	const auto sMap = read("data3.txt", true);
	const auto s4 = getS4Jacobian();
	const auto s5 = atan(1.0/s4 - s4*s4);
	__EQUAL_JACOBIAN__();
}
BOOST_AUTO_TEST_CASE(jacobian_chain_3)
{	using namespace jacobian;
	using namespace boost::mp11;
	const auto sMap = read("data3.txt", true);
	const auto s4 = getS4Jacobian();
	const auto s41 = s4.convert2Independent(mp_size_t<4>());
	const auto s51 = atan(1.0/s41 - s41*s41);
	const auto s5 = s51.chainRule(s4, mp_size_t<4>());
	__EQUAL_JACOBIAN__();
}
template<typename T>
static auto sqr(const T&_d)
{	return _d*_d;
}
BOOST_AUTO_TEST_CASE(taylor_4)
{	using namespace taylor;
	using namespace boost::mp11;
	const auto sMap = read("data4.txt");
	const auto s4 = getS4Taylor(false, mp_size_t<4>());
	const auto s5 = exp(-sqr(atan(1.0/s4 - s4*s4)));
	__EQUAL_TAYLOR__();
}
BOOST_AUTO_TEST_CASE(jacobian_4)
{	using namespace taylor;
	using namespace boost::mp11;
	const auto sMap = read("data4.txt", true);
	const auto s4 = getS4Jacobian();
	const auto s5 = exp(-sqr(atan(1.0/s4 - s4*s4)));
	__EQUAL_JACOBIAN__();
}
BOOST_AUTO_TEST_CASE(taylor_5)
{	using namespace taylor;
	using namespace boost::mp11;
	const auto sMap = read("data5.txt");
	const auto s4 = getS4Taylor(false, mp_size_t<5>());
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
