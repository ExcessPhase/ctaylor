#pragma once
#include <cmath>
#include <array>
#include <utility>
#include <boost/math/special_functions/polygamma.hpp>
namespace taylor
{
namespace implementation
{
template<typename>
struct divide_by_n_p_1_impl;
template<std::size_t ...POS>
struct divide_by_n_p_1_impl<std::index_sequence<POS...> >
{	static constexpr const std::array<double, sizeof...(POS)> value = {1.0/double(POS + 1)...};
};
template<std::size_t ...POS>
constexpr const std::array<double, sizeof...(POS)> divide_by_n_p_1_impl<std::index_sequence<POS...> >::value;
#if 0
template<>
struct divide_by_n_p_1_impl<std::index_sequence<> >
{	static constexpr const std::array<double, 0> value = {};
};
#endif
template<std::size_t SIZE>
struct divide_by_n_p_1
{	typedef divide_by_n_p_1_impl<std::make_index_sequence<SIZE> > type;
};
template<std::size_t SIZE>
std::array<double, SIZE> exp(double _d)
{	std::array<double, SIZE> s;
	double d = std::exp(_d);
	s[0] = d;
	auto &r = divide_by_n_p_1<SIZE - 1>::type::value;
	for (std::size_t i = 1; i < SIZE; ++i)
	{	d *= r[i - 1];
		s[i] = d;
	}
	return s;
}
template<std::size_t SIZE>
std::array<double, SIZE> expm1(double _d)
{	std::array<double, SIZE> s;
	double d = std::exp(_d);
	s[0] = d - 1.0;
	auto &r = divide_by_n_p_1<SIZE - 1>::type::value;
	for (std::size_t i = 1; i < SIZE; ++i)
	{	d *= r[i - 1];
		s[i] = d;
	}
	return s;
}
static const auto s_dLog2 = std::log(2.0);
template<std::size_t SIZE>
std::array<double, SIZE> exp2(double _d)
{	std::array<double, SIZE> s;
	double d = std::exp2(_d);
	s[0] = d;
	auto &r = divide_by_n_p_1<SIZE - 1>::type::value;
	std::transform(
		r.cbegin(),
		r.cend(),
		s.begin(),
		std::next(s.begin()),
		[](const double _dR, const double _dD)
		{       return _dR*_dD*s_dLog2;
		}
	);
	return s;
}
template<typename>
struct n_p_1_divided_by_n_p_2_impl;
template<std::size_t ...POS>
struct n_p_1_divided_by_n_p_2_impl<std::index_sequence<POS...> >
{	static constexpr const std::array<double, sizeof...(POS)> value = {double(POS + 1)/double(POS + 2)...};
};
template<std::size_t ...POS>
constexpr const std::array<double, sizeof...(POS)> n_p_1_divided_by_n_p_2_impl<std::index_sequence<POS...> >::value;
template<std::size_t SIZE>
struct n_p_1_divided_by_n_p_2
{	typedef n_p_1_divided_by_n_p_2_impl<std::make_index_sequence<SIZE> > type;
};
template<std::size_t SIZE>
std::array<double, SIZE> log(const double _d)
{	std::array<double, SIZE> s;
	const double d1 = 1.0/_d;
	double d = d1;
	s[0] = std::log(_d);
	s[1] = d;
	auto &r = n_p_1_divided_by_n_p_2<SIZE - 2>::type::value;
	for (std::size_t i = 2; i < SIZE; ++i)
	{	d *= -d1*r[i - 2];
		s[i] = d;
	}
	return s;
}
template<std::size_t SIZE>
std::array<double, SIZE> log1p(const double _d)
{	std::array<double, SIZE> s;
	const double d1 = 1.0/(1.0 + _d);
	double d = d1;
	s[0] = std::log1p(_d);
	s[1] = d;
	auto &r = n_p_1_divided_by_n_p_2<SIZE - 2>::type::value;
	for (std::size_t i = 2; i < SIZE; ++i)
	{	d *= -d1*r[i - 2];
		s[i] = d;
	}
	return s;
}
static const auto s_dInvLog10 = 1.0/std::log(10.0);
template<std::size_t SIZE>
std::array<double, SIZE> log10(const double _d)
{	std::array<double, SIZE> s;
	const double d1 = 1.0/_d;
	double d = d1*s_dInvLog10;
	s[0] = std::log10(_d);
	s[1] = d;
	auto &r = n_p_1_divided_by_n_p_2<SIZE - 2>::type::value;
	for (std::size_t i = 2; i < SIZE; ++i)
	{	d *= -d1*r[i - 2];
		s[i] = d;
	}
	return s;
}
static const auto s_dInvLog2 = 1.0/std::log(2.0);
template<std::size_t SIZE>
std::array<double, SIZE> log2(const double _d)
{	std::array<double, SIZE> s;
	const double d1 = 1.0/_d;
	double d = d1*s_dInvLog2;
	s[0] = std::log2(_d);
	s[1] = d;
	auto &r = n_p_1_divided_by_n_p_2<SIZE - 2>::type::value;
	for (std::size_t i = 2; i < SIZE; ++i)
	{	d *= -d1*r[i - 2];
		s[i] = d;
	}
	return s;
}
template<std::size_t SIZE>
std::array<double, SIZE> sin(const double _d)
{	std::array<double, SIZE> s;
	const double ds = std::sin(_d);
	const double dc = std::cos(_d);
	double d = 1.0;
	//0=sin
	//1=cos
	//2=-sin
	//3=-cos
	//4=0
#if 1
	s[0] = ds;
	auto &r = divide_by_n_p_1<SIZE - 1>::type::value;
	for (std::size_t i = 1; i < SIZE; ++i)
	{	d *= r[i - 1];
		s[i] = (i & 1 ? i & 2 ? -dc : dc : i & 2 ? -ds : ds)*d;
	}
#else
	for (std::size_t i = 0; i < SIZE; ++i)
	{	s[i] = (i & 1 ? i & 2 ? -dc : dc : i & 2 ? -ds : ds)*d;
		d /= i + 1;
	}
#endif
	return s;
}
template<std::size_t SIZE>
std::array<double, SIZE> cos(const double _d)
{	std::array<double, SIZE> s;
	const double ds = std::sin(_d);
	const double dc = std::cos(_d);
	double d = 1.0;
	//0=cos
	//1=-sin
	//2=-cos
	//3=sin
	//4=0
#if 1
	s[0] = dc;
	auto &r = divide_by_n_p_1<SIZE - 1>::type::value;
	for (std::size_t i = 1; i < SIZE; ++i)
	{	d *= r[i - 1];
		s[i] = (i & 1 ? i & 2 ? ds : -ds : i & 2 ? -dc: dc)*d;
	}
#else
	for (std::size_t i = 0; i < SIZE; ++i)
	{	s[i] = (i & 1 ? i & 2 ? ds : -ds : i & 2 ? -dc: dc)*d;
		d /= i + 1;
	}
#endif
	return s;
}
template<std::size_t SIZE>
std::array<double, SIZE> sinh(const double _d)
{	std::array<double, SIZE> s;
	const double ds = std::sinh(_d);
	const double dc = std::cosh(_d);
	double d = 1.0;
#if 1
	s[0] = ds;
	auto &r = divide_by_n_p_1<SIZE - 1>::type::value;
	for (std::size_t i = 1; i < SIZE; ++i)
	{	d *= r[i - 1];
		s[i] = (i & 1 ? dc : ds)*d;
	}
#else
	for (std::size_t i = 0; i < SIZE; ++i)
	{	s[i] = (i & 1 ? dc : ds)*d;
		d /= i + 1;
	}
#endif
	return s;
}
template<std::size_t SIZE>
std::array<double, SIZE> cosh(const double _d)
{	std::array<double, SIZE> s;
	const double ds = std::sinh(_d);
	const double dc = std::cosh(_d);
	double d = 1.0;
#if 1
	s[0] = dc;
	auto &r = divide_by_n_p_1<SIZE - 1>::type::value;
	for (std::size_t i = 1; i < SIZE; ++i)
	{	d *= r[i - 1];
		s[i] = (i & 1 ? ds : dc)*d;
	}
#else
	for (std::size_t i = 0; i < SIZE; ++i)
	{	s[i] = (i & 1 ? ds : dc)*d;
		d /= i + 1;
	}
#endif
	return s;
}
template<std::size_t SIZE>
std::array<double, SIZE> sqrt(const double _d)
{	std::array<double, SIZE> s;
	const double d0 = std::sqrt(_d);
	const double d1 = 1.0/_d;
	double dPow = 0.5;
	double d = d0;
	//0=x^0.5
	//1=0.5*x^-0.5
	//2=-0.25*x^-1.5
	//3=
#if 1
	s[0] = d;
	auto &r = divide_by_n_p_1<SIZE - 1>::type::value;
	for (std::size_t i = 1; i < SIZE; ++i)
	{	d *= d1*dPow*r[i - 1];	// sqrt(_d)/_d*0.5, sqrt(_d)/_d*0.5/_d*-0.5
		dPow -= 1.0;	// -0.5
		s[i] = d;	// sqrt(_d),
	}
#else
	for (std::size_t i = 0; i < SIZE; ++i)
	{	s[i] = d;	// sqrt(_d),
		d *= d1*dPow/(i + 1);	// sqrt(_d)/_d*0.5, sqrt(_d)/_d*0.5/_d*-0.5
		dPow -= 1.0;	// -0.5
	}
#endif
	return s;
}
template<std::size_t SIZE>
std::array<double, SIZE> cbrt(const double _d)
{	std::array<double, SIZE> s;
	const double d0 = std::cbrt(_d);
	const double d1 = 1.0/_d;
	double dPow = 1.0/3.0;
	double d = d0;
	//x^1.3
	//pow*previous/x
	//pow -= 1.0
#if 1
	s[0] = d;
	auto &r = divide_by_n_p_1<SIZE - 1>::type::value;
	for (std::size_t i = 1; i < SIZE; ++i)
	{	d *= d1*dPow*r[i - 1];
		dPow -= 1.0;
		s[i] = d;
	}
#else
	for (std::size_t i = 0; i < SIZE; ++i)
	{	s[i] = d;
		d *= d1*dPow/(i + 1);
		dPow -= 1.0;
	}
#endif
	return s;
}
template<std::size_t SIZE>
std::array<double, SIZE> inverse(const double _d)
{	std::array<double, SIZE> s;
	const double ds = 1.0/_d;
	double d = ds;
	//f0=x^-1/0!
	//f1=-x^-2/1!
	//f2=2x^-3/2!
	//f3=-6x^-4/3!
#if 1
	s[0] = d;
	for (std::size_t i = 1; i < SIZE; ++i)
	{	d *= ds;
		s[i] = i & 1 ? -d : d;
	}
#else
	for (std::size_t i = 0; i < SIZE; ++i)
	{	s[i] = i & 1 ? -d : d;
		d *= ds;
	}
#endif
	return s;
}
template<std::size_t SIZE>
std::array<double, SIZE> polygamma(const int _i, const double _d)
{	std::array<double, SIZE> s;
	auto d = 1.0;
	auto &r = divide_by_n_p_1<SIZE - 1>::type::value;
	for (std::size_t i = 0; i < SIZE; ++i)
	{	if (i > 1)
			d *= r[i-1];
		s[i] = boost::math::polygamma(_i + i, _d)*d;
	}
	return s;
}
template<std::size_t SIZE>
std::array<double, SIZE> lgamma(const double _d)
{	std::array<double, SIZE> s;
	s[0] = std::lgamma(_d);
	const auto sPG = polygamma<SIZE-1>(0, _d);
	auto &r = divide_by_n_p_1<SIZE - 1>::type::value;
	for (std::size_t i = 1; i < SIZE; ++i)
		s[i] = r[i - 1]*sPG[i - 1];
	return s;
}
}
}
