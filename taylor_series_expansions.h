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
template<std::size_t SIZE>
struct divide_by_n_p_1
{	typedef divide_by_n_p_1_impl<std::make_index_sequence<SIZE> > type;
};
template<std::size_t N>
struct inverseFactorial;
template<typename>
struct inverseFactorialArrayImpl;
template<std::size_t ...ARGS>
struct inverseFactorialArrayImpl<std::index_sequence<ARGS...> >
{	static constexpr const std::array<double, sizeof...(ARGS)> value = {inverseFactorial<ARGS>::value...};
};
template<std::size_t ...ARGS>
constexpr const std::array<double, sizeof...(ARGS)> inverseFactorialArrayImpl<std::index_sequence<ARGS...> >::value;
template<std::size_t N>
struct inverseFactorialArray
{	typedef inverseFactorialArrayImpl<std::make_index_sequence<N> > type;
};
template<std::size_t N>
struct inverseFactorial
{	static constexpr const double value = inverseFactorial<N-1>::value/N;
};
template<>
struct inverseFactorial<0>
{	static constexpr const double value = 1.0;
};
template<std::size_t SIZE>
std::array<double, SIZE> exp(double _d)
{	std::array<double, SIZE> s;
	const double d = std::exp(_d);
	s[0] = d;
	auto &r = divide_by_n_p_1<SIZE - 1>::type::value;
	std::transform(
		r.cbegin(),
		r.cend(),
		s.cbegin(),
		std::next(s.begin()),
		[](const double _dR, const double _dD)
		{	return _dR*_dD;
		}
	);
	return s;
}
template<std::size_t SIZE>
std::array<double, SIZE> expm1(double _d)
{	std::array<double, SIZE> s;
	const double d = std::exp(_d);
	auto &r = divide_by_n_p_1<SIZE - 1>::type::value;
	s[0] = d;
	std::transform(
		r.cbegin(),
		r.cend(),
		s.cbegin(),
		std::next(s.begin()),
		[](const double _dR, const double _dD)
		{	return _dR*_dD;
		}
	);
	s[0] -= 1.0;
	return s;
}
static const auto s_dLog2 = std::log(2.0);
template<std::size_t SIZE>
std::array<double, SIZE> exp2(double _d)
{	std::array<double, SIZE> s;
	const double d = std::exp2(_d);
	s[0] = d;
	auto &r = divide_by_n_p_1<SIZE - 1>::type::value;
	std::transform(
		r.cbegin(),
		r.cend(),
		s.cbegin(),
		std::next(s.begin()),
		[](const double _dR, const double _dD)
		{	return _dR*_dD*s_dLog2;
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
	//double d = d1;
	s[0] = std::log(_d);
	s[1] = d1;
	auto &r = n_p_1_divided_by_n_p_2<SIZE - 2>::type::value;
	std::transform(
		r.cbegin(),
		r.cend(),
		std::next(s.cbegin()),
		std::next(s.begin(), 2),
		[d1](const double _dR, const double _dD)
		{	return -d1*_dR*_dD;
		}
	);
	return s;
}
template<std::size_t SIZE>
std::array<double, SIZE> log1p(const double _d)
{	std::array<double, SIZE> s;
	const double d1 = 1.0/(1.0 + _d);
	s[0] = std::log1p(_d);
	s[1] = d1;
	auto &r = n_p_1_divided_by_n_p_2<SIZE - 2>::type::value;
	std::transform(
		r.cbegin(),
		r.cend(),
		std::next(s.cbegin()),
		std::next(s.begin(), 2),
		[d1](const double _dR, const double _dD)
		{	return -d1*_dR*_dD;
		}
	);
	return s;
}
static const auto s_dInvLog10 = 1.0/std::log(10.0);
template<std::size_t SIZE>
std::array<double, SIZE> log10(const double _d)
{	std::array<double, SIZE> s;
	const double d1 = 1.0/_d;
	s[0] = std::log10(_d);
	s[1] = d1*s_dInvLog10;
	auto &r = n_p_1_divided_by_n_p_2<SIZE - 2>::type::value;
	std::transform(
		r.cbegin(),
		r.cend(),
		std::next(s.cbegin()),
		std::next(s.begin(), 2),
		[d1](const double _dR, const double _dD)
		{	return -d1*_dR*_dD;
		}
	);
	return s;
}
static const auto s_dInvLog2 = 1.0/std::log(2.0);
template<std::size_t SIZE>
std::array<double, SIZE> log2(const double _d)
{	std::array<double, SIZE> s;
	const double d1 = 1.0/_d;
	s[0] = std::log2(_d);
	s[1] = d1*s_dInvLog2;
	auto &r = n_p_1_divided_by_n_p_2<SIZE - 2>::type::value;
	std::transform(
		r.cbegin(),
		r.cend(),
		std::next(s.cbegin()),
		std::next(s.begin(), 2),
		[d1](const double _dR, const double _dD)
		{	return -d1*_dR*_dD;
		}
	);
	return s;
}
template<std::size_t SIZE>
std::array<double, SIZE> sin(const double _d)
{	std::array<double, SIZE> s;
	const double ds = std::sin(_d);
	const double dc = std::cos(_d);
	//0=sin
	//1=cos
	//2=-sin
	//3=-cos
	//4=0
	//auto &r = divide_by_n_p_1<SIZE - 1>::type::value;
	auto &r = inverseFactorialArray<SIZE>::type::value;
	std::transform(
		r.cbegin(),
		r.cend(),
		s.begin(),
		[&](const double &_dR)
		{	const auto i = &_dR - r.data();
			return (i & 1 ? i & 2 ? -dc : dc : i & 2 ? -ds : ds)*_dR;
		}
	);
	return s;
}
template<std::size_t SIZE>
std::array<double, SIZE> cos(const double _d)
{	std::array<double, SIZE> s;
	const double ds = std::sin(_d);
	const double dc = std::cos(_d);
	auto &r = inverseFactorialArray<SIZE>::type::value;
	std::transform(
		r.cbegin(),
		r.cend(),
		s.begin(),
		[&](const double &_dR)
		{	const auto i = &_dR - r.data();
			return (i & 1 ? i & 2 ? ds : -ds : i & 2 ? -dc : dc)*_dR;
		}
	);
	return s;
}
template<std::size_t SIZE>
std::array<double, SIZE> sinh(const double _d)
{	std::array<double, SIZE> s;
	const double ds = std::sinh(_d);
	const double dc = std::cosh(_d);
	auto &r = inverseFactorialArray<SIZE>::type::value;
	std::transform(
		r.cbegin(),
		r.cend(),
		s.begin(),
		[&](const double &_dR)
		{	const auto i = &_dR - r.data();
			return (i & 1 ? dc : ds)*_dR;
		}
	);
	return s;
}
template<std::size_t SIZE>
std::array<double, SIZE> cosh(const double _d)
{	std::array<double, SIZE> s;
	const double ds = std::sinh(_d);
	const double dc = std::cosh(_d);
	auto &r = inverseFactorialArray<SIZE>::type::value;
	std::transform(
		r.cbegin(),
		r.cend(),
		s.begin(),
		[&](const double &_dR)
		{	const auto i = &_dR - r.data();
			return (i & 1 ? ds : dc)*_dR;
		}
	);
	return s;
}
template<std::size_t SIZE>
std::array<double, SIZE> sqrt(const double _d)
{	std::array<double, SIZE> s;
	const double d0 = std::sqrt(_d);
	const double d1 = 1.0/_d;
	auto &r = divide_by_n_p_1<SIZE - 1>::type::value;
	s[0] = d0;
	std::transform(
		r.cbegin(),
		r.cend(),
		s.cbegin(),
		std::next(s.begin()),
		[d1, &r](const double &_dR, const double _dP)
		{	return _dP*d1*(0.5 - (&_dR - r.data()))*_dR;
		}
	);
	return s;
}
template<std::size_t SIZE>
std::array<double, SIZE> cbrt(const double _d)
{	std::array<double, SIZE> s;
	const double d0 = std::cbrt(_d);
	const double d1 = 1.0/_d;
	auto &r = divide_by_n_p_1<SIZE - 1>::type::value;
	s[0] = d0;
	std::transform(
		r.cbegin(),
		r.cend(),
		s.cbegin(),
		std::next(s.begin()),
		[d1, &r](const double &_dR, const double _dP)
		{	constexpr const auto dOneThird = 1.0/3.0;
			return _dP*d1*(dOneThird - (&_dR - r.data()))*_dR;
		}
	);
	return s;
}
template<std::size_t SIZE>
std::array<double, SIZE> inverse(const double _d)
{	std::array<double, SIZE> s;
	const double ds = -1.0/_d;
	//f0=x^-1/0!
	//f1=-x^-2/1!
	//f2=2x^-3/2!
	//f3=-6x^-4/3!
#if 1
	s[0] = -ds;
	std::transform(
		s.cbegin(),
		std::prev(s.cend()),
		std::next(s.begin()),
		[ds](const double _d)
		{	return ds*_d;
		}
	);
#else
	double d = -1.0;
	for (std::size_t i = 0; i < SIZE; ++i)
	{	d *= ds;
		s[i] = d;
	}
#endif
	return s;
}
template<std::size_t SIZE>
std::array<double, SIZE> polygamma(const int _i, const double _d)
{	std::array<double, SIZE> s;
	auto &r = inverseFactorialArray<SIZE>::type::value;
	std::transform(
		r.cbegin(),
		r.cend(),
		s.begin(),
		[&r, _i, _d](const double &_dR)
		{	return _dR*boost::math::polygamma(_i + int(&_dR - r.data()), _d);
		}
	);
	return s;
}
template<std::size_t SIZE>
std::array<double, SIZE> lgamma(const double _d)
{	std::array<double, SIZE> s;
	s[0] = std::lgamma(_d);
	const auto sPG = polygamma<SIZE-1>(0, _d);
	auto &r = divide_by_n_p_1<SIZE - 1>::type::value;
	std::transform(
		r.cbegin(),
		r.cend(),
		sPG.cbegin(),
		std::next(s.begin()),
		[](const auto _dR, const double _dPG)
		{	return _dR*_dPG;
		}
	);
	return s;
}
}
}
