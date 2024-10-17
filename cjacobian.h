#pragma once
#include <iostream>
#include <type_traits>
#include <array>
#include <limits>
#include <functional>
#include <boost/mp11.hpp>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include "merge_sorted_sets.h"
namespace jacobian
{
namespace implementation
{
using namespace boost::mp11;


template<typename, typename>
struct createStdArrayImpl;
template<typename VECTOR, std::size_t ...I>
struct createStdArrayImpl<VECTOR, std::index_sequence<I...> >
{	static constexpr const std::array<std::size_t, sizeof...(I)> value = {mp_at_c<VECTOR, I>::value...};
};
template<typename VECTOR, std::size_t ...I>
constexpr const std::array<std::size_t, sizeof...(I)> createStdArrayImpl<VECTOR, std::index_sequence<I...> >::value;
template<typename VECTOR>
struct createStdArray
{	typedef createStdArrayImpl<VECTOR, std::make_index_sequence<mp_size<VECTOR>::value> > type;
};
template<typename, typename>
struct createStdArrayImpl2;
template<typename VECTOR, std::size_t ...I>
struct createStdArrayImpl2<VECTOR, std::index_sequence<I...> >
{	static constexpr const std::array<std::pair<std::size_t, std::size_t>, sizeof...(I)> value = {
		std::make_pair(mp_first<mp_at_c<VECTOR, I> >::value, mp_second<mp_at_c<VECTOR, I> >::value)...
	};
};
template<typename VECTOR, std::size_t ...I>
constexpr const std::array<std::pair<std::size_t, std::size_t>, sizeof...(I)> createStdArrayImpl2<VECTOR, std::index_sequence<I...> >::value;
template<typename VECTOR>
struct createStdArray2
{	typedef createStdArrayImpl2<VECTOR, std::make_index_sequence<mp_size<VECTOR>::value> > type;
};
template<typename TARGET, typename SOURCE>
struct createIndicies
{	template<typename STATE, typename TARGET_VALUE>
	using findSourcePosition = mp_push_back<
		STATE,
		mp_find<
			SOURCE,
			TARGET_VALUE
		>
	>;
	typedef mp_fold<
		TARGET,
		mp_list<>,
		findSourcePosition
	> type;
};
template<typename TARGET, typename SOURCE0, typename SOURCE1>
struct createIndicies2
{	template<typename STATE, typename TARGET_VALUE>
	using findSourcePosition = mp_push_back<
		STATE,
		mp_list<
			mp_find<
				SOURCE0,
				TARGET_VALUE
			>,
			mp_find<
				SOURCE1,
				TARGET_VALUE
			>
		>
	>;
	typedef mp_fold<
		TARGET,
		mp_list<>,
		findSourcePosition
	> type;
};
template<typename T0, typename T1>
struct merge
{	typedef typename taylor::merge_sorted_sets<
		mp_less,
		T0,
		T1
	>::type type;
};
template<typename VECTOR>
struct cjacobian
{	static constexpr std::size_t SIZE = mp_size<VECTOR>::value + 1;
	typedef std::array<double, SIZE> ARRAY;
	ARRAY m_s;
	cjacobian(void) = default;
	cjacobian(const cjacobian&) = default;
	cjacobian(cjacobian&&) = default;
	cjacobian&operator=(const cjacobian&) = default;
	cjacobian&operator=(cjacobian&&) = default;
	cjacobian(const double _d)
		:m_s({})
	{	m_s.back() = _d;
	}
	cjacobian(const double _d, const bool)
		:m_s({1.0, _d})
	{	static_assert(SIZE == 2, "size must be exactly 2!");
	}
	friend double value(const cjacobian&_r)
	{	return _r.m_s.back();
	}
	template<typename T1>
	cjacobian(const cjacobian<T1>&_r)
		:m_s(convert<T1>(_r.m_s))
	{
	}
	template<typename T1>
	static ARRAY convert(const typename cjacobian<T1>::ARRAY&_r)
	{	ARRAY s;
		typedef typename createIndicies<VECTOR, T1>::type INDICIES;
		auto &r = createStdArray<INDICIES>::type::value;
		std::transform(
			r.begin(),
			r.end(),
			s.begin(),
			[&](const std::size_t _i)
			{	if (_i == cjacobian<T1>::SIZE - 1)
					return 0.0;
				else
					return _r[_i];
			}
		);
		s.back() = _r.back();
		return s;
	}
	template<typename T1>
	cjacobian&operator=(const cjacobian<T1>&_r)
	{	m_s = convert<T1>(_r.m_s);
		return *this;
	}
	template<typename T1>
	auto operator*(const cjacobian<T1>&_r) const
	{	typedef typename merge<VECTOR, T1>::type MERGED;
		typedef typename createIndicies2<MERGED, VECTOR, T1>::type INDICIES;
		auto &r = createStdArray2<INDICIES>::type::value;
		cjacobian<MERGED> s;
		std::transform(
			r.begin(),
			r.end(),
			s.m_s.begin(),
			[&](const std::pair<std::size_t, std::size_t>&_rI)
			{	if (_rI.first == SIZE - 1)
					return _r.m_s[_rI.second]*m_s.back();
				else
				if (_rI.second == cjacobian<T1>::SIZE - 1)
					return m_s[_rI.first]*_r.m_s.back();
				else
					return m_s[_rI.first]*_r.m_s.back() + _r.m_s[_rI.second]*m_s.back();
			}
		);
		s.m_s.back() = m_s.back() * _r.m_s.back();
		return s;
	}
	template<typename T1>
	auto operator/(const cjacobian<T1>&_r) const
	{	typedef typename merge<VECTOR, T1>::type MERGED;
		typedef typename createIndicies2<MERGED, VECTOR, T1>::type INDICIES;
		auto &r = createStdArray2<INDICIES>::type::value;
		cjacobian<MERGED> s;
		const auto dInv = 1.0/_r.m_s.back();
		const auto dInv2 = dInv*dInv*m_s.back();
		std::transform(
			r.begin(),
			r.end(),
			s.m_s.begin(),
			[&](const std::pair<std::size_t, std::size_t>&_rI)
			{	if (_rI.first == SIZE - 1)
					return -dInv2*_r.m_s[_rI.second];
				else
				if (_rI.second == cjacobian<T1>::SIZE - 1)
					return dInv*m_s[_rI.first];
				else
					return dInv*m_s[_rI.first] - dInv2*_r.m_s[_rI.second];
			}
		);
		s.m_s.back() = m_s.back() * dInv;
		return s;
	}
	template<typename T1>
	auto operator+(const cjacobian<T1>&_r) const
	{	typedef typename merge<VECTOR, T1>::type MERGED;
		typedef typename createIndicies2<MERGED, VECTOR, T1>::type INDICIES;
		auto &r = createStdArray2<INDICIES>::type::value;
		cjacobian<MERGED> s;
		std::transform(
			r.begin(),
			r.end(),
			s.m_s.begin(),
			[&](const std::pair<std::size_t, std::size_t>&_rI)
			{	if (_rI.first == SIZE - 1)
					return _r.m_s[_rI.second];
				else
				if (_rI.second == cjacobian<T1>::SIZE - 1)
					return m_s[_rI.first];
				else
					return m_s[_rI.first] + _r.m_s[_rI.second];
			}
		);
		s.m_s.back() = m_s.back() + _r.m_s.back();
		return s;
	}
	template<typename T1>
	auto operator-(const cjacobian<T1>&_r) const
	{	typedef typename merge<VECTOR, T1>::type MERGED;
		typedef typename createIndicies2<MERGED, VECTOR, T1>::type INDICIES;
		auto &r = createStdArray2<INDICIES>::type::value;
		cjacobian<MERGED> s;
		std::transform(
			r.begin(),
			r.end(),
			s.m_s.begin(),
			[&](const std::pair<std::size_t, std::size_t>&_rI)
			{	if (_rI.first == SIZE - 1)
					return -_r.m_s[_rI.second];
				else
				if (_rI.second == cjacobian<T1>::SIZE - 1)
					return m_s[_rI.first];
				else
					return m_s[_rI.first] - _r.m_s[_rI.second];
			}
		);
		s.m_s.back() = m_s.back() - _r.m_s.back();
		return s;
	}
	cjacobian operator+(const double _d) const
	{	cjacobian s(*this);
		s.m_s.back() += _d;
		return s;
	}
	cjacobian operator-(const double _d) const
	{	cjacobian s(*this);
		s.m_s.back() -= _d;
		return s;
	}
	cjacobian operator*(const double _d) const
	{	cjacobian s;
		for (std::size_t i = 0; i < SIZE; ++i)
			s.m_s[i] = m_s[i]*_d;
		return s;
	}
	cjacobian operator/(const double _d) const
	{	cjacobian s;
		const auto dInv = 1.0/_d;
		for (std::size_t i = 0; i < SIZE; ++i)
			s.m_s[i] = m_s[i]*dInv;
		return s;
	}
	cjacobian operator+(const cjacobian&_r) const
	{	cjacobian s;
		for (std::size_t i = 0; i < SIZE; ++i)
			s.m_s[i] = m_s[i] + _r.m_s[i];
		return s;
	}
	cjacobian operator-(const cjacobian&_r) const
	{	cjacobian s;
		for (std::size_t i = 0; i < SIZE; ++i)
			s.m_s[i] = m_s[i] - _r.m_s[i];
		return s;
	}
	cjacobian operator*(const cjacobian&_r) const
	{	cjacobian s;
		for (std::size_t i = 0; i < SIZE - 1; ++i)
			s.m_s[i] = m_s[i]*_r.m_s.back() + _r.m_s[i]*m_s.back();
		s.m_s.back() = m_s.back()*_r.m_s.back();
		return s;
	}
	cjacobian operator/(const cjacobian&_r) const
	{	cjacobian s;
		const auto dInv = 1.0/_r.m_s.back();
		const auto dInv2 = dInv*dInv*m_s.back();
		for (std::size_t i = 0; i < SIZE - 1; ++i)
			s.m_s[i] = dInv*m_s[i] - dInv2*_r.m_s[i];
		s.m_s.back() = m_s.back()*dInv;
		return s;
	}
	cjacobian operator-(void) const
	{	cjacobian s;
		for (std::size_t i = 0; i < SIZE; ++i)
			s.m_s[i] = -m_s[i];
		return s;
	}
	friend cjacobian operator+(const double _d, const cjacobian&_r)
	{	cjacobian s(_r);
		s.m_s.back() += _d;
		return s;
	}
	friend cjacobian operator-(const double _d, const cjacobian&_r)
	{	cjacobian s(-_r);
		s.m_s.back() += _d;
		return s;
	}
	friend cjacobian operator*(const double _d, const cjacobian&_r)
	{	cjacobian s;
		for (std::size_t i = 0; i < SIZE; ++i)
			s.m_s[i] = _d*_r.m_s[i];
		return s;
	}
	friend cjacobian operator/(const double _d, const cjacobian&_r)
	{	cjacobian s;
		const auto dInv = 1.0/_r.m_s.back();
		const auto dInv2 = dInv*dInv*_d;
		for (std::size_t i = 0; i < SIZE - 1; ++i)
			s.m_s[i] = -dInv2*_r.m_s[i];
		s.m_s.back() = _d*dInv;
		return s;
	}
	template<typename T1>
	bool operator<(const cjacobian<T1>&_r) const
	{	return m_s.back() < value(_r);
	}
	template<typename T1>
	bool operator>(const cjacobian<T1>&_r) const
	{	return m_s.back() > value(_r);
	}
	template<typename T1>
	bool operator<=(const cjacobian<T1>&_r) const
	{	return m_s.back() <= value(_r);
	}
	template<typename T1>
	bool operator>=(const cjacobian<T1>&_r) const
	{	return m_s.back() >= value(_r);
	}
	template<typename T1>
	bool operator==(const cjacobian<T1>&_r) const
	{	return m_s.back() == value(_r);
	}
	template<typename T1>
	bool operator!=(const cjacobian<T1>&_r) const
	{	return m_s.back() != value(_r);
	}

	bool operator<(const double _r) const
	{	return m_s.back() < _r;
	}
	bool operator>(const double _r) const
	{	return m_s.back() > _r;
	}
	bool operator<=(const double _r) const
	{	return m_s.back() <= _r;
	}
	bool operator>=(const double _r) const
	{	return m_s.back() >= _r;
	}
	bool operator==(const double _r) const
	{	return m_s.back() == _r;
	}
	bool operator!=(const double _r) const
	{	return m_s.back() != _r;
	}

	friend bool operator<(const double _d, const cjacobian&_r)
	{	return _d < value(_r);
	}
	friend bool operator>(const double _d, const cjacobian&_r)
	{	return _d > value(_r);
	}
	friend bool operator<=(const double _d, const cjacobian&_r)
	{	return _d <= value(_r);
	}
	friend bool operator>=(const double _d, const cjacobian&_r)
	{	return _d >= value(_r);
	}
	friend bool operator==(const double _d, const cjacobian&_r)
	{	return _d == value(_r);
	}
	friend bool operator!=(const double _d, const cjacobian&_r)
	{	return _d != value(_r);
	}

	friend std::ostream &operator<<(std::ostream&_rS, const cjacobian&_r)
	{	_rS << "(";
		for (std::size_t i = 0; i < SIZE; ++i)
			_rS << _r.m_s[i] << ",";
		_rS << ")";
		return _rS;
	}
	template<typename T1>
	friend auto pow(const cjacobian&_r0, const cjacobian<T1>&_r1)
	{	return exp(log(_r0)*_r1);
	}
	friend auto pow(const cjacobian&_r0, const double _d1)
	{	const auto d1 = std::pow(value(_r0), _d1 - 1.0);
		const auto d = d1*_d1;
		cjacobian s;
		for (std::size_t i = 0; i < SIZE - 1; ++i)
			s.m_s[i] = _r0.m_s[i]*d;
		s.m_s.back() = d1*value(_r0);
		return s;
	}
	friend auto pow(const double _d0, const cjacobian&_r1)
	{	return exp(log(_d0)*_r1);
	}
	template<std::pair<double, double>(*P)(const double)>
	static cjacobian nonlinear(const cjacobian&_r)
	{	cjacobian s;
		const auto sPair = P(value(_r));
		for (std::size_t i = 0; i < SIZE - 1; ++i)
			s.m_s[i] = sPair.second*_r.m_s[i];
		s.m_s.back() = sPair.first;
		return s;
	}
	static auto sin(const double _d)
	{	return std::make_pair(std::sin(_d), std::cos(_d));
	}
	static auto cos(const double _d)
	{	return std::make_pair(std::cos(_d), -std::sin(_d));
	}
	static auto tan(const double _d)
	{	const auto d = std::tan(_d);
		return std::make_pair(d, 1.0 + d*d);
	}
	static auto asin(const double _d)
	{	return std::make_pair(std::asin(_d), 1.0/std::sqrt(1.0 - _d*_d));
	}
	static auto acos(const double _d)
	{	return std::make_pair(std::acos(_d), -1.0/std::sqrt(1.0 - _d*_d));
	}
	static auto atan(const double _d)
	{	return std::make_pair(std::atan(_d), 1.0/(1.0 + _d*_d));
	}

	static auto sinh(const double _d)
	{	return std::make_pair(std::sinh(_d), std::cosh(_d));
	}
	static auto cosh(const double _d)
	{	return std::make_pair(std::cosh(_d), std::sinh(_d));
	}
	static auto tanh(const double _d)
	{	const auto d = std::tanh(_d);
		return std::make_pair(d, 1.0 - d*d);
	}
	static auto asinh(const double _d)
	{	return std::make_pair(std::asinh(_d), 1.0/std::sqrt(1.0 + _d*_d));
	}
	static auto acosh(const double _d)
	{	return std::make_pair(std::acosh(_d), 1.0/std::sqrt(_d*_d - 1));
	}
	static auto atanh(const double _d)
	{	return std::make_pair(std::atan(_d), 1.0/(1.0 - _d*_d));
	}

	static auto exp(const double _d)
	{	const auto d = std::exp(_d);
		return std::make_pair(d, d);
	}
	static auto log(const double _d)
	{	return std::make_pair(std::log(_d), 1.0/_d);
	}
	static const double s_dLog10;
	static auto log10(const double _d)
	{	return std::make_pair(std::log10(_d), 1.0/(s_dLog10*_d));
	}
	static auto sqrt(const double _d)
	{	const auto d = std::sqrt(_d);
		return std::make_pair(d, 0.5/d);
	}
	friend auto sqr(const cjacobian&_r)
	{	return _r*_r;
	}
#define __create__(sin)\
	friend cjacobian sin(const cjacobian&_r)\
	{	return nonlinear<sin>(_r);\
	}
	__create__(sqrt)
	__create__(exp)
	__create__(log)
	__create__(log10)
	__create__(sin)
	__create__(cos)
	__create__(tan)
	__create__(asin)
	__create__(acos)
	__create__(atan)
	__create__(sinh)
	__create__(cosh)
	__create__(tanh)
	__create__(asinh)
	__create__(acosh)
	__create__(atanh)
#undef __create__
};
template<typename VECTOR>
const double cjacobian<VECTOR>::s_dLog10 = std::log(10.0);
template<typename, typename>
struct common_type;
template<typename T>
struct common_type<cjacobian<T>, double>
{	typedef cjacobian<T> type;
};
template<typename T>
struct common_type<double, cjacobian<T> >
{	typedef cjacobian<T> type;
};
template<typename T0, typename T1>
struct common_type<cjacobian<T0>, cjacobian<T1> >
{	typedef cjacobian<typename merge<T0, T1>::type> type;
};
template<typename T0, typename T1, typename ...R0, typename ...R1>
struct common_type<
	std::tuple<T0, R0...>,
	std::tuple<T1, R1...>
>
{	typedef std::tuple<typename common_type<T0, T1>::type> FIRST;
	typedef typename common_type<std::tuple<R0...>, std::tuple<R1...> >::type REST;
	typedef decltype(std::tuple_cat(std::declval<FIRST>(), std::declval<REST>())) type;
};
template<>
struct common_type<
	std::tuple<>,
	std::tuple<>
>
{	typedef std::tuple<> type;
};
template<typename T, typename F>
typename common_type<
	typename std::decay<decltype(std::declval<T>()())>::type,
	typename std::decay<decltype(std::declval<F>()())>::type
>::type if_(
	const bool _b,
	T &&_rT,
	F&&_rF
)
{	if (_b)
		return std::forward<T>(_rT)();
	else
		return std::forward<F>(_rF)();
}
}
using implementation::cjacobian;
using implementation::if_;
}
