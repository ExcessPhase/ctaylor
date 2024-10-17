#pragma once
#include <iostream>
#include <type_traits>
#include <array>
#include <limits>
#include <functional>
#include <boost/mp11.hpp>
#include <cstdlib>
#include <algorithm>
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
	double getValue(void) const
	{	return m_s.back();
	}
	template<typename T1>
	cjacobian(const cjacobian<T1>&_r)
		:m_s(convert<T1>(_r.m_s))
	{
	}
	template<typename T1>
	static ARRAY convert(const typename cjacobian<T1>::ARRAY&_r)
	{	ARRAY s;
//template<typename TARGET, typename SOURCE>
//struct createIndicies
		typedef typename createIndicies<VECTOR, T1>::type INDICIES;
//template<typename VECTOR>
//struct createStdArray
		auto &r = createStdArray<INDICIES>::type::value;
		std::transform(
			r.begin(),
			r.end(),
			s.begin(),
			[&](const std::size_t _i)
			{	if (_i == cjacobian<T1>::SIZE)
					return 0.0;
				else
					return _r[_i];
			}
		);
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
	friend std::ostream &operator<<(std::ostream&_rS, const cjacobian&_r)
	{	_rS << "(";
		for (std::size_t i = 0; i < SIZE; ++i)
			_rS << _r.m_s[i] << ",";
		_rS << ")";
		return _rS;
	}
};
}
using implementation::cjacobian;
}
