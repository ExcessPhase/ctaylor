/// published under MIT license
/// Author: Peter Foelsche
/// October-17th 2024
/// Austin, TX, USA
/// email:	peter_foelsche@outlook.com
/// A sparse, dual number implementation for calculating the 0th and 1th order of derivatives
/// Refer to VBIC95.cpp for a example usage (__JACOBIAN__ defined).
/// Requires C++14 and boost
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
#include <iterator>
#include <boost/iterator/permutation_iterator.hpp>
#include "merge_sorted_sets.h"
#include "initializer_list.h"
#include <boost/math/special_functions/polygamma.hpp>
namespace jacobian
{
namespace implementation
{
using namespace boost::mp11;


	/// returns a type which is sufficiently large to contain integer from 0...SIZE
	/// SIZE itself is used as an invalid ID and thus must also be in the range
template<typename SIZE>
struct getTypeFromSize
{	typedef typename std::conditional<
		(SIZE::value <= std::numeric_limits<unsigned int>::max()),
		typename std::conditional<
			(SIZE::value <= std::numeric_limits<unsigned short>::max()),
			typename std::conditional<
				(SIZE::value <= std::numeric_limits<unsigned char>::max()),
				unsigned char,
				unsigned short
			>::type,
			unsigned int
		>::type,
		std::size_t
	>::type type;
};
	/// create an initialized constexpr std::array
	/// needed to create the index_sequence used by createStdArrayImpl
template<typename VECTOR, typename SIZE>
struct createStdArray
{	typedef typename getTypeFromSize<SIZE>::type TYPE;
	typedef foelsche::init_list::convertToStdInitializerList<VECTOR, TYPE> type;
};
	/// used to call createStdArrayImpl2 and create the index_sequence
template<typename VECTOR, typename SIZE>
struct createStdArray2
{	typedef typename getTypeFromSize<SIZE>::type TYPE;
	typedef std::initializer_list<TYPE> IL;
	typedef foelsche::init_list::convertToStdInitializerList<
		mp_transform<mp_first, VECTOR>,
		TYPE
	> first;
	typedef foelsche::init_list::convertToStdInitializerList<
		mp_transform<mp_second, VECTOR>,
		TYPE
	> second;
	typedef std::pair<IL, IL> PAIR;
	static constexpr const PAIR value = {first::value, second::value};
};
template<typename VECTOR, typename SIZE>
constexpr const typename createStdArray2<VECTOR, SIZE>::PAIR createStdArray2<VECTOR, SIZE>::value;
	/// create a vector for every entry in TARGET by attempting to find the corresponding entry in SOURCE
	/// or size-of-SOURCE for if not found
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
	/// similar to createIndicies but with two sources
	/// accordingly the created vector is a vector of pairs
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
	/// for merging two meta vectors of IDs identifying independent variables
template<typename T0, typename T1>
struct merge
{	typedef typename taylor::merge_sorted_sets<
		mp_less,
		T0,
		T1
	>::type type;
};
	/// this is the dual number class
template<typename VECTOR>
struct cjacobian
{		/// SIZE of the vector is number of independent variables + 1 (0th derivative)
	static constexpr std::size_t SIZE = mp_size<VECTOR>::value + 1;
	typedef std::array<double, SIZE> ARRAY;
		/// here are all the derivatives stored. 0th derivative (value) at the end
	ARRAY m_s;
	cjacobian(void) = default;
	cjacobian(const cjacobian&) = default;
	cjacobian(cjacobian&&) = default;
	cjacobian&operator=(const cjacobian&) = default;
	cjacobian&operator=(cjacobian&&) = default;
		/// constructing from a plain value
		/// all derivatives are zero
	cjacobian(const double _d)
		:m_s({})
	{	m_s.back() = _d;
	}
		/// create an independent variable
		/// derivative is 1.0
	cjacobian(const double _d, const bool)
		:m_s({1.0, _d})
	{	static_assert(SIZE == 2, "size must be exactly 2!");
	}
		/// retrieve the 0th derivative (value)
	friend double value(const cjacobian&_r)
	{	return _r.m_s.back();
	}
		/// retrieve the derivative vs. the independent variable with index I
		/// does not compile if there is no such derivative
	template<std::size_t I>
	double getDer(const mp_size_t<I>&) const
	{	static_assert(mp_find<VECTOR, mp_size_t<I> >::value + 1 < SIZE, "derivative pattern not found!");
		return m_s[mp_find<VECTOR, mp_size_t<I> >::value];
	}
		/// a conversion constructor
		/// it is impossible to forget derivatives
	template<typename T1>
	cjacobian(const cjacobian<T1>&_r)
		:m_s(convert<T1>(_r.m_s))
	{
	}
		/// for converting between different cjacobian types
		/// target must contain all the source content and more
	template<typename T1>
	static ARRAY convert(const typename cjacobian<T1>::ARRAY&_r)
	{	ARRAY s;
		typedef typename createIndicies<VECTOR, T1>::type INDICIES;
		typedef mp_plus<mp_size<T1>, mp_size_t<1> > SIZE_T;
		auto &r = createStdArray<INDICIES, SIZE_T>::type::value;
		typedef typename getTypeFromSize<SIZE_T>::type TYPE;
		std::transform(
			r.begin(),
			r.end(),
			s.begin(),
			[&](const TYPE _i)
			{	if (_i == cjacobian<T1>::SIZE - 1)
					return 0.0;
				else
					return _r[_i];
			}
		);
		s.back() = _r.back();
		return s;
	}
		/// a conversion constructor
		/// impossible to forget derivatives
	template<typename T1>
	cjacobian&operator=(const cjacobian<T1>&_r)
	{	m_s = convert<T1>(_r.m_s);
		return *this;
	}
		/// assignment operator from a double
		/// all derivatives will be zero
	cjacobian&operator=(const double _r)
	{	std::fill(m_s.begin(), std::prev(m_s.end()), 0.0);
		m_s.back() = _r;
		return *this;
	}
		/// create a new independent variable for chainrule to reduce the number of carried derivatives
	template<std::size_t ENUM>
	cjacobian<mp_list<mp_size_t<ENUM> > > convert2Independent(const mp_size_t<ENUM>&) const
	{	return {value(*this), true};
	}
		/// substitutes one derivative by the ones passed in the first argument
		/// might have to be called multiple times
		/// the first argument must have been one on which convert2Independent() was called.
		/// ENUM must be identical to the ENUM passed to convert2Independent()
	template<typename T1, std::size_t ENUM>
	auto chainRule(const cjacobian<T1>&_r, const mp_size_t<ENUM>&) const
	{	return chainRule(_r, mp_size_t<ENUM>(), mp_set_contains<VECTOR, mp_size_t<ENUM> >());
	}
	template<typename T1, std::size_t ENUM>
	auto chainRule(const cjacobian<T1>&, const mp_size_t<ENUM>&, const mp_false&) const
	{	return *this;
	}
	template<typename T1, std::size_t ENUM>
	auto chainRule(const cjacobian<T1>&_r, const mp_size_t<ENUM>&, const mp_true&) const
	{	typedef mp_find<VECTOR, mp_size_t<ENUM> > START;
		typedef mp_size_t<START::value + 1> END;
		typedef mp_erase<VECTOR, START, END> REMOVED;
		typedef cjacobian<REMOVED> REMOVEDJ;
		REMOVEDJ s;
		std::copy(
			m_s.cbegin(),
			m_s.cbegin() + START::value,
			s.m_s.begin()
		);
		std::copy(
			m_s.cbegin() + START::value + 1,
			m_s.cend(),
			s.m_s.begin() + START::value
		);
		return s + getDer(mp_size_t<ENUM>())*(_r - value(_r));
	}
	template<typename T1>
	cjacobian&operator+=(const cjacobian<T1>&_r)
	{	typedef typename createIndicies<T1, VECTOR>::type INDICIES;
		auto &r = createStdArray<INDICIES, mp_size<VECTOR> >::type::value;
		typedef typename getTypeFromSize<mp_size<INDICIES> >::type TYPE;
		std::transform(
			_r.m_s.begin(),
			std::prev(_r.m_s.end()),
			boost::make_permutation_iterator(m_s.cbegin(), r.begin()),
			boost::make_permutation_iterator(m_s.begin(), r.begin()),
			[](const double _d0, const double _d1)
			{	return _d1 + _d0;
			}
		);
		m_s.back() += _r.m_s.back();
		return *this;
	}
	template<typename T1>
	cjacobian&operator-=(const cjacobian<T1>&_r)
	{	typedef typename createIndicies<T1, VECTOR>::type INDICIES;
		auto &r = createStdArray<INDICIES, mp_size<VECTOR> >::type::value;
		typedef typename getTypeFromSize<mp_size<INDICIES> >::type TYPE;
		std::transform(
			_r.m_s.begin(),
			std::prev(_r.m_s.end()),
			boost::make_permutation_iterator(m_s.cbegin(), r.begin()),
			boost::make_permutation_iterator(m_s.begin(), r.begin()),
			[](const double _d0, const double _d1)
			{	return _d1 - _d0;
			}
		);
		m_s.back() -= _r.m_s.back();
		return *this;
	}
		/// LHS and RHS are the same type
	cjacobian&operator+=(const cjacobian &_r)
	{	for (std::size_t i = 0; i < SIZE; ++i)
			m_s[i] += _r.m_s[i];
		return *this;
	}
		/// LHS and RHS are the same type
	cjacobian&operator-=(const cjacobian &_r)
	{	for (std::size_t i = 0; i < SIZE; ++i)
			m_s[i] -= _r.m_s[i];
		return *this;
	}
	cjacobian&operator*=(const double _d)
	{	for (std::size_t i = 0; i < SIZE; ++i)
			m_s[i] *= _d;
		return *this;
	}
	cjacobian&operator/=(const double _d)
	{	const auto d1 = 1.0/_d;
		for (std::size_t i = 0; i < SIZE; ++i)
			m_s[i] *= d1;
		return *this;
	}
	cjacobian&operator+=(const double _d)
	{	m_s.back() += _d;
		return *this;
	}
	cjacobian&operator-=(const double _d)
	{	m_s.back() -= _d;
		return *this;
	}
	template<typename T1>
	auto operator*(const cjacobian<T1>&_r) const
	{	typedef typename merge<VECTOR, T1>::type MERGED;
		typedef typename createIndicies2<MERGED, VECTOR, T1>::type INDICIES;
		typedef mp_plus<
			mp_max<
				mp_size<VECTOR>,
				mp_size<T1>
			>,
			mp_size_t<1>
		> SIZE_T;
		auto &r = createStdArray2<INDICIES, SIZE_T>::value;
		typedef typename getTypeFromSize<SIZE_T>::type TYPE;
		cjacobian<MERGED> s;
		std::transform(
			r.first.begin(),
			r.first.end(),
			r.second.begin(),
			s.m_s.begin(),
			[&](const TYPE _i0, const TYPE _i1)
			{	if (_i0 == SIZE - 1)
					return _r.m_s[_i1]*m_s.back();
				else
				if (_i1 == cjacobian<T1>::SIZE - 1)
					return m_s[_i0]*_r.m_s.back();
				else
					return m_s[_i0]*_r.m_s.back() + _r.m_s[_i1]*m_s.back();
			}
		);
		s.m_s.back() = m_s.back() * _r.m_s.back();
		return s;
	}
	template<typename T1>
	auto operator/(const cjacobian<T1>&_r) const
	{	typedef typename merge<VECTOR, T1>::type MERGED;
		typedef typename createIndicies2<MERGED, VECTOR, T1>::type INDICIES;
		typedef mp_plus<
			mp_max<
				mp_size<VECTOR>,
				mp_size<T1>
			>,
			mp_size_t<1>
		> SIZE_T;
		auto &r = createStdArray2<INDICIES, SIZE_T>::value;
		typedef typename getTypeFromSize<SIZE_T>::type TYPE;
		cjacobian<MERGED> s;
		const auto dInv = 1.0/_r.m_s.back();
		const auto dInv2 = dInv*dInv*m_s.back();
		std::transform(
			r.first.begin(),
			r.first.end(),
			r.second.begin(),
			s.m_s.begin(),
			[&](const TYPE _i0, const TYPE _i1)
			{	if (_i0 == SIZE - 1)
					return -dInv2*_r.m_s[_i1];
				else
				if (_i1 == cjacobian<T1>::SIZE - 1)
					return dInv*m_s[_i0];
				else
					return dInv*m_s[_i0] - dInv2*_r.m_s[_i1];
			}
		);
		s.m_s.back() = m_s.back() * dInv;
		return s;
	}
	template<typename T1>
	auto operator+(const cjacobian<T1>&_r) const
	{	typedef typename merge<VECTOR, T1>::type MERGED;
		typedef typename createIndicies2<MERGED, VECTOR, T1>::type INDICIES;
		typedef mp_plus<
			mp_max<
				mp_size<VECTOR>,
				mp_size<T1>
			>,
			mp_size_t<1>
		> SIZE_T;
		auto &r = createStdArray2<INDICIES, SIZE_T>::value;
		typedef typename getTypeFromSize<SIZE_T>::type TYPE;
		cjacobian<MERGED> s;
		std::transform(
			r.first.begin(),
			r.first.end(),
			r.second.begin(),
			s.m_s.begin(),
			[&](const TYPE _i0, const TYPE _i1)
			{	if (_i0 == SIZE - 1)
					return _r.m_s[_i1];
				else
				if (_i1 == cjacobian<T1>::SIZE - 1)
					return m_s[_i0];
				else
					return m_s[_i0] + _r.m_s[_i1];
			}
		);
		s.m_s.back() = m_s.back() + _r.m_s.back();
		return s;
	}
	template<typename T1>
	auto operator-(const cjacobian<T1>&_r) const
	{	typedef typename merge<VECTOR, T1>::type MERGED;
		typedef typename createIndicies2<MERGED, VECTOR, T1>::type INDICIES;
		typedef mp_plus<
			mp_max<
				mp_size<VECTOR>,
				mp_size<T1>
			>,
			mp_size_t<1>
		> SIZE_T;
		auto &r = createStdArray2<INDICIES, SIZE_T>::value;
		typedef typename getTypeFromSize<SIZE_T>::type TYPE;
		cjacobian<MERGED> s;
		std::transform(
			r.first.begin(),
			r.first.end(),
			r.second.begin(),
			s.m_s.begin(),
			[&](const TYPE _i0, const TYPE _i1)
			{	if (_i0 == SIZE - 1)
					return -_r.m_s[_i1];
				else
				if (_i1 == cjacobian<T1>::SIZE - 1)
					return m_s[_i0];
				else
					return m_s[_i0] - _r.m_s[_i1];
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
	typedef std::pair<double, double> doublePair;
	template<doublePair(*P)(const double)>
	static cjacobian nonlinear(const cjacobian&_r)
	{	cjacobian s;
		const auto sPair = P(value(_r));
		for (std::size_t i = 0; i < SIZE - 1; ++i)
			s.m_s[i] = sPair.second*_r.m_s[i];
		s.m_s.back() = sPair.first;
		return s;
	}
	static doublePair sin_(const double _d)
	{	return std::make_pair(std::sin(_d), std::cos(_d));
	}
	static doublePair cos_(const double _d)
	{	return std::make_pair(std::cos(_d), -std::sin(_d));
	}
	static doublePair tan_(const double _d)
	{	const auto d = std::tan(_d);
		return std::make_pair(d, 1.0 + d*d);
	}
	static doublePair asin_(const double _d)
	{	return std::make_pair(std::asin(_d), 1.0/std::sqrt(1.0 - _d*_d));
	}
	static doublePair acos_(const double _d)
	{	return std::make_pair(std::acos(_d), -1.0/std::sqrt(1.0 - _d*_d));
	}
	static doublePair atan_(const double _d)
	{	return std::make_pair(std::atan(_d), 1.0/(1.0 + _d*_d));
	}
	static doublePair sinh_(const double _d)
	{	return std::make_pair(std::sinh(_d), std::cosh(_d));
	}
	static doublePair cosh_(const double _d)
	{	return std::make_pair(std::cosh(_d), std::sinh(_d));
	}
	static doublePair tanh_(const double _d)
	{	const auto d = std::tanh(_d);
		return std::make_pair(d, 1.0 - d*d);
	}
	static doublePair asinh_(const double _d)
	{	return std::make_pair(std::asinh(_d), 1.0/std::sqrt(1.0 + _d*_d));
	}
	static doublePair acosh_(const double _d)
	{	return std::make_pair(std::acosh(_d), 1.0/std::sqrt(_d*_d - 1));
	}
	static doublePair atanh_(const double _d)
	{	return std::make_pair(std::atanh(_d), 1.0/(1.0 - _d*_d));
	}
	static doublePair exp_(const double _d)
	{	const auto d = std::exp(_d);
		return std::make_pair(d, d);
	}
	static doublePair log_(const double _d)
	{	return std::make_pair(std::log(_d), 1.0/_d);
	}
	static doublePair tgamma_(const double _d)
	{	const auto d = std::tgamma(_d);
		return std::make_pair(d, d*boost::math::polygamma(0, _d));
	}
#if defined(__GNUC__) && !defined(__clang__)
	static constexpr double s_dLog10 = std::log(10.0);
#else
	static const double s_dLog10;
#endif
	static doublePair log10_(const double _d)
	{	return std::make_pair(std::log10(_d), 1.0/(s_dLog10*_d));
	}
	static doublePair sqrt_(const double _d)
	{	const auto d = std::sqrt(_d);
		return std::make_pair(d, 0.5/d);
	}
	static doublePair cbrt_(const double _d)
	{	const double d0 = std::cbrt(_d);
		return std::make_pair(
			d0,
			1.0/(3.0*d0*d0)
		);
	}
	friend auto sqr(const cjacobian&_r)
	{	return _r*_r;
	}
	friend bool isnan(const cjacobian&_r)
	{	return std::any_of(
			_r.m_s.begin(),
			_r.m_s.end(),
			static_cast<bool(*)(double)>(&std::isnan)
		);
	}
	friend bool isfinite(const cjacobian&_r)
	{	return std::all_of(
			_r.m_s.begin(),
			_r.m_s.end(),
			static_cast<bool(*)(double)>(&std::isfinite)
		);
	}
	template<typename T0, typename T1>
	static auto atan2_(const T0 &_rY, const T1&_rX) -> decltype(atan(_rY/_rX))
	{	if (_rX > 0.0)
			return atan(_rY/_rX);
		else
		if (_rX < 0.0 && _rY >= 0.0)
			return atan(_rY/_rX) + M_PI;
		else
		if (_rX < 0.0 && _rY < 0.0)
			return atan(_rY/_rX) - M_PI;
		else
		if (_rX == 0.0 && _rY > 0.0)
			return M_PI/2.0;
		else
		if (_rX == 0.0 && _rY < 0.0)
			return -M_PI/2.0;
		else
			return atan(_rY/_rX);
	}
	template<typename T1>
	friend auto atan2(
		const cjacobian&_rY,
		const cjacobian<T1>&_rX
	)
	{	return atan2_(_rY, _rX);
	}
	friend auto atan2(const cjacobian&_rY, const double _dX)
	{	return atan2_(_rY, _dX);
	}
	friend auto atan2(const double _dY, const cjacobian&_rX)
	{	return atan2_(_dY, _rX);
	}
	friend auto max(const cjacobian&_r0, const double _r1)
	{	return if_(
			_r0 > _r1,
			[&](void)
			{	return _r0;
			},
			[&](void)
			{	return _r1;
			}
		);
	}
	friend auto max(const double _r0, const cjacobian&_r1)
	{	return if_(
			_r0 > _r1,
			[&](void)
			{	return _r0;
			},
			[&](void)
			{	return _r1;
			}
		);
	}
		/// if this does not compile, check the order of ENUM-order pairs
	template<typename T1>
	friend auto max(const cjacobian&_r0, const cjacobian<T1>&_r1)
	{	return if_(
			_r0 > _r1,
			[&](void)
			{	return _r0;
			},
			[&](void)
			{	return _r1;
			}
		);
	}
	friend auto min(const cjacobian&_r0, const double _r1)
	{	return if_(
			_r0 < _r1,
			[&](void)
			{	return _r0;
			},
			[&](void)
			{	return _r1;
			}
		);
	}
	friend auto min(const double _r0, const cjacobian&_r1)
	{	return if_(
			_r0 < _r1,
			[&](void)
			{	return _r0;
			},
			[&](void)
			{	return _r1;
			}
		);
	}
		/// if this does not compile, check the order of ENUM-order pairs
	template<typename T1>
	friend auto min(const cjacobian&_r0, const cjacobian<T1>&_r1)
	{	return if_(
			_r0 < _r1,
			[&](void)
			{	return _r0;
			},
			[&](void)
			{	return _r1;
			}
		);
	}
	friend auto fmod(const cjacobian&_r0, const double _r1)
	{	return _r0 - static_cast<int>(value(_r0)/_r1)*_r1;
	}
	friend auto fmod(const double _r0, const cjacobian&_r1)
	{	return _r0 - static_cast<int>(_r0/value(_r1))*_r1;
	}
	template<typename T1>
	friend auto fmod(const cjacobian&_r0, const cjacobian<T1>&_r1)
	{	return _r0 - static_cast<int>(value(_r0)/value(_r1))*_r1;
	}
	static doublePair erf_(const double _d)
	{	return std::make_pair(
			std::erf(_d),
			s_dTwoOverSqrtPi*std::exp(-_d*_d)
		);
	}
	static doublePair erfc_(const double _d)
	{	return std::make_pair(
			std::erfc(_d),
			-s_dTwoOverSqrtPi*std::exp(-_d*_d)
		);
	}
	static doublePair polygamma_(const int _i, const double _d)
	{	return std::make_pair(
			boost::math::polygamma(_i, _d),
			boost::math::polygamma(_i + 1, _d)
		);
	}
	friend cjacobian polygamma(const int _i, const cjacobian&_r)\
	{	cjacobian s;
		const auto sPair = polygamma_(_i, value(_r));
		for (std::size_t i = 0; i < SIZE - 1; ++i)
			s.m_s[i] = sPair.second*_r.m_s[i];
		s.m_s.back() = sPair.first;
		return s;
	}
#define __create__(sin)\
	friend cjacobian sin(const cjacobian&_r)\
	{	return nonlinear<sin##_>(_r);\
	}
	__create__(tgamma)
	__create__(cbrt)
	__create__(erfc)
	__create__(erf)
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
	template<typename T1>
	friend auto hypot(
		const cjacobian&_r0,
		const cjacobian<T1>&_r1
	)
	{	return sqrt(sqr(_r0) + sqr(_r1));
	}
	friend auto hypot(const cjacobian&_rX, const double _dY)
	{	return sqrt(sqr(_rX) + _dY*_dY);
	}
	friend auto hypot(const double _dX, const cjacobian&_rY)
	{	return sqrt(sqr(_rY) + _dX*_dX);
	}
#if defined(__GNUC__) && !defined(__clang__)
	static constexpr double s_dTwoOverSqrtPi = 2.0/std::sqrt(M_PI);
#else
	static const double s_dTwoOverSqrtPi;
#endif
};
#if !defined(__GNUC__) || defined(__clang__)
template<typename T>
const double cjacobian<T>::s_dTwoOverSqrtPi = 2.0/std::sqrt(M_PI);
template<typename VECTOR>
const double cjacobian<VECTOR>::s_dLog10 = std::log(10.0);
#endif
template<typename A, typename B>
struct common_type
{	typedef typename std::common_type<A, B>::type type;
};
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
