/// published under MIT license
/// Author: Peter Foelsche
/// October-1th 2024
/// Austin, TX, USA
/// email:	peter_foelsche@outlook.com
/// A sparse, dual number implementation for calculating not just the 1th order of derivatives
/// Refer to ctaylor.cpp for a example usage
/// Compile time under Visual C++ 2022 tends to be much longer than using g++
/// Requires C++14
#pragma once
#include <iostream>
#include <type_traits>
#include <array>
#include <limits>
#include <functional>
#include <boost/mp11.hpp>
#include <cstdlib>
#include <algorithm>
#include <numeric>
#include "merge_sorted_sets.h"
#include "taylor_series_expansions.h"

namespace taylor
{
namespace implementation
{
template<typename SIZE>
struct getTypeFromSize;
	/// meta function for merging different result types
template<typename, typename>
struct common_type;
	/// function for merginig different result tyoes
template<typename T, typename F>
typename common_type<
	typename std::decay<decltype(std::declval<T>()())>::type,
	typename std::decay<decltype(std::declval<F>()())>::type
>::type if_(
	const bool _b,
	T &&_rT,
	F&&_rF
);
using namespace boost::mp11;
	/// to be passed to mp_for_each with a vector/set argument to ctaylor
	/// for debgging purposes
struct output
{	std::ostream&m_r;
	output(std::ostream&_r)
		:m_r(_r)
	{
	}
	template<std::size_t I>
	void operator()(const mp_size_t<I>&) const
	{	m_r << I << ",";
	}
	void operator()(const mp_list<>&) const
	{	m_r << "(), ";
	}
	void operator()(const mp_list<>&, const mp_true&) const
	{
	}
	template<typename FIRST, typename ...REST>
	void operator()(const mp_list<FIRST, REST...>&, const mp_true&) const
	{	(*this)(FIRST());
		(*this)(mp_list<REST...>(), mp_true());
	}
	template<typename FIRST, typename ...REST>
	void operator()(const mp_list<FIRST, REST...>&) const
	{	m_r << "(";
		(*this)(FIRST());
		(*this)(mp_list<REST...>(), mp_true());
		m_r << ")";
	}
};
template<typename ...REST>
std::ostream&operator<<(std::ostream&_r, const mp_list<REST...>&)
{	_r << "(";
	mp_for_each<
		mp_list<REST...>
	>(output(_r));
	_r << ")";
	return _r;
}
template<std::size_t I>
struct factorial
{	static constexpr double value = I*factorial<I - 1>::value;
};
template<>
struct factorial<0>
{	static constexpr double value = 1.0;
};
template<typename>
struct accumulatedFactorial;
template<>
struct accumulatedFactorial<mp_list<> >
{	static constexpr double value = 1.0;
};
template<std::size_t ENUM, std::size_t ORDER>
struct accumulatedFactorial<
	mp_list<
		pair<
			mp_size_t<ENUM>,
			mp_size_t<ORDER>
		>
	>
>
{	static constexpr double value = factorial<ORDER>::value;
};
template<typename FIRST, typename ...REST>
struct accumulatedFactorial<
	mp_list<
		FIRST,
		REST...
	>
>
{	static constexpr double value = accumulatedFactorial<mp_list<FIRST> >::value*accumulatedFactorial<mp_list<REST...> >::value;
};
template<typename T>
struct containsValue;
template<typename T>
struct containsValue2;
/// creates a type ready for passing to ctaylor
/// with the value as the first element
/// and the 1th order derivative as the second element
/// ENUM indicating the independent variable
template<std::size_t ENUM>
using makeIndependent = mp_list<
	mp_list<>,
	mp_list<
		pair<
			mp_size_t<ENUM>,
			mp_size_t<1>	/// the order
		>
	>
>;
/// for determining the order -- accumulating the second part of the pair
template<typename SUM, typename PAIR>
using add_second=std::integral_constant<
	std::size_t,
	SUM::value + PAIR::second_type::value
>;
	/// for determining the order
template<typename LIST>
struct order
{	static constexpr const auto value = mp_fold<LIST, mp_size_t<0>, add_second>::value;
};
template<typename T0, typename T1>
struct lexicographical_compare;
template<typename T>
struct lexicographical_compare<T, T>
{	typedef mp_false type;
};
template<>
struct lexicographical_compare<mp_list<>, mp_list<> >
{	typedef mp_false type;
};
template<typename ...T>
struct lexicographical_compare<mp_list<T...>, mp_list<> >
{	typedef mp_false type;
};
template<typename ...T>
struct lexicographical_compare<mp_list<>, mp_list<T...> >
{	typedef mp_true type;
};
/// doing the comparison starting from the rear
/// makes certain, that multiplying a taylor polynomial
/// with an element of anothher polynomial yields an already sorted polynomial
template<typename ...R0, typename ...R1>
struct lexicographical_compare<mp_list<R0...>, mp_list<R1...> >
{	typedef typename std::conditional<
		(mp_back<mp_list<R0...> >::first_type::value < mp_back<mp_list<R1...> >::first_type::value),
		mp_identity<mp_true>,
		typename std::conditional<
			(mp_back<mp_list<R0...> >::first_type::value > mp_back<mp_list<R1...> >::first_type::value),
			mp_identity<mp_false>,
			typename std::conditional<
				(mp_back<mp_list<R0...> >::second_type::value < mp_back<mp_list<R1...> >::second_type::value),
				mp_identity<mp_true>,
				typename std::conditional<
					(mp_back<mp_list<R0...> >::second_type::value > mp_back<mp_list<R1...> >::second_type::value),
					mp_identity<mp_false>,
					lexicographical_compare<mp_pop_back<mp_list<R0...> >, mp_pop_back<mp_list<R1...> > >
				>::type
			>::type
		>::type
	>::type::type type;
};
template<typename T0, typename T1>
struct compareListOfPairs
{	typedef typename std::conditional<
		(order<T0>::value < order<T1>::value),
		mp_identity<mp_true>,
		typename std::conditional<
			(order<T0>::value > order<T1>::value),
			mp_identity<mp_false>,
			lexicographical_compare<T0, T1>
		>::type
	>::type::type type;
};
template<typename T>
struct compareListOfPairs<T, T>
{	typedef mp_false type;
};
template<typename T0, typename T1>
struct compareListOfPairs2
{	typedef typename compareListOfPairs<typename T0::first_type, typename T1::first_type>::type type;
};
#if 0
template<typename LIST>
struct listOfListsIsSorted;

template<>
struct listOfListsIsSorted<mp_list<> >
{	typedef mp_true type;
};
template<typename T>
struct listOfListsIsSorted<mp_list<T> >
{	typedef mp_true type;
};
template<typename T0, typename T1, typename ...REST>
struct listOfListsIsSorted<mp_list<T0, T1, REST...> >
{	typedef mp_and<
		typename compareListOfPairs<T0, T1>::type,
		typename listOfListsIsSorted<mp_pop_front<mp_list<T0, T1, REST...> > >::type
	> type;
};
#endif
/// find positions of elements in SOURCE in TARET
template<typename TARGET, typename SOURCE, typename SIZE, bool CHECK=true>
struct findPositions
{
#ifndef NDEBUG
	static_assert(!CHECK || mp_size<TARGET>::value >= mp_size<SOURCE>::value, "size of target must be larger than size of source!");
	static_assert(mp_is_set<TARGET>::value, "TARGET must be a set!");
	static_assert(mp_is_set<SOURCE>::value, "SOURCE must be a set!");
	static_assert(!CHECK || std::is_same<TARGET, mp_set_union<TARGET, SOURCE> >::value, "TARGET must contain all elements in SOURCE");
#endif
	typedef typename getTypeFromSize<SIZE>::type TYPE;
#if 0
	template<typename STATE, typename SOURCE_ELEMENT>
	using checkPosition = mp_list<
		typename std::conditional<
			(mp_first<STATE>::value == std::numeric_limits<TYPE>::max()),
			std::conditional<
				std::is_same<
					SOURCE_ELEMENT,
					mp_third<STATE>
				>::value,
				mp_second<STATE>,
				mp_first<STATE>
			>,
			mp_identity<mp_first<STATE> >
		>::type::type,
		mp_size_t<mp_second<STATE>::value + 1>,
		mp_third<STATE>
	>;

	template<typename STATE, typename TARGET_ELEMENT>
	using findPosition = mp_push_back<
		STATE,
		mp_list<
			mp_size<STATE>,
			mp_first<
				mp_fold<
					SOURCE,
					mp_list<
						mp_size_t<std::numeric_limits<TYPE>::max()>, // the result
						mp_size_t<0>,	// the next position
						TARGET_ELEMENT
					>,
					checkPosition
				>
			>
		>
	>;
	typedef mp_fold<
		TARGET,
		mp_list<>,
		findPosition
	> type;
#else
	template<typename RESULT, typename T>
	using findElement = mp_push_back<
		RESULT,
		mp_list<
			mp_size<RESULT>,
			typename std::conditional<
				(mp_find<SOURCE, T>::value == mp_size<SOURCE>::value),
				mp_size_t<std::numeric_limits<TYPE>::max()>,
				mp_find<SOURCE, T>
			>::type
		>
	>;
	typedef mp_fold<
		TARGET,
		mp_list<>,
		findElement
	> type;
#endif
};
template<typename L0, typename L1>
using combine = mp_list<
	mp_first<L0>,
	mp_second<L0>,
	mp_second<L1>
>;
/// find positions of elements in SOURCE in TARET
template<typename TARGET, typename SOURCE0, typename SOURCE1>
struct findPositions2
{
#ifndef NDEBUG
	static_assert((mp_size<TARGET>::value >= mp_size<SOURCE0>::value), "size of target must be larger than size of source!");
	static_assert(mp_is_set<TARGET>::value, "TARGET must be a set!");
	static_assert(mp_is_set<SOURCE0>::value, "SOURCE must be a set!");
	static_assert(mp_is_set<SOURCE1>::value, "SOURCE must be a set!");
	static_assert(std::is_same<TARGET, mp_set_union<TARGET, SOURCE0, SOURCE1> >::value, "TARGET must contain all elements in SOURCE");
#endif
	typedef mp_plus<mp_max<mp_size<SOURCE0>, mp_size<SOURCE1> >, mp_size_t<1> > SIZE;
	typedef mp_transform<
		combine,
		typename findPositions<TARGET, SOURCE0, SIZE>::type,
		typename findPositions<TARGET, SOURCE1, SIZE>::type
	> type;
};
template<typename A, typename B>
struct combineTwoPairs;
template<typename A, typename ...B, typename ...C>
struct combineTwoPairs<
	pair<A, mp_list<B...> >,
	pair<A, mp_list<C...> >
>
{	typedef pair<
		A,
		mp_append<
			mp_list<B...>,
			mp_list<C...>
		>
	> type;
};
/// new merge<>
/// COMPARE= compareListOfPairs2
/// MERGE=combineTwoPairs
/// CONTAINS_VALUE=containsValue2
/// merge two sets of list_of_list
template<
	typename T0,
	typename T1,
	template<typename, typename> class COMPARE=compareListOfPairs,
	template<typename, typename> class MERGE=combineTwo,
	template<typename> class CONTAINS_VALUE=containsValue
>
struct merge
{
#ifndef NDEBUG
	static_assert(mp_is_set<T0>::value, "must be a set!");
	static_assert(mp_is_set<T1>::value, "must be a set!");
#endif
	typedef typename merge_sorted_sets<
		COMPARE,
		T0,
		T1,
		MERGE
	>::type type;
#ifndef NDEBUG
	static_assert(
		CONTAINS_VALUE<type>::type::value == mp_or<
			typename CONTAINS_VALUE<T0>::type,
			typename CONTAINS_VALUE<T1>::type
		>::value, "value in merge result!");
#endif
};
template<
	typename T,
	template<typename, typename> class COMPARE,
	template<typename, typename> class MERGE,
	template<typename> class CONTAINS_VALUE
>
struct merge<T, mp_list<>, COMPARE, MERGE, CONTAINS_VALUE>
{	static_assert(mp_is_set<T>::value, "must be a set!");
	typedef T type;
};
template<
	typename T,
	template<typename, typename> class COMPARE,
	template<typename, typename> class MERGE,
	template<typename> class CONTAINS_VALUE
>
struct merge<mp_list<>, T, COMPARE, MERGE, CONTAINS_VALUE>
{	static_assert(mp_is_set<T>::value, "must be a set!");
	typedef T type;
};
template<
	typename T,
	template<typename, typename> class COMPARE,
	template<typename, typename> class MERGE,
	template<typename> class CONTAINS_VALUE
>
struct merge<T, T, COMPARE, MERGE, CONTAINS_VALUE>
{	typedef T type;
};
template<
	template<typename, typename> class COMPARE,
	template<typename, typename> class MERGE,
	template<typename> class CONTAINS_VALUE
>
struct merge<mp_list<>, mp_list<>, COMPARE, MERGE, CONTAINS_VALUE>
{	typedef mp_list<> type;
};

template<typename, typename, typename>
struct convertToStdInitializerListImpl;
template<typename LIST, std::size_t ...INDICES, typename TYPE>
struct convertToStdInitializerListImpl<LIST, std::index_sequence<INDICES...>, TYPE>
{	static constexpr const std::initializer_list<std::pair<TYPE, TYPE> > value =
	{	std::make_pair(
			TYPE(mp_first<mp_at_c<LIST, INDICES> >::value),
			TYPE(mp_second<mp_at_c<LIST, INDICES> >::value)
		)...
	};
};
template<typename LIST, std::size_t ...INDICES, typename TYPE>
constexpr const std::initializer_list<std::pair<TYPE, TYPE> >
convertToStdInitializerListImpl<LIST, std::index_sequence<INDICES...>, TYPE>::value;
template<typename LIST, typename TYPE>
struct convertToStdInitializerList
{	typedef convertToStdInitializerListImpl<LIST, std::make_index_sequence<mp_size<LIST>::value>, TYPE> type;
};
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
template<typename, typename, typename>
struct convertToStdArray3Impl;
template<typename LIST, std::size_t ...INDICES, typename SIZE>
struct convertToStdArray3Impl<LIST, std::index_sequence<INDICES...>, SIZE>
{	typedef typename getTypeFromSize<SIZE>::type TYPE;
	static constexpr const std::array<
		std::initializer_list<std::pair<TYPE, TYPE> >,
		mp_size<LIST>::value
	> value = {
		convertToStdInitializerList<mp_at_c<LIST, INDICES>, TYPE>::type::value...
	};
};
template<typename LIST, std::size_t ...INDICES, typename SIZE>
constexpr const std::array<
	std::initializer_list<std::pair<typename getTypeFromSize<SIZE>::type, typename getTypeFromSize<SIZE>::type> >,
	mp_size<LIST>::value
>
convertToStdArray3Impl<LIST, std::index_sequence<INDICES...>, SIZE>::value;
template<typename LIST, typename SIZE>
struct convertToStdArray3
{	typedef convertToStdArray3Impl<
		LIST,
		std::make_index_sequence<mp_size<LIST>::value>,
		SIZE
	> type;
};
	/// convert meta ARRAY into std::array
template<typename, typename, typename>
struct convertToStdArrayImpl;
template<typename LIST, std::size_t ...INDICES, typename SIZE>
struct convertToStdArrayImpl<LIST, std::index_sequence<INDICES...>, SIZE>
{	typedef typename getTypeFromSize<SIZE>::type TYPE;
	static constexpr const std::array<
		std::pair<TYPE, TYPE>,
		mp_size<LIST>::value
	> value = {std::make_pair(TYPE(mp_second<mp_at_c<LIST, INDICES> >::value), TYPE(mp_third<mp_at_c<LIST, INDICES> >::value))...};
};
template<typename LIST, std::size_t ...INDICES, typename SIZE>
const std::array<
	std::pair<typename getTypeFromSize<SIZE>::type, typename getTypeFromSize<SIZE>::type>,
	mp_size<LIST>::value
> convertToStdArrayImpl<LIST, std::index_sequence<INDICES...>, SIZE>::value;

template<typename LIST, typename SIZE>
struct convertToStdArray
{	typedef convertToStdArrayImpl<LIST, std::make_index_sequence<mp_size<LIST>::value>, SIZE> type;
};

template<typename, typename, typename>
struct convertToStdArray2Impl;
template<typename LIST, std::size_t ...INDICES, typename SIZE>
struct convertToStdArray2Impl<LIST, std::index_sequence<INDICES...>, SIZE>
{	typedef typename getTypeFromSize<SIZE>::type TYPE;
	static constexpr const std::array<
		TYPE,
		mp_size<LIST>::value
	> value = {TYPE(mp_second<mp_at_c<LIST, INDICES> >::value)...};
};
template<typename LIST, std::size_t ...INDICES, typename SIZE>
constexpr const std::array<
	typename getTypeFromSize<SIZE>::type,
	mp_size<LIST>::value
> convertToStdArray2Impl<LIST, std::index_sequence<INDICES...>, SIZE>::value;
template<typename LIST, typename SIZE>
struct convertToStdArray2
{	typedef convertToStdArray2Impl<LIST, std::make_index_sequence<mp_size<LIST>::value>, SIZE> type;
};
	/// for creating the type result of multiplying one element of a ctaylor array with another
	/// second and third arguments are an element of the first template argument of ctaylor
	/// calls itself recursively
template<typename, typename, typename>
struct multiply_1_1_R;
template<typename RESULT>
struct multiply_1_1_R<RESULT, mp_list<>, mp_list<> >
{	typedef RESULT type;
};
template<typename RESULT, typename T, typename ...REST>
struct multiply_1_1_R<RESULT, mp_list<T, REST...>, mp_list<> >
{	typedef typename multiply_1_1_R<
		mp_push_back<
			RESULT,
			T
		>,
		mp_list<REST...>,
		mp_list<>
	>::type type;
};
template<typename RESULT, typename T, typename ...REST>
struct multiply_1_1_R<RESULT, mp_list<>, mp_list<T, REST...> >
{	typedef typename multiply_1_1_R<
		mp_push_back<RESULT, T>,
		mp_list<>,
		mp_list<REST...>
	>::type type;
};
template<typename RESULT, typename T0, typename ...R0, typename T1, typename ...R1>
struct multiply_1_1_R<RESULT, mp_list<T0, R0...>, mp_list<T1, R1...> >
{	typedef typename std::conditional<
		(T0::first_type::value < T1::first_type::value),
		multiply_1_1_R<
			mp_push_back<RESULT, T0>,
			mp_list<R0...>,
			mp_list<T1, R1...>
		>,
		typename std::conditional<
			(T0::first_type::value > T1::first_type::value),
			multiply_1_1_R<
				mp_push_back<RESULT, T1>,
				mp_list<T0, R0...>,
				mp_list<R1...>
			>,
			multiply_1_1_R<
				mp_push_back<
					RESULT,
					pair<
						typename T0::first_type,
						mp_size_t<T0::second_type::value + T1::second_type::value>
					>
				>,
				mp_list<R0...>,
				mp_list<R1...>
			>
		>::type
	>::type::type type;
};
	/// type is only accessed in case of it is needed
template<typename STATE, typename T0E>
struct multiply_1_1_E
{	typedef mp_push_back<
		mp_first<STATE>,
		pair<
			typename multiply_1_1_R<
				mp_list<>,
				mp_first<T0E>,
				mp_first<mp_second<STATE> >// T1E
			>::type,
			mp_list<
				mp_list<
					mp_second<T0E>,
					mp_second<mp_second<STATE> >
				>
			>
		>
	> type;
};
	/// invokes multiply_1_1_R
	/// only if the resuling order would be smaller or equal MAX
template<typename STATE, typename T0E>
using multiply_1_1 = mp_list<
	typename std::conditional<
		(order<mp_first<T0E> >::value + order<mp_first<mp_second<STATE> > >::value <= mp_third<STATE>::value),
		multiply_1_1_E<STATE, T0E>,
		mp_identity<mp_first<STATE> >
	>::type::type,
	mp_second<STATE>,//T1E
	mp_third<STATE>//MAX
>;
	/// multiplies all elements of a template argument to ctaylor with one element
template<typename STATE, typename T1E>
using multiply_2_1 = mp_list<
	mp_first<STATE>,//T0
	typename merge<
		mp_first<
			mp_fold<
				mp_first<STATE>,	// T0
				mp_list<
					mp_list<>,
					T1E,
					mp_third<STATE>//MAX
				>,
				multiply_1_1
			>
		>,
		mp_second<STATE>,
		compareListOfPairs2,
		combineTwoPairs,
		containsValue2
	>::type,
	mp_third<STATE>//MAX
>;
template<typename A, typename B>
using make_pair = pair<A, B>;
// Transform function to pair each type with its index
template<typename L>
using add_index = mp_transform<
	make_pair,
	L,
	mp_iota_c<mp_size<L>::value>
>;
	/// multiplies two template arguments to ctaylor with each other
template<typename T0, typename T1, std::size_t MAX>
using multiply_2_2 = mp_second<
	mp_fold<
		add_index<T1>,
		mp_list<
			add_index<T0>,
			mp_list<>,
			mp_size_t<MAX>
		>,
		multiply_2_1
	>
>;
template<typename T0, typename T1>
struct TypeDisplayer
{	static_assert(
		std::is_same<
			T0,
			T1
		>::value,
		"types are not identical!"
	);
};
	/// the class
	/// first template argument is a vector of a vector of pairs of independent variable enum and order
	/// all vectors must be sorted
	/// usually first entry is mp_list<> indicating 0th derivative or value
	/// MAX indicates maximum order of derivatives calculated
	/// MAX should be minimally 1 for it to work
	/// MAX should be minimally 2 for application of this class to make sense otherwise jacobian.h ought to be used
template<typename T, std::size_t MAX>
struct ctaylor
{	typedef T SET;
	static constexpr std::size_t SIZE = mp_size<T>::value;
	//static_assert(SIZE > 0, "size must be at least one!");
	static_assert(mp_is_set<T>::value, "must be a set!");
	typedef std::array<double, SIZE> ARRAY;
	ARRAY m_s;
	ctaylor(void) = default;
	ctaylor(ctaylor&&) = default;
	ctaylor(const ctaylor&) = default;
	ctaylor&operator=(const ctaylor&) = default;
	ctaylor&operator=(ctaylor&&) = default;
		/// copying all elements from the source array to the destination array
		/// the destination array is necessarily larger
		/// and some elements will have to be initialized with zero
		/// used by copy constructor and assignment operator
	template<typename T1, bool CHECK = true>
	static ARRAY convert(const typename ctaylor<T1, MAX>::ARRAY&_r, const mp_bool<CHECK>& = mp_bool<CHECK>())
	{	typedef mp_plus<mp_size<T1>, mp_size_t<1> > SIZE;
		typedef typename findPositions<T, T1, SIZE, CHECK>::type SOURCE_POSITIONS;
		ARRAY s;
		typedef typename getTypeFromSize<SIZE>::type TYPE;
		auto &rT = convertToStdArray2<SOURCE_POSITIONS, SIZE>::type::value;
		std::transform(
			rT.begin(),
			rT.end(),
			s.begin(),
			[&](const std::size_t _i)
			{	if (_i == std::numeric_limits<TYPE>::max())
					return 0.0;
				else
					return _r[_i];
			}
		);
		return s;
	}
	template<
		typename U=T,
		typename std::enable_if<
			(mp_size<mp_first<U> >::value == 0),
			int
		>::type = 0
	>
	ctaylor &operator=(const double _d)
	{	m_s[0] = _d;
		std::fill(m_s.begin() + 1, m_s.end(), 0.0);
		return *this;
	}
	template<typename T1>
	ctaylor &operator=(const ctaylor<T1, MAX>&_r)
	{	static_assert(ctaylor<T1, MAX>::SIZE < SIZE, "RHS size must be smaller!");
		m_s = convert<T1, true>(_r.m_s);
		return *this;
	}
	template<
		typename U=T,
		typename std::enable_if<
			(mp_size<mp_first<U> >::value == 0),
			int
		>::type = 0
	>
	ctaylor &operator+=(const double _d)
	{	m_s[0] += _d;
		return *this;
	}
	template<
		typename U=T,
		typename std::enable_if<
			(mp_size<mp_first<U> >::value == 0),
			int
		>::type = 0
	>
	ctaylor &operator-=(const double _d)
	{	m_s[0] -= _d;
		return *this;
	}
	ctaylor &operator*=(const double _d)
	{	for (std::size_t i = 0; i < SIZE; ++i)
			m_s[i] *= _d;
		return *this;
	}
	ctaylor &operator/=(const double _d)
	{	const auto d1 = 1.0/_d;
		for (std::size_t i = 0; i < SIZE; ++i)
			m_s[i] *= d1;
		return *this;
	}
	ctaylor &operator+=(const ctaylor&_r)
	{	for (std::size_t i = 0; i < SIZE; ++i)
			m_s[i] += _r.m_s[i];
		return *this;
	}
	ctaylor &operator-=(const ctaylor&_r)
	{	for (std::size_t i = 0; i < SIZE; ++i)
			m_s[i] -= _r.m_s[i];
		return *this;
	}
	template<typename T1>
	ctaylor &operator+=(const ctaylor<T1, MAX>&_r)
	{	typedef mp_plus<mp_size<T1>, mp_size_t<1> > SIZE;
		typedef typename findPositions<T, T1, SIZE, true>::type SOURCE_POSITIONS;
		typedef typename getTypeFromSize<SIZE>::type TYPE;
		auto &rT = convertToStdArray2<SOURCE_POSITIONS, SIZE>::type::value;
		for (std::size_t i = 0; i < ctaylor::SIZE; ++i)
		{	const auto iT = rT[i];
			if (iT != std::numeric_limits<TYPE>::max())
				m_s[i] += _r.m_s[iT];
		}
		return *this;
	}
	template<typename T1>
	ctaylor &operator-=(const ctaylor<T1, MAX>&_r)
	{	typedef mp_plus<mp_size<T1>, mp_size_t<1> > SIZE;
		typedef typename findPositions<T, T1, SIZE, true>::type SOURCE_POSITIONS;
		typedef typename getTypeFromSize<SIZE>::type TYPE;
		auto &rT = convertToStdArray2<SOURCE_POSITIONS, SIZE>::type::value;
		for (std::size_t i = 0; i < ctaylor::SIZE; ++i)
		{	const auto iT = rT[i];
			if (iT != std::numeric_limits<TYPE>::max())
				m_s[i] -= _r.m_s[iT];
		}
		return *this;
	}
	template<
		typename U=T,
		typename std::enable_if<
			(mp_size<U>::value == 2	//{{}, {{enum, order}}}
				&& mp_size<mp_first<U> >::value == 0
				&& mp_size<mp_second<U> >::value == 1
				&& mp_size<mp_first<mp_second<U> > >::value == 2
				&& mp_second<mp_first<mp_second<U> > >::value == 1
			),
			int
		>::type = 0
	>
	ctaylor(const double _d, const bool)
		:m_s({_d, 1.0})
	{
	}
	ctaylor(const double _d)
		:m_s({_d})
	{
	}
	template<typename T1, bool CHECK = true>
	ctaylor(const ctaylor<T1, MAX>&_r, const mp_bool<CHECK>& = mp_bool<CHECK>())
		:m_s(convert<T1, CHECK>(_r.m_s))
	{	static_assert(!CHECK || ctaylor<T1, MAX>::SIZE < SIZE, "RHS size must be smaller!");
	}
	ctaylor operator+(const ctaylor&_r) const
	{	ctaylor s;
		for (std::size_t i = 0; i < SIZE; ++i)
			s.m_s[i] = m_s[i] + _r.m_s[i];
		return s;
	}
	ctaylor operator-(const ctaylor&_r) const
	{	ctaylor s;
		for (std::size_t i = 0; i < SIZE; ++i)
			s.m_s[i] = m_s[i] - _r.m_s[i];
		return s;
	}
	template<typename T1>
	auto operator+(const ctaylor<T1, MAX>&_r) const
	{	typedef typename merge<T, T1>::type TT;
		typedef typename findPositions2<TT, T, T1>::type SOURCE_POSITIONS;
		ctaylor<TT, MAX> s;
		typedef mp_plus<
			mp_max<
				mp_size<T>,
				mp_size<T1>
			>,
			mp_size_t<1>
		> SIZE;
		auto &rT = convertToStdArray<
			SOURCE_POSITIONS,
			SIZE
		>::type::value;
		typedef typename getTypeFromSize<SIZE>::type TYPE;
		std::transform(
			rT.begin(),
			rT.end(),
			s.m_s.begin(),
			[&](const std::pair<TYPE, TYPE>&_rI)
			{	return _rI.first != std::numeric_limits<TYPE>::max()
				? (_rI.second != std::numeric_limits<TYPE>::max()
					? m_s[_rI.first] + _r.m_s[_rI.second]
					: m_s[_rI.first]
				)
				: (_rI.second != std::numeric_limits<TYPE>::max()
					? _r.m_s[_rI.second]
					: 0.0
				);
			}
		);
		return s;
	}
	template<typename T1>
	auto operator-(const ctaylor<T1, MAX>&_r) const
	{	typedef typename merge<T, T1>::type TT;
		typedef typename findPositions2<TT, T, T1>::type SOURCE_POSITIONS;
		ctaylor<TT, MAX> s;
		typedef mp_plus<
			mp_max<
				mp_size<T>,
				mp_size<T1>
			>,
			mp_size_t<1>
		> SIZE;
		auto &rT = convertToStdArray<SOURCE_POSITIONS, SIZE>::type::value;
		typedef typename getTypeFromSize<SIZE>::type TYPE;
		std::transform(
			rT.begin(),
			rT.end(),
			s.m_s.begin(),
			[&](const std::pair<TYPE, TYPE>&_rI)
			{	return _rI.first != std::numeric_limits<TYPE>::max()
				? (_rI.second != std::numeric_limits<TYPE>::max()
					? m_s[_rI.first] - _r.m_s[_rI.second]
					: m_s[_rI.first]
				)
				: (_rI.second != std::numeric_limits<TYPE>::max()
					? -_r.m_s[_rI.second]
					: 0.0
				);
			}
		);
		return s;
	}
	auto operator-(void) const
	{	ctaylor s;
		for (std::size_t i = 0; i < SIZE; ++i)
			s.m_s[i] = -m_s[i];
		return s;
	}
	friend auto operator-(const double _d, const ctaylor&_r)
	{	ctaylor s;
		for (std::size_t i = 0; i < SIZE; ++i)
			s.m_s[i] = -_r.m_s[i];
		s.m_s[0] += _d;
		return s;
	}
	template<
		typename U=T,
		typename std::enable_if<
			(mp_size<mp_first<U> >::value == 0),
			int
		>::type = 0
	>
	ctaylor operator-(const double _d) const
	{	static_assert(mp_size<mp_first<T> >::value == 0, "must have value!");
		ctaylor s(*this);
		s.m_s[0] -= _d;
		return s;
	}
	ctaylor operator*(const double _d) const
	{	ctaylor s;
		for (std::size_t i = 0; i < SIZE; ++i)
			s.m_s[i] = m_s[i]*_d;
		return s;
	}
	friend ctaylor operator*(const double _d, const ctaylor&_r)
	{	return _r*_d;
	}
	template<typename T1>
	auto operator*(const ctaylor<T1, MAX>&_r) const
	{	typedef std::integral_constant<bool, (mp_size<T1>::value <= mp_size<T>::value)> BOOL;
		return multiply(_r, BOOL());
	}
	template<typename T1>
	auto multiply(const ctaylor<T1, MAX>&_r, const std::false_type&) const
	{	return _r**this;
	}
	template<typename T1>
	auto multiply(const ctaylor<T1, MAX>&_r, const std::true_type&) const
	{	typedef multiply_2_2<T, T1, MAX> CALCULATED_PAIRS;
		typedef mp_transform<
			mp_first,
			CALCULATED_PAIRS
		> CALCULATED;
		typedef mp_transform<
			mp_second,
			CALCULATED_PAIRS
		> POSITIONS;
		ctaylor<CALCULATED, MAX> s;
		auto &r = convertToStdArray3<POSITIONS, mp_max<mp_size<T>, mp_size<T1> > >::type::value;
		typedef typename getTypeFromSize<mp_max<mp_size<T>, mp_size<T1> > >::type TYPE;
		std::transform(
			r.begin(),
			r.end(),
			s.m_s.begin(),
			[&](const std::initializer_list<std::pair<TYPE, TYPE> >&_rIL)
			{	return std::accumulate(
					_rIL.begin(),
					_rIL.end(),
					double(),
					[&](const double _d, const std::pair<TYPE, TYPE> &_rP)
					{	return _d + m_s[_rP.first]*_r.m_s[_rP.second];
					}
				);
			}
		);
		return s;
	}
	template<
		typename U=T,
		typename std::enable_if<
			(mp_size<mp_first<U> >::value == 0),
			int
		>::type = 0
	>
	ctaylor operator+(const double _d) const
	{	ctaylor s(*this);
		s.m_s[0] += _d;
		return s;
	}
	template<
		typename U=T,
		typename std::enable_if<
			(mp_size<mp_first<U> >::value != 0),
			int
		>::type = 0
	>
	auto operator+(const double _d) const
	{	ctaylor<typename merge<T, mp_list<mp_list<> > >::type, MAX> s(*this);
		s.m_s[0] += _d;
		return s;
	}
	friend auto operator+(const double _d, const ctaylor&_r)
	{	return _r + _d;
	}
	template<std::size_t DIM, std::size_t POSM>
	auto apply(const std::array<double, DIM>&_r, const mp_size_t<POSM>&) const
	{	return _r[DIM - POSM] + apply(_r, mp_size_t<POSM-1>())**this;
	}
	template<std::size_t DIM>
	auto apply(const std::array<double, DIM>&_r, const mp_size_t<1>&) const
	{	return _r[DIM - 1];
	}
	template<
		typename T1,
		typename U=T,
		typename std::enable_if<
			(mp_size<mp_first<U> >::value == 0)
				&& (mp_size<mp_first<T1> >::value == 0),
			int
		>::type = 0
	>
	friend auto pow(
		const ctaylor&_r0,
		const ctaylor<T1, MAX>&_r1
	)
	{	return exp(log(_r0)*_r1);
	}
	template<
		typename U=T,
		typename std::enable_if<
			(mp_size<mp_first<U> >::value == 0),
			int
		>::type = 0
	>
	friend auto pow(
		const double _d0,
		const ctaylor&_r1
	)
	{	return exp(std::log(_d0)*_r1);
	}
	template<
		typename U=T,
		typename std::enable_if<
			(mp_size<mp_first<U> >::value == 0),
			int
		>::type = 0
	>
	friend auto pow(
		const ctaylor&_r0,
		const double _d1
	)
#if 0
	{	return exp(log(_r0)*_d1);
	}
#else
	{	std::array<double, MAX + 1> s;
		const auto d = 1.0/_r0.m_s[0];
		s[0] = std::pow(_r0.m_s[0], _d1);	//x^n
		auto &r = divide_by_n_p_1<MAX>::type::value;
		for (std::size_t i = 1; i < MAX + 1; ++i)
			s[i] = s[i - 1]*d*(_d1 - (i - 1))*r[i - 1];//n*x^(n - 1)
		return dropValue(_r0).apply(s, mp_size_t<MAX + 1>());
	}
#endif
	template<
		typename U=T,
		typename std::enable_if<
			(mp_size<mp_first<U> >::value == 0),
			int
		>::type = 0
	>
	friend auto abs(const ctaylor&_r)
	{	if (_r < 0.0)
			return -_r;
		else
			return _r;
	}
#if defined(__GNUC__) && !defined(__clang__)
	static constexpr double dTwoOverSqrtPi = 2.0/std::sqrt(M_PI);
#else
	static const double dTwoOverSqrtPi;
#endif
	template<
		std::size_t MAXM=MAX,
		typename std::enable_if<
			(MAXM == 1),
			int
		>::type = 0
	>
	friend auto erfc(const ctaylor&_r, const mp_size_t<MAXM>&)
	{	const auto d0 = std::erfc(value(_r));
		const auto d1 = -dTwoOverSqrtPi*std::exp(-value(_r)*value(_r));
		std::array<double, MAX + 1> s;
		s[0] = d0;
		s[1] = d1;
		return dropValue(_r).apply(s, mp_size_t<MAX + 1>());
	}
	template<
		std::size_t MAXM=MAX,
		typename std::enable_if<
			(MAXM > 1),
			int
		>::type = 0
	>	
	friend auto erfc(const ctaylor&_r, const mp_size_t<MAXM>&)
	{	const auto d0 = std::erfc(value(_r));
		const auto s1 = -dTwoOverSqrtPi*exp(-sqr(ctaylor<makeIndependent<0>, MAX - 1>(value(_r), true)));
		auto &r = divide_by_n_p_1<MAX>::type::value;
		std::array<double, MAX + 1> s;
		s[0] = d0;
		for (std::size_t i = 1; i < MAX + 1; ++i)
			s[i] = s1.m_s[i - 1]*r[i - 1];
		return dropValue(_r).apply(s, mp_size_t<MAX + 1>());
	}
	template<
		typename U=T,
		typename std::enable_if<
			(mp_size<mp_first<U> >::value == 0),
			int
		>::type = 0
	>
	friend auto erfc(const ctaylor&_r)
	{	return erfc(_r, mp_size_t<MAX>());
	}

	template<
		std::size_t MAXM=MAX,
		typename std::enable_if<
			(MAXM == 1),
			int
		>::type = 0
	>	
	friend auto erf(const ctaylor&_r, const mp_size_t<MAXM>&)
	{	const auto d0 = std::erf(value(_r));
		const auto d1 = dTwoOverSqrtPi*std::exp(-value(_r)*value(_r));
		std::array<double, MAX + 1> s;
		s[0] = d0;
		s[1] = d1;
		return dropValue(_r).apply(s, mp_size_t<MAX + 1>());
	}
	template<
		std::size_t MAXM=MAX,
		typename std::enable_if<
			(MAXM > 1),
			int
		>::type = 0
	>	
	friend auto erf(const ctaylor&_r, const mp_size_t<MAXM>&)
	{	const auto d0 = std::erf(value(_r));
		const auto s1 = dTwoOverSqrtPi*exp(-sqr(ctaylor<makeIndependent<0>, MAX - 1>(value(_r), true)));
		auto &r = divide_by_n_p_1<MAX>::type::value;
		std::array<double, MAX + 1> s;
		s[0] = d0;
		for (std::size_t i = 1; i < MAX + 1; ++i)
			s[i] = s1.m_s[i - 1]*r[i - 1];
		return dropValue(_r).apply(s, mp_size_t<MAX + 1>());
	}
	template<
		typename U=T,
		typename std::enable_if<
			(mp_size<mp_first<U> >::value == 0),
			int
		>::type = 0
	>
	friend auto erf(const ctaylor&_r)
	{	return erf(_r, mp_size_t<MAX>());
	}

	template<
		typename U=T,
		typename std::enable_if<
			(mp_size<mp_first<U> >::value == 0),
			int
		>::type = 0
	>
	friend auto value(const ctaylor&_r)
	{	return _r.m_s[0];
	}
	template<
		typename U=T,
		typename std::enable_if<
			(mp_size<mp_first<U> >::value == 0),
			int
		>::type = 0
	>
	friend auto tan(const ctaylor&_r)
	{	return tan(_r, mp_size_t<MAX>());
	}
	template<
		std::size_t MAXM=MAX,
		typename std::enable_if<
			(MAXM == 1),
			int
		>::type = 0
	>	
	friend auto tan(const ctaylor&_r, const mp_size_t<MAXM>&)
	{	const auto sTan = std::tan(value(_r));
		const auto s1 = 1.0 + sTan*sTan;
		std::array<double, MAX + 1> s;
		s[0] = sTan;
		s[1] = s1;
		return dropValue(_r).apply(s, mp_size_t<MAX + 1>());
		
	}
	template<
		std::size_t MAXM=MAX,
		typename std::enable_if<
			(MAXM > 1),
			int
		>::type = 0
	>	
	friend auto tan(const ctaylor&_r, const mp_size_t<MAXM>&)
	{	const auto sTan = tan(ctaylor<makeIndependent<0>, MAX - 1>(value(_r), true));
		const auto s1 = 1.0 + sqr(sTan);
		auto &r = divide_by_n_p_1<MAX>::type::value;
		std::array<double, MAX + 1> s;
		s[0] = sTan.m_s[0];
		for (std::size_t i = 1; i < MAX + 1; ++i)
			s[i] = s1.m_s[i - 1]*r[i - 1];
		return dropValue(_r).apply(s, mp_size_t<MAX + 1>());
	}
	template<
		typename U=T,
		typename std::enable_if<
			(mp_size<mp_first<U> >::value == 0),
			int
		>::type = 0
	>
	friend auto tanh(const ctaylor&_r)
	{	return tanh(_r, mp_size_t<MAX>());
	}
	template<
		std::size_t MAXM=MAX,
		typename std::enable_if<
			(MAXM == 1),
			int
		>::type = 0
	>	
	friend auto tanh(const ctaylor&_r, const mp_size_t<MAXM>&)
	{	const auto sTanh = std::tanh(value(_r));
		const auto s1 = 1.0 - sTanh*sTanh;
		std::array<double, MAX + 1> s;
		s[0] = sTanh;
		s[1] = s1;
		return dropValue(_r).apply(s, mp_size_t<MAX + 1>());
		
	}
	template<
		std::size_t MAXM=MAX,
		typename std::enable_if<
			(MAXM > 1),
			int
		>::type = 0
	>	
	friend auto tanh(const ctaylor&_r, const mp_size_t<MAXM>&)
	{	const auto sTanh = tanh(ctaylor<makeIndependent<0>, MAX - 1>(value(_r), true));
		const auto s1 = 1.0 - sqr(sTanh);
		auto &r = divide_by_n_p_1<MAX>::type::value;
		std::array<double, MAX + 1> s;
		s[0] = sTanh.m_s[0];
		for (std::size_t i = 1; i < MAX + 1; ++i)
			s[i] = s1.m_s[i - 1]*r[i - 1];
		return dropValue(_r).apply(s, mp_size_t<MAX + 1>());
	}
	template<
		typename U=T,
		typename std::enable_if<
			(mp_size<mp_first<U> >::value == 0),
			int
		>::type = 0
	>
	friend auto asin(const ctaylor&_r)
	{	const auto d0 = std::asin(value(_r));
		const auto s1 = 1.0/sqrt(1.0 - sqr(ctaylor<makeIndependent<0>, MAX - 1>(value(_r), true)));
		auto &r = divide_by_n_p_1<MAX>::type::value;
		std::array<double, MAX + 1> s;
		s[0] = d0;
		for (std::size_t i = 1; i < MAX + 1; ++i)
			s[i] = s1.m_s[i - 1]*r[i - 1];
		return dropValue(_r).apply(s, mp_size_t<MAX + 1>());
	}
	template<
		typename U=T,
		typename std::enable_if<
			(mp_size<mp_first<U> >::value == 0),
			int
		>::type = 0
	>
	friend auto acos(const ctaylor&_r)
	{	const auto d0 = std::acos(value(_r));
		const auto s1 = -1.0/sqrt(1.0 - sqr(ctaylor<makeIndependent<0>, MAX - 1>(value(_r), true)));
		auto &r = divide_by_n_p_1<MAX>::type::value;
		std::array<double, MAX + 1> s;
		s[0] = d0;
		for (std::size_t i = 1; i < MAX + 1; ++i)
			s[i] = s1.m_s[i - 1]*r[i - 1];
		return dropValue(_r).apply(s, mp_size_t<MAX + 1>());
	}
	template<
		typename U=T,
		typename std::enable_if<
			(mp_size<mp_first<U> >::value == 0),
			int
		>::type = 0
	>
	friend auto atan(const ctaylor&_r)
	{	const auto d0 = std::atan(value(_r));
		const auto s1 = 1.0/(1.0 + sqr(ctaylor<makeIndependent<0>, MAX - 1>(value(_r), true)));
		auto &r = divide_by_n_p_1<MAX>::type::value;
		std::array<double, MAX + 1> s;
		s[0] = d0;
		for (std::size_t i = 1; i < MAX + 1; ++i)
			s[i] = s1.m_s[i - 1]*r[i - 1];
		return dropValue(_r).apply(s, mp_size_t<MAX + 1>());
	}
	template<
		typename U=T,
		typename std::enable_if<
			(mp_size<mp_first<U> >::value == 0),
			int
		>::type = 0
	>
	friend auto asinh(const ctaylor&_r)
	{	const auto d0 = std::asinh(value(_r));
		const auto s1 = 1.0/sqrt(1.0 + sqr(ctaylor<makeIndependent<0>, MAX - 1>(value(_r), true)));
		auto &r = divide_by_n_p_1<MAX>::type::value;
		std::array<double, MAX + 1> s;
		s[0] = d0;
		for (std::size_t i = 1; i < MAX + 1; ++i)
			s[i] = s1.m_s[i - 1]*r[i - 1];
		return dropValue(_r).apply(s, mp_size_t<MAX + 1>());
	}
	template<
		typename U=T,
		typename std::enable_if<
			(mp_size<mp_first<U> >::value == 0),
			int
		>::type = 0
	>
	friend auto acosh(const ctaylor&_r)
	{	const auto d0 = std::acosh(value(_r));
		const auto s1 = 1.0/sqrt(sqr(ctaylor<makeIndependent<0>, MAX - 1>(value(_r), true)) - 1.0);
		auto &r = divide_by_n_p_1<MAX>::type::value;
		std::array<double, MAX + 1> s;
		s[0] = d0;
		for (std::size_t i = 1; i < MAX + 1; ++i)
			s[i] = s1.m_s[i - 1]*r[i - 1];
		return dropValue(_r).apply(s, mp_size_t<MAX + 1>());
	}
	template<
		typename U=T,
		typename std::enable_if<
			(mp_size<mp_first<U> >::value == 0),
			int
		>::type = 0
	>
	friend auto atanh(const ctaylor&_r)
	{	const auto d0 = std::atanh(value(_r));
		const auto s1 = 1.0/(1.0 - sqr(ctaylor<makeIndependent<0>, MAX - 1>(value(_r), true)));
		auto &r = divide_by_n_p_1<MAX>::type::value;
		std::array<double, MAX + 1> s;
		s[0] = d0;
		for (std::size_t i = 1; i < MAX + 1; ++i)
			s[i] = s1.m_s[i - 1]*r[i - 1];
		return dropValue(_r).apply(s, mp_size_t<MAX + 1>());
	}
	template<typename T1>
	auto operator/(const ctaylor<T1, MAX>&_r) const
	{	return *this*dropValue(_r).apply(
			inverse<MAX + 1>(value(_r)),
			mp_size_t<MAX + 1>()
		);
	}
	auto operator/(const double _d) const
	{	return *this*(1.0/_d);
	}
	friend auto operator/(const double _d, const ctaylor&_r)
	{	return _d*dropValue(_r).apply(
			inverse<MAX + 1>(value(_r)),
			mp_size_t<MAX + 1>()
		);
	}
#define __CREATE_NONLINEAR__(sin)\
	template<\
		typename U=T,\
		typename std::enable_if<\
			(mp_size<mp_first<U> >::value == 0),\
			int\
		>::type = 0\
	>\
	friend auto sin(const ctaylor&_r)\
	{	return dropValue(_r).apply(\
			sin<MAX + 1>(value(_r)),\
			mp_size_t<MAX + 1>()\
		);\
	}
	__CREATE_NONLINEAR__(exp)
	__CREATE_NONLINEAR__(log)
	__CREATE_NONLINEAR__(log10)
	__CREATE_NONLINEAR__(sqrt)
	__CREATE_NONLINEAR__(cbrt)
	__CREATE_NONLINEAR__(sin)
	__CREATE_NONLINEAR__(cos)
	__CREATE_NONLINEAR__(sinh)
	__CREATE_NONLINEAR__(cosh)
	struct output
	{	std::ostream&m_rS;
		const ctaylor&m_rT;
		output(
			std::ostream&_rS,
			const ctaylor&_rT
		)
			:m_rS(_rS),
			m_rT(_rT)
		{
		}
		template<typename ...LIST>
		void operator()(const mp_list<LIST...>&) const
		{	m_rS << m_rT.m_s[mp_find<SET, mp_list<LIST...> >::value];
			mp_for_each<mp_list<LIST...> >(*this);
			m_rS << "\n";
		}
		template<std::size_t I, std::size_t O>
		void operator()(const pair<mp_size_t<I>, mp_size_t<O> >&) const
		{	m_rS << "*x" << I << "^" << O;
		}
		template<std::size_t I>
		void operator()(const pair<mp_size_t<I>, mp_size_t<1> >&) const
		{	m_rS << "*x" << I;
		}
	};
	friend std::ostream &operator<<(std::ostream&_rS, const ctaylor&_r)
	{	_rS << "ctaylor(";
		mp_for_each<SET>(output(_rS, _r));
		//for (const auto d : _r.m_s)
			//_rS << d << ",";
		return _rS << ")";
	}
	template<
		typename U=T,
		typename std::enable_if<
			(mp_size<mp_first<U> >::value == 0),
			int
		>::type = 0
	>
	friend auto dropValue(const ctaylor&_r)
	{	return ctaylor<mp_pop_front<T>, MAX>(_r, mp_false());
	}
	friend auto sqr(const ctaylor&_r)
	{	return _r*_r;
	}
	template<
		typename T1,
		typename U=T,
		typename std::enable_if<
			(mp_size<mp_first<U> >::value == 0)
				&& (mp_size<mp_first<T1> >::value == 0),
			int
		>::type = 0
	>
	friend auto hypot(
		const ctaylor&_r0,
		const ctaylor<T1, MAX>&_r1
	)
	{	return sqrt(sqr(_r0) + sqr(_r1));
	}
	friend auto hypot(const ctaylor&_rX, const double _dY)
	{	return sqrt(sqr(_rX) + _dY*_dY);
	}
	friend auto hypot(const double _dX, const ctaylor&_rY)
	{	return sqrt(sqr(_rY) + _dX*_dX);
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
	template<
		typename T1,
		typename U=T,
		typename std::enable_if<
			(mp_size<mp_first<U> >::value == 0)
				&& (mp_size<mp_first<T1> >::value == 0),
			int
		>::type = 0
	>
	friend auto atan2(
		const ctaylor&_rY,
		const ctaylor<T1, MAX>&_rX
	)
	{	return atan2_(_rY, _rX);
	}
	friend auto atan2(const ctaylor&_rY, const double _dX)
	{	return atan2_(_rY, _dX);
	}
	friend auto atan2(const double _dY, const ctaylor&_rX)
	{	return atan2_(_dY, _rX);
	}
	friend auto max(const ctaylor&_r0, const double _r1)
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
	friend auto max(const double _r0, const ctaylor&_r1)
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
	friend auto max(const ctaylor&_r0, const ctaylor<T1, MAX>&_r1)
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
	friend auto min(const ctaylor&_r0, const double _r1)
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
	friend auto min(const double _r0, const ctaylor&_r1)
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
	friend auto min(const ctaylor&_r0, const ctaylor<T1, MAX>&_r1)
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
	friend auto fmod(const ctaylor&_r0, const double _r1)
	{	return _r0 - static_cast<int>(value(_r0)/_r1)*_r1;
	}
	friend auto fmod(const double _r0, const ctaylor&_r1)
	{	return _r0 - static_cast<int>(_r0/value(_r1))*_r1;
	}
	template<typename T1>
	friend auto fmod(const ctaylor&_r0, const ctaylor<T1, MAX>&_r1)
	{	return _r0 - static_cast<int>(value(_r0)/value(_r1))*_r1;
	}
		/// if this does not compile, check the order of ENUM-order pairs
	template<typename LIST_OF_PAIRS>
	auto getDer(const LIST_OF_PAIRS&) const
	{	static_assert(mp_find<T, LIST_OF_PAIRS>::value < SIZE, "derivative pattern not found! (wrong sorting?)");
		return m_s[mp_find<T, LIST_OF_PAIRS>::value]*accumulatedFactorial<LIST_OF_PAIRS>::value;
	}

	template<typename T1>
	bool operator==(const ctaylor<T1, MAX>&_r) const
	{	return value(*this) == value(_r);
	}
	bool operator==(const double _d) const
	{	return value(*this) == _d;
	}
	friend bool operator==(const double _d, const ctaylor&_r)
	{	return _d == value(_r);
	}

	template<typename T1>
	bool operator!=(const ctaylor<T1, MAX>&_r) const
	{	return value(*this) != value(_r);
	}
	bool operator!=(const double _d) const
	{	return value(*this) != _d;
	}
	friend bool operator!=(const double _d, const ctaylor&_r)
	{	return _d != value(_r);
	}

	template<typename T1>
	bool operator<(const ctaylor<T1, MAX>&_r) const
	{	return value(*this) < value(_r);
	}
	bool operator<(const double _d) const
	{	return value(*this) < _d;
	}
	friend bool operator<(const double _d, const ctaylor&_r)
	{	return _d < value(_r);
	}

	template<typename T1>
	bool operator>(const ctaylor<T1, MAX>&_r) const
	{	return value(*this) > value(_r);
	}
	bool operator>(const double _d) const
	{	return value(*this) > _d;
	}
	friend bool operator>(const double _d, const ctaylor&_r)
	{	return _d > value(_r);
	}

	template<typename T1>
	bool operator<=(const ctaylor<T1, MAX>&_r) const
	{	return value(*this) <= value(_r);
	}
	bool operator<=(const double _d) const
	{	return value(*this) <= _d;
	}
	friend bool operator<=(const double _d, const ctaylor&_r)
	{	return _d <= value(_r);
	}

	template<typename T1>
	bool operator>=(const ctaylor<T1, MAX>&_r) const
	{	return value(*this) >= value(_r);
	}
	bool operator>=(const double _d) const
	{	return value(*this) >= _d;
	}
	friend bool operator>=(const double _d, const ctaylor&_r)
	{	return _d >= value(_r);
	}
	friend bool isnan(const ctaylor&_r)
	{	return std::any_of(
			_r.m_s.begin(),
			_r.m_s.end(),
			static_cast<bool(*)(double)>(&std::isnan)
		);
	}
	friend bool isfinite(const ctaylor&_r)
	{	return std::all_of(
			_r.m_s.begin(),
			_r.m_s.end(),
			static_cast<bool(*)(double)>(&std::isfinite)
		);
	}
};
#if !defined(__GNUC__) || defined(__clang__)
template<typename T, std::size_t MAX>
const double ctaylor<T, MAX>::dTwoOverSqrtPi = 2.0/std::sqrt(M_PI);
#endif
template<typename T>
struct containsValue2
{	typedef mp_empty<
		mp_first<
			mp_first<T>
		>
	> type;
};
template<typename T>
struct containsValue
{	typedef mp_empty<
		mp_first<T>
	> type;
};
template<>
struct containsValue<mp_list<> >
{	typedef mp_false type;
};
template<
	typename NEW,
	typename T1,
	typename LIST_LHS_PAIR_LIST_AT_POS
>
struct check
{	static_assert(
		containsValue<NEW>::type::value == mp_and<
			typename containsValue<T1>::type,
			typename containsValue<LIST_LHS_PAIR_LIST_AT_POS>::type
		>::value,
		"problem"
	);
};
template<typename, typename>
struct common_type;
template<typename T, std::size_t MAX>
struct common_type<ctaylor<T, MAX>, double>
{	typedef ctaylor<T, MAX> type;
};
template<typename T, std::size_t MAX>
struct common_type<double, ctaylor<T, MAX> >
{	typedef ctaylor<T, MAX> type;
};
template<typename T0, typename T1, std::size_t MAX>
struct common_type<ctaylor<T0, MAX>, ctaylor<T1, MAX> >
{	typedef ctaylor<typename merge<T0, T1>::type, MAX> type;
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
using implementation::ctaylor;
using implementation::makeIndependent;
using implementation::if_;
using namespace implementation;
}
