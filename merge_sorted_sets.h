#pragma once
#include <boost/mp11.hpp>
#include <type_traits>

namespace taylor
{
namespace implementation
{
//namespace mp = boost::mp11;
using namespace boost::mp11;

template<template<typename, typename> class F, typename Set1, typename Set2>
struct merge_sorted_sets;

template<template<typename, typename> class F>
struct merge_sorted_sets<F, mp_list<>, mp_list<> >
{	using type = mp_list<>;
};

template<template<typename, typename> class F, typename... Ts>
struct merge_sorted_sets<F, mp_list<Ts...>, mp_list<> >
{	using type = mp_list<Ts...>;
};

template<template<typename, typename> class F, typename... Ts>
struct merge_sorted_sets<F, mp_list<>, mp_list<Ts...> >
{	using type = mp_list<Ts...>;
};

template<template<typename, typename> class F, typename T1, typename... Ts1, typename T2, typename... Ts2>
struct merge_sorted_sets<F, mp_list<T1, Ts1...>, mp_list<T2, Ts2...> >
{
#if 1
	typedef typename std::conditional<
		F<T1, T2>::type::value,
		mp_list<
			merge_sorted_sets<
				F,
				mp_list<Ts1...>,
				mp_list<T2, Ts2...>
			>,
			T1
		>,
		typename std::conditional<
			F<T2, T1>::type::value,
			mp_list<
				merge_sorted_sets<
					F,
					mp_list<T1, Ts1...>,
					mp_list<Ts2...>
				>,
				T2
			>,
			mp_list<
				merge_sorted_sets<
					F,
					mp_list<Ts1...>,
					mp_list<Ts2...>
				>,
				T2
			>
		>::type
	>::type tmp;
	typedef mp_push_front<
		typename mp_first<tmp>::type,
		mp_second<tmp>
	> type;
#else
	typedef typename std::conditional<
		F<T1, T2>::type::value,
		mp_push_front<
			typename merge_sorted_sets<
				F,
				mp_list<Ts1...>,
				mp_list<T2, Ts2...>
			>::type,
			T1
		>,
		typename std::conditional<
			F<T2, T1>::type::value,
			mp_push_front<
				typename merge_sorted_sets<
					F,
					mp_list<T1, Ts1...>,
					mp_list<Ts2...>
				>::type,
				T2
			>,
			mp_push_front<
				typename merge_sorted_sets<
					F,
					mp_list<Ts1...>,
					mp_list<Ts2...>
				>::type,
				T2
			>
		>::type
	>::type type;
#endif
};
}
using implementation::merge_sorted_sets;
}
