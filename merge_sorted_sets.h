#pragma once
#include <boost/mp11.hpp>
#include <type_traits>

namespace taylor
{
namespace implementation
{
//namespace mp = boost::mp11;
using namespace boost::mp11;

template<typename, typename>
struct combineTwo;
template<typename A>
struct combineTwo<A, A>
{	typedef A type;
};
template<template<typename, typename> class F, typename Set1, typename Set2, template<typename, typename> class MERGE=combineTwo>
struct merge_sorted_sets;

template<template<typename, typename> class F, template<typename, typename> class MERGE>
struct merge_sorted_sets<F, mp_list<>, mp_list<>, MERGE>
{	using type = mp_list<>;
};

template<template<typename, typename> class F, typename... Ts, template<typename, typename> class MERGE>
struct merge_sorted_sets<F, mp_list<Ts...>, mp_list<>, MERGE>
{	using type = mp_list<Ts...>;
};

template<template<typename, typename> class F, typename... Ts, template<typename, typename> class MERGE>
struct merge_sorted_sets<F, mp_list<>, mp_list<Ts...>, MERGE>
{	using type = mp_list<Ts...>;
};

template<template<typename, typename> class F, typename T1, typename... Ts1, typename T2, typename... Ts2, template<typename, typename> class MERGE>
struct merge_sorted_sets<F, mp_list<T1, Ts1...>, mp_list<T2, Ts2...>, MERGE>
{
	struct defer_true
	{	typedef mp_list<
			merge_sorted_sets<
				F,
				mp_list<Ts1...>,
				mp_list<T2, Ts2...>,
				MERGE
			>,
			mp_identity<T1>
		> type;
	};
	struct defer_false
	{	typedef mp_if<
			typename F<T2, T1>::type,
			mp_list<
				merge_sorted_sets<
					F,
					mp_list<T1, Ts1...>,
					mp_list<Ts2...>,
					MERGE
				>,
				mp_identity<T2>
			>,
			mp_list<
				merge_sorted_sets<
					F,
					mp_list<Ts1...>,
					mp_list<Ts2...>,
					MERGE
				>,
				MERGE<T2, T1>
			>
		> type;
	};
	typedef typename mp_if<
		typename F<T1, T2>::type,
		defer_true,
		defer_false
	>::type tmp;
	typedef mp_push_front<
		typename mp_first<tmp>::type,
		typename mp_second<tmp>::type
	> type;
};
}
using implementation::merge_sorted_sets;
}
