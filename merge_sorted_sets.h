#pragma once
#include <boost/mp11.hpp>
#include <type_traits>
#include <utility>

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
template<template<typename, typename> class F, typename ...T, template<typename, typename> class MERGE>
struct merge_sorted_sets<F, mp_list<T...>, mp_list<T...>, MERGE>
{	using type = mp_list<T...>;
};

template<template<typename, typename> class F, typename... Ts, template<typename, typename> class MERGE>
struct merge_sorted_sets<F, mp_list<Ts...>, mp_list<>, MERGE>
{	using type = mp_list<Ts...>;
};

template<template<typename, typename> class F, typename... Ts, template<typename, typename> class MERGE>
struct merge_sorted_sets<F, mp_list<>, mp_list<Ts...>, MERGE>
{	using type = mp_list<Ts...>;
};

template<template<typename, typename> class F, typename... Ts1, typename... Ts2, template<typename, typename> class MERGE>
struct merge_sorted_sets<F, mp_list<Ts1...>, mp_list<Ts2...>, MERGE>
{	typedef mp_list<Ts1...> TS1;
	typedef mp_list<Ts2...> TS2;
	typedef mp_front<TS1> T1;
	typedef mp_front<TS2> T2;
	struct defer_true
	{	typedef std::pair<
			merge_sorted_sets<
				F,
				mp_pop_front<TS1>,
				TS2,
				MERGE
			>,
			mp_identity<T1>
		> type;
	};
	struct defer_false
	{	typedef mp_if<
			typename F<T2, T1>::type,
			std::pair<
				merge_sorted_sets<
					F,
					TS1,
					mp_pop_front<TS2>,
					MERGE
				>,
				mp_identity<T2>
			>,
			std::pair<
				merge_sorted_sets<
					F,
					mp_pop_front<TS1>,
					mp_pop_front<TS2>,
					MERGE
				>,
				MERGE<
					mp_front<
						TS2
					>,
					mp_front<
						TS1
					>
				>
			>
		> type;
	};
	typedef typename mp_if<
		typename F<T1, T2>::type,
		defer_true,
		defer_false
	>::type tmp;
	using type = mp_push_front<
		typename tmp::first_type::type,
		typename tmp::second_type::type
	>;
};
}
using implementation::merge_sorted_sets;
}
