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
template<
	template<typename, typename> class F,
	typename SETS,
	template<typename, typename> class MERGE=combineTwo
>
struct merge_sorted_set_of_sets
{	template<typename T0, typename T1>
	using compare=typename F<T0, T1>::type;
	template<typename LIST>
	using getMin = mp_front<
		LIST
	>;
	typedef mp_remove_if<
		SETS,
		mp_empty
	> SETS1;
	typedef mp_transform<
		getMin,
		SETS1
	> MINS;
	typedef mp_min_element<
		MINS,
		compare
	> MIN;
	template<typename T>
	using isMin = mp_and<
		mp_not<compare<MIN, T> >,
		mp_not<compare<T, MIN> >
	>;
	template<typename LIST>
	using removeNotMin = mp_filter<isMin, LIST>;
	typedef mp_transform<
		removeNotMin,
		SETS1
	> SETS_NOT_MIN_REMOVED;
	template<typename T>
	using not_empty = mp_not<mp_empty<T> >;
	typedef mp_filter<not_empty, SETS_NOT_MIN_REMOVED> SETS_NOT_MIN_REMOVED_3;
	template<typename A, typename B>
	using merge=typename MERGE<A, mp_front<B> >::type;
	typedef mp_fold<
		SETS_NOT_MIN_REMOVED_3,
		MIN,
		merge
	> MINC;
	template<typename T>
	using isNotMin = mp_or<
		compare<MIN, T>,
		compare<T, MIN>
	>;
	template<typename LIST>
	using removeMin = mp_filter<isNotMin, LIST>;
	typedef mp_transform<
		removeMin,
		SETS1
	> SETS_MIN_REMOVED;
	typedef mp_remove_if<
		SETS_MIN_REMOVED,
		mp_empty
	> SETS_MIN_REMOVED_2;
	typedef typename std::conditional<
		mp_empty<SETS_MIN_REMOVED_2>::value,
		mp_identity<mp_list<> >,
		merge_sorted_set_of_sets<F, SETS_MIN_REMOVED_2, MERGE>
	>::type::type REST;
	typedef mp_push_front<REST, MINC> type;
};
template<
	template<typename, typename> class F,
	typename Set1,
	typename Set2,
	template<typename, typename> class MERGE=combineTwo
>
struct merge_sorted_sets
{	typedef typename merge_sorted_set_of_sets<F, mp_list<Set1, Set2>, MERGE>::type type;
};
}
using implementation::merge_sorted_sets;
}
