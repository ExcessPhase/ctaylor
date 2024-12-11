#include <boost/mp11/utility.hpp>
#include <boost/mp11.hpp>
#include <type_traits>
#include <iostream>

	using namespace boost::mp11;
struct output1
{	std::ostream&m_r;
	output1(std::ostream&_r)
		:m_r(_r)
	{
	}
	template<std::size_t I, bool PRINT>
	void operator()(const mp_size_t<I>&, const std::integral_constant<bool, PRINT>&) const
	{	m_r << I << ",";
	}
	template<bool PRINT>
	void operator()(const mp_list<>&, const std::integral_constant<bool, PRINT>&) const
	{	if (PRINT)
			m_r << "()";
	}
	template<typename FIRST, typename ...ARGS, bool PRINT>
	void operator()(const mp_list<FIRST, ARGS...>&, const std::integral_constant<bool, PRINT>&) const
	{	if (PRINT)
			m_r << "(";
		(*this)(FIRST(), std::true_type());
		(*this)(mp_list<ARGS...>(), std::false_type());
		if (PRINT)
			m_r << ")";
	}

	template<std::size_t I>
	void operator()(const mp_size_t<I>&) const
	{	(*this)(mp_size_t<I>(), std::true_type());
	}
	void operator()(const mp_list<>&) const
	{	(*this)(mp_list<>(), std::true_type());
	}
	template<typename FIRST, typename ...ARGS>
	void operator()(const mp_list<FIRST, ARGS...>&) const
	{	(*this)(mp_list<FIRST, ARGS...>(), std::true_type());
	}
};
#include "ctaylor.h"
int main(int, char**)
{
	using namespace taylor;
	using namespace taylor::implementation;
	using namespace boost::mp11;
	
	typedef mp_list<
		mp_size_t<0>,
		mp_size_t<2>
	> T0;
	typedef mp_list<
		mp_size_t<1>,
		mp_size_t<3>
	> T1;
	typedef merge_sorted_set_of_sets<
		mp_less,
		mp_list<T0, T1>,
		combineTwo
	> MERGE;
	mp_for_each<MERGE::type>(output1(std::cerr));
	std::cerr << "\n";
	//mp_for_each<MERGE::type>(output1(std::cerr));
	//std::cerr << "\n";
#if 0
	typedef multiply_2_2<T3, T3, MAX> T4;
	typedef mp_transform<
		mp_first,
		T4
	> T5;
	typedef multiply_2_2<T5, T5, MAX> T6;
	typedef mp_transform<
		mp_first,
		T6
	> T7;
#endif
#if 0
	const auto s0 = T0(1.2, false);
	const auto s1 = T1(1.3, false);
	const auto s2 = s0*s1;
	std::cerr << value(s0)*value(s1) << "=" << s2 << "\n";
	const auto s3 = s2*s2;
	std::cerr << s3 << "\n";
	const auto s4 = s3*s3;
	std::cerr << s4 << "\n";
#endif
}
