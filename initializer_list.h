#pragma once
#include <initializer_list>
#include <boost/mp11.hpp>
#include <limits>
#include <utility>

namespace foelsche
{
namespace init_list
{
using namespace boost::mp11;


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
struct convertToStdInitializerListImpl;


template<typename LIST, std::size_t ...INDICES, typename TYPE>
struct convertToStdInitializerListImpl<LIST, std::index_sequence<INDICES...>, TYPE>
{	static constexpr const std::initializer_list<TYPE> value =
	{	TYPE(mp_at_c<LIST, INDICES>::value)...
	};
};


template<typename LIST, typename TYPE>
struct convertToStdInitializerList
{	typedef convertToStdInitializerListImpl<
		LIST,
		std::make_index_sequence<mp_size<LIST>::value>,
		TYPE
	> type;
};


template<typename LIST_A, typename LIST_B, typename TYPE>
struct convertToPair
{	static_assert(mp_size<LIST_A>::value == mp_size<LIST_B>::value, "size must be identical!");
	typedef std::initializer_list<TYPE> IL;
	static constexpr const std::pair<IL, IL> value =
	{	convertToStdInitializerList<LIST_A, TYPE>::type::value,
		convertToStdInitializerList<LIST_B, TYPE>::type::value
	};
};

}
}