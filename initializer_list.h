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


template<typename LIST, typename TYPE>
struct convertToStdInitializerList;


template<typename ...ITEMS, typename TYPE>
struct convertToStdInitializerList<mp_list<ITEMS...>, TYPE>
{	static constexpr const std::initializer_list<TYPE> value =
	{	ITEMS::value...
	};
};
template<typename ...ITEMS, typename TYPE>
constexpr const std::initializer_list<TYPE> convertToStdInitializerList<mp_list<ITEMS...>, TYPE>::value;


template<typename LIST_A, typename LIST_B, typename TYPE>
struct convertToPair
{	static_assert(mp_size<LIST_A>::value == mp_size<LIST_B>::value, "size must be identical!");
	typedef std::initializer_list<TYPE> IL;
	static constexpr const std::pair<IL, IL> value =
	{	convertToStdInitializerList<LIST_A, TYPE>::value,
		convertToStdInitializerList<LIST_B, TYPE>::value
	};
};

}
}
