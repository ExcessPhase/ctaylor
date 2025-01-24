#include "test.h"
#define BOOST_TEST_MODULE mytests
#include <boost/test/included/unit_test.hpp>


bool compare::operator()(const ind2Expr&_r0, const ind2Expr&_r1) const
{	const auto i0 = order(_r0);
	const auto i1 = order(_r1);
	if (i0 < i1)
		return true;
	else
	if (i1 < i0)
		return false;
	else
		return std::lexicographical_compare(
			_r0.crbegin(),
			_r0.crend(),
			_r1.crbegin(),
			_r1.crend()
		);
}
