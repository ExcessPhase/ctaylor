#pragma once
#include <map>
#include <vector>
//#include <boost/multiprecision/cpp_dec_float.hpp>
namespace vbic95
{
enum class enumNodes;
}
namespace lufac
{
typedef std::map<std::size_t, double> index2Double;
typedef std::map<std::size_t, index2Double> index2Index2Double;
typedef std::map<std::size_t, index2Index2Double> index2Index2Index2Double;
index2Index2Double factor(const index2Index2Double&_r, index2Double&);
index2Double solve(const index2Index2Double&_r, const index2Double&_rY);
}
