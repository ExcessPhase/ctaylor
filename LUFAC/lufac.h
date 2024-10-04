#pragma once
#include <map>
#include <vector>
namespace vbic95
{
enum class enumNodes;
}
namespace lufac
{
typedef std::map<std::size_t, double> index2Double;
typedef std::map<std::size_t, index2Double> index2Index2Double;
std::pair<std::vector<std::size_t>, index2Index2Double> factor(const index2Index2Double&_r);
index2Double solve(const index2Index2Double&_r, const index2Double&_rY, const std::vector<std::size_t> &_rI2E);
}
