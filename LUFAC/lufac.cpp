#include "lufac.h"
#include <set>
#include <iterator>
#include <iostream>
#include <stdexcept>
#include <numeric>
#include <algorithm>
namespace lufac
{
double findMaxInColumn(const index2Index2Double::const_iterator &_pRow, const index2Index2Double&_rM)
{	const auto iRow = _pRow->first;
	return std::abs(std::max_element(
		_pRow,
		_rM.end(),
		[iRow](const index2Index2Double::value_type&_r0, const index2Index2Double::value_type&_r1)
		{	const auto p1 = _r1.second.find(iRow);
			if (p1 != _r1.second.end())
			{	const auto p0 = _r0.second.find(iRow);
				if (p0 != _r0.second.end())
					return abs(p0->second) < abs(p1->second);
				else
					return true;
			}
			else
				return false;
		}
	)->second.at(iRow));
}
index2Index2Double factor(const index2Index2Double&_r, index2Double&_rV)
{	const auto iDim = _r.size();
#if 0
	const auto sCounts = [&](void)
	{	std::vector<std::size_t> sCC(iDim), sRC(iDim);
		for (const auto &rR : _r)
		{	sRC[rR.first] = rR.second.size();
			for (const auto &rC : rR.second)
				++sCC[rC.first];
		}
		return std::make_pair(sRC, sCC);
	}();
	const auto sN2RC = [&](void)
	{	std::map<std::size_t, std::set<std::size_t> > sN2C, sN2R;
		for (std::size_t i = 0; i < iDim; ++i)
		{	sN2C[sCounts.second[i]].insert(i);
			sN2R[sCounts.first[i]].insert(i);
		}
		return std::make_pair(sN2R, sN2C);
	}();
	const auto sI2ER = [&](void)
	{	std::vector<std::size_t> sI2ER;
		sI2ER.reserve(iDim);
		for (auto p = sN2RC.first.rbegin(), pEnd = sN2RC.first.rend();
			p != pEnd;
			++p)
			for (const auto i : p->second)
				sI2ER.push_back(i);
		return sI2ER;
	}();
	const auto sI2EC = [&](void)
	{	std::vector<std::size_t> sI2EC;
		sI2EC.reserve(iDim);
		for (auto p = sN2RC.second.rbegin(), pEnd = sN2RC.second.rend();
			p != pEnd;
			++p)
			for (const auto i : p->second)
				sI2EC.push_back(i);
		return sI2EC;
	}();
	const auto sE2IC = [&](void)
	{	std::vector<std::size_t> sE2IC(iDim);
		for (std::size_t i = 0; i < iDim; ++i)
			sE2IC[sI2EC[i]] = i;
		return sE2IC;
	}();
	auto sM = [&](void)
	{	
		index2Index2Double s;
		for (auto iR = 0; iR < iDim; ++iR)
		{	auto &rRNew = s[iR];
			for (const auto &rCOld : _r.at(sI2ER[iR]))
				rRNew[sE2IC[rCOld.first]] = rCOld.second;
		}
		return s;
	}();
#else
	auto sM = _r;
#endif
	//std::vector<std::size_t> sI2E(iDim);
	//std::iota(sI2E.begin(), sI2E.end(), std::size_t());
	for (auto pRow = sM.begin(); pRow != sM.end(); ++pRow)
	{	auto pPivot = pRow->second.find(pRow->first);
		const auto dMax = findMaxInColumn(pRow, sM);
		if (dMax == 0.0 || !std::isfinite(dMax) || std::isnan(dMax))
			throw std::logic_error("matrix is singular");
		while (true)
			if(pPivot != pRow->second.end() && std::abs(pPivot->second) >= dMax)
			{	const auto dPivot = pPivot->second = 1.0/pPivot->second;
				for (auto pRowTarget = std::next(pRow); pRowTarget != sM.end(); ++pRowTarget)
				{	const auto pFind = pRowTarget->second.find(pRow->first);
					if (pFind != pRowTarget->second.end())
					{	const auto d = pFind->second*dPivot;
						for (auto pS = std::next(pPivot);
							pS != pRow->second.end();
							++pS
						)
							pRowTarget->second[pS->first] -= d*pS->second;
					}
				}
				break;
			}
			else
			{	bool b = false;
				for (auto p = std::next(pRow); p != sM.end(); ++p)
				{	pPivot = p->second.find(pRow->first);
					if (pPivot == p->second.end())
						continue;
					else
					if (std::abs(pPivot->second) < dMax)
						continue;
					else
					{	std::swap(p->second, pRow->second);
						std::swap(_rV[pRow->first], _rV[p->first]);
						//std::swap(sI2E[pRow->first], sI2E[p->first]);
						b = true;
						break;
					}
				}
				if (!b)
					throw std::logic_error("Did not find a pivot!");
			}
	}
	return sM;
}
index2Double solve(const index2Index2Double&_rM, const index2Double&_rY)
{	index2Double s(_rY);
	const auto iDim = _rY.size();
	for (std::size_t iRow = 0; iRow < iDim; ++iRow)
	{	const auto dInvPivot = _rM.at(iRow).at(iRow);
		auto &rRow = s[iRow];
		for (std::size_t iDestRow = iRow + 1; iDestRow < iDim; ++iDestRow)
		{	const auto pFind = _rM.at(iDestRow).find(iRow);
			if (pFind != _rM.at(iDestRow).end())
				s[iDestRow] -= dInvPivot*pFind->second*rRow;
		}
	}
	for (std::ptrdiff_t iRow = iDim - 1; iRow >= 0; --iRow)
	{	const auto &rRow = _rM.at(iRow);
		auto &rY = s[iRow];
		for (auto p = rRow.rbegin(); p != rRow.rend() && p->first > iRow; ++p)
			rY -= p->second*s[p->first];
		rY *= rRow.at(iRow);
	}
	return s;
}	
}
