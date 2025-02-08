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
taylorMap read(const char *const _p, const bool _bOnlyFirstOrder)
{	taylorMap sMap;
	std::ifstream sFile(_p);
	if (!sFile)
		throw std::runtime_error(std::string("Cannot open file \"") + _p + "\"!");
	std::string sLine;
		/// a real number followed by an series of zero or more X[0-9]*
	const std::regex sRegEx(R"(^([+-]?(?:\d+(?:\.\d*)?)(?:[eE][+-]?\d+)?)((?:\*X\d+(?:\^\d+)?)*)\s*$)");
		/// a single *X[0-9]*
	const std::regex sVarRegEx(R"((?:\*X(\d+)(?:\^(\d+))?)\s*)");
		/// only a series of X0*X1 with an implicit 1.0
	const std::regex sRegEx2(R"(^((?:X\d+(?:\^\d+)?)(?:\*X\d+(?:\^\d+)?)*)\s*$)");
		/// the leading X0 without a '*'
	const std::regex sVarRegEx2(R"((?:X(\d+)(?:\^(\d+))?)\s*)");

	//const auto sPath = std::filesystem::current_path();
	//std::cerr << sPath << "\n";
	while (std::getline(sFile, sLine))
	{	std::smatch sMatch;
		if (std::regex_match(sLine, sMatch, sRegEx))
		{ // Extract the numerical coefficient
			const double d = std::stod(sMatch[1]);
			ind2Expr sI2E;
			const std::string sVars = sMatch[2]; // match.position(5) gives position after the coefficient
			std::smatch sVarsMatch;
			auto pStart = sVars.cbegin();
			while (pStart != sVars.cend())
				if (std::regex_search(pStart, sVars.cend(), sVarsMatch, sVarRegEx))
				{	//const std::string &sVar = sVarsMatch[1];
					const auto iID = std::atoi(std::string(sVarsMatch[1]).c_str());
					const std::string sExp = sVarsMatch[2];
					const auto iExp = sExp.empty() ? 1 : std::atoi(sExp.c_str());
					sI2E.emplace(iID, iExp);
					pStart = sVarsMatch.suffix().first;
					//sVars = sVars.substr(sVarsMatch.str(0).size());
				}
				else
				{	std::cerr << "unmatched line: " << std::string(pStart, sVars.cend()) << std::endl;
					break;
				}
			if (!_bOnlyFirstOrder || order(sI2E) < 2)
				sMap.emplace(sI2E, d);
		}
		else
		if (std::regex_match(sLine, sMatch, sRegEx2))
		{ // Extract the numerical coefficient
			ind2Expr sI2E;
			const std::string sVars = sMatch[1]; // match.position(5) gives position after the coefficient
			std::smatch sVarsMatch;
			auto pStart = sVars.cbegin();
			while (pStart != sVars.cend())
				if (std::regex_search(pStart, sVars.cend(), sVarsMatch, pStart == sVars.cbegin() ? sVarRegEx2 : sVarRegEx))
				{	//const std::string &sVar = sVarsMatch[1];
					const auto iID = std::atoi(std::string(sVarsMatch[1]).c_str());
					const std::string sExp = sVarsMatch[2];
					const auto iExp = sExp.empty() ? 1 : std::atoi(sExp.c_str());
					sI2E.emplace(iID, iExp);
					pStart = sVarsMatch.suffix().first;
					//sVars = sVars.substr(sVarsMatch.str(0).size());
				}
				else
				{	std::cerr << "unmatched line: " << std::string(pStart, sVars.cend()) << std::endl;
					break;
				}
			if (!_bOnlyFirstOrder || order(sI2E) < 2)
				sMap.emplace(sI2E, 1.0);
		}
		else
			std::cerr << "unmatched line: " << sLine << std::endl;
	}
	return sMap;
}
