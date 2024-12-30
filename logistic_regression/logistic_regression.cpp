#include "../cjacobian.h"
#include <cmath>
#include <numeric>
#include <vector>
#include <boost/iterator/zip_iterator.hpp>
#include <cassert>
#include <utility>

std::pair<double, double> logistic(
	const std::vector<double> &_rX,
	const std::vector<double>&_rBeta
)
{	using namespace jacobian;
	using namespace boost::mp11;
	assert(_rX.size() == _rBeta.size());
	const auto s = 1.0/(1.0 + exp(
		-cjacobian<mp_list<mp_size_t<0> > >(
			std::accumulate(
				boost::make_zip_iterator(boost::make_tuple(_rX.cbegin(), _rBeta.cbegin())),
				boost::make_zip_iterator(boost::make_tuple(_rX.cend(), _rBeta.cend())),
				0.0,
				[&](const double _dPrev, const boost::tuple<const double&, const double&>&_r)
				{	return _dPrev + _r.get<0>()*_r.get<1>();
				}
			),
			true
		)
	));
	return std::make_pair(value(s), s.getDer(mp_size_t<0>()));
}
std::pair<double, double> entropy(const double _dY, const std::pair<double, double>& _rP)
{	using namespace jacobian;
	using namespace boost::mp11;

	const auto sP = cjacobian<mp_list<mp_size_t<0> > >(_rP.first, true);
	
	const auto s = - _dY*log(sP) - (1.0 - _dY)*log(1.0 - sP);
	return std::make_pair(value(s), s.getDer(mp_size_t<0>())*_rP.second);
}
auto entropy(
	const double _dY,
	const std::vector<double> &_rX,
	const std::vector<double>&_rBeta
)
{	return entropy(_dY, logistic(_rX, _rBeta));
}
#include <iostream>

int main()
{	static constexpr std::size_t SIZE = 100;
	const auto sX = [&](void)
	{	std::vector<double> s;
		s.reserve(SIZE);
		for (std::size_t i = 0; i < SIZE; ++i)
			s.push_back((i + 1) / 1000.0);
		return s;
	}();
	const auto s = entropy(
		0.5,
		sX,
		[&](void)
		{	std::vector<double> s;
			s.reserve(SIZE);
			for (std::size_t i = 0; i < SIZE; ++i)
				s.push_back(1.0 + i / 100.0);
			return s;
		}()
	);
	std::cout << "value=" << s.first << "\n";
	for (const auto d : sX)
		std::cout << "der=" << d*s.second << "\n";
}
