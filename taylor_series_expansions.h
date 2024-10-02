#pragma once
#include <cmath>
#include <array>
namespace taylor
{
template<std::size_t SIZE>
std::array<double, SIZE> exp(double _d)
{	std::array<double, SIZE> s;
	double d = std::exp(_d);
	for (std::size_t i = 0; i < SIZE; ++i)
	{	s[i] = d;
		d /= i + 1;
	}
	return s;
}
template<std::size_t SIZE>
std::array<double, SIZE> log(const double _d)
{	std::array<double, SIZE> s;
	const double d1 = 1.0/_d;
	double d = d1;
	s[0] = std::log(_d);
	for (std::size_t i = 1; i < SIZE; ++i)
	{	s[i] = d;
		d *= -d1*i/(i + 1);
	}
	return s;
}
static const auto s_dLog10 = 1.0/std::log(10.0);
template<std::size_t SIZE>
std::array<double, SIZE> log10(const double _d)
{	std::array<double, SIZE> s;
	const double d1 = 1.0/_d;
	double d = d1*s_dLog10;
	s[0] = std::log10(_d);
	for (std::size_t i = 1; i < SIZE; ++i)
	{	s[i] = d;
		d *= -d1*i/(i + 1);
	}
	return s;
}
template<std::size_t SIZE>
std::array<double, SIZE> sin(const double _d)
{	std::array<double, SIZE> s;
	const double ds = std::sin(_d);
	const double dc = std::cos(_d);
	double d = 1.0;
	//0=sin
	//1=cos
	//2=-sin
	//3=-cos
	//4=0
	for (std::size_t i = 0; i < SIZE; ++i)
	{	s[i] = (i & 1 ? i & 2 ? -dc : dc : i & 2 ? -ds : ds)*d;
		d /= i + 1;
	}
	return s;
}
template<std::size_t SIZE>
std::array<double, SIZE> cos(const double _d)
{	std::array<double, SIZE> s;
	const double ds = std::sin(_d);
	const double dc = std::cos(_d);
	double d = 1.0;
	//0=cos
	//1=-sin
	//2=-cos
	//3=sin
	//4=0
	for (std::size_t i = 0; i < SIZE; ++i)
	{	s[i] = (i & 1 ? i & 2 ? ds : -ds : i & 2 ? -dc: dc)*d;
		d /= i + 1;
	}
	return s;
}
template<std::size_t SIZE>
std::array<double, SIZE> sinh(const double _d)
{	std::array<double, SIZE> s;
	const double ds = std::sinh(_d);
	const double dc = std::cosh(_d);
	double d = 1.0;
	for (std::size_t i = 0; i < SIZE; ++i)
	{	s[i] = (i & 1 ? dc : ds)*d;
		d /= i + 1;
	}
	return s;
}
template<std::size_t SIZE>
std::array<double, SIZE> cosh(const double _d)
{	std::array<double, SIZE> s;
	const double ds = std::sinh(_d);
	const double dc = std::cosh(_d);
	double d = 1.0;
	for (std::size_t i = 0; i < SIZE; ++i)
	{	s[i] = (i & 1 ? ds : dc)*d;
		d /= i + 1;
	}
	return s;
}
template<std::size_t SIZE>
std::array<double, SIZE> sqrt(const double _d)
{	std::array<double, SIZE> s;
	const double d0 = std::sqrt(_d);
	const double d1 = 1.0/_d;
	double dPow = 0.5;
	double d = d0;
	//0=x^0.5
	//1=0.5*x^-0.5
	//2=-0.25*x^-1.5
	//3=
	for (std::size_t i = 0; i < SIZE; ++i)
	{	s[i] = d;	// sqrt(_d),
		d *= d1*dPow/(i + 1);	// sqrt(_d)/_d*0.5, sqrt(_d)/_d*0.5/_d*-0.5
		dPow -= 1.0;	// -0.5
	}
	return s;
}
template<std::size_t SIZE>
std::array<double, SIZE> inverse(const double _d)
{	std::array<double, SIZE> s;
	const double ds = 1.0/_d;
	double d = ds;
	//f0=x^-1/0!
	//f1=-x^-2/1!
	//f2=2x^-3/2!
	//f3=-6x^-4/3!
	for (std::size_t i = 0; i < SIZE; ++i)
	{	s[i] = i & 1 ? -d : d;
		d *= ds;
	}
	return s;
}
}
#if 0
((n & 1) ? ((n & 2) ? -cosh_x0 : cosh_x0) : ((n & 2) ? -sinh_x0 : sinh_x0)) * term
#endif
