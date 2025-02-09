#include "test.h"
#include <boost/test/unit_test.hpp>
#include <boost/math/constants/constants.hpp>

using namespace boost::math::constants;
template <typename X>
static auto Phi(X const& x) {
  return 0.5 * erfc(-one_div_root_two<double>() * x);
}

template <typename Price, typename Sigma, typename Tau, typename Rate>
static auto black_scholes_option_price_put(
                                                            double K,
                                                            Price const& S,
                                                            Sigma const& sigma,
                                                            Tau const& tau,
                                                            Rate const& r) {
  using namespace std;
  auto const d1 = (log(S / K) + (r + sigma * sigma / 2.0) * tau) / (sigma * sqrt(tau));
  auto const d2 = (log(S / K) + (r - sigma * sigma / 2.0) * tau) / (sigma * sqrt(tau));
      return exp(-r * tau) * K * Phi(-d2) - S * Phi(-d1);
}


BOOST_AUTO_TEST_CASE(taylor_35)
{	using namespace taylor;
	using namespace boost::mp11;
	const auto sMap = read("data35.txt");
	double const K = 100.0;  // Strike price.
	constexpr std::size_t MAX = 2;
	auto const S = ctaylor<makeIndependent<0>, MAX>(105, true);
	auto const sigma = ctaylor<makeIndependent<1>, MAX>(5, true);
	auto const tau = ctaylor<makeIndependent<2>, MAX>(30.0 / 365, true);
	auto const r = ctaylor<makeIndependent<3>, MAX>(1.25 / 100, true);
	const auto s5 = black_scholes_option_price_put(K, S, sigma, tau, r);
	__EQUAL_TAYLOR__();
}
BOOST_AUTO_TEST_CASE(jacobian_35)
{	using namespace jacobian;
	using namespace boost::mp11;
	const auto sMap = read("data35.txt", true);
	double const K = 100.0;  // Strike price.
	constexpr std::size_t MAX = 2;
	auto const S = cjacobian<mp_list<mp_size_t<0> > >(105, true);
	auto const sigma = cjacobian<mp_list<mp_size_t<1> > >(5, true);
	auto const tau = cjacobian<mp_list<mp_size_t<2> > >(30.0 / 365, true);
	auto const r = cjacobian<mp_list<mp_size_t<3> > >(1.25 / 100, true);
	const auto s5 = black_scholes_option_price_put(K, S, sigma, tau, r);
	__EQUAL_JACOBIAN__();
}
