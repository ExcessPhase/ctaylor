#include "pch.h"
#include "CppUnitTest.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace ctaylortest
{
	TEST_CLASS(ctaylortest)
	{
	public:
		
		TEST_METHOD(TestMethod1)
		{	using namespace taylor;
			using namespace boost::mp11;

			constexpr const std::size_t MAX = 2;
			const auto s0 = ctaylor<makeIndependent<0>, MAX>(1.2, false);
				/// create an independent variable for x1 (this is what the unused boolean is for)
			const auto s1 = ctaylor<makeIndependent<1>, MAX>(1.3, false);
			const auto s2 = ctaylor<makeIndependent<2>, MAX>(1.4, false);
			const auto s3 = ctaylor<makeIndependent<3>, MAX>(1.5, false);
				/// some calculation
			const auto s4 = -s0 + s1 - s2 + s1*s2 - s0*s1 + s2*s3;
			//const auto s41 = s4.convert2Independent(mp_size_t<4>());
			const auto s5 = fmod(s4*s4, 1.0 - s4*s4);
			//const auto s52 = fmod(s41*s41, 1.0 - s41*s41);
			//const auto s51 = s52.chainRule(s4, mp_size_t<4>());
			constexpr const double s5_value = 0.0112;
			constexpr const double s5_x0 = 39.008;
			constexpr const double s5_x1 = -20.352;
			constexpr const double s5_x2 = -30.528;
			constexpr const double s5_x3 = -23.744;
			constexpr const double s5_x0_x0 = -42.32*taylor::implementation::factorial<2>::value;
			constexpr const double s5_x0_x1 = 61.12;
			constexpr const double s5_x0_x2 = 66.24;
			constexpr const double s5_x0_x3 = 51.52;

			constexpr const double s5_x1_x1 = -11.52*taylor::implementation::factorial<2>::value;
			constexpr const double s5_x1_x2 = -51.52;
			constexpr const double s5_x1_x3 = -26.88;

			constexpr const double s5_x2_x2 = -25.92*taylor::implementation::factorial<2>::value;
			constexpr const double s5_x2_x3 = -57.28;

			constexpr const double s5_x3_x3 = -15.68*taylor::implementation::factorial<2>::value;
			constexpr const double epsilon = 1e-9;
			static_assert(decltype(s5)::SIZE == 15, "SIZE");
			Assert::IsTrue(std::fabs(value(s5) - s5_value) < std::fabs(s5_value*epsilon), L"s5 value test");
			Assert::IsTrue(std::fabs(s5.getDer(mp_list<pair<mp_size_t<0>, mp_size_t<1> > >()) - s5_x0) < std::fabs(s5_x0*epsilon), L"s5_x0 test");
			Assert::IsTrue(std::fabs(s5.getDer(mp_list<pair<mp_size_t<1>, mp_size_t<1> > >()) - s5_x1) < std::fabs(s5_x1*epsilon), L"s5_x1 test");
			Assert::IsTrue(std::fabs(s5.getDer(mp_list<pair<mp_size_t<2>, mp_size_t<1> > >()) - s5_x2) < std::fabs(s5_x2*epsilon), L"s5_x2 test");
			Assert::IsTrue(std::fabs(s5.getDer(mp_list<pair<mp_size_t<3>, mp_size_t<1> > >()) - s5_x3) < std::fabs(s5_x3*epsilon), L"s5_x3 test");
			Assert::IsTrue(std::fabs(s5.getDer(mp_list<pair<mp_size_t<0>, mp_size_t<2> > >()) - s5_x0_x0) < std::fabs(s5_x0_x0*epsilon), L"s5_x0_x0 test");
			Assert::IsTrue(std::fabs(s5.getDer(mp_list<pair<mp_size_t<0>, mp_size_t<1> >, pair<mp_size_t<1>, mp_size_t<1> > >()) - s5_x0_x1) < std::fabs(s5_x0_x1*epsilon), L"s5_x0_x1 test");
			Assert::IsTrue(std::fabs(s5.getDer(mp_list<pair<mp_size_t<0>, mp_size_t<1> >, pair<mp_size_t<2>, mp_size_t<1> > >()) - s5_x0_x2) < std::fabs(s5_x0_x2*epsilon), L"s5_x0_x2 test");
			Assert::IsTrue(std::fabs(s5.getDer(mp_list<pair<mp_size_t<0>, mp_size_t<1> >, pair<mp_size_t<3>, mp_size_t<1> > >()) - s5_x0_x3) < std::fabs(s5_x0_x3*epsilon), L"s5_x0_x3 test");

			Assert::IsTrue(std::fabs(s5.getDer(mp_list<pair<mp_size_t<1>, mp_size_t<2> > >()) - s5_x1_x1) < std::fabs(s5_x1_x1*epsilon), L"s5_x1_x1 test");
			Assert::IsTrue(std::fabs(s5.getDer(mp_list<pair<mp_size_t<1>, mp_size_t<1> >, pair<mp_size_t<2>, mp_size_t<1> > >()) - s5_x1_x2) < std::fabs(s5_x1_x2*epsilon), L"s5_x1_x2 test");
			Assert::IsTrue(std::fabs(s5.getDer(mp_list<pair<mp_size_t<1>, mp_size_t<1> >, pair<mp_size_t<3>, mp_size_t<1> > >()) - s5_x1_x3) < std::fabs(s5_x1_x3*epsilon), L"s5_x1_x3 test");

			Assert::IsTrue(std::fabs(s5.getDer(mp_list<pair<mp_size_t<2>, mp_size_t<2> > >()) - s5_x2_x2) < std::fabs(s5_x2_x2*epsilon), L"s5_x2_x2 test");
			Assert::IsTrue(std::fabs(s5.getDer(mp_list<pair<mp_size_t<2>, mp_size_t<1> >, pair<mp_size_t<3>, mp_size_t<1> > >()) - s5_x2_x3) < std::fabs(s5_x2_x3*epsilon), L"s5_x2_x3 test");

			Assert::IsTrue(std::fabs(s5.getDer(mp_list<pair<mp_size_t<3>, mp_size_t<2> > >()) - s5_x3_x3) < std::fabs(s5_x3_x3*epsilon), L"s5_x3_x3 test");
#if 0
(-15.68*X3^2)-57.28*X2*X3-26.88*X1*X3+51.52*X0*X3-23.744*X3-25.92*X2^2
-51.52*X1*X2+66.24*X0*X2-30.528*X2-11.52*X1^2+61.12*X0*X1
-20.352*X1-42.32*X0^2+39.008*X0+0.0112
#endif
		}
	};
}
