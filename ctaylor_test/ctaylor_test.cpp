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
			constexpr const double s5_x0_x0 = -42.32;
			constexpr const double s5_x0_x1 = 61.12;
			constexpr const double s5_x0_x2 = 66.24;
			constexpr const double s5_x0_x3 = 51.52;

			constexpr const double s5_x1_x1 = -11.52;
			constexpr const double s5_x1_x2 = -51.52;
			constexpr const double s5_x1_x3 = -26.88;

			constexpr const double s5_x2_x2 = -25.92;
			constexpr const double s5_x2_x3 = -57.28;

			constexpr const double s5_x3_x3 = -15.68;
#define __TEST1__(ID, s5_x0) Assert::IsTrue(std::fabs(s5.getDer(mp_list<pair<mp_size_t<ID>, mp_size_t<1> > >()) - s5_x0) < std::fabs(s5_x0*epsilon), #s5_x0 L" test")
#define __TEST2__(ID, s5_x0_x0) Assert::IsTrue(std::fabs(s5.getDer(mp_list<pair<mp_size_t<ID>, mp_size_t<2> > >()) - s5_x0_x0*taylor::implementation::factorial<2>::value) < std::fabs(s5_x0_x0*taylor::implementation::factorial<2>::value*epsilon), #s5_x0_x0 L" test")
#define __TEST11__(ID0, ID1, s5_x0_x1) Assert::IsTrue(std::fabs(s5.getDer(mp_list<pair<mp_size_t<ID0>, mp_size_t<1> >, pair<mp_size_t<ID1>, mp_size_t<1> > >()) - s5_x0_x1) < std::fabs(s5_x0_x1*epsilon), #s5_x0_x1 L" test")


#define __TEST__(__SIZE__)\
			constexpr const double epsilon = 1e-12;\
			static_assert(decltype(s5)::SIZE == __SIZE__, "SIZE");\
			Assert::IsTrue(std::fabs(value(s5) - s5_value) < std::fabs(s5_value*epsilon), L"s5 value test");\
			__TEST1__(0, s5_x0);\
			__TEST1__(1, s5_x1);\
			__TEST1__(2, s5_x2);\
			__TEST1__(3, s5_x3);\
			__TEST2__(0, s5_x0_x0);\
			__TEST2__(1, s5_x1_x1);\
			__TEST2__(2, s5_x2_x2);\
			__TEST2__(3, s5_x3_x3);\
			__TEST11__(0, 1, s5_x0_x1);\
			__TEST11__(0, 2, s5_x0_x2);\
			__TEST11__(0, 3, s5_x0_x3);\
			__TEST11__(1, 2, s5_x1_x2);\
			__TEST11__(1, 3, s5_x1_x3);\
			__TEST11__(2, 3, s5_x2_x3);
__TEST__(15)
#if 0
(-15.68*X3^2)-57.28*X2*X3-26.88*X1*X3+51.52*X0*X3-23.744*X3-25.92*X2^2
-51.52*X1*X2+66.24*X0*X2-30.528*X2-11.52*X1^2+61.12*X0*X1
-20.352*X1-42.32*X0^2+39.008*X0+0.0112
#endif
		}
		TEST_METHOD(TestMethod2)
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
			const auto s5 = hypot(s4*s4, 1.0 - s4*s4);
			//const auto s52 = fmod(s41*s41, 1.0 - s41*s41);
			//const auto s51 = s52.chainRule(s4, mp_size_t<4>());
			constexpr const double s5_value = 1.130377777559343;
			constexpr const double s5_x0 = -5.379924588689765;
			constexpr const double s5_x1 = +2.806917176707705;
			constexpr const double s5_x2 = +4.210375765061556;
			constexpr const double s5_x3 = +3.274736706158988;
			constexpr const double s5_x0_x0 = +14.06721719559887;
			constexpr const double s5_x0_x1 = -17.01793298179438;
			constexpr const double s5_x0_x2 = -22.01825300180694;
			constexpr const double s5_x0_x3 = -17.12530789029429;

			constexpr const double s5_x1_x1 = +3.829261391618601;
			constexpr const double s5_x1_x2 = +13.82688182211221;
			constexpr const double s5_x1_x3 = +8.93494324711007;

			constexpr const double s5_x2_x2 = +8.61583813114185;
			constexpr const double s5_x2_x3 = +15.74151251792152;

			constexpr const double s5_x3_x3 = 5.212050227480872;
__TEST__(15)
#if 0
5.212050227480872*X3^2
+15.74151251792152*X2*X3
+8.93494324711007*X1*X3
-17.12530789029429*X0*X3
+3.274736706158988*X3
+8.61583813114185*X2^2
+13.82688182211221*X1*X2
-22.01825300180694*X0*X2
+4.210375765061556*X2
+3.829261391618601*X1^2
-17.01793298179438*X0*X1
+2.806917176707705*X1
+14.06721719559887*X0^2
-5.379924588689765*X0
+1.130377777559343
#endif
		}
		TEST_METHOD(TestMethod3)
		{	using namespace taylor;
			using namespace boost::mp11;

			constexpr const std::size_t MAX = 3;
			const auto s0 = ctaylor<makeIndependent<0>, MAX>(1.2, false);
				/// create an independent variable for x1 (this is what the unused boolean is for)
			const auto s1 = ctaylor<makeIndependent<1>, MAX>(1.3, false);
			const auto s2 = ctaylor<makeIndependent<2>, MAX>(1.4, false);
			const auto s3 = ctaylor<makeIndependent<3>, MAX>(1.5, false);
				/// some calculation
			const auto s4 = -s0 + s1 - s2 + s1*s2 - s0*s1 + s2*s3;
			const auto s5 = atan(1/s4-s4*s4);

			constexpr const double s5_x3_x3_x3 = 18.79819770369593;//*X3^3
			constexpr const double s5_x2_x3_x3 = +76.3608003892592;//*X2*X3^2
			constexpr const double s5_x1_x3_x3 = +48.33822266664668;//*X1*X3^2
			constexpr const double s5_x0_x3_x3 = -92.64826011107282;//*X0*X3^2

			constexpr const double s5_x2_x2_x3 = +98.17817192904754;//*X2^2*X3
			constexpr const double s5_x1_x2_x3 = +131.4547244371999;//*X1*X2*X3
			constexpr const double s5_x0_x2_x3 = -244.5690779251623;//*X0*X2*X3
			constexpr const double s5_x1_x1_x3 = +41.43276228569716;//*X1^2*X3
			constexpr const double s5_x0_x1_x3 = -162.6790551511282;//*X0*X1*X3
			constexpr const double s5_x0_x0_x3 = +152.2078558967624;//*X0^2*X3
			constexpr const double s5_x2_x2_x2 = +39.95302077549369;//*X2^3
			constexpr const double s5_x1_x2_x2 = +84.86049833721629;//*X1*X2^2
			constexpr const double s5_x0_x2_x2 = -153.1532463060591;//*X0*X2^2
			constexpr const double s5_x1_x1_x2 = +56.57366555814418;//*X1^2*X2
			constexpr const double s5_x0_x1_x2 = -215.4894799767113;//*X0*X1*X2
			constexpr const double s5_x0_x0_x2 = +195.6958147244089;//*X0^2*X2
			constexpr const double s5_x1_x1_x1 = +11.83793208162775;//*X1^3
			constexpr const double s5_x0_x1_x1 = -71.37108066017889;//*X0*X1^2
			constexpr const double s5_x0_x0_x1 = +136.7945712653429;//*X0^2*X1
			constexpr const double s5_x0_x0_x0 = -83.3519210863223;//*X0^3
			
			constexpr const double s5_x3_x3 = +2.697426472502412;//*X3^2
			constexpr const double s5_x2_x3 = +4.020913604699608;//*X2*X3
			constexpr const double s5_x1_x3 = +4.624159667146992;//*X1*X3
			constexpr const double s5_x0_x3 = -8.862972695365068;//*X0*X3
			constexpr const double s5_x2_x2 = +4.459011107606028;//*X2^2
			constexpr const double s5_x1_x2 = +3.030022247453823;//*X1*X2
			constexpr const double s5_x0_x2 = -11.39525060832651;//*X0*X2
			constexpr const double s5_x1_x1 = +1.981782714491568;//*X1^2
			constexpr const double s5_x0_x1 = -4.681507842863463;//*X0*X1
			constexpr const double s5_x0_x0 = +7.280298999764164;//*X0^2
			
			constexpr const double s5_x3 = -4.081456254429233;//*X3
			constexpr const double s5_x2 = -5.247586612837586;//*X2
			constexpr const double s5_x1 = -3.498391075225057;//*X1
			constexpr const double s5_x0 = +6.705249560848027;//*X0
			constexpr const double s5_value = -0.178290309737224;//

#define __TEST_3__\
			__TEST__(35)\
			__TEST3__(0, s5_x0_x0_x0);\
			__TEST3__(1, s5_x1_x1_x1);\
			__TEST3__(2, s5_x2_x2_x2);\
			__TEST3__(3, s5_x3_x3_x3);\
			__TEST21__(0, 1, s5_x0_x0_x1);\
			__TEST21__(0, 2, s5_x0_x0_x2);\
			__TEST21__(0, 3, s5_x0_x0_x3);\
			__TEST21__(1, 2, s5_x1_x1_x2);\
			__TEST21__(1, 3, s5_x1_x1_x3);\
			__TEST21__(2, 3, s5_x2_x2_x3);\
			__TEST12__(0, 1, s5_x0_x1_x1);\
			__TEST12__(0, 2, s5_x0_x2_x2);\
			__TEST12__(0, 3, s5_x0_x3_x3);\
			__TEST12__(1, 2, s5_x1_x2_x2);\
			__TEST12__(1, 3, s5_x1_x3_x3);\
			__TEST12__(2, 3, s5_x2_x3_x3);\
			__TEST111__(0, 1, 2, s5_x0_x1_x2);\
			__TEST111__(1, 2, 3, s5_x1_x2_x3);\
			__TEST111__(0, 2, 3, s5_x0_x2_x3);\

#define __TEST3__(ID, s5_x0_x0_x0) Assert::IsTrue(std::fabs(s5.getDer(mp_list<pair<mp_size_t<ID>, mp_size_t<3> > >()) - s5_x0_x0_x0*taylor::implementation::factorial<3>::value) < std::fabs(s5_x0_x0_x0*taylor::implementation::factorial<3>::value*epsilon), #s5_x0_x0_x0 L" test")
#define __TEST111__(ID0, ID1, ID2, s5_x0_x1_x2) Assert::IsTrue(std::fabs(s5.getDer(mp_list<pair<mp_size_t<ID0>, mp_size_t<1> >, pair<mp_size_t<ID1>, mp_size_t<1> >, pair<mp_size_t<ID2>, mp_size_t<1> > >()) - s5_x0_x1_x2) < std::fabs(s5_x0_x1_x2*epsilon), #s5_x0_x1_x2 L" test")
#define __TEST21__(ID0, ID1, s5_x0_x0_x1) Assert::IsTrue(std::fabs(s5.getDer(mp_list<pair<mp_size_t<ID0>, mp_size_t<2> >, pair<mp_size_t<ID1>, mp_size_t<1> > >()) - s5_x0_x0_x1*taylor::implementation::factorial<2>::value) < std::fabs(s5_x0_x0_x1*taylor::implementation::factorial<2>::value*epsilon), #s5_x0_x0_x1 L" test")
#define __TEST12__(ID0, ID1, s5_x0_x0_x1) Assert::IsTrue(std::fabs(s5.getDer(mp_list<pair<mp_size_t<ID0>, mp_size_t<1> >, pair<mp_size_t<ID1>, mp_size_t<2> > >()) - s5_x0_x0_x1*taylor::implementation::factorial<2>::value) < std::fabs(s5_x0_x0_x1*taylor::implementation::factorial<2>::value*epsilon), #s5_x0_x0_x1 L" test")
__TEST_3__
#if 0
//atan(1/x4-x4*x4)
18.79819770369593*X3^3
+76.3608003892592*X2*X3^2
+48.33822266664668*X1*X3^2
-92.64826011107282*X0*X3^2
+2.697426472502412*X3^2
+98.17817192904754*X2^2*X3
+131.4547244371999*X1*X2*X3
-244.5690779251623*X0*X2*X3
+4.020913604699608*X2*X3
+41.43276228569716*X1^2*X3
-162.6790551511282*X0*X1*X3
+4.624159667146992*X1*X3
+152.2078558967624*X0^2*X3
-8.862972695365068*X0*X3
-4.081456254429233*X3
+39.95302077549369*X2^3
+84.86049833721629*X1*X2^2
-153.1532463060591*X0*X2^2
+4.459011107606028*X2^2
+56.57366555814418*X1^2*X2
-215.4894799767113*X0*X1*X2
+3.030022247453823*X1*X2
+195.6958147244089*X0^2*X2
-11.39525060832651*X0*X2
-5.247586612837586*X2
+11.83793208162775*X1^3
-71.37108066017889*X0*X1^2
+1.981782714491568*X1^2
+136.7945712653429*X0^2*X1
-4.681507842863463*X0*X1
-3.498391075225057*X1
-83.3519210863223*X0^3
+7.280298999764164*X0^2
+6.705249560848027*X0
-0.178290309737224
#endif
		}
	};
}
