/// \file TFormulaGradientTests.h
///
/// \brief The file contain unit tests which test the clad-based gradient
///        computations.
///
/// \author Vassil Vassilev <vvasilev@cern.ch>
///
/// \date Oct, 2018
///
/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <Math/MinimizerOptions.h>
#include <TFormula.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TH1.h>

#include "gtest/gtest.h"
#include "gmock/gmock.h"

// Copied from TFileMergerTests.cxx.
// FIXME: Factor out in a new testing library in ROOT.
namespace {
using testing::StartsWith;
using testing::StrEq;
using testing::internal::GetCapturedStderr;
using testing::internal::CaptureStderr;
using testing::internal::RE;
class ExpectedDiagRAII {
public:
   enum ExpectedDiagKind {
      EDK_NoDiag = 0,
      EDK_Info,
      EDK_Warning,
      EDK_Error
   };
private:
   ExpectedDiagKind fDiagKind;
   std::string fExpectedRoutine;
   std::string fExpectedDiag;
   void pop()
   {
      // Diagnostics in ROOT have the format:
      // Error|Warning|Info|...| in <Routine>: free text
      std::string Seen = GetCapturedStderr();

      // Try to reconstruct the precise expected string.
      std::string Expected;
      switch(fDiagKind) {
      default:
         assert (0 && "Unsupported diag kind.");
         break;
      case EDK_NoDiag:
         EXPECT_THAT(Seen, StrEq(""));
         return;
      case EDK_Error:
         Expected = "Error";
         break;
      case EDK_Warning:
         Expected = "Warning";
         break;
      case EDK_Info:
         Expected = "Info";
         break;
      }

      // Check if the Diag kind matches what we saw.
      EXPECT_THAT(Seen, StartsWith(Expected));

      Expected += " in ";
      Expected += "<" + fExpectedRoutine + ">: ";

      // Check if the routine matches what we saw.
      EXPECT_THAT(Seen, StartsWith(Expected));

      Expected += fExpectedDiag;

      // The captured stderr also includes new lines.
      Expected += "\n";

      EXPECT_THAT(Seen, StrEq(Expected));
   }

public:
   ExpectedDiagRAII(ExpectedDiagKind DiagKind): fDiagKind(DiagKind) {
      assert(DiagKind == ExpectedDiagRAII::EDK_NoDiag);
      CaptureStderr();
   }

   ExpectedDiagRAII(ExpectedDiagKind DiagKind, std::string InRoutine,
                    std::string E)
      : fDiagKind(DiagKind), fExpectedRoutine(InRoutine), fExpectedDiag(E) {
      CaptureStderr();
   }
   ~ExpectedDiagRAII() { pop(); }
};
}

#define ROOT_EXPECT_ERROR(expression, where, expected_diag )            \
   {                                                                    \
      ExpectedDiagRAII EE(ExpectedDiagRAII::EDK_Error, where,           \
                          expected_diag);                               \
      expression;                                                       \
   }

#define ROOT_EXPECT_WARNING(expression, where, expected_diag)           \
   {                                                                    \
      ExpectedDiagRAII EE(ExpectedDiagRAII::EDK_Warning, where,         \
                          expected_diag);                               \
      expression;                                                       \
   }

#define ROOT_EXPECT_INFO(expression, where, expected_diag)              \
   {                                                                    \
      ExpectedDiagRAII EE(ExpectedDiagRAII::EDK_Info, where,            \
                          expected_diag);                               \
      expression;                                                       \
   }

#define ROOT_EXPECT_NODIAG(expression)                                  \
   {                                                                    \
      ExpectedDiagRAII EE(ExpectedDiagRAII::EDK_NoDiag);                \
      expression;                                                       \
   }


TEST(TFormulaGradientPar, Sanity)
{
   TFormula f("f", "x*std::sin([0]) - y*std::cos([1])");
   double p[] = {30, 60};
   f.SetParameters(p);
   double x[] = {1, 2};
   TFormula::GradientStorage result(2);
   f.GradientPar(x, result);

   ASSERT_FLOAT_EQ(x[0] * std::cos(30), result[0]);
   ASSERT_FLOAT_EQ(-x[1] * -std::sin(60), result[1]);
}

TEST(TFormulaGradientPar, ResultUpsize)
{
   TFormula f("f", "std::sin([1]) - std::cos([0])");
   double p[] = {60, 30};
   f.SetParameters(p);
   TFormula::GradientStorage result;
   double x[] = {2, 1};

   ASSERT_TRUE(0 == result.size());
   ROOT_EXPECT_WARNING(f.GradientPar(x, result),
                       "TFormula::GradientPar",
                       "The size of gradient result is 0 but 2 is required. Resizing."
                       );

   ASSERT_FLOAT_EQ(std::cos(30), result[1]);
   ASSERT_FLOAT_EQ(std::sin(60), result[0]);
   ASSERT_TRUE(2 == result.size());
}

TEST(TFormulaGradientPar, ResultDownsize)
{
   TFormula f("f", "std::sin([0])");
   double p[] = {60};
   f.SetParameters(p);
   TFormula::GradientStorage result(2);
   double x[] = {1};

   ASSERT_TRUE(2 == result.size());

   ROOT_EXPECT_NODIAG(f.GradientPar(x, result));

   ASSERT_FLOAT_EQ(std::cos(60), result[0]);
   ASSERT_TRUE(2 == result.size());
}

TEST(TFormulaGradientPar, GausCrossCheck)
{
   auto h = new TF1("f1", "gaus");
   // auto h = new TF1("f1", "landau"); -- inheritently does not work. See DIFLAN
   //crystalball, breitwigner, cheb3, bigaus,
   //auto h = new TF1("f1", "");
   double p[] = {3, 1, 2};
   h->SetParameters(p);
   double x[] = {0};
   TFormula::GradientStorage result_clad(3);
   h->GetFormula()->GradientPar(x, result_clad);

   TFormula::GradientStorage result_num(3);
   h->GradientPar(x, result_num.data());

   ASSERT_FLOAT_EQ(result_num[0], result_clad[0]);
   ASSERT_FLOAT_EQ(result_num[1], result_clad[1]);
   ASSERT_FLOAT_EQ(result_num[2], result_clad[2]);
}

// FIXME: See vgvassilev/clad#105
static void InitBreitWigner() {
   static bool FirstCall = true;
   if (FirstCall) {
      std::string clad_breitwigner_pdf_grad = R"code(
namespace custom_derivatives {
   void breitwigner_pdf_grad(double x, double gamma, double x0, double *_result) {
      double _t0 = 1 / (3.1415926535897931 * ((x - x0) * (x - x0) + (gamma / 2.) * (gamma / 2.)));
      double _t1 = _t0 / 2.;
      _result[1UL] += _t1;
      double _t2 = _t0 * -gamma / (2. * 2.);
      double _t3 = 1 * -(gamma / 2.) / ((3.1415926535897931 * ((x - x0) * (x - x0) + (gamma / 2.) * (gamma / 2.))) * (3.1415926535897931 * ((x - x0) * (x - x0) + (gamma / 2.) * (gamma / 2.))));
      double _t4 = _t3 * ((x - x0) * (x - x0) + (gamma / 2.) * (gamma / 2.));
      double _t5 = 3.1415926535897931 * _t3;
      double _t6 = _t5 * (x - x0);
      _result[0UL] += _t6;
      _result[2UL] += -_t6;
      double _t7 = (x - x0) * _t5;
      _result[0UL] += _t7;
      _result[2UL] += -_t7;
      double _t8 = _t5 * (gamma / 2.);
      double _t9 = _t8 / 2.;
      _result[1UL] += _t9;
      double _t10 = _t8 * -gamma / (2. * 2.);
      double _t11 = (gamma / 2.) * _t5;
      double _t12 = _t11 / 2.;
      _result[1UL] += _t12;
      double _t13 = _t11 * -gamma / (2. * 2.);
   }
}
)code";
   gInterpreter->Declare(clad_breitwigner_pdf_grad.c_str());
   FirstCall = false;
   }
}

TEST(TFormulaGradientPar, BreitWignerCrossCheck)
{
   InitBreitWigner();
   auto h = new TF1("f1", "breitwigner");
   double p[] = {3, 1, 2.1};
   h->SetParameters(p);
   double x[] = {0};
   TFormula::GradientStorage result_clad(3);
   TFormula* formula = h->GetFormula();
   formula->GradientPar(x, result_clad);
   TFormula::GradientStorage result_num(3);
   h->GradientPar(x, result_num.data());

   ASSERT_FLOAT_EQ(result_num[0], result_clad[0]);
   ASSERT_FLOAT_EQ(result_num[1], result_clad[1]);
   ASSERT_FLOAT_EQ(result_num[2], result_clad[2]);
}

TEST(TFormulaGradientPar, BreitWignerCrossCheckAccuracyDemo)
{
   InitBreitWigner();
   auto h = new TF1("f1", "breitwigner");
   double p[] = {3, 1, 2};
   h->SetParameters(p);
   double x[] = {0};
   TFormula::GradientStorage result_clad(3);
   TFormula* formula = h->GetFormula();
   formula->GradientPar(x, result_clad);
   TFormula::GradientStorage result_num(3);
   h->GradientPar(x, result_num.data());

   // This is a classical example why clad is better.
   // The gradient with respect to gamma leads to a cancellation when gamma is
   // set to 2. This is not a problem for clad yielding the correct result of 0
   ASSERT_FLOAT_EQ(0, result_clad[2]);

   // However, that is not the case for the numerical approach where we give
   // a small but non-zero result.
   EXPECT_NEAR(0, result_num[2], /*abs_error*/1e-13);
}

// FIXME: Add more: crystalball, cheb3, bigaus?

TEST(TFormulaGradientPar, GetGradFormula)
{
   TFormula f("f", "gaus");
   double p[] = {3, 1, 2};
   f.SetParameters(p);
   ASSERT_TRUE(f.GenerateGradientPar());
   std::string s = f.GetGradientFormula().Data();
   ASSERT_THAT(s, testing::ContainsRegex("void TFormula____id[0-9]*_grad"));
}

// TEST(TFormulaGradientPar, BiGausCrossCheck)
// {
//    auto h = new TF1("f1", "bigaus");
//    // auto h = new TF1("f1", "landau"); -- inheritently does not work. See DIFLAN
//    //crystalball, cheb3, bigaus,
//    //auto h = new TF1("f1", "");
//    double p[] = {1, 2, 3, 4, 5};
//    h->SetParameters(p);
//    double x[] = {0, 1};
//    TFormula::GradientStorage result_clad(5);
//    h->GetFormula()->GradientPar(x, result_clad);

//    TFormula::GradientStorage result_num(5);
//    h->GradientPar(x, result_num.data());

//    ASSERT_FLOAT_EQ(result_num[0], result_clad[0]);
//    ASSERT_FLOAT_EQ(result_num[1], result_clad[1]);
//    ASSERT_FLOAT_EQ(result_num[2], result_clad[2]);
//    ASSERT_FLOAT_EQ(result_num[3], result_clad[3]);
//    ASSERT_FLOAT_EQ(result_num[4], result_clad[4]);
// }

// // PERFORMANCE:
// // Compare the numerical vs AD performance of GradientPar
// TEST(TFormulaGradientPar, GausTFormulaPerf)
// {
//    auto h = new TF1("f1", "gaus");
//    double p[] = {3, 1, 2};
//    h->SetParameters(p);
//    double x[] = {0};
//    TFormula::GradientStorage result_clad(3);

//    // Clad
//    TFormula *f = h->GetFormula();
//    // Note calling GenerateGradientPar is super expensive.
//    f->GenerateGradientPar();
//    auto start = std::chrono::system_clock::now();
//    for (unsigned i = 0; i<10000; ++i)
//       f->GradientPar(x, result_clad.data());
//    auto end = std::chrono::system_clock::now();
//    auto elapsed = end - start;

//    std::cout << "Clad:" << elapsed.count() << '\n';

//    // Numerical
//    TFormula::GradientStorage result_num(3);
//    auto start2 = std::chrono::system_clock::now();
//    for (unsigned i = 0; i<10000; ++i)
//       h->GradientPar(x, result_num.data());
//    auto end2 = std::chrono::system_clock::now();
//    auto elapsed2 = end2 - start2;
//    std::cout << "Num:" << elapsed2.count() << '\n';

//    // For i = 10000, we get: Clad ~622 vs Num ~6431
//    // However, if we include the f->GenerateGradientPar() in the counting we
//    // get: Clad ~69275 vs Num ~7714
// }

TEST(TFormulaGradientPar, GausPerf)
{
   auto f1 = new TF1("f1", "gaus");
   auto h1 = new TH1D("h1", "h1", 100000, -5, 5);
   double p1[] = {1, 0, 1.5};
   f1->SetParameters(p1);
   h1->FillRandom("f1", 100000);

   double p2[] = {100, 1, 3};

   f1->SetParameters(p2);


   auto start = std::chrono::system_clock::now();
   auto r1 = h1->Fit(f1, "S"); // numerical
   auto end = std::chrono::system_clock::now();
   auto elapsed = end - start;
   std::cout << "Num:" << elapsed.count() << '\n';

   f1->SetParameters(p2);

   auto start0 = std::chrono::system_clock::now();
   f1->GetFormula()->GenerateGradientPar(); // makes it om par
   auto end0 = std::chrono::system_clock::now();
   std::chrono::duration<double> elapsed0 = end0 - start0;
   std::cout << "Clad Gen:" << elapsed0.count() << 's\n';

   auto start1 = std::chrono::system_clock::now();
   auto r2 = h1->Fit(f1, "S G"); // clad
   auto end1 = std::chrono::system_clock::now();
   auto elapsed1 = end1 - start1;
   std::cout << "Clad:" << elapsed1.count() << '\n';

   ASSERT_FLOAT_EQ(r1->MinFcnValue(), r2->MinFcnValue());
}

// TEST(TFormulaGradientPar, GausPerfNum)
// {
//    auto f1 = new TF1("f1", "gaus");
//    auto h1 = new TH1D("h1", "h1", 1000000, -5, 5);
//    double p1[] = {1, 0, 1.5};
//    f1->SetParameters(p1);
//    h1->FillRandom("f1", 1000000);

//    double p2[] = {100, 1, 3};

//    f1->SetParameters(p2);
//    ROOT::Math::MinimizerOptions::SetDefaultStrategy(0);
//    auto start = std::chrono::system_clock::now();
//    auto r1 = h1->Fit(f1, "S"); // numerical
//    auto end = std::chrono::system_clock::now();
//    auto elapsed = end - start;
//    std::cout << "Num:" << elapsed.count() << '\n';
// }

// TEST(TFormulaGradientPar, GausPerfClad)
// {
//    auto f2 = new TF1("f2", "gaus");
//    auto h2 = new TH1D("h2", "h2", 1000000, -5, 5);
//    double p1[] = {1, 0, 1.5};
//    f2->SetParameters(p1);
//    h2->FillRandom("f2", 1000000);

//    double p2[] = {100, 1, 3};

//    f2->SetParameters(p2);
//    ROOT::Math::MinimizerOptions::SetDefaultStrategy(0);
//    auto start = std::chrono::system_clock::now();
//    // FIXME: Does not call our code.
//    auto r2 = h2->Fit(f2, "S G"); // clad
//    auto end = std::chrono::system_clock::now();
//    auto elapsed = end - start;
//    std::cout << "Clad" << elapsed.count() << '\n';
// }

// // For GausPerf Clad vs Num: ~7790376 vs ~2304988
