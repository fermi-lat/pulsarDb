/** \file test_pulsarDb.cxx
    \brief Test code for pulsarDb package.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "pulsarDb/AbsoluteTime.h"
#include "pulsarDb/GlastTime.h"
#include "pulsarDb/PulsarDb.h"
#include "pulsarDb/TimeSystemTime.h"
#include "pulsarDb/TimingModel.h"

#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"
#include "st_facilities/Env.h"

#include "tip/FileSummary.h"
#include "tip/IFileSvc.h"
#include "tip/Table.h"

using namespace tip;

/** \class ErrorMsg
    \brief Trivial helper class for formatting output.
*/
class ErrorMsg {
  public:
    ErrorMsg(const std::string & method) { std::cerr << "PulsarDbTest::" << method << ": "; s_status = 1; }

    static int getStatus() { return s_status; }

  private:
    static int s_status;
};

const ErrorMsg & operator <<(const ErrorMsg & msg, std::ostream & (*op)(std::ostream &)) { (*op)(std::cerr); return msg; }

template <typename T>
const ErrorMsg & operator <<(const ErrorMsg & msg, const T & t) { std::cerr << t; return msg; }

int ErrorMsg::s_status = 0;

/** \class PulsarDbTest
    \brief Test application.
*/
class PulsarDbTest : public st_app::StApp {
  public:
    virtual ~PulsarDbTest() throw() {}

    /// Do all tests.
    virtual void run();

    /// Test recognition of pulsars by alternate names given in the ALTERNATIVE_NAMES extension.
    virtual void testAlternateName();

    /// Test error cases to make sure errors are detected properly.
    virtual void testBadInterval();

    /// Test method which chooses the best ephemeris from several which could be used.
    virtual void testChooser();

    /// Test finding pulsars using their PSR J name.
    virtual void testExplicitName();

    /// Test filtering expressions.
    virtual void testExpression();

    /// Test filtering which doesn't actually narrow the selection.
    virtual void testNoOp();

    /// Test filtering based on a time range (finding ephemerides which overlaps the range.)
    virtual void testTime();

    /// Test AbsoluteTime class.
    virtual void testAbsoluteTime();

    /// Test PulsarEph classes.
    virtual void testPulsarEph();

    /// Test TimingModel class.
    virtual void testTimingModel();

  private:
    std::string m_in_file;
    std::string m_out_file;
};

using namespace pulsarDb;

void PulsarDbTest::run() {

  // Find data directory for this app.
  std::string data_dir = st_facilities::Env::getDataDir("pulsarDb");

  // Find test file.
  m_in_file = st_facilities::Env::appendFileName(data_dir, "groD4-dc2v4.fits");

  // Output file.
  m_out_file = "spud.fits";

  // Successful tests.
  testNoOp();
  testExplicitName();
  testAlternateName();
  testTime();
  testAbsoluteTime();
  testPulsarEph();
  testTimingModel();
  testExpression();
  testChooser();

  // Failures.
  testBadInterval();

  if (0 != ErrorMsg::getStatus()) throw std::runtime_error("PulsarDbTest::run: test failed.");
}

void PulsarDbTest::testAlternateName() {
  std::string method_name = "testAlternateName";

  // Create pulsar database object.
  PulsarDb database(m_in_file);
  // Guess who?
  database.filterName("CRab");

  // Confirm that the correct number of ephemerides were found.
  int num_eph = database.getNumEph();
  if (36 != num_eph)
    ErrorMsg(method_name) << "there are " << num_eph << " ephemerides for the crab, not 36" << std::endl;

  // Filter on the crab's PSR B name (no-op).
  database.filterName("PsR b0531+21");
  num_eph = database.getNumEph();
  if (36 != num_eph)
    ErrorMsg(method_name) << "there are " << num_eph << " ephemerides for the crab, not 36" << std::endl;

  // Write this output to form basis for comparing future tests.
  database.save("crab_db.fits");
}

void PulsarDbTest::testBadInterval() {
  std::string method_name = "testBadInterval";

  // Create pulsar database object.
  PulsarDb database(m_in_file);

  try {
    // Invalid interval, with start time later than stop time.
    database.filterInterval(54500., 54499.);
    ErrorMsg(method_name) << "filterInterval(54500., 54499.) did not throw an exception" << std::endl;
  } catch (const std::exception & x) {
    // This is fine.
  }

}

void PulsarDbTest::testChooser() {
  std::string method_name = "testChooser";

  // Create pulsar database object.
  PulsarDb database(m_in_file);

  std::string pulsar_name = "PSR J0139+5814";

  // Select pulsar we want to use.
  database.filterName(pulsar_name);

  // Confirm that the correct number of ephemerides were found.
  int num_eph = database.getNumEph();
  if (8 != num_eph)
    ErrorMsg(method_name) << "there are " << num_eph << " ephemerides for " << pulsar_name << ", not 8" << std::endl;

  // Write this output to form basis for comparing future tests.
  database.save("chooser_db.fits");

  double pick_time = 54012.5;

  PulsarEphCont eph_cont;

  database.getEph(eph_cont);

  // Test one with no tiebreaking needed.
  const PulsarEph * chosen = &eph_cont.chooseEph(pick_time);
  if (54262. != chosen->epoch())
    ErrorMsg(method_name) << "for time " << pick_time << ", chooser chose ephemeris with EPOCH == " << chosen->epoch() << std::endl;

  // Test one with tiebreaking.
  pick_time = 53545.5;
  chosen = &eph_cont.chooseEph(pick_time);
  if (53891. != chosen->epoch())
    ErrorMsg(method_name) << "for time " << pick_time << ", chooser chose ephemeris with EPOCH == " << chosen->epoch() << std::endl;

  // Test one which is too early.
  pick_time = 53544.5;
  try {
    chosen = &eph_cont.chooseEph(pick_time);
    ErrorMsg(method_name) << "for time " << pick_time << ", chooser chose ephemeris with EPOCH == " << chosen->epoch() << std::endl;
  } catch (const std::runtime_error &) {
    // This is to be expected.
  }

  // Test one which is too late.
  pick_time = 55579.5;
  try {
    chosen = &eph_cont.chooseEph(pick_time);
    ErrorMsg(method_name) << "for time " << pick_time << ", chooser chose ephemeris with EPOCH == " << chosen->epoch() << std::endl;
  } catch (const std::runtime_error &) {
    // This is to be expected.
  }

  // Try one which is too late, but without being strict about validity.
  pick_time = 55579.5;
  try {
    chosen = &eph_cont.chooseEph(pick_time, false);
  } catch (const std::runtime_error &) {
    ErrorMsg(method_name) << "for time " << pick_time << ", chooser did not find an ephemeris even with strict_validity == false" <<
      std::endl;
  }

  // Make a selection which will result in an empty container of ephemerides.
  database.filterName("Aunt Gertrude");
  if (0 != database.getNumEph()) {
    ErrorMsg(method_name) << "What? Aunt Gertrude is a PULSAR???" << std::endl;
  }
  
  database.getEph(eph_cont);

  // Try to choose an ephemeris from the empty set.
  pick_time = 55579.5;
  try {
    chosen = &eph_cont.chooseEph(pick_time);
    ErrorMsg(method_name) << "chooser chose ephemeris from an empty set of candidates." << std::endl;
  } catch (const std::runtime_error &) {
    // This is to be expected.
  }

}

void PulsarDbTest::testExplicitName() {
  std::string method_name = "testExplicitName";

  // Create pulsar database object.
  PulsarDb database(m_in_file);

  // Filter a pulsar known to be present.
  database.filterName("PSr j0323+3944");

  // Confirm that the correct number of ephemerides were found.
  int num_eph = database.getNumEph();
  if (2 != num_eph)
    ErrorMsg(method_name) << "there are " << num_eph << " ephemerides for PSR J0323+3944, not 2" << std::endl;
}

void PulsarDbTest::testExpression() {
  std::string method_name = "testExpression";

  // Create pulsar database object.
  PulsarDb database(m_in_file);

  // Filter using more complex criteria.
  database.filter("F2 != 0.");

  // Confirm that the correct number of ephemerides were found.
  int num_eph = database.getNumEph();
  if (837 != num_eph)
    ErrorMsg(method_name) << "found " << num_eph << " ephemerides with F2 != 0., not 837" << std::endl;

  // Test saving this for basis of comparing future test output.
  database.save("f2not0_db.fits");
}

void PulsarDbTest::testNoOp() {
  std::string method_name = "testExpression";

  // Create pulsar database object.
  PulsarDb database(m_in_file);

  // Make sure the test data has expected size.
  int num_eph = database.getNumEph();
  if (1141 != num_eph)
    ErrorMsg(method_name) << "there are initially " << num_eph << " ephemerides, not 1141" << std::endl;

  // Perform filtering which doesn't exclude any ephemerides.
  database.filterInterval(0., 1.e6);

  // Make sure the test data still has expected size.
  num_eph = database.getNumEph();
  if (1141 != num_eph)
    ErrorMsg(method_name) << "PulsarDbTest::run, after no-op filterInterval there are " << num_eph << " ephemerides, not 1141" << std::endl;

  // Another no-op is to use "any".
  database.filterName("aNy");

  // Make sure the test data still has expected size.
  num_eph = database.getNumEph();
  if (1141 != num_eph)
    ErrorMsg(method_name) << "PulsarDbTest::run, after no-op filterName there are " << num_eph << " ephemerides, not 1141" << std::endl;

  // Another no-op is to use a silly expression
  database.filter("#row<1142");

  // Make sure the test data still has expected size.
  num_eph = database.getNumEph();
  if (1141 != num_eph)
    ErrorMsg(method_name) << "PulsarDbTest::run, after no-op filter there are " << num_eph << " ephemerides, not 1141" << std::endl;

  // Another no-op is to use a blank string.
  database.filter("\t  \t");

  // Make sure the test data still has expected size.
  num_eph = database.getNumEph();
  if (1141 != num_eph)
    ErrorMsg(method_name) << "PulsarDbTest::run, after no-op filter there are " << num_eph << " ephemerides, not 1141" << std::endl;
}

void PulsarDbTest::testTime() {
  std::string method_name = "testTime";

  // Create pulsar database object.
  PulsarDb database(m_in_file);

  // Give a time filter.
  database.filterInterval(53400., 53800.);

  // Make sure the test data has expected size.
  int num_eph = database.getNumEph();
  if (406 != num_eph)
    ErrorMsg(method_name) << "after filterInterval(53400., 53800.) there are " << num_eph << " ephemerides, not 406" << std::endl;

}

void PulsarDbTest::testAbsoluteTime() {
  std::string method_name = "testAbsoluteTime";

  std::cerr.precision(24);

  TtSystemTime tt = 32.184L / 86400.;
  TaiSystemTime tai(tt);
  if (0.L != tai.mjd())
    ErrorMsg(method_name) << "after TaiSystemTime tai(tt), tai.mjd() returned " << tai.mjd() << ", not 0.L." << std::endl;

  tai = TaiSystemTime(-tt.mjd());
  if (-tt.mjd() != tai.mjd())
    ErrorMsg(method_name) << "after tai = TaiSystemTime(-tt.mjd()), tai.mjd() returned " << tai.mjd() << ", not " << -tt.mjd() <<
      std::endl;

  tt = tai;
  if (0.L != tt.mjd())
    ErrorMsg(method_name) << "after tt = tai, tt.mjd() returned " << tt.mjd() << ", not 0.L." << std::endl;

  GlastTtSystemTime gtt;
  tt = gtt;
  if (54101.L != tt.mjd())
    ErrorMsg(method_name) << "after tt = gtt, tt.mjd() returned " << tt.mjd() << ", not 54101.L." << std::endl;

  tai = gtt;
  if (54101.L - 32.184L / 86400.L != tai.mjd())
    ErrorMsg(method_name) << "after tai = gtt, tai.mjd() returned " << tai.mjd() << ", not " << 54101.L - 32.184L / 86400.L <<
      std::endl;

  gtt = GlastTtSystemTime(100.);
  gtt = tt;
  if (0. != gtt.elapsed())
    ErrorMsg(method_name) << "after gtt = tt, gtt.elapsed() returned " << gtt.elapsed() << ", not 0.L." << std::endl;

  gtt = GlastTtSystemTime(100.);
  gtt = tai;
  if (0. != gtt.elapsed())
    ErrorMsg(method_name) << "after gtt = tai, gtt.elapsed() returned " << gtt.elapsed() << ", not 0.L." << std::endl;
}

void PulsarDbTest::testPulsarEph() {
  std::string method_name = "testPulsarEph";

  double epsilon = 1.e-8;

  // Create a frequency ephemeris.
  FrequencyEph f_eph(0., 1., 123.456789, 0.11, 1.125e-2, -2.25e-4, 6.75e-6);

  // Create a period ephemeris.
  // This is a set of values known to be the inverses of the frequency coefficients above.
  PeriodEph p_eph(0., 1., 123.456789, 0.11, 88.8888888888888888888889, 1.777777777777777777777778, 0.0177777777777777777777778);

  // They should agree completely.
  if (fabs(f_eph.epoch() / p_eph.epoch() - 1.) > epsilon)
    ErrorMsg(method_name) << "FrequencyEph and PeriodEph give different values for epoch" << std::endl;

  if (fabs(f_eph.phi0() / p_eph.phi0() - 1.) > epsilon)
    ErrorMsg(method_name) << "FrequencyEph and PeriodEph give different values for phi0" << std::endl;

  if (fabs(f_eph.f0() / p_eph.f0() - 1.) > epsilon)
    ErrorMsg(method_name) << "FrequencyEph and PeriodEph give different values for f0" << std::endl;

  if (fabs(f_eph.f1() / p_eph.f1() - 1.) > epsilon)
    ErrorMsg(method_name) << "FrequencyEph and PeriodEph give different values for f1" << std::endl;

  if (fabs(f_eph.f2() / p_eph.f2() - 1.) > epsilon)
    ErrorMsg(method_name) << "FrequencyEph and PeriodEph give different values for f2" << std::endl;
}

void PulsarDbTest::testTimingModel() {
  std::string method_name = "testPulsarEph";

  double epsilon = 1.e-8;

  // Create a frequency ephemeris.
  FrequencyEph f_eph(0., 1., 123.456789, 0.11, 1.125e-2, -2.25e-4, 6.75e-6);

  TimingModel model;

  double phase = model.calcPhase(f_eph, 223.456789);

  // Result determined independently.
  if (fabs(phase/.235 - 1.) > epsilon)
    ErrorMsg(method_name) << "TimingModel::calcPhase produced phase == " << phase << " not .235" << std::endl;

  // Change ephemeris to produce a noticeable effect.
  FrequencyEph f_eph2(0., 1., 123.4567891234567, .11, 1.125e-2, -2.25e-4, 13.5e-6);
  double ev_time = 223.4567891234567;

  // Test frequency computation.
  FrequencyEph f_eph3 = model.calcEphemeris(f_eph2, ev_time);
  double correct_f0 = 5.625e-2;
  if (fabs(f_eph3.f0() - correct_f0) > epsilon) {
    ErrorMsg(method_name) << "TimingModel::calcEphemeris produced f0 == " << f_eph3.f0() << " not " << correct_f0 << std::endl;
  }
  
  double correct_f1 = 11.25e-4;
  if (fabs(f_eph3.f1() - correct_f1) > epsilon) {
    ErrorMsg(method_name) << "TimingModel::calcEphemeris produced f1 == " << f_eph3.f1() << " not " << correct_f1 << std::endl;
  }
  
  double correct_f2 = 13.5e-6;
  if (fabs(f_eph3.f2() - correct_f2) > epsilon) {
    ErrorMsg(method_name) << "TimingModel::calcEphemeris produced f2 == " << f_eph3.f2() << " not " << correct_f2 << std::endl;
  }

  double pdot_t = model.calcPdotCorr(f_eph2, ev_time);
  double correct_t = 323.4567891234567;

  // For this test, time difference between these two values must be << 1.e-6. (1 microsecond.)
  epsilon = 1.e-6;
  if (fabs(pdot_t - correct_t) > epsilon) {
    ErrorMsg(method_name) << "TimingModel::calcPdotCorr produced pdot-corrected time == " << pdot_t << " not " <<
      correct_t << std::endl;
  }

}

st_app::StAppFactory<PulsarDbTest> g_factory("gtpulsardb");
