/** \file test_pulsarDb.cxx
    \brief Test code for pulsarDb package.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#include <cstdio>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "pulsarDb/AbsoluteTime.h"
#include "pulsarDb/CanonicalTime.h"
#include "pulsarDb/EphChooser.h"
#include "pulsarDb/GlastTime.h"
#include "pulsarDb/OrbitalEph.h"
#include "pulsarDb/PulsarDb.h"
#include "pulsarDb/PulsarEph.h"
#include "pulsarDb/TextPulsarDb.h"
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

    /// Test appending ephemerides to an existing file.
    virtual void testAppend();

    /// Test TextPulsarDb class.
    virtual void testTextPulsarDb();

    /// Test converter class, period->freq.
    virtual void testPeriodConverter();

    /// Test orbital ephemerides classes.
    virtual void testOrbitalEph();

    /// Test duration class.
    virtual void testDuration();

  private:
    std::string m_data_dir;
    std::string m_in_file;
    std::string m_out_file;
    std::string m_tpl_file;
};

using namespace pulsarDb;

void PulsarDbTest::run() {

  // Find data directory for this app.
  m_data_dir = st_facilities::Env::getDataDir("pulsarDb");

  // Find test file.
  m_in_file = st_facilities::Env::appendFileName(m_data_dir, "groD4-dc2v4.fits");

  // Output file.
  m_out_file = "spud.fits";

  // Find template file.
  m_tpl_file = st_facilities::Env::appendFileName(m_data_dir, "PulsarEph.tpl");

  // Successful tests.
  testNoOp();
  testExplicitName();
  testAlternateName();
  testTime();
  testAbsoluteTime();
  testPulsarEph();
  testTimingModel();
  testAppend();
  testExpression();
  testChooser();
  testTextPulsarDb();
  testPeriodConverter();
  testOrbitalEph();
  testDuration();

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
  remove("crab_db.fits");
  database.save("crab_db.fits", m_tpl_file);
}

void PulsarDbTest::testBadInterval() {
  std::string method_name = "testBadInterval";

  // Create pulsar database object.
  PulsarDb database(m_in_file);

  try {
    // Invalid interval, with start time later than stop time.
    database.filterInterval(54500., 54499.);
    ErrorMsg(method_name) << "filterInterval(54500., 54499.) did not throw an exception" << std::endl;
  } catch (const std::exception &) {
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
  remove("chooser_db.fits");
  database.save("chooser_db.fits", m_tpl_file);

  TdbTime pick_time = 54012.5;

  PulsarEphCont eph_cont;

  database.getEph(eph_cont);

  EphChooser chooser;

  // Test one with no tiebreaking needed.
  const PulsarEph * chosen = &chooser.choose(eph_cont, pick_time);
  if (TdbTime(54262.) != chosen->epoch())
    ErrorMsg(method_name) << "for time " << pick_time << ", chooser chose ephemeris with EPOCH == " << chosen->epoch() << std::endl;

  // Test one with tiebreaking.
  pick_time.setMjd(53545.5);
  chosen = &chooser.choose(eph_cont, pick_time);
  if (TdbTime(53891.) != chosen->epoch())
    ErrorMsg(method_name) << "for time " << pick_time << ", chooser chose ephemeris with EPOCH == " << chosen->epoch() << std::endl;

  // Test one which is too early.
  pick_time.setMjd(53544.5);
  try {
    chosen = &chooser.choose(eph_cont, pick_time);
    ErrorMsg(method_name) << "for time " << pick_time << ", chooser chose ephemeris with EPOCH == " << chosen->epoch() << std::endl;
  } catch (const std::runtime_error &) {
    // This is to be expected.
  }

  // Test one which is too late.
  pick_time.setMjd(55579.5);
  try {
    chosen = &chooser.choose(eph_cont, pick_time);
    ErrorMsg(method_name) << "for time " << pick_time << ", chooser chose ephemeris with EPOCH == " << chosen->epoch() << std::endl;
  } catch (const std::runtime_error &) {
    // This is to be expected.
  }

  // Try one which is too late, but without being strict about validity.
  pick_time.setMjd(55579.5);
  try {
    chosen = &(SloppyEphChooser().choose(eph_cont, pick_time));
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
  pick_time.setMjd(55579.5);
  try {
    chosen = &chooser.choose(eph_cont, pick_time);
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
  remove("f2not0_db.fits");
  database.save("f2not0_db.fits", m_tpl_file);
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

  TtTime tt = 32.184L / 86400.;
  TaiTime tai(tt);
  if (0.L != tai.mjd())
    ErrorMsg(method_name) << "after TaiTime tai(tt), tai.mjd() returned " << tai.mjd() << ", not 0.L." << std::endl;

  tai = TaiTime(-tt.mjd());
  if (-tt.mjd() != tai.mjd())
    ErrorMsg(method_name) << "after tai = TaiTime(-tt.mjd()), tai.mjd() returned " << tai.mjd() << ", not " << -tt.mjd() <<
      std::endl;

  tt = tai;
  if (0.L != tt.mjd())
    ErrorMsg(method_name) << "after tt = tai, tt.mjd() returned " << tt.mjd() << ", not 0.L." << std::endl;

  GlastTtTime gtt;
  tt = gtt;
  if (54101.L != tt.mjd())
    ErrorMsg(method_name) << "after tt = gtt, tt.mjd() returned " << tt.mjd() << ", not 54101.L." << std::endl;

  tai = gtt;
  if (54101.L - 32.184L / 86400.L != tai.mjd())
    ErrorMsg(method_name) << "after tai = gtt, tai.mjd() returned " << tai.mjd() << ", not " << 54101.L - 32.184L / 86400.L <<
      std::endl;

  gtt = GlastTtTime(100.);
  gtt = tt;
  if (0. != gtt.elapsed())
    ErrorMsg(method_name) << "after gtt = tt, gtt.elapsed() returned " << gtt.elapsed() << ", not 0.L." << std::endl;

  gtt = GlastTtTime(100.);
  gtt = tai;
  if (0. != gtt.elapsed())
    ErrorMsg(method_name) << "after gtt = tai, gtt.elapsed() returned " << gtt.elapsed() << ", not 0.L." << std::endl;

  GlastTtTime gtt1(100.);
  GlastTtTime gtt2(200.);
  if (gtt1.elapsed() == gtt2.elapsed())
    ErrorMsg(method_name) << "After initializing gtt1 and gtt2 with different values, elapsed() returned same value." << std::endl;

  AbsoluteTime * abs_ref1(&gtt1);
  AbsoluteTime * abs_ref2(&gtt2);
  *abs_ref1 = *abs_ref2;
  if (gtt1.elapsed() != gtt2.elapsed())
    ErrorMsg(method_name) << "After *abs_ref1 = *abs_ref2, gtt1.elapsed() returned " << gtt1.elapsed() << ", not " <<
      gtt2.elapsed() << " as expected." << std::endl;

  TaiTime tai1(GlastTtTime(100.));
  long double expected = 54101.L - 32.184L / 86400.L + 100.L / 86400.L;
  if (expected != tai1.mjd())
    ErrorMsg(method_name) << "After creating tai1 from GlastTtTime(100.), tai1.mjd() returned " << tai1.mjd() << ", not " <<
      expected << ", as expected." << std::endl;

  TaiTime tai2(GlastTtTime(200.));
  expected = 54101.L - 32.184L / 86400.L + 200.L / 86400.L;
  if (expected != tai2.mjd())
    ErrorMsg(method_name) << "After creating tai2 from GlastTtTime(200.), tai2.mjd() returned " << tai2.mjd() << ", not " <<
      expected << ", as expected." << std::endl;

  abs_ref1 = &tai1;
  // Note: abs_ref2 still -> gtt2.
  *abs_ref1 = *abs_ref2;
  if (tai1.mjd() != tai2.mjd())
    ErrorMsg(method_name) << "After *abs_ref1 = *abs_ref2, tai1.mjd() returned " << tai1.mjd() << ", not " <<
      tai2.mjd() << " as expected." << std::endl;
}

void PulsarDbTest::testPulsarEph() {
  std::string method_name = "testPulsarEph";

  double epsilon = 1.e-8;

  // Create a frequency ephemeris.
  FrequencyEph f_eph(GlastTtTime(0.), GlastTtTime(1.), GlastTtTime(123.456789), 0.11, 1.125e-2, -2.25e-4, 6.75e-6);

  // Create a period ephemeris.
  // This is a set of values known to be the inverses of the frequency coefficients above.
  PeriodEph p_eph(GlastTtTime(0.), GlastTtTime(1.), GlastTtTime(123.456789), 0.11, 88.8888888888888888888889,
    1.777777777777777777777778, 0.0177777777777777777777778);

  // Create a database ephemeris.
  TimingModel model;
  DatabaseEph db_eph(GlastTtTime(0.), GlastTtTime(1.), GlastTtTime(123.456789), GlastTtTime(223.456789),
    1.125e-2, -2.25e-4, 6.75e-6);

  // All three should agree completely.
  // First, compare frequency & period.
  const long double nano_sec = 1.e-9;
  if (Duration(nano_sec, UnitSec) < (f_eph.epoch() - p_eph.epoch()))
    ErrorMsg(method_name) << "FrequencyEph and PeriodEph give different values for epoch" << std::endl;

  if (fabs(f_eph.phi0() / p_eph.phi0() - 1.) > epsilon)
    ErrorMsg(method_name) << "FrequencyEph and PeriodEph give different values for phi0" << std::endl;

  if (fabs(f_eph.f0() / p_eph.f0() - 1.) > epsilon)
    ErrorMsg(method_name) << "FrequencyEph and PeriodEph give different values for f0" << std::endl;

  if (fabs(f_eph.f1() / p_eph.f1() - 1.) > epsilon)
    ErrorMsg(method_name) << "FrequencyEph and PeriodEph give different values for f1" << std::endl;

  if (fabs(f_eph.f2() / p_eph.f2() - 1.) > epsilon)
    ErrorMsg(method_name) << "FrequencyEph and PeriodEph give different values for f2" << std::endl;

  // Next, compare frequency & database.
  if (Duration(nano_sec, UnitSec) < (f_eph.epoch() - db_eph.epoch()))
    ErrorMsg(method_name) << "FrequencyEph and DatabaseEph give different values for epoch" << std::endl;

  // Correct phi0 is 1. - (1.125e-2 * 100 - 2.25e-4*10000 / 2. + 6.75e-6 * 1000000 / 6.) = -.125 -> .875
  if (fabs(db_eph.phi0() / .875 - 1.) > epsilon)
    ErrorMsg(method_name) << "FrequencyEph and DatabaseEph give different values for phi0" << std::endl;

  if (fabs(f_eph.f0() / db_eph.f0() - 1.) > epsilon)
    ErrorMsg(method_name) << "FrequencyEph and DatabaseEph give different values for f0" << std::endl;

  if (fabs(f_eph.f1() / db_eph.f1() - 1.) > epsilon)
    ErrorMsg(method_name) << "FrequencyEph and DatabaseEph give different values for f1" << std::endl;

  if (fabs(f_eph.f2() / db_eph.f2() - 1.) > epsilon)
    ErrorMsg(method_name) << "FrequencyEph and DatabaseEph give different values for f2" << std::endl;

}

void PulsarDbTest::testTimingModel() {
  std::string method_name = "testTimingModel";

  double epsilon = 1.e-8;

  // Create a frequency ephemeris.
  FrequencyEph f_eph(GlastTtTime(0.), GlastTtTime(1.), GlastTtTime(123.456789), 0.11, 1.125e-2, -2.25e-4, 6.75e-6);

  TimingModel model;

  double phase = model.calcPhase(f_eph, GlastTtTime(223.456789));

  // Result determined independently.
  if (fabs(phase/.235 - 1.) > epsilon)
    ErrorMsg(method_name) << "TimingModel::calcPhase produced phase == " << phase << " not .235" << std::endl;

  // Change ephemeris to produce a noticeable effect.
  FrequencyEph f_eph2(GlastTtTime(0.), GlastTtTime(1.), GlastTtTime(123.4567891234567), .11, 1.125e-2, -2.25e-4, 13.5e-6);
  GlastTtTime ev_time = 223.4567891234567;

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

  model.correctPdot(f_eph2, ev_time);
  double pdot_t = ev_time.elapsed();
  double correct_t = 323.4567891234567;

  // For this test, time difference between these two values must be << 1.e-6. (1 microsecond.)
  epsilon = 1.e-6;
  if (fabs(pdot_t - correct_t) > epsilon) {
    ErrorMsg(method_name) << "TimingModel::calcPdotCorr produced pdot-corrected time == " << pdot_t << " not " <<
      correct_t << std::endl;
  }

}

void PulsarDbTest::testAppend() {
  std::string method_name = "testAppend";

  PulsarDb database(st_facilities::Env::appendFileName(m_data_dir, "groD4-dc2v4.fits"));
  remove("groD4-dc2v4_twice.fits");
  database.save("groD4-dc2v4_twice.fits", m_tpl_file);
  database.save("groD4-dc2v4_twice.fits", m_tpl_file);
}

void PulsarDbTest::testTextPulsarDb() {
  std::string method_name = "testTextPulsarDb";

  // Ingest one of each type of table.
  TextPulsarDb text_psrdb_spin(st_facilities::Env::appendFileName(m_data_dir, "psrdb_spin.txt"), m_tpl_file);
  TextPulsarDb text_psrdb_binary(st_facilities::Env::appendFileName(m_data_dir, "psrdb_binary.txt"), m_tpl_file);
  TextPulsarDb text_psrdb_obs(st_facilities::Env::appendFileName(m_data_dir, "psrdb_obs.txt"), m_tpl_file);
  TextPulsarDb text_psrdb_name(st_facilities::Env::appendFileName(m_data_dir, "psrdb_name.txt"), m_tpl_file);

  // Save all tables into one FITS file.
  std::string filename("psrdb_all.fits");
  remove(filename.c_str());
  text_psrdb_binary.save(filename, m_tpl_file);
  text_psrdb_spin.save(filename, m_tpl_file);
  text_psrdb_obs.save(filename, m_tpl_file);
  text_psrdb_name.save(filename, m_tpl_file);

  // Copy an existing FITS database.
  PulsarDb fits_psrdb("crab_db.fits");
  filename = "psrdb_append.fits";
  remove(filename.c_str());
  fits_psrdb.save(filename, m_tpl_file);

  // Append all tables into the FITS database.
  text_psrdb_name.save(filename, m_tpl_file);
  text_psrdb_obs.save(filename, m_tpl_file);
  text_psrdb_binary.save(filename, m_tpl_file);
  text_psrdb_spin.save(filename, m_tpl_file);
}

void PulsarDbTest::testPeriodConverter() {
  std::string method_name = "testPeriodConverter";

  double p0[] = { 0.033, 0.089, 0.237, 0.197, 0.102};
  double p1[] = { 422.e-15, 2.0e-15, 11.4e-15, 5.8e-15, 92.2e-15 };
  double p2[] = { 1E-19, 2E-20, 2E-20, 32E-20, 1E-23 };
  double epoch[] = {51900., 51900., 51900., 51900., 51900. };
  double phi0[] = { 0.02, .03, 0.21, 0.12, 0.11 };

  for (size_t ii = 0; ii < 5; ++ii) {
    GlastTdbTime dummy(0.);
    std::cout.precision(15);
    PeriodEph period(dummy, dummy, dummy, 0., p0[ii], p1[ii], p2[ii]);
    std::cout << period.f0() << " " << period.f1() << " " << period.f2() << " ";
    std::cout << int(epoch[ii] - 1) << " " << 1. - p0[ii] * phi0[ii] / 86400. << " " << epoch[ii] << " 0. ";
    std::cout << "54100 54200 GLAT F \"JPL DE200\"" << std::endl;
  }
}

void PulsarDbTest::testOrbitalEph() {
  std::string method_name = "testOrbitalEph";

  // Binary parameters: (PB, PBDOT, A1, XDOT, ECC, ECCDOT, OM, OMDOT, T0, GAMMA, SHAPIRO_R, SHAPIRO_S)
  double binary_par[] = {
    27906.980897, -2.43e-12, 2.3417598, 0.0, 0.61713101, 0.0, 220.142729, 4.22662, 45888.1172487, 0.004295, 0.0, 0.0
  };

  OrbitalEph eph1(binary_par);

  // MJD's: { {original-MJD, modulated-MJD}, ... }
  long double mjd_test_values[][2] = {
    { 4.59881172486599971307e+04L, 4.59881172823346967320e+04L },
    { 4.59881519708822161192e+04L, 4.59881520057480382526e+04L },
    { 4.59881866931044423836e+04L, 4.59881867233097334768e+04L },
    { 4.59882214153266613721e+04L, 4.59882214303721880313e+04L },
    { 4.59882561375488876365e+04L, 4.59882561254377882811e+04L },
    { 4.59882908597711066250e+04L, 4.59882908532530402717e+04L },
    { 4.59883255819933328894e+04L, 4.59883255877721779576e+04L },
    { 4.59883603042155518779e+04L, 4.59883603213319299883e+04L },
    { 4.59883950264377781423e+04L, 4.59883950526724878252e+04L },
    { 4.59884297486599971307e+04L, 4.59884297811300504257e+04L },
    { 4.59884644708822161192e+04L, 4.59884645058955188297e+04L }
  };

  // Permitted difference is 100 nanoseconds.
  long double delta = 100.L * 1.e-9 / 86400.;

  TimingModel model;
  for (size_t ii = 0; ii != sizeof(mjd_test_values)/sizeof(long double[2]); ++ii) {
    TdbTime tdb_mjd(mjd_test_values[ii][0]);
    model.modulateBinary(eph1, tdb_mjd);
    if (fabs(tdb_mjd.mjd() - mjd_test_values[ii][1]) > delta)
      ErrorMsg(method_name) << "Binary modulation of " << mjd_test_values[ii][0] << " was computed to be " << tdb_mjd.mjd() <<
        ", not " << mjd_test_values[ii][1] << ", as expected." << std::endl;
  }

  for (size_t ii = 0; ii != sizeof(mjd_test_values)/sizeof(long double[2]); ++ii) {
    TdbTime tdb_mjd(mjd_test_values[ii][1]);
    model.demodulateBinary(eph1, tdb_mjd);
    if (fabs(tdb_mjd.mjd() - mjd_test_values[ii][0]) > delta)
      ErrorMsg(method_name) << "Binary demodulation of " << mjd_test_values[ii][1] << " was computed to be " << tdb_mjd.mjd() <<
        ", not " << mjd_test_values[ii][0] << ", as expected." << std::endl;
  }

} 

void PulsarDbTest::testDuration() {
  std::string method_name = "testDuration";

  Duration six_days(6., UnitDay);
  if (6. != six_days.day())
    ErrorMsg(method_name) << "After Duration six_days(6., UnitDay), six_days.day() returned " << six_days.day() <<
      ", not 6. as expected." << std::endl;
  if (6. * 86400. != six_days.sec())
    ErrorMsg(method_name) << "After Duration six_days(6., UnitDay), six_days.sec() returned " << six_days.sec() <<
      ", not " << 6. * 86400. << " as expected." << std::endl;

  Duration six_sec(6., UnitSec);
  if (6. / 86400. != six_sec.day())
    ErrorMsg(method_name) << "After Duration six_sec(6., UnitDay), six_sec.day() returned " << six_sec.day() <<
      ", not " << 6. / 86400. << " as expected." << std::endl;
  if (6. != six_sec.sec())
    ErrorMsg(method_name) << "After Duration six_sec(6., UnitDay), six_sec.sec() returned " << six_sec.sec() <<
      ", not " << 6. << " as expected." << std::endl;
} 

st_app::StAppFactory<PulsarDbTest> g_factory("gtpulsardb");
