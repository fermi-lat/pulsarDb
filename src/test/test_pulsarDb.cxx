/** \file test_pulsarDb.cxx
    \brief Test code for pulsarDb package.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "pulsarDb/EphChooser.h"
#include "pulsarDb/EphComputer.h"
#include "pulsarDb/OrbitalEph.h"
#include "pulsarDb/PulsarDb.h"
#include "pulsarDb/PulsarEph.h"
#include "pulsarDb/TextPulsarDb.h"
#include "pulsarDb/TimingModel.h"

#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"
#include "st_facilities/Env.h"

#include "timeSystem/AbsoluteTime.h"
#include "timeSystem/ElapsedTime.h"
#include "timeSystem/TimeConstant.h"
#include "timeSystem/TimeRep.h"

#include "tip/FileSummary.h"
#include "tip/IFileSvc.h"
#include "tip/Table.h"

using namespace timeSystem;
using namespace tip;

static long s_orig_mjdref = 54101;
static long s_current_mjdref = 51910;
//static long double s_current_mjdref = s_orig_mjdref;
static long s_mjd_offset = s_current_mjdref - s_orig_mjdref;

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

    /// Test EphComputer class.
    virtual void testEphComputer();

    /// Test Eph getter in pulsarDb class.
    virtual void testEphGetter();

  protected:
    void testEquality(const std::string & context, const pulsarDb::PulsarEph & eph1, const pulsarDb::PulsarEph & eph2) const;

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
  testPulsarEph();
  testTimingModel();
  testAppend();
  testExpression();
  testChooser();
  testTextPulsarDb();
  testPeriodConverter();
  testOrbitalEph();
  testEphComputer();
  testEphGetter();

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
    database.filterInterval(54500.L, 54499.L);
    ErrorMsg(method_name) << "filterInterval(" << 54500.L << ", " <<
      54499.L << ") did not throw an exception" << std::endl;
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

  MjdRep mjd_tdb("TDB", 54012, .5);
  AbsoluteTime pick_time = mjd_tdb.getTime();
//  TdbTime pick_time = 54012.5L;

  PulsarEphCont eph_cont;

  database.getEph(eph_cont);

  EphChooser chooser;

  // Test one with no tiebreaking needed.
  const PulsarEph * chosen = &chooser.choose(eph_cont, pick_time);
  mjd_tdb.setValue(54262, 0.);
  ElapsedTime tolerance("TDB", Duration(0, 1.e-9)); // 1 nanosecond.
  if (!mjd_tdb.getTime().equivalentTo(chosen->epoch(), tolerance))
    ErrorMsg(method_name) << "for time " << pick_time << ", chooser chose ephemeris with EPOCH == " << chosen->epoch() <<
      ", not " << mjd_tdb.getTime() << " as expected." << std::endl;
//  if (TdbTime(54262.L) != chosen->epoch())
//    ErrorMsg(method_name) << "for time " << pick_time << ", chooser chose ephemeris with EPOCH == " << chosen->epoch() <<
//      ", not " << TdbTime(54262.L) << " as expected." << std::endl;

  // Test one with tiebreaking.
  mjd_tdb.setValue(53545, .5);
  pick_time = mjd_tdb.getTime();
//  pick_time.setMjd(Duration(53545, .5 * SecPerDay()));
  chosen = &chooser.choose(eph_cont, pick_time);
  mjd_tdb.setValue(53891, 0.);
  if (!mjd_tdb.getTime().equivalentTo(chosen->epoch(), tolerance))
    ErrorMsg(method_name) << "for time " << pick_time << ", chooser chose ephemeris with EPOCH == " << chosen->epoch() <<
      ", not " << mjd_tdb.getTime() << " as expected." << std::endl;
//  if (TdbTime(53891.) != chosen->epoch())
//    ErrorMsg(method_name) << "for time " << pick_time << ", chooser chose ephemeris with EPOCH == " << chosen->epoch() <<
//      ", not " << TdbTime(53891.L) << " as expected." << std::endl;

  // Test one which is too early.
  mjd_tdb.setValue(53544, .5);
  pick_time = mjd_tdb.getTime();
//  pick_time.setMjd(Duration(53544, .5 * SecPerDay()));
  try {
    chosen = &chooser.choose(eph_cont, pick_time);
    ErrorMsg(method_name) << "for time " << pick_time << ", chooser chose ephemeris with EPOCH == " << chosen->epoch() << std::endl;
  } catch (const std::runtime_error &) {
    // This is to be expected.
  }

  // Test one which is too late.
  mjd_tdb.setValue(55579, .5);
  pick_time = mjd_tdb.getTime();
//  pick_time.setMjd(Duration(55579, .5 * SecPerDay()));
  try {
    chosen = &chooser.choose(eph_cont, pick_time);
    ErrorMsg(method_name) << "for time " << pick_time << ", chooser chose ephemeris with EPOCH == " << chosen->epoch() << std::endl;
  } catch (const std::runtime_error &) {
    // This is to be expected.
  }

  // Try one which is too late, but without being strict about validity.
  mjd_tdb.setValue(55579, .5);
  pick_time = mjd_tdb.getTime();
//  pick_time.setMjd(Duration(55579, .5 * SecPerDay()));
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
  mjd_tdb.setValue(55579, .5);
  pick_time = mjd_tdb.getTime();
//  pick_time.setMjd(Duration(55579, .5 * SecPerDay()));
  try {
    chosen = &chooser.choose(eph_cont, pick_time);
    ErrorMsg(method_name) << "chooser chose ephemeris from an empty set of candidates." << std::endl;
  } catch (const std::runtime_error &) {
    // This is to be expected.
  }

  // Test choosing OrbitalEph.
  // Get new independent access to database, to keep independent from the tests above.
  PulsarDb database2(m_in_file);

  OrbitalEphCont orbital_cont;
  database2.filterName("PSR J1834-0010");
  database2.getEph(orbital_cont);
  mjd_tdb.setValue(52500, 0.);
  pick_time = mjd_tdb.getTime();
//  pick_time.setMjd(Duration(52500, 0.));
  try {
    const OrbitalEph & orbital_eph = chooser.choose(orbital_cont, pick_time);
    long double expected_mjd = 5.206084100795000e+04;
    if (orbital_eph[T0] != expected_mjd) {
      ErrorMsg(method_name) << "chooser chose orbital ephemeris with time " << orbital_eph.t0() << ", not " <<
        expected_mjd << std::endl;
    }
  } catch (const std::runtime_error & x) {
    ErrorMsg(method_name) << "for time " << pick_time << ", chooser had trouble choosing orbital eph: " << x.what() << std::endl;
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

void PulsarDbTest::testPulsarEph() {
  std::string method_name = "testPulsarEph";

  MetRep glast_tt("TT", 51910, 0., 0.);
  glast_tt.setValue(0.);
  AbsoluteTime since = glast_tt.getTime();
  glast_tt.setValue(1.);
  AbsoluteTime until = glast_tt.getTime();
  glast_tt.setValue(123.456789);
  AbsoluteTime epoch = glast_tt.getTime();

  // Create a frequency ephemeris.
  FrequencyEph f_eph(since, until, epoch, 0.875, 1.125e-2, -2.25e-4, 6.75e-6);
//  FrequencyEph f_eph(GlastTtTime(0.), GlastTtTime(1.), GlastTtTime(123.456789), 0.875, 1.125e-2, -2.25e-4, 6.75e-6);

  // Create a period ephemeris.
  // This is a set of values known to be the inverses of the frequency coefficients above.
  PeriodEph p_eph(since, until, epoch, 0.875, 88.8888888888888888888889,
    1.777777777777777777777778, 0.0177777777777777777777778);

  // First, compare frequency & period.
  testEquality("FrequencyEph and PeriodEph", f_eph, p_eph);
}

void PulsarDbTest::testTimingModel() {
  std::string method_name = "testTimingModel";

  MetRep glast_tt("TT", 51910, 0., 0.);
  glast_tt.setValue(0.);
  AbsoluteTime since = glast_tt.getTime();
  glast_tt.setValue(1.);
  AbsoluteTime until = glast_tt.getTime();
  glast_tt.setValue(123.456789);
  AbsoluteTime epoch = glast_tt.getTime();

  long double epsilon = 1.e-8;

  // Create a frequency ephemeris.
  FrequencyEph f_eph(since, until, epoch, 0.11, 1.125e-2, -2.25e-4, 6.75e-6);

  TimingModel model;

  glast_tt.setValue(223.456789);
  AbsoluteTime pick_time = glast_tt.getTime();
  long double phase = model.calcPulsePhase(f_eph, pick_time);
//  long double phase = model.calcPulsePhase(f_eph, GlastTtTime(223.456789));

  // Result determined independently.
  if (fabs(phase/.235 - 1.) > epsilon)
    ErrorMsg(method_name) << "TimingModel::calcPulsePhase produced phase == " << phase << " not .235" << std::endl;
 
  // Change ephemeris to produce a noticeable effect.
  glast_tt.setValue(123.4567891234567);
  epoch = glast_tt.getTime();
  FrequencyEph f_eph2(since, until, epoch, .11, 1.125e-2, -2.25e-4, 13.5e-6);
//  FrequencyEph f_eph2(GlastTtTime(0.), GlastTtTime(1.), GlastTtTime(123.4567891234567), .11, 1.125e-2, -2.25e-4, 13.5e-6);
  glast_tt.setValue(223.4567891234567);
  AbsoluteTime ev_time = glast_tt.getTime();
//  GlastTtTime ev_time = 223.4567891234567;

  // Test frequency computation.
  FrequencyEph f_eph3 = model.calcEphemeris(f_eph2, ev_time);
  long double correct_f0 = 5.625e-2;
  if (fabs(f_eph3.f0() - correct_f0) > epsilon) {
    ErrorMsg(method_name) << "TimingModel::calcEphemeris produced f0 == " << f_eph3.f0() << " not " << correct_f0 << std::endl;
  }
  
  long double correct_f1 = 11.25e-4;
  if (fabs(f_eph3.f1() - correct_f1) > epsilon) {
    ErrorMsg(method_name) << "TimingModel::calcEphemeris produced f1 == " << f_eph3.f1() << " not " << correct_f1 << std::endl;
  }
  
  long double correct_f2 = 13.5e-6;
  if (fabs(f_eph3.f2() - correct_f2) > epsilon) {
    ErrorMsg(method_name) << "TimingModel::calcEphemeris produced f2 == " << f_eph3.f2() << " not " << correct_f2 << std::endl;
  }

  // correct_epoch == ev_time;
  // Note: epsilon now is a time difference of 10ns
//  if (fabs((f_eph3.epoch() - ev_time).sec()) > epsilon) {
  ElapsedTime tolerance("TT", Duration(0, 1.e-9)); // 1 nanosecond.
  if (!f_eph3.epoch().equivalentTo(ev_time, tolerance)) {
    ErrorMsg(method_name) << "TimingModel::calcEphemeris produced epoch == " << f_eph3.epoch() << " not " << ev_time << std::endl;
  }

  model.cancelPdot(f_eph2, ev_time);
  ev_time.getTime(glast_tt);
  long double pdot_t = glast_tt.getValue();
//  long double pdot_t = ev_time.elapsed();
  long double correct_t = 323.4567891234567;

  // For this test, time difference between these two values must be << 1.e-6. (1 microsecond.)
  epsilon = 1.e-6;
  if (fabs(pdot_t - correct_t) > epsilon) {
    ErrorMsg(method_name) << "TimingModel::calcPdotCorr produced pdot-corrected time == " << pdot_t << " not " <<
      correct_t << std::endl;
  }

  MetRep glast_tdb("TDB", 51910, 0., 123.456789);
  AbsoluteTime t0 = glast_tdb.getTime();
  OrbitalEph o_eph(1000., .2, 0., 0., 0., 0., 0., 0., t0, 0., 0., 0.);
//  OrbitalEph o_eph(1000., .2, 0., 0., 0., 0., 0., 0., GlastTdbTime(123.456789), 0., 0., 0.);
  glast_tdb.setValue(223.456789);
  ev_time = glast_tdb.getTime();
  phase = model.calcOrbitalPhase(o_eph, ev_time);
//  phase = model.calcOrbitalPhase(o_eph, GlastTdbTime(223.456789));

  // Result determined independently.
  if (fabs(phase/.099 - 1.) > epsilon)
    ErrorMsg(method_name) << "TimingModel::calcOrbitalPhase produced phase == " << phase << " not .099" << std::endl;
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

  long double p0[] = { 0.033, 0.089, 0.237, 0.197, 0.102 };
  long double p1[] = { 422.e-15, 2.0e-15, 11.4e-15, 5.8e-15, 92.2e-15 };
  long double p2[] = { 1E-19, 2E-20, 2E-20, 32E-20, 1E-23 };
  long double epoch[] = {51900., 51900., 51900., 51900., 51900. };
  long double phi0[] = { 0.02, .03, 0.21, 0.12, 0.11 };

  for (size_t ii = 0; ii < 5; ++ii) {
    AbsoluteTime dummy("TDB", Duration(0, 0.), Duration(0, 0.));
//    GlastTdbTime dummy(0.);
    std::cout.precision(15);
    PeriodEph period(dummy, dummy, dummy, 0., p0[ii], p1[ii], p2[ii]);
    std::cout << period.f0() << " " << period.f1() << " " << period.f2() << " ";
    std::cout << int(epoch[ii] - 1) << " " << 1. - p0[ii] * phi0[ii] / 86400. << " " << epoch[ii] << " 0. ";
    std::cout << int(54100.L) << " " << int(54200.L) << " GLAT F \"JPL DE200\"" << std::endl;
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
  Duration mjd_test_values[][2] = {
    { Duration(45988, SecPerDay() * .00001172486599971307e+04),
      Duration(45988, SecPerDay() * .00001172823346967320e+04) },
    { Duration(45988, SecPerDay() * .00001519708822161192e+04),
      Duration(45988, SecPerDay() * .00001520057480382526e+04) },
    { Duration(45988, SecPerDay() * .00001866931044423836e+04),
      Duration(45988, SecPerDay() * .00001867233097334768e+04) },
    { Duration(45988, SecPerDay() * .00002214153266613721e+04),
      Duration(45988, SecPerDay() * .00002214303721880313e+04) },
    { Duration(45988, SecPerDay() * .00002561375488876365e+04),
      Duration(45988, SecPerDay() * .00002561254377882811e+04) },
    { Duration(45988, SecPerDay() * .00002908597711066250e+04),
      Duration(45988, SecPerDay() * .00002908532530402717e+04) },
    { Duration(45988, SecPerDay() * .00003255819933328894e+04),
      Duration(45988, SecPerDay() * .00003255877721779576e+04) },
    { Duration(45988, SecPerDay() * .00003603042155518779e+04),
      Duration(45988, SecPerDay() * .00003603213319299883e+04) },
    { Duration(45988, SecPerDay() * .00003950264377781423e+04),
      Duration(45988, SecPerDay() * .00003950526724878252e+04) },
    { Duration(45988, SecPerDay() * .00004297486599971307e+04),
      Duration(45988, SecPerDay() * .00004297811300504257e+04) },
    { Duration(45988, SecPerDay() * .00004644708822161192e+04),
      Duration(45988, SecPerDay() * .00004645058955188297e+04) }
  };

  // Permitted difference is 10 microseconds. This is not really good enough, because long double is different on Windows.
  //long double delta = 10000.L * 1.e-9 / 86400.;
  //long double delta = 100.L * 1.e-9 / 86400.;
  long double delta = 100.L * 1.e-9;
  ElapsedTime tolerance("TDB", Duration(0, delta));

  TimingModel model;
  std::cerr.precision(24);
  for (size_t ii = 0; ii != sizeof(mjd_test_values)/sizeof(Duration[2]); ++ii) {
//    TdbTime tdb_mjd(mjd_test_values[ii][0]);
    AbsoluteTime tdb_mjd("TDB", mjd_test_values[ii][0], Duration(0, 0.));
    AbsoluteTime expected_tdb_mjd("TDB", mjd_test_values[ii][1], Duration(0, 0.));
    model.modulateBinary(eph1, tdb_mjd);
//    if (std::fabs((tdb_mjd.getMjd() - mjd_test_values[ii][1]).sec()) > delta) {
    if (!tdb_mjd.equivalentTo(expected_tdb_mjd, tolerance)) {
//      ErrorMsg(method_name) << "Difference is " << std::fabs((tdb_mjd.getMjd() - mjd_test_values[ii][1]).day()) <<
//        " days." << std::endl;
      ErrorMsg(method_name) << "Binary modulation of " << mjd_test_values[ii][0] << " was computed to be " << tdb_mjd <<
        ", not " << mjd_test_values[ii][1] << ", as expected." << std::endl;
    }
  }

  for (size_t ii = 0; ii != sizeof(mjd_test_values)/sizeof(Duration[2]); ++ii) {
//    TdbTime tdb_mjd(mjd_test_values[ii][1]);
    AbsoluteTime tdb_mjd("TDB", mjd_test_values[ii][1], Duration(0, 0.));
    AbsoluteTime expected_tdb_mjd("TDB", mjd_test_values[ii][0], Duration(0, 0.));
    model.demodulateBinary(eph1, tdb_mjd);
//    if (std::fabs((tdb_mjd.getMjd() - mjd_test_values[ii][0]).sec()) > delta) {
    if (!tdb_mjd.equivalentTo(expected_tdb_mjd, tolerance)) {
//      ErrorMsg(method_name) << "Difference is " << std::fabs((tdb_mjd.getMjd() - mjd_test_values[ii][0]).day()) <<
//        " days." << std::endl;
      ErrorMsg(method_name) << "Binary demodulation of " << mjd_test_values[ii][1] << " was computed to be " << tdb_mjd <<
        ", not " << mjd_test_values[ii][0] << ", as expected." << std::endl;
    }
  }
} 

void PulsarDbTest::testEphComputer() {
  std::string method_name = "testEphComputer";

  // Set up pieces needed for computation.
  TimingModel model;
  EphChooser chooser;

  // Get access to database.
  PulsarDb database(m_in_file);

  // Filter a pulsar known to be present.
  database.filterName("PSr j0323+3944");

  // First perform computations without the computer.
  PulsarEphCont eph_cont;
  database.getEph(eph_cont);

  MetRep glast_tdb("TDB", 54101, 0., 100.);
  AbsoluteTime expected_gtdb = glast_tdb.getTime();
//  GlastTdbTime expected_gtdb(100. - s_mjd_offset * 86400.L);
  double expected_elapsed = glast_tdb.getValue();
  const PulsarEph & eph(chooser.choose(eph_cont, expected_gtdb));
  model.cancelPdot(eph, expected_gtdb);
  long double expected_pulse_phase = model.calcPulsePhase(eph, expected_gtdb);
  FrequencyEph expected_eph(model.calcEphemeris(eph, expected_gtdb));
 
  // Repeat computations using the EphComputer class, and compare results.
  // Create the computer.
  EphComputer computer;

  // Load the data from the database.
  computer.load(database);

  // Test cancelPdot, and compare result to previous result.

  glast_tdb.setValue(100.);
  AbsoluteTime gtdb = glast_tdb.getTime();
//  GlastTdbTime gtdb(100. - s_mjd_offset * 86400.L);
  computer.cancelPdot(gtdb);
//  if (expected_gtdb.elapsed() != gtdb.elapsed())
  if (expected_elapsed != glast_tdb.getValue())
    ErrorMsg(method_name) << "EphComputer::cancelPdot returned elapsed time " << glast_tdb.getValue() << ", not " <<
      expected_elapsed << ", as expected." << std::endl;

  // Test calcPulsePhase, and compare result to previous result.
  long double pulse_phase = computer.calcPulsePhase(expected_gtdb);
  if (expected_pulse_phase != pulse_phase)
    ErrorMsg(method_name) << "EphComputer::calcPulsePhase returned phase " << pulse_phase << ", not " <<
      expected_pulse_phase << ", as expected." << std::endl;

  // Test calcPulsarEph, and compare result to previous result.
  FrequencyEph freq_eph = computer.calcPulsarEph(expected_gtdb);
  testEquality("EphComputer::calcPulsarEph and TimingModel::calcEphemeris", freq_eph, expected_eph);

  // Test binary modulation/demodulation.
  // Get new independent access to database, to keep independent from the tests above.
  PulsarDb database2(m_in_file);

  // Select a particular pulsar.
  database2.filterName("PSR J1834-0010");

  // First perform computations without the computer.
  OrbitalEphCont orbital_eph_cont;
  database2.getEph(orbital_eph_cont);

  const OrbitalEph & orbital_eph(chooser.choose(orbital_eph_cont, gtdb));
  expected_gtdb = gtdb;
  model.modulateBinary(orbital_eph, expected_gtdb);
  
  computer.load(database2);

  computer.modulateBinary(gtdb);
  expected_gtdb.getTime(glast_tdb);
  expected_elapsed = glast_tdb.getValue();
  gtdb.getTime(glast_tdb);
//  if (expected_gtdb.elapsed() != gtdb.elapsed())
  if (expected_elapsed != glast_tdb.getValue())
    ErrorMsg(method_name) << "After EphComputer::modulateBinary, elapsed time was " << glast_tdb.getValue() << ", not " <<
      expected_elapsed << ", as expected." << std::endl;

  expected_gtdb = gtdb;
  model.demodulateBinary(orbital_eph, expected_gtdb);
  
  computer.demodulateBinary(gtdb);
  expected_gtdb.getTime(glast_tdb);
  expected_elapsed = glast_tdb.getValue();
  gtdb.getTime(glast_tdb);
//  if (expected_gtdb.elapsed() != gtdb.elapsed())
  if (expected_elapsed != glast_tdb.getValue())
    ErrorMsg(method_name) << "After EphComputer::demodulateBinary, elapsed time was " << glast_tdb.getValue() << ", not " <<
      expected_elapsed << ", as expected." << std::endl;
}

void PulsarDbTest::testEphGetter() {
  std::string method_name = "testEphGetter";

  PulsarDb database(m_in_file);
  PulsarEphCont pulsar_eph_cont;
  database.getEph(pulsar_eph_cont);
  if (pulsar_eph_cont.size() != size_t(database.getNumEph()))
    ErrorMsg(method_name) << "PulsarDb::getEph(PulsarEphCont &) got " << pulsar_eph_cont.size() << " ephemerides, not " <<
      database.getNumEph() << ", as expected." << std::endl;

  OrbitalEphCont orbital_eph_cont;
  database.getEph(orbital_eph_cont);
  size_t expected_orbital = 19;
  if (orbital_eph_cont.size() != expected_orbital) 
    ErrorMsg(method_name) << "PulsarDb::getEph(OrbitalEphCont &) got " << orbital_eph_cont.size() << " ephemerides, not " <<
      expected_orbital << ", as expected." << std::endl;
}

void PulsarDbTest::testEquality(const std::string & context, const PulsarEph & eph1, const PulsarEph & eph2) const {
  std::string method_name = "testEquality";

  long double epsilon = 1.e-8;
  const long double nano_sec = 1.e-9;
  ElapsedTime tolerance("TT", Duration(0, nano_sec));

//  if (Duration(0, nano_sec) < (eph1.epoch() - eph2.epoch()))
  if (!eph1.epoch().equivalentTo(eph2.epoch(), tolerance))
    ErrorMsg(method_name) << context << " give different values for epoch" << std::endl;

  char * field[] = { "phi0", "f0", "f1" , "f2" };
  double value1[] = { eph1.phi0(), eph1.f0(), eph1.f1(), eph1.f2() };
  double value2[] = { eph2.phi0(), eph2.f0(), eph2.f1(), eph2.f2() };
  for (int ii = 0; ii != sizeof(value1) / sizeof(double); ++ii) {
    if (0. == value1[ii] || 0. == value2[ii]) {
      if (fabs(value1[ii] + value2[ii]) > std::numeric_limits<double>::epsilon())
        ErrorMsg(method_name) << context << " give absolutely different values for " << field[ii] << std::endl;
    } else if (fabs(value1[ii] / value2[ii] - 1.) > epsilon) {
      ErrorMsg(method_name) << context << " give fractionally different values for " << field[ii] << std::endl;
    }
  }
}

st_app::StAppFactory<PulsarDbTest> g_factory("gtpulsardb");
