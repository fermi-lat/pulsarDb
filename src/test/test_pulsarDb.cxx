/** \file test_pulsarDb.cxx
    \brief Test code for pulsarDb package.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "facilities/commonUtilities.h"

#include "pulsarDb/EphChooser.h"
#include "pulsarDb/EphComputer.h"
#include "pulsarDb/FrequencyEph.h"
#include "pulsarDb/OrbitalEph.h"
#include "pulsarDb/PdotCanceler.h"
#include "pulsarDb/PeriodEph.h"
#include "pulsarDb/PulsarDb.h"
#include "pulsarDb/PulsarEph.h"
#include "pulsarDb/SimpleDdEph.h"

#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "timeSystem/AbsoluteTime.h"
#include "timeSystem/ElapsedTime.h"
#include "timeSystem/TimeConstant.h"
#include "timeSystem/TimeRep.h"

#include "tip/IFileSvc.h"

using namespace timeSystem;
using namespace pulsarDb;

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

class EphRoutingInfo;

/** \class PulsarDbTest
    \brief Test application.
*/
class PulsarDbTest : public st_app::StApp {
  public:
    virtual ~PulsarDbTest() throw() {}

    /// Do all tests.
    virtual void run();

    /// Test filtering which doesn't actually narrow the selection.
    virtual void testNoOp();

    /// Test finding pulsars using their PSR J name.
    virtual void testExplicitName();

    /// Test recognition of pulsars by alternate names given in the ALTERNATIVE_NAMES extension.
    virtual void testAlternateName();

    /// Test filtering based on a time range (finding ephemerides which overlaps the range.)
    virtual void testTime();

    /// Test error cases to make sure errors are detected properly.
    virtual void testBadInterval();

    /// Test filtering based on a solar system ephemeris name.
    virtual void testSolarEph();

    /// Test filtering expressions.
    virtual void testExpression();

    /// Test appending ephemerides to an existing file.
    virtual void testAppend();

    /// Test ephemerides database in text format.
    virtual void testTextPulsarDb();

    /// Test FrequencyEph class (a subclass of PulsarEph class).
    virtual void testFrequencyEph();

    /// Test PeriodEph classes (a subclass of PulsarEph class).
    virtual void testPeriodEph();

    /// Test SimpleDdEph class (a subclass of OrbitalEph class).
    virtual void testSimpleDdEph();

    /// Test PdotCanceler class.
    virtual void testPdotCanceler();

    /// Test method which chooses the best ephemeris from several which could be used.
    virtual void testChooser();

    /// Test EphComputer class.
    virtual void testEphComputer();

    /// Test Eph getter in pulsarDb class.
    virtual void testEphGetter();

    /// Test support for multiple ephemeris models.
    virtual void testMultipleEphModel();

  private:
    typedef std::list<std::pair<std::string, std::string> > ExtInfoCont;

    std::string m_data_dir;
    std::string m_in_file;
    std::string m_tpl_file;
    std::string m_creator;

    /// Helper method for testMultipleEphModel, to test loading ephemerides from FITS database files.
    void testLoadingFits(const std::string & method_name, PulsarDb & database, const std::string & tpl_file,
      bool load_original, bool expected_to_fail);

    /// Helper method for testMultipleEphModel, to test loading ephemerides from text database files.
    void testLoadingText(const std::string & method_name, PulsarDb & database, const ExtInfoCont & ext_info_cont,
      bool load_original, bool expected_to_fail);

    /// Helper method for testMultipleEphModel, check ephemerides returned by PulsarDb::getEph method.
    void checkEphRouting(const std::string & method_name, const PulsarDb & database,
      const std::map<std::string, EphRoutingInfo> & expected_route_dict) const;
};

void PulsarDbTest::run() {

  // Find data directory for this app.
  m_data_dir = facilities::commonUtilities::getDataPath("pulsarDb");

  // Find test file.
  m_in_file = facilities::commonUtilities::joinPath(m_data_dir, "groD4-dc2v4r2.fits");

  // Find template file.
  m_tpl_file = facilities::commonUtilities::joinPath(m_data_dir, "PulsarDb.tpl");

  // Set a default value for CREATOR header keyword.
  m_creator = "test_pulsarDb";

  // Test filtering of database entries.
  testNoOp();
  testExplicitName();
  testAlternateName();
  testTime();
  testBadInterval();
  testSolarEph();
  testExpression();

  // Test database manipulation.
  testAppend();
  testTextPulsarDb(); // TODO: This test seems to interfere test MultipleEphModel. See below for details.

  // Test ephemeris computation.
  testFrequencyEph();
  testPeriodEph();
  testSimpleDdEph();
  testPdotCanceler();

  // Test ephemeris manipulation.
  testChooser();
  testEphComputer(); // TODO: This test seems to interfere test MultipleEphModel. See below for details.
  testEphGetter(); // TODO: This test seems to interfere test MultipleEphModel. See below for details.
  testMultipleEphModel();
  // TODO: Nail down a mysterious behavior of PulsarDb::loadFits/Text method.
  // Symptoms: testTextPulsarDb, testEphComputer, and testEphGetter seem to interfere testMultipleEphModel.
  //           testMultipleEphModel produces an error if filterRow("#row>0") is not performed on all tables
  //           in the memory FITS file at the end of PulsarDb::loadFits/Text methods, AND if one of the
  //           interferring tests are performed. Otherwise, it doesn't produce an error. Why?
#endif

  if (0 != ErrorMsg::getStatus()) throw std::runtime_error("PulsarDbTest::run: test failed.");
}

void PulsarDbTest::testNoOp() {
  std::string method_name = "testNoOp";

  // Create pulsar database object.
  PulsarDb database(m_tpl_file);
  database.load(m_in_file);

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

void PulsarDbTest::testExplicitName() {
  std::string method_name = "testExplicitName";

  // Create pulsar database object.
  PulsarDb database(m_tpl_file);
  database.load(m_in_file);

  // Filter a pulsar known to be present.
  database.filterName("PSr j0323+3944");

  // Confirm that the correct number of ephemerides were found.
  int num_eph = database.getNumEph();
  if (2 != num_eph)
    ErrorMsg(method_name) << "there are " << num_eph << " ephemerides for PSR J0323+3944, not 2" << std::endl;
}

void PulsarDbTest::testAlternateName() {
  std::string method_name = "testAlternateName";

  // Create pulsar database object.
  PulsarDb database(m_tpl_file);
  database.load(m_in_file);

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
  database.save("crab_db.fits", m_creator);
}

void PulsarDbTest::testTime() {
  std::string method_name = "testTime";

  // Create pulsar database object.
  PulsarDb database(m_tpl_file);
  database.load(m_in_file);

  // Give a time filter.
  database.filterInterval(53400., 53800.);

  // Make sure the test data has expected size.
  int num_eph = database.getNumEph();
  if (406 != num_eph)
    ErrorMsg(method_name) << "after filterInterval(53400., 53800.) there are " << num_eph << " ephemerides, not 406" << std::endl;

}

void PulsarDbTest::testBadInterval() {
  std::string method_name = "testBadInterval";

  // Create pulsar database object.
  PulsarDb database(m_tpl_file);
  database.load(m_in_file);

  try {
    // Invalid interval, with start time later than stop time.
    database.filterInterval(54500., 54499.);
    ErrorMsg(method_name) << "filterInterval(" << 54500. << ", " <<
      54499. << ") did not throw an exception" << std::endl;
  } catch (const std::exception &) {
    // This is fine.
  }

}

void PulsarDbTest::testSolarEph() {
  std::string method_name = "testSolarEph";

  // Create pulsar database object.
  PulsarDb database(m_tpl_file);
  database.load(m_in_file);

  // Give a time filter.
  database.filterSolarEph("JPL DE405");

  // Make sure the test data has expected size (SPIN_PARAMETERS extension).
  int num_eph = database.getNumEph();
  if (5 != num_eph)
    ErrorMsg(method_name) << "after filterSolarEph(\"JPL DE405\") there are " << num_eph << " spin ephemerides, not 5" << std::endl;

  // Make sure the test data has expected size (ORBITAL_PARAMETERS extension).
  num_eph = database.getNumEph(false);
  if (4 != num_eph)
    ErrorMsg(method_name) << "after filterSolarEph(\"JPL DE405\") there are " << num_eph << " orbital ephemerides, not 4" << std::endl;
}

void PulsarDbTest::testExpression() {
  std::string method_name = "testExpression";

  // Create pulsar database object.
  PulsarDb database(m_tpl_file);
  database.load(m_in_file);

  // Filter using more complex criteria.
  database.filter("F2 != 0.");

  // Confirm that the correct number of ephemerides were found.
  int num_eph = database.getNumEph();
  if (837 != num_eph)
    ErrorMsg(method_name) << "found " << num_eph << " ephemerides with F2 != 0., not 837" << std::endl;

  // Test saving this for basis of comparing future test output.
  remove("f2not0_db.fits");
  database.save("f2not0_db.fits", m_creator);
}

void PulsarDbTest::testAppend() {
  std::string method_name = "testAppend";
  PulsarDb database(m_tpl_file);
  database.load(m_in_file);
  database.load(m_in_file);
  remove("twice_db.fits");
  database.save("twice_db.fits", m_creator);
}

void PulsarDbTest::testTextPulsarDb() {
  std::string method_name = "testTextPulsarDb";

  // Ingest one of each type of table.
  PulsarDb database(m_tpl_file);
  database.load(facilities::commonUtilities::joinPath(m_data_dir, "psrdb_spin.txt"));
  database.load(facilities::commonUtilities::joinPath(m_data_dir, "psrdb_binary.txt"));
  database.load(facilities::commonUtilities::joinPath(m_data_dir, "psrdb_obs.txt"));
  database.load(facilities::commonUtilities::joinPath(m_data_dir, "psrdb_name.txt"));

  // Save all tables into one FITS file.
  std::string filename1("psrdb_all.fits");
  remove(filename1.c_str());
  database.save(filename1, m_creator);

  // Copy an existing FITS database.
  PulsarDb fits_psrdb(m_tpl_file);
  fits_psrdb.load("crab_db.fits");

  // Append all tables into the FITS database.
  fits_psrdb.load(filename1);

  // Save all tables into another FITS file.
  std::string filename2 = "psrdb_append.fits";
  remove(filename2.c_str());
  fits_psrdb.save(filename2, m_creator);
}

void PulsarDbTest::testFrequencyEph() {
  std::string method_name = "testFrequencyEph";

  MetRep glast_tt("TT", 51910, 0., 0.);
  glast_tt.setValue(0.);
  AbsoluteTime since(glast_tt);
  glast_tt.setValue(1.);
  AbsoluteTime until(glast_tt);
  glast_tt.setValue(123.456789);
  AbsoluteTime epoch(glast_tt);

  double epsilon = 1.e-8;

  // Create a frequency ephemeris.
  FrequencyEph f_eph("TDB", since, until, epoch, 22., 45., 0.11, 1.125e-2, -2.25e-4, 6.75e-6);

  glast_tt.setValue(223.456789);
  AbsoluteTime pick_time(glast_tt);

  // Test pulse phase computations.
  double phase = f_eph.calcPulsePhase(pick_time);
  if (fabs(phase/.235 - 1.) > epsilon)
    ErrorMsg(method_name) << "FrequencyEph::calcPulsePhase produced phase == " << phase << " not .235" << std::endl;
 
  // Test pulse phase computations, with a non-zero global phase offset.
  phase = f_eph.calcPulsePhase(pick_time, 0.1234);
  if (fabs(phase/.3584 - 1.) > epsilon)
    ErrorMsg(method_name) << "FrequencyEph::calcPulsePhase produced phase == " << phase << " not .3584" << std::endl;
 
  // Change ephemeris to produce a noticeable effect.
  glast_tt.setValue(123.4567891234567);
  epoch = glast_tt;
  FrequencyEph f_eph2("TDB", since, until, epoch, 22., 45., .11, 1.125e-2, -2.25e-4, 13.5e-6);
  glast_tt.setValue(223.4567891234567);
  AbsoluteTime ev_time(glast_tt);

  // Test coordinates.
  std::pair<double, double> computed_ra_dec = f_eph2.calcSkyPosition(ev_time);
  double computed_ra = computed_ra_dec.first;
  double correct_ra = 22.;
  if (fabs(computed_ra - correct_ra) > epsilon) {
    ErrorMsg(method_name) << "FrequencyEph::calcSkyPosition produced ra == " << computed_ra << " not " << correct_ra << std::endl;
  }
  
  double computed_dec = computed_ra_dec.second;
  double correct_dec = 45.;
  if (fabs(computed_dec - correct_dec) > epsilon) {
    ErrorMsg(method_name) << "FrequencyEph::calcSkyPosition produced dec == " << computed_dec << " not " << correct_dec << std::endl;
  }
  
  // Test phase computation.
  double computed_phi0 = f_eph2.calcPulsePhase(ev_time);
  if (fabs(computed_phi0/.36 - 1.) > epsilon)
    ErrorMsg(method_name) << "FrequencyEph::calcPulsePhase produced phi0 == " << computed_phi0 << " not .36" << std::endl;
 
  // Test frequency computation.
  double computed_f0 = f_eph2.calcFrequency(ev_time, 0);
  double correct_f0 = 5.625e-2;
  if (fabs(computed_f0 - correct_f0) > epsilon) {
    ErrorMsg(method_name) << "FrequencyEph::calcFrequency produced f0 == " << computed_f0 << " not " << correct_f0 << std::endl;
  }
  
  double computed_f1 = f_eph2.calcFrequency(ev_time, 1);
  double correct_f1 = 11.25e-4;
  if (fabs(computed_f1 - correct_f1) > epsilon) {
    ErrorMsg(method_name) << "FrequencyEph::calcFrequency produced f1 == " << computed_f1 << " not " << correct_f1 << std::endl;
  }
  
  double computed_f2 = f_eph2.calcFrequency(ev_time, 2);
  double correct_f2 = 13.5e-6;
  if (fabs(computed_f2 - correct_f2) > epsilon) {
    ErrorMsg(method_name) << "FrequencyEph::calcFrequency produced f2 == " << computed_f2 << " not " << correct_f2 << std::endl;
  }

  // The following test is no longer needed or possible because dt is private. It is implicitly
  // tested in phase (and other) calculations.
#if 0
  // Create a frequency ephemeris with unit time 5 s, to test PulsarEph::dt method with one argument.
  FrequencyEph f_eph4("TDB", since, until, epoch, 22., 45., 0.11, 1.125e-2, -2.25e-4, 6.75e-6, 5.);
  double delta_t = f_eph4.dt(ev_time);
  if (fabs(delta_t/20. - 1.) > epsilon)
    ErrorMsg(method_name) << "PulsarEph::dt() produced delta_t == " << delta_t << ", not 20. as expected." << std::endl;

  // Test PulsarEph::dt method with two arguments.
  glast_tdb.setValue(23.456789);
  AbsoluteTime time_origin(glast_tdb);
  delta_t = f_eph4.dt(ev_time, time_origin);
  if (fabs(delta_t/40. - 1.) > epsilon)
    ErrorMsg(method_name) << "PulsarEph::dt() produced delta_t == " << delta_t << ", not 20. as expected." << std::endl;
#endif
}

class SimplePeriodEph {
  public:
    SimplePeriodEph(double p0, double p1, double p2): m_p0(p0), m_p1(p1), m_p2(p2) {}

    double calcFrequency(double dt, double step, int order) {
      if (order < 0) throw std::runtime_error("Bad test parameter given: a negative derivative order");
      if (order == 0) return 1. / (m_p0 + m_p1*dt + m_p2/2.*dt*dt);
      else return (calcFrequency(dt+step/2., step, order-1) - calcFrequency(dt-step/2., step, order-1)) / step;
    }

  private:
    double m_p0, m_p1, m_p2;
};

void PulsarDbTest::testPeriodEph() {
  std::string method_name = "testPeriodEph";

  MetRep glast_tdb("TDB", 51910, 0., 0.);
  glast_tdb.setValue(0.);
  AbsoluteTime since(glast_tdb);
  glast_tdb.setValue(1.);
  AbsoluteTime until(glast_tdb);
  glast_tdb.setValue(123.456789);
  AbsoluteTime epoch(glast_tdb);

  // Create a frequency ephemeris.
  FrequencyEph f_eph("TDB", since, until, epoch, 22., 45., 0.875, 1.125e-2, -2.25e-4, 6.75e-6);

  // Create a period ephemeris.
  // This is a set of values known to be the inverses of the frequency coefficients above.
  PeriodEph p_eph("TDB", since, until, epoch, 22., 45., 0.875, 88.8888888888888888888889,
    1.777777777777777777777778, 0.0177777777777777777777778);

  // Compare frequency & period.
  double epsilon = 1.e-8;
  const double nano_sec = 1.e-9;
  ElapsedTime tolerance("TDB", Duration(0, nano_sec));

  if (!f_eph.getEpoch().equivalentTo(p_eph.getEpoch(), tolerance))
    ErrorMsg(method_name) << "FrequencyEph and PeriodEph give different values for epoch" << std::endl;

  char * field[] = { "ra", "dec", "phi0", "f0", "f1" , "f2" };
  std::pair<double, double> f_ra_dec = f_eph.calcSkyPosition(epoch);
  double value1[] = { f_ra_dec.first, f_ra_dec.second, f_eph.calcPulsePhase(epoch), f_eph.calcFrequency(epoch, 0),
                      f_eph.calcFrequency(epoch, 1), f_eph.calcFrequency(epoch, 2) };
  std::pair<double, double> p_ra_dec = f_eph.calcSkyPosition(epoch);
  double value2[] = { p_ra_dec.first, p_ra_dec.second, p_eph.calcPulsePhase(epoch), p_eph.calcFrequency(epoch, 0),
                      p_eph.calcFrequency(epoch, 1), p_eph.calcFrequency(epoch, 2) };
  for (int ii = 0; ii != sizeof(value1) / sizeof(double); ++ii) {
    if (0. == value1[ii] || 0. == value2[ii]) {
      if (fabs(value1[ii] + value2[ii]) > std::numeric_limits<double>::epsilon())
        ErrorMsg(method_name) << "FrequencyEph and PeriodEph give absolutely different values for " << field[ii] <<
          " (" << value1[ii] << " and " << value2[ii] << ")" << std::endl;
    } else if (fabs(value1[ii] / value2[ii] - 1.) > epsilon) {
      ErrorMsg(method_name) << "FrequencyEph and PeriodEph give fractionally different values for " << field[ii] <<
        " (" << value1[ii] << " and " << value2[ii] << ")" << std::endl;
    }
  }

  // Test frequency computation away from the reference epoch.
  double time_since_epoch = 1000.;
  double step_size = 100.;
  glast_tdb = epoch;
  double dbl_epoch = glast_tdb.getValue();
  glast_tdb.setValue(dbl_epoch + time_since_epoch);
  AbsoluteTime abs_time(glast_tdb);
  double ra = 22.;
  double dec = 45.;
  double phi0 = 0.875;
  double p0 = 1.23456789;
  double p1 = 9.87654321e-5;
  double p2 = 1.357902468e-10;
  PeriodEph p_eph1("TDB", since, until, epoch, ra, dec, phi0, p0, p1, p2);
  SimplePeriodEph s_eph1(p0, p1, p2);

  epsilon = 1.e-3; // Note: Need a loose tolerance because SimplePeriodEph::calcFrequency is not that precise.
  for (int order = 0; order < 5; ++order) {
    double result = p_eph1.calcFrequency(abs_time, order);
    double expected = s_eph1.calcFrequency(time_since_epoch, step_size, order);
    bool test_failed = true;
    if (0. == result || 0. == expected) test_failed = fabs(result + expected) > std::numeric_limits<double>::epsilon();
    else test_failed = std::fabs(result / expected - 1.) > epsilon;
    if (test_failed) {
      ErrorMsg(method_name) << "PeriodEph::calcFrequency(abs_time, " << order << ") returned " << result <<
        ", not " << expected << " as expected" << std::endl;
    }
  }
}

void PulsarDbTest::testSimpleDdEph() {
  std::string method_name = "testSimpleDdEph";

  MetRep glast_tdb("TDB", 51910, 0., 123.456789);
  AbsoluteTime t0(glast_tdb);
  SimpleDdEph o_eph("TDB", 1000., .2, 0., 0., 0., 0., 0., 0., t0, 0., 0., 0.);
  glast_tdb.setValue(223.456789);
  AbsoluteTime ev_time(glast_tdb);

  double epsilon = 1.e-8;

  // Test orbital phase computations.
  double phase = o_eph.calcOrbitalPhase(ev_time);
  if (fabs(phase/.099 - 1.) > epsilon)
    ErrorMsg(method_name) << "SimpleDdEph::calcOrbitalPhase produced phase == " << phase << " not .099" << std::endl;

  // Test orbital phase computations, with a non-zero global phase offset.
  phase = o_eph.calcOrbitalPhase(ev_time, 0.1234);
  if (fabs(phase/.2224 - 1.) > epsilon)
    ErrorMsg(method_name) << "SimpleDdEph::calcOrbitalPhase produced phase == " << phase << " not .2224" << std::endl;

  // The following test is no longer needed or possible because dt is private. It is implicitly
  // tested in phase (and other) calculations.
#if 0
  // Create an orbital ephemeris with unit time 5 s, to test SimpleDdEph::dt method.
  SimpleDdEph o_eph2("TDB", 1000., .2, 0., 0., 0., 0., 0., 0., t0, 0., 0., 0., 5.);
  delta_t = o_eph2.dt(ev_time);
  if (fabs(delta_t/20. - 1.) > epsilon)
    ErrorMsg(method_name) << "SimpleDdEph::dt() produced delta_t == " << delta_t << ", not 20. as expected." << std::endl;
#endif

  // Test binary modulation and demodulation.
  // Binary parameters: (PB, PBDOT, A1, XDOT, ECC, ECCDOT, OM, OMDOT, T0, GAMMA, SHAPIRO_R, SHAPIRO_S)
  // = (27906.980897, -2.43e-12, 2.3417598, 0.0, 0.61713101, 0.0, 220.142729, 4.22662, 45888.1172487, 0.004295, 0.0, 0.0)
  AbsoluteTime abs_t0("TDB", Duration(45888, SecPerDay() * .1172487), Duration(0, 0.));
  SimpleDdEph eph1("TDB", 27906.980897, -2.43e-12, 2.3417598, 0.0, 0.61713101, 0.0, 220.142729, 4.22662,
                   abs_t0, 0.004295, 0.0, 0.0);

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

  // Permitted difference is 100 ns.
  double delta = 100. * 1.e-9;
  ElapsedTime tolerance("TDB", Duration(0, delta));

  std::cerr.precision(24);
  for (size_t ii = 0; ii != sizeof(mjd_test_values)/sizeof(Duration[2]); ++ii) {
    AbsoluteTime tdb_mjd("TDB", mjd_test_values[ii][0], Duration(0, 0.));
    AbsoluteTime expected_tdb_mjd("TDB", mjd_test_values[ii][1], Duration(0, 0.));
    eph1.modulateBinary(tdb_mjd);
    if (!tdb_mjd.equivalentTo(expected_tdb_mjd, tolerance)) {
      ErrorMsg(method_name) << "Binary modulation of " << mjd_test_values[ii][0] << " was computed to be " << tdb_mjd <<
        ", not " << mjd_test_values[ii][1] << ", as expected." << std::endl;
    }
  }

  for (size_t ii = 0; ii != sizeof(mjd_test_values)/sizeof(Duration[2]); ++ii) {
    AbsoluteTime tdb_mjd("TDB", mjd_test_values[ii][1], Duration(0, 0.));
    AbsoluteTime expected_tdb_mjd("TDB", mjd_test_values[ii][0], Duration(0, 0.));
    eph1.demodulateBinary(tdb_mjd);
    if (!tdb_mjd.equivalentTo(expected_tdb_mjd, tolerance)) {
      ErrorMsg(method_name) << "Binary demodulation of " << mjd_test_values[ii][1] << " was computed to be " << tdb_mjd <<
        ", not " << mjd_test_values[ii][0] << ", as expected." << std::endl;
    }
  }
} 

void PulsarDbTest::testPdotCanceler() {
  std::string method_name = "testPdotCanceler";

  // Set test parameters.
  MetRep glast_tt("TT", 51910, 0., 0.);
  glast_tt.setValue(123.4567891234567);
  AbsoluteTime origin(glast_tt);
  glast_tt.setValue(223.4567891234567);
  AbsoluteTime ev_time1(glast_tt);
  AbsoluteTime ev_time2(glast_tt);
  double f0 = 1.125e-2;
  double f1 = -2.25e-4;
  double f2 = 13.5e-6;
  double correct_t = 323.4567891234567;
  double epsilon = 1.e-6; // 1 micro-second.

  // Test PdotCanceler created from literal numbers.
  std::vector<double> fdot_ratio(2);
  fdot_ratio[0] = f1 / f0;
  fdot_ratio[1] = f2 / f0;
  PdotCanceler canceler1("TDB", origin, fdot_ratio);

  canceler1.cancelPdot(ev_time1);
  glast_tt = ev_time1;
  double pdot_t = glast_tt.getValue();
  if (fabs(pdot_t - correct_t) > epsilon) {
    ErrorMsg(method_name) << "After constructed from literal numbers, PdotCanceler::cancelPdot produced pdot-corrected time == "
      << pdot_t << " not " << correct_t << std::endl;
  }

  // Test PdotCanceler created from a PulsarEph object.
  glast_tt.setValue(0.);
  AbsoluteTime since(glast_tt);
  glast_tt.setValue(1.);
  AbsoluteTime until(glast_tt);
  FrequencyEph f_eph("TDB", since, until, origin, 22., 45., .11, f0, f1, f2);
  PdotCanceler canceler2(origin, f_eph, 2);

  canceler2.cancelPdot(ev_time2);
  glast_tt = ev_time2;
  pdot_t = glast_tt.getValue();
  if (fabs(pdot_t - correct_t) > epsilon) {
    ErrorMsg(method_name) << "After constructed from PulsarEph, PdotCanceler::cancelPdot produced pdot-corrected time == "
      << pdot_t << " not " << correct_t << std::endl;
  }
}

void PulsarDbTest::testChooser() {
  std::string method_name = "testChooser";

  // Create pulsar database object.
  PulsarDb database(m_tpl_file);
  database.load(m_in_file);

  // Register PulsarEph and OrbitalEph subclasses for various ephemeris models.
  database.registerPulsarEph<FrequencyEph>("FREQ");
  database.registerOrbitalEph<SimpleDdEph>("DD");

  std::string pulsar_name = "PSR J0139+5814";

  // Select pulsar we want to use.
  database.filterName(pulsar_name);

  // Confirm that the correct number of ephemerides were found.
  int num_eph = database.getNumEph();
  if (8 != num_eph)
    ErrorMsg(method_name) << "there are " << num_eph << " ephemerides for " << pulsar_name << ", not 8" << std::endl;

  // Write this output to form basis for comparing future tests.
  remove("chooser_db.fits");
  database.save("chooser_db.fits", m_creator);

  MjdRep mjd_tdb("TDB", 54012, .5);
  AbsoluteTime pick_time(mjd_tdb);

  PulsarEphCont eph_cont;

  database.getEph(eph_cont);

  StrictEphChooser chooser;

  // Test one with no tiebreaking needed.
  const PulsarEph * chosen = &chooser.choose(eph_cont, pick_time);
  mjd_tdb.setValue(54262, 0.);
  ElapsedTime tolerance("TDB", Duration(0, 1.e-9)); // 1 nanosecond.
  if (!AbsoluteTime(mjd_tdb).equivalentTo(chosen->getEpoch(), tolerance))
    ErrorMsg(method_name) << "for time " << pick_time << ", chooser chose ephemeris with EPOCH == " << chosen->getEpoch() <<
      ", not " << AbsoluteTime(mjd_tdb) << " as expected." << std::endl;

  // Test one with tiebreaking.
  mjd_tdb.setValue(53545, .5);
  pick_time = mjd_tdb;
  chosen = &chooser.choose(eph_cont, pick_time);
  mjd_tdb.setValue(53891, 0.);
  if (!AbsoluteTime(mjd_tdb).equivalentTo(chosen->getEpoch(), tolerance))
    ErrorMsg(method_name) << "for time " << pick_time << ", chooser chose ephemeris with EPOCH == " << chosen->getEpoch() <<
      ", not " << AbsoluteTime(mjd_tdb) << " as expected." << std::endl;

  // Test one which is too early.
  mjd_tdb.setValue(53544, .5);
  pick_time = mjd_tdb;
  try {
    chosen = &chooser.choose(eph_cont, pick_time);
    ErrorMsg(method_name) << "for time " << pick_time << ", chooser chose ephemeris with EPOCH == " << chosen->getEpoch() << std::endl;
  } catch (const std::runtime_error &) {
    // This is to be expected.
  }

  // Test one which is too late.
  mjd_tdb.setValue(55579, .5);
  pick_time = mjd_tdb;
  try {
    chosen = &chooser.choose(eph_cont, pick_time);
    ErrorMsg(method_name) << "for time " << pick_time << ", chooser chose ephemeris with EPOCH == " << chosen->getEpoch() << std::endl;
  } catch (const std::runtime_error &) {
    // This is to be expected.
  }

  // Try one which is too late, but without being strict about validity.
  mjd_tdb.setValue(55579, .5);
  pick_time = mjd_tdb;
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
  pick_time = mjd_tdb;
  try {
    chosen = &chooser.choose(eph_cont, pick_time);
    ErrorMsg(method_name) << "chooser chose ephemeris from an empty set of candidates." << std::endl;
  } catch (const std::runtime_error &) {
    // This is to be expected.
  }

  // Test choosing OrbitalEph.
  // Get new independent access to database, to keep independent from the tests above.
  PulsarDb database2(m_tpl_file);
  database2.load(m_in_file);

  // Register PulsarEph and OrbitalEph subclasses for various ephemeris models.
  database2.registerPulsarEph<FrequencyEph>("FREQ");
  database2.registerOrbitalEph<SimpleDdEph>("DD");

  OrbitalEphCont orbital_cont;
  database2.filterName("PSR J1834-0010");
  database2.getEph(orbital_cont);
  mjd_tdb.setValue(52500, 0.);
  pick_time = mjd_tdb;
  try {
    const OrbitalEph & orbital_eph = chooser.choose(orbital_cont, pick_time);
    mjd_tdb.setValue(52060, .84100795);
    AbsoluteTime expected_t0(mjd_tdb);
    ElapsedTime tolerance("TT", Duration(0, 1.e-9)); // 1 nanosecond.
    if (!orbital_eph.t0().equivalentTo(expected_t0, tolerance)) {
      ErrorMsg(method_name) << "chooser chose orbital ephemeris with time " << orbital_eph.t0() << ", not " <<
        expected_t0 << std::endl;
    }
  } catch (const std::runtime_error & x) {
    ErrorMsg(method_name) << "for time " << pick_time << ", chooser had trouble choosing orbital eph: " << x.what() << std::endl;
  }

  // Test tiebreaking with different tolerances.
  eph_cont.clear();
  Duration origin(51910, 0.);
  AbsoluteTime valid_since("TDB", origin, Duration(0, 101.));
  AbsoluteTime valid_until("TDB", origin, Duration(0, 200.));
  AbsoluteTime epoch("TDB", origin, Duration(0, 150.));
  AbsoluteTime t("TDB", origin, Duration(0, 120.));
  eph_cont.push_back(new FrequencyEph("TT", valid_since, valid_until, epoch, 22., 45., 0., 1., 0., 0.));
  valid_since = AbsoluteTime("TDB", origin, Duration(0, 100.));
  eph_cont.push_back(new FrequencyEph("TDB", valid_since, valid_until, epoch, 22., 45., 0., 1., 0., 0.));

  StrictEphChooser strict_chooser(ElapsedTime("TDB", Duration(0, .99)));
  chosen = &strict_chooser.choose(eph_cont, t);
  if ("TT" != chosen->getSystem().getName())
    ErrorMsg(method_name) << "for time " << t << " with tolerance .99, chooser chose eph with " << chosen->getSystem().getName() <<
      ", not TT as expected" << std::endl;

  strict_chooser = StrictEphChooser(ElapsedTime("TDB", Duration(0, 1.01)));
  chosen = &strict_chooser.choose(eph_cont, t);
  if ("TDB" != chosen->getSystem().getName())
    ErrorMsg(method_name) << "for time " << t << " with tolerance 1.01, chooser chose eph with " << chosen->getSystem().getName() <<
      ", not TDB as expected" << std::endl;

  // Clean up.
  for (PulsarEphCont::reverse_iterator itor = eph_cont.rbegin(); itor != eph_cont.rend(); ++itor) delete *itor;

  // Test sloppy chooser around two disjoint ephemerides.
  eph_cont.clear();
  epoch = MjdRep("TDB", 51910, 0.);
  eph_cont.push_back(new FrequencyEph("TT", AbsoluteTime(MjdRep("TDB", 51910, 0.)), AbsoluteTime(MjdRep("TDB", 51920, 0.)), epoch, 22., 45., 0., 1., 0., 0.));
  eph_cont.push_back(new FrequencyEph("TDB", AbsoluteTime(MjdRep("TDB", 51930, 0.)), AbsoluteTime(MjdRep("TDB", 51940, 0.)), epoch, 22., 45., 0., 1., 0., 0.));

  SloppyEphChooser sloppy_chooser;
  t = MjdRep("TDB", 51905, 0.);
  chosen = &sloppy_chooser.choose(eph_cont, t);
  if ("TT" != chosen->getSystem().getName())
    ErrorMsg(method_name) << "for time before either ephemeris: " << t << ", chooser chose eph with " <<
      chosen->getSystem().getName() << ", not TT as expected" << std::endl;

  t = MjdRep("TDB", 51915, 0.);
  chosen = &sloppy_chooser.choose(eph_cont, t);
  if ("TT" != chosen->getSystem().getName())
    ErrorMsg(method_name) << "for time during first ephemeris: " << t << ", chooser chose eph with " <<
      chosen->getSystem().getName() << ", not TT as expected" << std::endl;

  t = MjdRep("TDB", 51921, 0.);
  chosen = &sloppy_chooser.choose(eph_cont, t);
  if ("TT" != chosen->getSystem().getName())
    ErrorMsg(method_name) << "for time shortly after first ephemeris: " << t << ", chooser chose eph with " <<
      chosen->getSystem().getName() << ", not TT as expected" << std::endl;

  t = MjdRep("TDB", 51929, 0.);
  chosen = &sloppy_chooser.choose(eph_cont, t);
  if ("TDB" != chosen->getSystem().getName())
    ErrorMsg(method_name) << "for time shortly before second ephemeris: " << t << ", chooser chose eph with " <<
      chosen->getSystem().getName() << ", not TDB as expected" << std::endl;

  t = MjdRep("TDB", 51935, 0.);
  chosen = &sloppy_chooser.choose(eph_cont, t);
  if ("TDB" != chosen->getSystem().getName())
    ErrorMsg(method_name) << "for time during second ephemeris: " << t << ", chooser chose eph with " <<
      chosen->getSystem().getName() << ", not TDB as expected" << std::endl;

  t = MjdRep("TDB", 51945, 0.);
  chosen = &sloppy_chooser.choose(eph_cont, t);
  if ("TDB" != chosen->getSystem().getName())
    ErrorMsg(method_name) << "for time after second ephemeris: " << t << ", chooser chose eph with " <<
      chosen->getSystem().getName() << ", not TDB as expected" << std::endl;

  // Test choice from prehistory, say 100000 years before origin of MJD.
  t = MjdRep("TDB", -36525000, 0.);
  try {
    chosen = &sloppy_chooser.choose(eph_cont, t);
    if ("TT" != chosen->getSystem().getName())
      ErrorMsg(method_name) << "for time a long time before the first ephemeris: " << t << ", chooser chose eph with " <<
        chosen->getSystem().getName() << ", not TT as expected" << std::endl;
  } catch (const std::exception & x) {
    ErrorMsg(method_name) << "for time a long time before the first ephemeris: " << t << ", chooser threw exception: " <<
      std::endl << x.what() << std::endl;
  }

  // Test choice from far future, say 100000 years after origin of MJD.
  t = MjdRep("TDB", 36525000, 0.);
  try {
    chosen = &sloppy_chooser.choose(eph_cont, t);
    if ("TDB" != chosen->getSystem().getName())
      ErrorMsg(method_name) << "for time a long time after the second ephemeris: " << t << ", chooser chose eph with " <<
        chosen->getSystem().getName() << ", not TDB as expected" << std::endl;
  } catch (const std::exception & x) {
    ErrorMsg(method_name) << "for time a long time after the second ephemeris: " << t << ", chooser threw exception: " <<
      std::endl << x.what() << std::endl;
  }

  // Clean up.
  for (PulsarEphCont::reverse_iterator itor = eph_cont.rbegin(); itor != eph_cont.rend(); ++itor) delete *itor;
}

void PulsarDbTest::testEphComputer() {
  std::string method_name = "testEphComputer";

  // Set up pieces needed for computation.
  StrictEphChooser chooser;

  // Get access to database.
  PulsarDb database(m_tpl_file);
  database.load(m_in_file);

  // Register PulsarEph and OrbitalEph subclasses for various ephemeris models.
  database.registerPulsarEph<FrequencyEph>("FREQ");
  database.registerOrbitalEph<SimpleDdEph>("DD");

  // Filter a pulsar known to be present.
  database.filterName("PSr j0323+3944");

  // First perform computations without the computer.
  PulsarEphCont eph_cont;
  database.getEph(eph_cont);

  MetRep glast_tdb("TDB", 54101, 0., 100.);
  AbsoluteTime expected_gtdb(glast_tdb);
  const PulsarEph & eph(chooser.choose(eph_cont, expected_gtdb));
  PdotCanceler canceler(eph.getEpoch(), eph, 2);
  canceler.cancelPdot(expected_gtdb);
  glast_tdb = expected_gtdb;
  double expected_elapsed = glast_tdb.getValue();

  // Repeat computations using the EphComputer class, and compare results.
  // Create the computer.
  EphComputer computer;

  // Load the data from the database.
  computer.load(database);

  // Test cancelPdot, and compare result to previous result.
  glast_tdb.setValue(100.);
  AbsoluteTime gtdb(glast_tdb);
  std::vector<double> fdot_ratio(2, 0.);
  double f0 = eph.calcFrequency(eph.getEpoch(), 0);
  fdot_ratio[0] = eph.calcFrequency(eph.getEpoch(), 1) / f0;
  fdot_ratio[1] = eph.calcFrequency(eph.getEpoch(), 2) / f0;
  computer.setPdotCancelParameter(eph.getSystem().getName(), eph.getEpoch(), fdot_ratio);
  computer.cancelPdot(gtdb);
  glast_tdb = gtdb;
  if (expected_elapsed != glast_tdb.getValue())
    ErrorMsg(method_name) << "Given time system name, time origin, and fdot ratios, " <<
      "EphComputer::cancelPdot returned elapsed time " << glast_tdb.getValue() << ", not " <<
      expected_elapsed << ", as expected." << std::endl;

  // Test cancelPdot after setting pdot parameters in a different way, and compare result to previous result.
  glast_tdb.setValue(100.);
  gtdb = glast_tdb;
  computer.setPdotCancelParameter(eph.getEpoch(), eph, 2);
  computer.cancelPdot(gtdb);
  glast_tdb = gtdb;
  if (expected_elapsed != glast_tdb.getValue())
    ErrorMsg(method_name) << "Given time origin, pulsar ephemeris, and maximum derivative order, " <<
      "EphComputer::cancelPdot returned elapsed time " << glast_tdb.getValue() << ", not " <<
      expected_elapsed << ", as expected." << std::endl;

  // Test cancelPdot after setting pdot parameters in yet another way, and compare result to previous result.
  glast_tdb.setValue(100.);
  gtdb = glast_tdb;
  computer.setPdotCancelParameter(eph.getEpoch(), 2);
  computer.cancelPdot(gtdb);
  glast_tdb = gtdb;
  if (expected_elapsed != glast_tdb.getValue())
    ErrorMsg(method_name) << "Given time origin and maximum derivative order, "
      "EphComputer::cancelPdot returned elapsed time " << glast_tdb.getValue() << ", not " <<
      expected_elapsed << ", as expected." << std::endl;

  // Test calcPulsePhase, by comparing it with PulsarEph::calcPulsePhase.
  double expected_pulse_phase = eph.calcPulsePhase(expected_gtdb);
  double pulse_phase = computer.calcPulsePhase(expected_gtdb);
  if (expected_pulse_phase != pulse_phase)
    ErrorMsg(method_name) << "EphComputer::calcPulsePhase returned phase " << pulse_phase << ", not " <<
      expected_pulse_phase << ", as expected." << std::endl;

  // Test calcPulsePhase, by comparing it with PulsarEph::calcPulsePhase, with a non-zero global phase offset.
  expected_pulse_phase = eph.calcPulsePhase(expected_gtdb, 0.1234);
  pulse_phase = computer.calcPulsePhase(expected_gtdb, 0.1234);
  if (expected_pulse_phase != pulse_phase)
    ErrorMsg(method_name) << "EphComputer::calcPulsePhase returned phase " << pulse_phase << ", not " <<
      expected_pulse_phase << ", as expected." << std::endl;

  // Test calcSkyPosition, by comparing it with FrequencyEph::calcSkyPosition.
  std::pair<double, double> expected_ra_dec = eph.calcSkyPosition(expected_gtdb);
  std::pair<double, double> ra_dec = computer.calcSkyPosition(expected_gtdb);
  if (expected_ra_dec != ra_dec)
    ErrorMsg(method_name) << "EphComputer::calcSkyPosition returned (RA, Dec) = (" << ra_dec.first << ", " << ra_dec.second <<
      "), not (" << expected_ra_dec.first << ", " << expected_ra_dec.second << "), as expected." << std::endl;

  // Test binary modulation/demodulation.
  // Get new independent access to database, to keep independent from the tests above.
  PulsarDb database2(m_tpl_file);
  database2.load(m_in_file);

  // Register PulsarEph and OrbitalEph subclasses for various ephemeris models.
  database2.registerPulsarEph<FrequencyEph>("FREQ");
  database2.registerOrbitalEph<SimpleDdEph>("DD");

  // Select a particular pulsar.
  database2.filterName("PSR J1834-0010");
  OrbitalEphCont orbital_eph_cont;
  database2.getEph(orbital_eph_cont);
  const OrbitalEph & orbital_eph(chooser.choose(orbital_eph_cont, gtdb));
  expected_gtdb = gtdb;
  computer.load(database2);

  // Test calcOrbitalPhase, by comparing it with OrbitalEph::calcOrbitalPhase.
  double expected_orbital_phase = orbital_eph.calcOrbitalPhase(expected_gtdb);
  double orbital_phase = computer.calcOrbitalPhase(expected_gtdb);
  if (expected_orbital_phase != orbital_phase)
    ErrorMsg(method_name) << "EphComputer::calcOrbitalPhase returned phase " << orbital_phase << ", not " <<
      expected_orbital_phase << ", as expected." << std::endl;

  // Test calcOrbitalPhase, by comparing it with OrbitalEph::calcOrbitalPhase, with a non-zero global phase offset.
  expected_orbital_phase = orbital_eph.calcOrbitalPhase(expected_gtdb, 0.1234);
  orbital_phase = computer.calcOrbitalPhase(expected_gtdb, 0.1234);
  if (expected_orbital_phase != orbital_phase)
    ErrorMsg(method_name) << "EphComputer::calcOrbitalPhase returned phase " << orbital_phase << ", not " <<
      expected_orbital_phase << ", as expected." << std::endl;

  // Test binary modulation/demodulation.
  // First perform computations without the computer.
  orbital_eph.modulateBinary(expected_gtdb);

  computer.modulateBinary(gtdb);
  glast_tdb = expected_gtdb;
  expected_elapsed = glast_tdb.getValue();
  glast_tdb = gtdb;
  if (expected_elapsed != glast_tdb.getValue())
    ErrorMsg(method_name) << "After EphComputer::modulateBinary, elapsed time was " << glast_tdb.getValue() << ", not " <<
      expected_elapsed << ", as expected." << std::endl;

  expected_gtdb = gtdb;
  orbital_eph.demodulateBinary(expected_gtdb);
  
  computer.demodulateBinary(gtdb);
  glast_tdb = expected_gtdb;
  expected_elapsed = glast_tdb.getValue();
  glast_tdb = gtdb;
  if (expected_elapsed != glast_tdb.getValue())
    ErrorMsg(method_name) << "After EphComputer::demodulateBinary, elapsed time was " << glast_tdb.getValue() << ", not " <<
      expected_elapsed << ", as expected." << std::endl;

  // Test loadPulsarEph.
  // Get containers of spin and orbital ephemerides.
  PulsarEphCont & pulsar_eph_ref(computer.getPulsarEphCont());
  OrbitalEphCont & orbital_eph_ref(computer.getOrbitalEphCont());

  // Clear out the computer.
  pulsar_eph_ref.clear();
  orbital_eph_ref.clear();

  if (!pulsar_eph_ref.empty())
    ErrorMsg(method_name) << "After clearing computer's PulsarEphCont, there were " << pulsar_eph_ref.size() <<
      " ephemerides, not 0 as expected." << std::endl;

  if (!orbital_eph_ref.empty())
    ErrorMsg(method_name) << "After clearing computer's OrbitalEphCont, there were " << orbital_eph_ref.size() <<
      " ephemerides, not 0 as expected." << std::endl;

  // Get access to database.
  PulsarDb database3(m_tpl_file);
  database3.load(m_in_file);

  // Register PulsarEph and OrbitalEph subclasses for various ephemeris models.
  database3.registerPulsarEph<FrequencyEph>("FREQ");
  database3.registerOrbitalEph<SimpleDdEph>("DD");

  // Filter a pulsar known to be present.
  database3.filterName("PSR J1959+2048");

  // Load just the spin parameters from this database.
  computer.loadPulsarEph(database3);

  if (7 != pulsar_eph_ref.size())
    ErrorMsg(method_name) << "After loading spin pulsar ephemerides, there were " << pulsar_eph_ref.size() <<
      " ephemerides, not 7 as expected." << std::endl;

  if (!orbital_eph_ref.empty())
    ErrorMsg(method_name) << "After loading spin pulsar ephemerides, there were " << orbital_eph_ref.size() <<
      " ephemerides, not 0 as expected." << std::endl;

  // Now load just the orbital parameters from this database.
  computer.loadOrbitalEph(database3);

  if (7 != pulsar_eph_ref.size())
    ErrorMsg(method_name) << "After loading orbital pulsar ephemerides, there were " << pulsar_eph_ref.size() <<
      " ephemerides, not 7 as expected." << std::endl;

  if (7 != orbital_eph_ref.size())
    ErrorMsg(method_name) << "After loading orbital orbital ephemerides, there were " << orbital_eph_ref.size() <<
      " ephemerides, not 7 as expected." << std::endl;

}

void PulsarDbTest::testEphGetter() {
  std::string method_name = "testEphGetter";

  // Get access to database.
  PulsarDb database(m_tpl_file);
  database.load(m_in_file);

  // Register PulsarEph and OrbitalEph subclasses for various ephemeris models.
  database.registerPulsarEph<FrequencyEph>("FREQ");
  database.registerOrbitalEph<SimpleDdEph>("DD");

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

class EphRoutingInfo {
  public:
    EphRoutingInfo(): m_string_value(), m_ext_name(), m_model_name(), m_class_name() {}

    EphRoutingInfo(const std::string & string_value, const std::string & ext_name, const std::string & model_name,
      const std::string & class_name): m_string_value(string_value), m_ext_name(ext_name), m_model_name(model_name),
      m_class_name(class_name) {}

    EphRoutingInfo(const tip::Table::ConstRecord & record, const tip::Header & header, const std::string & base_name,
      int class_number): m_string_value(), m_ext_name(), m_model_name(), m_class_name() {
      // Get SRTING_VALUE and store it.
      try {
        record["STRING_VALUE"].get(m_string_value);
      } catch (const tip::TipException &) {
        m_string_value.clear();
      }

      // Get EXTNAME and store it.
      try {
        header["EXTNAME"].get(m_ext_name);
      } catch (const tip::TipException &) {
        m_ext_name.clear();
      }

      // Get EPHSTYLE and store it.
      try {
        header["EPHSTYLE"].get(m_model_name);
      } catch (const tip::TipException &) {
        m_model_name.clear();
      }

      // Construct a class name.
      std::ostringstream oss;
      oss << base_name << class_number;
      m_class_name = oss.str();
    }

    const std::string & getStringValue() const { return m_string_value; }
    const std::string & getExtensionName() const { return m_ext_name; }
    const std::string & getModelName() const { return m_model_name; }
    const std::string & getClassName() const { return m_class_name; }

  private:
    std::string m_string_value;
    std::string m_ext_name;
    std::string m_model_name;
    std::string m_class_name;
};

class BogusPulsarEphBase: public PulsarEph {
  public:
    virtual const EphRoutingInfo & getRoutingInfo() const {
      static const EphRoutingInfo s_bogus_info;
      return s_bogus_info;
    }

    virtual const TimeSystem & getSystem() const { return TimeSystem::getSystem("TDB"); }
    virtual const AbsoluteTime & getValidSince() const { return getBogusTime(); }
    virtual const AbsoluteTime & getValidUntil() const { return getBogusTime(); }
    virtual const AbsoluteTime & getEpoch() const { return getBogusTime(); }
    virtual PulsarEph * clone() const { return new BogusPulsarEphBase(*this); }
    virtual double calcPulsePhase(const timeSystem::AbsoluteTime & /* ev_time */, double /* phase_offset */ = 0.) const {
      return 0.;
    }
    virtual double calcFrequency(const AbsoluteTime & /* ev_time */, int /* derivative_order */ = 0) const { return 0.; }
    virtual std::pair<double, double> calcSkyPosition(const AbsoluteTime & /* ev_time */) const {
      return std::make_pair(0., 0.); 
    }

  protected:
    virtual void writeModelParameter(st_stream::OStream & /* os */) const {}

  private:
    inline const AbsoluteTime & getBogusTime() const {
      static const AbsoluteTime s_bogus_time("TDB", Duration(0, 0.), Duration(0, 0.));
      return s_bogus_time;
    }
};

template<int CLASSNUMBER>
class BogusPulsarEph: public BogusPulsarEphBase {
  public:
    BogusPulsarEph(const tip::Table::ConstRecord & record, const tip::Header & header):
      m_routing_info(record, header, "BogusPulsarEph", CLASSNUMBER) {}
    virtual ~BogusPulsarEph() {}
    virtual const EphRoutingInfo & getRoutingInfo() const { return m_routing_info; }

  private:
    EphRoutingInfo m_routing_info;
};

class BogusOrbitalEphBase: public OrbitalEph {
  public:
    BogusOrbitalEphBase(): OrbitalEph(timeSystem::ElapsedTime("TDB", Duration(0, 10.e-9)), 100) {}

    virtual const EphRoutingInfo & getRoutingInfo() const {
      static const EphRoutingInfo s_bogus_info;
      return s_bogus_info;
    }

    virtual const AbsoluteTime & t0() const {
      static const AbsoluteTime s_bogus_absolute_time("TDB", Duration(0, 0.), Duration(0, 0.));
      return s_bogus_absolute_time;
    }
    virtual double calcOrbitalPhase(const AbsoluteTime & /* ev_time */, double /* phase_offset */ = 0.) const {
      return 0.;
    }
    virtual ElapsedTime calcOrbitalDelay(const timeSystem::AbsoluteTime & /* ev_time */) const {
      static const ElapsedTime s_bogus_elapsed_time("TDB", Duration(0, 0.));
      return s_bogus_elapsed_time;
    }
    virtual OrbitalEph * clone() const { return new BogusOrbitalEphBase(*this); }

  protected:
    virtual void writeModelParameter(st_stream::OStream & /* os */) const {}
};

template<int CLASSNUMBER>
class BogusOrbitalEph: public BogusOrbitalEphBase {
  public:
    BogusOrbitalEph(const tip::Table::ConstRecord & record, const tip::Header & header):
      m_routing_info(record, header, "BogusOrbitalEph", CLASSNUMBER) {}
    virtual ~BogusOrbitalEph() {}
    virtual const EphRoutingInfo & getRoutingInfo() const { return m_routing_info; }

  private:
    EphRoutingInfo m_routing_info;
};

void PulsarDbTest::testMultipleEphModel() {
  std::string method_name = "testMultipleEphModel";

  std::auto_ptr<PulsarDb> database(0);

  // Test rejection of a template file w/o EPHSTYLE in SPIN_PARAMETERS extension.
  std::string tpl_file = facilities::commonUtilities::joinPath(m_data_dir, "test_pulsarDb_badspin.tpl");
  try {
    database.reset(new PulsarDb(tpl_file));
    ErrorMsg(method_name) << "PulsarDb::PulsarDb(\"" << tpl_file << "\") did not throw an exception" << std::endl;
  } catch (const std::exception &) {
    // This is fine.
  }

  // Test rejection of a template file w/o EPHSTYLE in ORBITAL_PARAMETERS extension.
  tpl_file = facilities::commonUtilities::joinPath(m_data_dir, "test_pulsarDb_badorbital.tpl");
  try {
    database.reset(new PulsarDb(tpl_file));
    ErrorMsg(method_name) << "PulsarDb::PulsarDb(\"" << tpl_file << "\") did not throw an exception" << std::endl;
  } catch (const std::exception &) {
    // This is fine.
  }

  // Test successful creation of a PulsarDb object with a correct template.
  tpl_file = facilities::commonUtilities::joinPath(m_data_dir, "test_pulsarDb.tpl");
  try {
    database.reset(new PulsarDb(tpl_file));
  } catch (const std::exception & x) {
    ErrorMsg(method_name) << "PulsarDb::PulsarDb(\"" << tpl_file << "\") threw exception: " << std::endl <<
      x.what() << std::endl;
  }

  // Test rejection of a wrong target extension for spin ephemerides in the original format.
  try {
    database.reset(new PulsarDb(tpl_file, 3, 4));
    ErrorMsg(method_name) << "PulsarDb::PulsarDb(\"" << tpl_file << "\", 3, 4) did not throw an exception" << std::endl;
  } catch (const std::exception &) {
    // This is fine.
  }

  // Test rejection of a wrong target extension for orbital ephemerides in the original format.
  try {
    database.reset(new PulsarDb(tpl_file, 2, 1));
    ErrorMsg(method_name) << "PulsarDb::PulsarDb(\"" << tpl_file << "\", 2, 1) did not throw an exception" << std::endl;
  } catch (const std::exception &) {
    // This is fine.
  }

  // Test successful creation of a PulsarDb object with a correct target extension for spin and orbital ephemerides
  // in the original format.
  try {
    database.reset(new PulsarDb(tpl_file, 2, 4));
  } catch (const std::exception & x) {
    ErrorMsg(method_name) << "PulsarDb::PulsarDb(\"" << tpl_file << "\", 2, 4) threw exception: " << std::endl <<
      x.what() << std::endl;
  }

  // Test loading ephemerides from FITS database files in the current format.
  database.reset(new PulsarDb(tpl_file));
  bool load_original = false;
  bool expected_to_fail = false;
  testLoadingFits(method_name + " (loading current FITS)", *database, tpl_file, load_original, expected_to_fail);

  // Test getting ephemerides that were loaded from FITS database files in the current format.
  database->registerPulsarEph<BogusPulsarEph<1> >("MODEL1");
  database->registerPulsarEph<BogusPulsarEph<2> >("MODEL2");
  database->registerOrbitalEph<BogusOrbitalEph<1> >("MODEL1");
  database->registerOrbitalEph<BogusOrbitalEph<2> >("MODEL2");
  std::map<std::string, EphRoutingInfo> expected_route_dict;
  expected_route_dict["VALUE1"] = EphRoutingInfo("VALUE1", "SPIN_PARAMETERS", "MODEL1", "BogusPulsarEph1");
  expected_route_dict["VALUE2"] = EphRoutingInfo("VALUE2", "SPIN_PARAMETERS", "MODEL2", "BogusPulsarEph2");
  expected_route_dict["VALUE3"] = EphRoutingInfo("VALUE3", "ORBITAL_PARAMETERS", "MODEL1", "BogusOrbitalEph1");
  expected_route_dict["VALUE4"] = EphRoutingInfo("VALUE4", "ORBITAL_PARAMETERS", "MODEL2", "BogusOrbitalEph2");
  checkEphRouting(method_name + " (checking current FITS)", *database, expected_route_dict);

  // Test loading ephemerides from FITS database files in the original format, with target extensions unspecified.
  database.reset(new PulsarDb(tpl_file));
  load_original = true;
  expected_to_fail = true;
  testLoadingFits(method_name + " (loading original FITS)", *database, tpl_file, load_original, expected_to_fail);

  // Test loading ephemerides from FITS database files in the original format, with target extensions specified.
  database.reset(new PulsarDb(tpl_file, 1, 3));
  load_original = true;
  expected_to_fail = false;
  testLoadingFits(method_name + " (loading original FITS)", *database, tpl_file, load_original, expected_to_fail);

  // Test getting ephemerides that were loaded from FITS database files in the original format.
  database->registerPulsarEph<BogusPulsarEph<1> >("MODEL1");
  database->registerPulsarEph<BogusPulsarEph<2> >("MODEL2");
  database->registerOrbitalEph<BogusOrbitalEph<1> >("MODEL1");
  database->registerOrbitalEph<BogusOrbitalEph<2> >("MODEL2");
  expected_route_dict.clear();
  expected_route_dict["VALUE1"] = EphRoutingInfo("VALUE1", "SPIN_PARAMETERS", "MODEL1", "BogusPulsarEph1");
  expected_route_dict["VALUE2"] = EphRoutingInfo("VALUE2", "SPIN_PARAMETERS", "MODEL1", "BogusPulsarEph1");
  expected_route_dict["VALUE3"] = EphRoutingInfo("VALUE3", "ORBITAL_PARAMETERS", "MODEL1", "BogusOrbitalEph1");
  expected_route_dict["VALUE4"] = EphRoutingInfo("VALUE4", "ORBITAL_PARAMETERS", "MODEL1", "BogusOrbitalEph1");
  checkEphRouting(method_name + " (checking original FITS)", *database, expected_route_dict);

  // Prepare information of extensions for loading ephemerides from text database files.
  std::list<std::pair<std::string, std::string> > ext_info_cont;
  ext_info_cont.push_back(std::make_pair("SPIN_PARAMETERS", "MODEL1"));
  ext_info_cont.push_back(std::make_pair("SPIN_PARAMETERS", "MODEL2"));
  ext_info_cont.push_back(std::make_pair("ORBITAL_PARAMETERS", "MODEL1"));
  ext_info_cont.push_back(std::make_pair("ORBITAL_PARAMETERS", "MODEL2"));

  // Test loading ephemerides from text database files in the current format.
  database.reset(new PulsarDb(tpl_file));
  load_original = false;
  expected_to_fail = false;
  testLoadingText(method_name + " (loading current TEXT)", *database, ext_info_cont, load_original, expected_to_fail);

  // Test getting ephemerides that were loaded from text database files.
  database->registerPulsarEph<BogusPulsarEph<1> >("MODEL1");
  database->registerPulsarEph<BogusPulsarEph<2> >("MODEL2");
  database->registerOrbitalEph<BogusOrbitalEph<1> >("MODEL1");
  database->registerOrbitalEph<BogusOrbitalEph<2> >("MODEL2");
  expected_route_dict.clear();
  expected_route_dict["VALUE1"] = EphRoutingInfo("VALUE1", "SPIN_PARAMETERS", "MODEL1", "BogusPulsarEph1");
  expected_route_dict["VALUE2"] = EphRoutingInfo("VALUE2", "SPIN_PARAMETERS", "MODEL2", "BogusPulsarEph2");
  expected_route_dict["VALUE3"] = EphRoutingInfo("VALUE3", "ORBITAL_PARAMETERS", "MODEL1", "BogusOrbitalEph1");
  expected_route_dict["VALUE4"] = EphRoutingInfo("VALUE4", "ORBITAL_PARAMETERS", "MODEL2", "BogusOrbitalEph2");
  checkEphRouting(method_name + " (checking current TEXT)", *database, expected_route_dict);

  // Test loading ephemerides from FITS database files in the original format, with target extensions unspecified.
  database.reset(new PulsarDb(tpl_file));
  load_original = true;
  expected_to_fail = true;
  testLoadingText(method_name + " (loading original TEXT)", *database, ext_info_cont, load_original, expected_to_fail);

  // Test loading ephemerides from FITS database files in the original format, with target extensions specified.
  database.reset(new PulsarDb(tpl_file, 1, 3));
  load_original = true;
  expected_to_fail = false;
  testLoadingText(method_name + " (loading original TEXT)", *database, ext_info_cont, load_original, expected_to_fail);

  // Test getting ephemerides that were loaded from FITS database files in the original format.
  database->registerPulsarEph<BogusPulsarEph<1> >("MODEL1");
  database->registerPulsarEph<BogusPulsarEph<2> >("MODEL2");
  database->registerOrbitalEph<BogusOrbitalEph<1> >("MODEL1");
  database->registerOrbitalEph<BogusOrbitalEph<2> >("MODEL2");
  expected_route_dict.clear();
  expected_route_dict["VALUE1"] = EphRoutingInfo("VALUE1", "SPIN_PARAMETERS", "MODEL1", "BogusPulsarEph1");
  expected_route_dict["VALUE2"] = EphRoutingInfo("VALUE2", "SPIN_PARAMETERS", "MODEL1", "BogusPulsarEph1");
  expected_route_dict["VALUE3"] = EphRoutingInfo("VALUE3", "ORBITAL_PARAMETERS", "MODEL1", "BogusOrbitalEph1");
  expected_route_dict["VALUE4"] = EphRoutingInfo("VALUE4", "ORBITAL_PARAMETERS", "MODEL1", "BogusOrbitalEph1");
  checkEphRouting(method_name + " (checking original TEXT)", *database, expected_route_dict);
}

void PulsarDbTest::testLoadingFits(const std::string & method_name, PulsarDb & database, const std::string & tpl_file,
  bool load_original, bool expected_to_fail) {
  std::string filename = "testdb.fits";

  // Create a FITS file to load ephemerides from.
  tip::IFileSvc & file_svc(tip::IFileSvc::instance());
  file_svc.createFile(filename, tpl_file, true);
  tip::FileSummary file_summary;
  file_svc.getFileSummary(filename, file_summary);
  for (std::size_t ext_number = 1; ext_number < file_summary.size(); ++ext_number) {
    // Open an extension.
    std::ostringstream oss;
    oss << ext_number;
    std::auto_ptr<tip::Table> table(file_svc.editTable(filename, oss.str()));

    // Put the integer number in STRING_VALUE column.
    tip::Table::Iterator record_itor = table->begin();
    tip::TableRecord & record(*record_itor);
    std::string string_value = "VALUE" + oss.str();
    record["STRING_VALUE"].set(string_value);

    // Erase EPHSTYLE header record, if the original format is requested.
    if (load_original) table->getHeader().erase("EPHSTYLE");
  }

  // Test loading ephemerides.
  try {
    // Load the text database.
    database.load(filename);
    if (expected_to_fail) {
      ErrorMsg(method_name) << "PulsarDb::load method did not throw exception for FITS file \"" << filename <<
        "\"" << std::endl;
    } else {
      // This is fine.
    }
  } catch (const std::exception & x) {
    if (expected_to_fail) {
      // This is fine.
    } else {
      ErrorMsg(method_name) << "PulsarDb::load method threw exception for FITS file \"" << filename <<
        "\": " << std::endl << x.what() << std::endl;
    }
  }
}

void PulsarDbTest::testLoadingText(const std::string & method_name, PulsarDb & database, const ExtInfoCont & ext_info_cont,
  bool load_original, bool expected_to_fail) {
  std::string filename("testdb.txt");

  int int_value = 1;
  for (ExtInfoCont::const_iterator itor = ext_info_cont.begin(); itor != ext_info_cont.end(); ++itor) {
    const std::string & ext_name(itor->first);
    const std::string & model_name(itor->second);
    std::ostringstream oss;
    oss << "VALUE" << int_value;
    const std::string string_value = oss.str();
  
    // Create a text database file.
    remove(filename.c_str());
    std::ofstream ofs(filename.c_str());
    ofs << ext_name << std::endl;
    if (!load_original) ofs << "EPHSTYLE = " << model_name << std::endl;
    ofs << "STRING_VALUE" << std::endl;
    ofs << string_value << std::endl;
    ofs.close();

    // Load the text database.
    try {
      database.load(filename);
      if (expected_to_fail) {
        ErrorMsg(method_name) << "PulsarDb::load method did not throw exception for text file \"" << filename <<
          "\" with EXTNAME=" << ext_name << ", EPHSTYLE=" << model_name << ", STRING_VALUE=" << string_value <<
          ": " << std::endl;
      } else {
        // This is fine.
      }
    } catch (const std::exception & x) {
      if (expected_to_fail) {
        // This is fine.
      } else {
        ErrorMsg(method_name) << "PulsarDb::load method threw exception for text file \"" << filename <<
          "\" with EXTNAME=" << ext_name << ", EPHSTYLE=" << model_name << ", STRING_VALUE=" << string_value <<
          ": " << std::endl << x.what() << std::endl;
      }
    }

    // Increment the column value to distinguish test database files.
    ++int_value;
  }
}

void PulsarDbTest::checkEphRouting(const std::string & method_name, const PulsarDb & database,
  const std::map<std::string, EphRoutingInfo> & expected_route_dict) const {
  // Get pulsar and orbital ephemerides out of database.
  PulsarEphCont pulsar_eph_cont;
  database.getEph(pulsar_eph_cont);
  OrbitalEphCont orbital_eph_cont;
  database.getEph(orbital_eph_cont);

  // Collect ephemeris routing information returned by PulsarDb::getEph method.
  std::list<EphRoutingInfo> returned_route_list;
  for (PulsarEphCont::const_iterator itor = pulsar_eph_cont.begin(); itor != pulsar_eph_cont.end(); ++itor) {
    BogusPulsarEphBase * eph(dynamic_cast<BogusPulsarEphBase *>(*itor));
    if (eph == 0) {
      ErrorMsg(method_name) << "PulsarDb::getEph(PulsarEphCont &) returned an object of an unregistered class" <<
        std::endl;
      continue;
    } else {
      returned_route_list.push_back(eph->getRoutingInfo());
    }
  }
  for (OrbitalEphCont::const_iterator itor = orbital_eph_cont.begin(); itor != orbital_eph_cont.end(); ++itor) {
    BogusOrbitalEphBase * eph(dynamic_cast<BogusOrbitalEphBase *>(*itor));
    if (eph == 0) {
      ErrorMsg(method_name) << "PulsarDb::getEph(OrbitalEphCont &) returned an object of an unregistered class" <<
        std::endl;
      continue;
    } else {
      returned_route_list.push_back(eph->getRoutingInfo());
    }
  }

  // Check whether all expected ephemerides are found only once in return ephemerides.
  for (std::map<std::string, EphRoutingInfo>::const_iterator expected_itor = expected_route_dict.begin();
    expected_itor != expected_route_dict.end(); ++expected_itor) {
    const std::string & string_value(expected_itor->first);
    std::size_t num_eph_found = 0;
    for (std::list<EphRoutingInfo>::const_iterator returned_itor = returned_route_list.begin();
      returned_itor != returned_route_list.end(); ++returned_itor) {
      if (string_value == returned_itor->getStringValue()) ++num_eph_found;
    }
    if (num_eph_found == 0) {
      ErrorMsg(method_name) << "Ephemeris with value \"" << string_value <<
        "\" was not returned by PulsarDb::getEph method" << std::endl;
    } else if (num_eph_found > 1) {
      ErrorMsg(method_name) << "Ephemeris with value \"" << string_value <<
        "\" was returned more than once by PulsarDb::getEph method" << std::endl;
    }
  }

  for (std::list<EphRoutingInfo>::const_iterator returned_itor = returned_route_list.begin();
    returned_itor != returned_route_list.end(); ++returned_itor) {
    // Get routing information of this ephemeris entry.
    std::string string_value = returned_itor->getStringValue();
    std::string ext_name = returned_itor->getExtensionName();
    std::string model_name = returned_itor->getModelName();
    std::string class_name = returned_itor->getClassName();

    // Look for this entry in the expected routing information dictionary.
    std::map<std::string, EphRoutingInfo>::const_iterator expected_itor = expected_route_dict.find(string_value);
    if (expected_itor == expected_route_dict.end()) {
      ErrorMsg(method_name) << "PulsarDb::getEph method returned an unexpected ephemeris data: " <<
        string_value << std::endl;

    } else {
      // Check an extension value of this ephemeris entry.
      std::string expected_ext_name = expected_itor->second.getExtensionName();
      if (ext_name != expected_ext_name) {
      ErrorMsg(method_name) << "Ephemeris with value \"" << string_value <<
        "\" was loaded into an extension with EXTNAME=" << ext_name << ", not " << expected_ext_name <<
        " as expected" << std::endl;
      }

      // Check EPHSTYLE keyword value of an extension that this ephemeris entry was coming through.
      std::string expected_model_name = expected_itor->second.getModelName();
      if (model_name != expected_model_name) {
        ErrorMsg(method_name) << "Ephemeris with value \"" << string_value <<
          "\" was loaded into an extension with EPHSTYLE=" << model_name << ", not " << expected_model_name <<
          " as expected" << std::endl;
      }

      // Check a class name that this ephemeris entry was passed to.
      std::string expected_class_name = expected_itor->second.getClassName();
      if (class_name != expected_class_name) {
        ErrorMsg(method_name) << "Ephemeris with value \"" << string_value <<
          "\" was passed to " << class_name << " class, not " << expected_class_name <<
          " as expected" << std::endl;
      }
    }
  }

  // Clear the contents of the pulsar ephemeris container.
  for (PulsarEphCont::iterator itor = pulsar_eph_cont.begin(); itor != pulsar_eph_cont.end(); ++itor) delete *itor;
  pulsar_eph_cont.clear();
}

st_app::StAppFactory<PulsarDbTest> g_factory("test_pulsarDb");
