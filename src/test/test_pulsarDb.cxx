/** \file test_pulsarDb.cxx
    \brief Test code for pulsarDb package.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "pulsarDb/GlastTime.h"
#include "pulsarDb/PulsarDb.h"

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

  // Test one with no tiebreaking needed.
  PulsarEph chosen = database.chooseEph(pick_time);
  if (54514 != chosen.m_until)
    ErrorMsg(method_name) << "for time " << pick_time << ", chooser chose ephemeris with VALID_UNTIL == " << chosen.m_until << std::endl;

  // Test one with tiebreaking.
  pick_time = 53545.5;
  chosen = database.chooseEph(pick_time);
  if (54238 != chosen.m_until)
    ErrorMsg(method_name) << "for time " << pick_time << ", chooser chose ephemeris with VALID_UNTIL == " << chosen.m_until << std::endl;

  // Test one which is too early.
  pick_time = 53544.5;
  try {
    chosen = database.chooseEph(pick_time);
    ErrorMsg(method_name) << "for time " << pick_time << ", chooser chose ephemeris with VALID_UNTIL == " << chosen.m_until << std::endl;
  } catch (const std::runtime_error &) {
    // This is to be expected.
  }

  // Test one which is too late.
  pick_time = 55579.5;
  try {
    chosen = database.chooseEph(pick_time);
    ErrorMsg(method_name) << "for time " << pick_time << ", chooser chose ephemeris with VALID_UNTIL == " << chosen.m_until << std::endl;
  } catch (const std::runtime_error &) {
    // This is to be expected.
  }

  // Extrapolate one which is too late.
  pick_time = 55579.5;
  try {
    chosen = database.chooseEph(pick_time, true);
  } catch (const std::runtime_error &) {
    ErrorMsg(method_name) << "for time " << pick_time << ", chooser could not extrapolate an ephemeris" << std::endl;
  }

  // Make a selection which will result in an empty container of ephemerides.
  database.filterName("Aunt Gertrude");
  if (0 != database.getNumEph()) {
    ErrorMsg(method_name) << "What? Aunt Gertrude is a PULSAR???" << std::endl;
  }
  
  // Try to choose an ephemeris from the empty set.
  pick_time = 55579.5;
  try {
    chosen = database.chooseEph(pick_time);
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

st_app::StAppFactory<PulsarDbTest> g_factory("gtpulsardb");
