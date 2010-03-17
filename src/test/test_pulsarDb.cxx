/** \file test_pulsarDb.cxx
    \brief Test code for pulsarDb package.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <list>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "pulsarDb/EphChooser.h"
#include "pulsarDb/EphComputer.h"
#include "pulsarDb/EphComputerApp.h"
#include "pulsarDb/EphStatus.h"
#include "pulsarDb/FrequencyEph.h"
#include "pulsarDb/OrbitalEph.h"
#include "pulsarDb/PdotCanceler.h"
#include "pulsarDb/PeriodEph.h"
#include "pulsarDb/PulsarDb.h"
#include "pulsarDb/PulsarDbApp.h"
#include "pulsarDb/PulsarEph.h"
#include "pulsarDb/SimpleDdEph.h"

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "st_stream/Stream.h"

#include "timeSystem/AbsoluteTime.h"
#include "timeSystem/ElapsedTime.h"
#include "timeSystem/CalendarFormat.h"
#include "timeSystem/MjdFormat.h"
#include "timeSystem/PulsarTestApp.h"
#include "timeSystem/TimeConstant.h"
#include "timeSystem/TimeInterval.h"

#include "tip/Extension.h"
#include "tip/FileSummary.h"
#include "tip/Header.h"
#include "tip/IFileSvc.h"
#include "tip/Table.h"
#include "tip/TipException.h"

static const std::string s_cvs_id("$Name:  $");

using namespace timeSystem;
using namespace pulsarDb;

class EphRoutingInfo;
class TextEph;

/** \class PulsarDbAppTester
    \brief Test PulsarDbApp application (gtpulsardb).
*/
class PulsarDbAppTester: public PulsarApplicationTester {
  public:
  /** \brief Construct a PulsarDbAppTester object.
      \param test_app Unit test appliction of pulsar tool package, under which this application tester is to run.
  */
  PulsarDbAppTester(PulsarTestApp & test_app);

  /// \brief Destruct this PulsarDbAppTester object.
  virtual ~PulsarDbAppTester() throw() {}

  /// \brief Returns an application object to be tested.
  virtual st_app::StApp * createApplication() const;

  /** \brief Return a logical true if the given header keyword is determined correct, and a logical false otherwise.
      \param keyword_name Name of the header keyword to be verified.
      \param out_keyword Header keyword taken from the output file to be verified.
      \param ref_keyword Header keyword taken from the reference file which out_keyword is checked against.
      \param error_stream Output stream for this method to put an error messages when verification fails.
  */
  virtual bool verify(const std::string & keyword_name, const tip::KeyRecord & out_keyword,
    const tip::KeyRecord & ref_keyword, std::ostream & error_stream) const;

  /** \brief Return a logical true if the given table cell is considered correct, and a logical false otherwise.
      \param column_name Name of the FITS column that the given table cell belongs to.
      \param out_cell Table cell taken from the output file to be verified.
      \param ref_cell Table cell taken from the reference file which out_cell is checked against.
      \param error_stream Output stream for this method to put an error message when verification fails.
  */
  virtual bool verify(const std::string & column_name, const tip::TableCell & out_cell, const tip::TableCell & ref_cell,
    std::ostream & error_stream) const;
};

PulsarDbAppTester::PulsarDbAppTester(PulsarTestApp & test_app): PulsarApplicationTester("gtpulsardb", test_app) {}

st_app::StApp * PulsarDbAppTester::createApplication() const {
  return new PulsarDbApp();
}

bool PulsarDbAppTester::verify(const std::string & /*keyword_name*/, const tip::KeyRecord & out_keyword,
  const tip::KeyRecord & ref_keyword, std::ostream & error_stream) const {
  // Extract keyword values as character strings.
  std::string out_value = out_keyword.getValue();
  std::string ref_value = ref_keyword.getValue();

  // Require an exact match as character strings.
  bool verified = (out_value == ref_value);
  if (!verified) error_stream << "Character string \"" << out_value << "\" not identical to \"" << ref_value << "\"";

  // Return the result.
  return verified;
}

bool PulsarDbAppTester::verify(const std::string & column_name, const tip::TableCell & out_cell,
  const tip::TableCell & ref_cell, std::ostream & error_stream) const {
  // Initialize return value.
  bool verified = false;

  // Compare columns.
  if ("PSRNAME" == column_name || "VALID_SINCE" == column_name || "VALID_UNTIL" == column_name || "EPOCH_INT" == column_name ||
    "TOAGEO_INT" == column_name || "TOABARY_INT" == column_name || "OBSERVER_CODE" == column_name || "BINARY_FLAG" == column_name ||
    "SOLAR_SYSTEM_EPHEMERIS" == column_name || "DESCRIPTION" == column_name || "OBSERVATORY" == column_name ||
    "CONTACT_PERSON" == column_name || "REFERENCE" == column_name || "ALTNAME" == column_name) {
    // Require an exact match as character strings.
    std::string out_value;
    std::string ref_value;
    out_cell.get(out_value);
    ref_cell.get(ref_value);
    verified = (out_value == ref_value);
    if (!verified) error_stream << "Character string \"" << out_value << "\" not identical to \"" << ref_value << "\"";

  } else {
    // Extract cell values as floating-point numbers.
    double out_value;
    double ref_value;
    out_cell.get(out_value);
    ref_cell.get(ref_value);
    error_stream.precision(std::numeric_limits<double>::digits10);

    if ("RA" == column_name || "DEC" == column_name) {
      // Require a match down to the 10th decimal point.
      verified = (std::fabs(out_value - ref_value) <= 1.e-10);
      if (!verified) {
        error_stream << "Coordinate " << out_value << " not equivalent to reference " << ref_value <<
          " with absolute tolerance of 1e-10 degrees.";
      }

    } else if ("EPOCH_FRAC" == column_name || "TOAGEO_FRAC" == column_name || "TOABARY_FRAC" == column_name) {
      // Require a match down to the 10th decimal point, which is approx. 10 microseconds.
      verified = (std::fabs(out_value - ref_value) <= 1.e-10);
      if (!verified) {
        error_stream << out_value << " days not equivalent to reference " << ref_value <<
          " days with absolute tolerance of 1e-10 days (8.64 microseconds).";
      }

    } else if ("EFFECTIVE_SINCE" == column_name || "EFFECTIVE_UNTIL" == column_name) {
      // Require a match down to the 8th decimal point, which is approx. 1 milliseconds.
      verified = (std::fabs(out_value - ref_value) <= 1.e-8);
      if (!verified) {
        error_stream << "MJD number " << out_value << " not equivalent to reference " << ref_value <<
          " with absolute tolerance of 1e-8 days (0.864 milliseconds).";
      }

    } else {
      // Compare as floating-point numbers.
      // NOTE: This requires out_value be exactly zero when ref_value is zero.
      double rel_tol = std::numeric_limits<double>::epsilon() * 1000.;
      verified = (std::fabs(out_value - ref_value) <= std::fabs(ref_value) * rel_tol);
      if (!verified) {
        error_stream << "Value " << out_value << " not equivalent to reference " << ref_value <<
          " with relative tolerance of " << rel_tol;
      }
    }
  }

  // Return the result.
  return verified;
}

/** \class EphComputerAppTester
    \brief Test EphComputerApp application (gtephem).
*/
class EphComputerAppTester: public PulsarApplicationTester {
  public:
  /** \brief Construct a EphComputerAppTester object.
      \param test_app Unit test appliction of pulsar tool package, under which this application tester is to run.
  */
  EphComputerAppTester(PulsarTestApp & test_app);

  /// \brief Destruct this EphComputerAppTester object.
  virtual ~EphComputerAppTester() throw() {}

  /// \brief Returns an application object to be tested.
  virtual st_app::StApp * createApplication() const;

  /** \brief Return a logical true if the given character string is considered correct, and a logical false otherwise.
      \param out_string Character string taken from the output file to be verified.
      \param ref_string Character string taken from the reference file which out_string is checked against.
      \param error_stream Output stream for this method to put an error message when verification fails.
  */
  virtual bool verify(const std::string & out_string, const std::string & ref_string, std::ostream & error_stream) const;
};

EphComputerAppTester::EphComputerAppTester(PulsarTestApp & test_app): PulsarApplicationTester("gtephem", test_app) {}

st_app::StApp * EphComputerAppTester::createApplication() const {
  return new EphComputerApp();
}

bool EphComputerAppTester::verify(const std::string & out_string, const std::string & ref_string, std::ostream & error_stream) const {
  // Initialize return value.
  bool verified = false;

  if (ref_string.find("seconds after") != std::string::npos) {
    // Require a match down to the 5th decimal point (10 microseconds for elapsed-time part).
    // NOTE: This block must come before "MJD"-block, because an absolute time string contains character string "MJD".
    verified = equivalent(out_string, ref_string, 1.e-5, 0.);
    if (!verified) {
      error_stream << "Absolute time not equivalent to reference with absolute tolerance of 1e-5 seconds." <<
        std::endl << "[OUT] " << out_string << std::endl << "[REF] " << ref_string;
    }

  } else if (ref_string.find("MJD") != std::string::npos) {
    // Require a match down to the 10th decimal point, which is approx. 10 microseconds.
    verified = equivalent(out_string, ref_string, 1.e-10, 0.);
    if (!verified) {
      error_stream << "MJD number not equivalent to reference with absolute tolerance of 1e-10 days (8.64 microseconds)." <<
        std::endl << "[OUT] " << out_string << std::endl << "[REF] " << ref_string;
    }

  } else if (ref_string.find("Reference Time") != std::string::npos) {
    // Require a match down to the 5th decimal point (10 microseconds for elapsed second part in ISO8601 format).
    // NOTE: This block must come after "MJD"-block to avoid checking Reference Time in MJD representation.
    verified = equivalent(out_string, ref_string, 1.e-5, 0.);
    if (!verified) {
      error_stream << "Absolute time not equivalent to reference with absolute tolerance of 1e-5 seconds." <<
        std::endl << "[OUT] " << out_string << std::endl << "[REF] " << ref_string;
    }

  } else if (ref_string.find("Right Ascension") != std::string::npos || ref_string.find("Declination") != std::string::npos ||
    ref_string.find("RA") != std::string::npos || ref_string.find("Dec") != std::string::npos) {
    // Require a match down to the 10th decimal point.
    verified = equivalent(out_string, ref_string, 1.e-10, 0.);
    if (!verified) {
      error_stream << "Coordinate not equivalent to reference with absolute tolerance of 1e-10 degrees." <<
        std::endl << "[OUT] " << out_string << std::endl << "[REF] " << ref_string;
    }

  } else if (ref_string.find("Pulse Phase") != std::string::npos || ref_string.find("Phi0") != std::string::npos) {
    // Require a match down to the 3rd decimal point.
    verified = equivalent(out_string, ref_string, 1.e-3, 0.);
    if (!verified) {
      error_stream << "Pulse phase not equivalent to reference with absolute tolerance of 1e-3." <<
        std::endl << "[OUT] " << out_string << std::endl << "[REF] " << ref_string;
    }

  } else {
    // Compare others as floating-point numbers with default tolerances.
    verified = equivalent(out_string, ref_string);
    if (!verified) {
      error_stream << "Line not equivalent to reference." <<
        std::endl << "[OUT] " << out_string << std::endl << "[REF] " << ref_string;
    }
  }

  // Return the result.
  return verified;
}

/** \class PulsarDbTestApp
    \brief Test pulsarDb package and applications in it.
*/
class PulsarDbTestApp : public PulsarTestApp {
  public:
    /// \brief Construct a PulsarDbTestApp object.
    PulsarDbTestApp();

    /// \brief Destruct this PulsarDbTestApp object.
    virtual ~PulsarDbTestApp() throw() {}

    /// \brief Do all tests.
    virtual void runTest();

    /// \brief Test filtering which doesn't actually narrow the selection.
    virtual void testNoOp();

    /// \brief Test finding pulsars using their PSR J name.
    virtual void testExplicitName();

    /// \brief Test recognition of pulsars by alternate names given in the ALTERNATIVE_NAMES extension.
    virtual void testAlternateName();

    /// \brief Test filtering based on a time range (finding ephemerides which overlaps the range.)
    virtual void testTime();

    /// \brief Test error cases to make sure errors are detected properly.
    virtual void testBadInterval();

    /// \brief Test filtering based on a solar system ephemeris name.
    virtual void testSolarEph();

    /// \brief Test filtering expressions.
    virtual void testExpression();

    /// \brief Test appending ephemerides to an existing file.
    virtual void testAppend();

    /// \brief Test ephemerides database in text format. Also test the getter for history records.
    virtual void testTextPulsarDb();

    /// \brief Test FrequencyEph class (a subclass of PulsarEph class).
    virtual void testFrequencyEph();

    /// \brief Test PeriodEph classes (a subclass of PulsarEph class).
    virtual void testPeriodEph();

    /// \brief Test SimpleDdEph class (a subclass of OrbitalEph class).
    virtual void testSimpleDdEph();

    /// \brief Test PdotCanceler class.
    virtual void testPdotCanceler();

    /// \brief Test method which chooses the best ephemeris from several which could be used.
    virtual void testChooser();

    /// \brief Test EphComputer class.
    virtual void testEphComputer();

    /// \brief Test Eph getter in PulsarDb class. Also test the getter for ephemeris remarks.
    virtual void testEphGetter();

    /// \brief Test support for multiple ephemeris models.
    virtual void testMultipleEphModel();

    /// \brief Test EphStatus class.
    virtual void testEphStatus();

    /// \brief Test binary-ness determination by PulsarDb class.
    virtual void testBinary();

    /// \brief Test PulsarDbApp class.
    virtual void testPulsarDbApp();

    /// \brief Test EphComputerApp class.
    virtual void testEphComputerApp();

  private:
    std::string m_in_file;
    std::string m_tpl_file;
    std::string m_creator;
    std::string m_author;

    /** \brief Helper method for testMultipleEphModel, to test loading ephemerides from FITS database files.
        \param test_subject Character string to identify a subject to be tested.
        \param database Pulsar ephemeris database object to load ephemeris to.
        \param tpl_file Name of the FITS template file to be used to create a FITS file to load ephemerides from.
        \param load_original Set to true to test loading ephemerides from a FITS file in the original format (w/o EPHSTYLE
               keyword). Set to false to test loading from a FITS file in the current format.
        \param expected_to_fail Set to true if this test is expected to fail. Set to false otherwise.
    */
    void testLoadingFits(const std::string & test_subject, PulsarDb & database, const std::string & tpl_file,
      bool load_original, bool expected_to_fail);

    /** \brief Helper method for testMultipleEphModel, to test loading ephemerides from text database files.
        \param test_subject Character string to identify a subject to be tested.
        \param database Pulsar ephemeris database object to load ephemeris to.
        \param ext_name Extension name of an input text ephemerides database to load ephemerides from.
        \param eph_style EPHSTYLE value of an input text ephemerides database to load ephemerides from.
        \param string_value Ephemeris data (STRING_VALUE column value) to load.
        \param load_original Set to true to test loading ephemerides from a FITS file in the original format (w/o EPHSTYLE
               keyword). Set to false to test loading from a FITS file in the current format.
        \param expected_to_fail Set to true if this test is expected to fail. Set to false otherwise.
    */
    void testLoadingText(const std::string & test_subject, PulsarDb & database, const std::string & ext_name,
      const std::string & eph_style, const std::string & string_value, bool load_original, bool expected_to_fail);

    /** \brief Helper method for testMultipleEphModel, to check ephemerides returned by PulsarDb::getEph method.
        \param test_subject Character string to identify a subject to be tested.
        \param database Pulsar ephemeris database object to load ephemeris to.
        \param expected_route_dict Dictionary of reference routing information for test results to be compared with.
    */
    void checkEphRouting(const std::string & test_subject, const PulsarDb & database,
      const std::map<std::string, EphRoutingInfo> & expected_route_dict);

    /** \brief Helper method to check text outputs returned by FormattedEph::write method.
        \param base_name Character string to be used to create temporary file names from.
        \param eph Ephemeris to be tested.
        \param text_eph TextEph object that holds character strings to be compared with a text output from write method of
               the given ephemeris object (eph).
    */
    void checkEphParameter(const std::string & base_name, const FormattedEph & eph, const TextEph & text_eph);
};

PulsarDbTestApp::PulsarDbTestApp(): PulsarTestApp("pulsarDb"), m_in_file(), m_tpl_file(), m_creator(), m_author() {
  setName("test_pulsarDb");
  setVersion(s_cvs_id);

  // Set test file name.
  m_in_file = prependDataPath("groD4-dc2v5.fits");

  // Set template file name.
  m_tpl_file = prependDataPath("PulsarDb.tpl");

  // Set a default value for CREATOR header keyword.
  m_creator = getName() + " " + getVersion();

  // Set a default value for AUTHOR header keyword.
  m_author = "Anonymous Tester";
}

void PulsarDbTestApp::runTest() {
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
  testTextPulsarDb();

  // Test ephemeris computation.
  testFrequencyEph();
  testPeriodEph();
  testSimpleDdEph();
  testPdotCanceler();

  // Test ephemeris manipulation.
  testChooser();
  testEphComputer();
  testEphGetter();
  testMultipleEphModel();
  testEphStatus();
  testBinary();

  // Test applications.
  testPulsarDbApp();
  testEphComputerApp();
}

void PulsarDbTestApp::testNoOp() {
  setMethod("testNoOp");
  PulsarDbAppTester tester(*this);

  // Create pulsar database object.
  PulsarDb database(m_tpl_file);
  database.load(m_in_file);

  // Make sure the test data has expected size.
  int num_eph = database.getNumEph();
  if (1141 != num_eph)
    err() << "there are initially " << num_eph << " ephemerides, not 1141" << std::endl;

  // Perform filtering which doesn't exclude any ephemerides.
  database.filterInterval(0., 1.e6);

  // Make sure the test data still has expected size.
  num_eph = database.getNumEph();
  if (1141 != num_eph)
    err() << "After no-op filterInterval there are " << num_eph << " ephemerides, not 1141" << std::endl;

  // Another no-op is to use "any".
  database.filterName("aNy");

  // Make sure the test data still has expected size.
  num_eph = database.getNumEph();
  if (1141 != num_eph)
    err() << "After no-op filterName there are " << num_eph << " ephemerides, not 1141" << std::endl;

  // Another no-op is to use a silly expression
  database.filterExpression("#row<1142");

  // Make sure the test data still has expected size.
  num_eph = database.getNumEph();
  if (1141 != num_eph)
    err() << "After no-op filter there are " << num_eph << " ephemerides, not 1141" << std::endl;

  // Another no-op is to use a blank string.
  database.filterExpression("\t  \t");

  // Make sure the test data still has expected size.
  num_eph = database.getNumEph();
  if (1141 != num_eph)
    err() << "After no-op filter there are " << num_eph << " ephemerides, not 1141" << std::endl;

  // Save the result for basis of comparing future test output.
  std::string outfile("noop_db.fits");
  remove(outfile.c_str());
  database.save(outfile, m_creator, m_author);

  // Check the result against its reference file in data/outref/ directory.
  tester.checkOutputFits(outfile, prependOutrefPath(outfile));
}

void PulsarDbTestApp::testExplicitName() {
  setMethod("testExplicitName");
  PulsarDbAppTester tester(*this);

  // Create pulsar database object.
  PulsarDb database(m_tpl_file);
  database.load(m_in_file);

  // Filter a pulsar known to be present.
  database.filterName("PSr j0323+3944");

  // Confirm that the correct number of ephemerides were found.
  int num_eph = database.getNumEph();
  if (2 != num_eph)
    err() << "there are " << num_eph << " ephemerides for PSR J0323+3944, not 2" << std::endl;

  // Save the result for basis of comparing future test output.
  std::string outfile("j0323_db.fits");
  remove(outfile.c_str());
  database.save(outfile, m_creator, m_author);

  // Check the result against its reference file in data/outref/ directory.
  tester.checkOutputFits(outfile, prependOutrefPath(outfile));
}

void PulsarDbTestApp::testAlternateName() {
  setMethod("testAlternateName");
  PulsarDbAppTester tester(*this);

  // Create pulsar database object.
  PulsarDb database(m_tpl_file);
  database.load(m_in_file);

  // Guess who?
  database.filterName("CRab");

  // Confirm that the correct number of ephemerides were found.
  int num_eph = database.getNumEph();
  if (36 != num_eph)
    err() << "there are " << num_eph << " ephemerides for the crab, not 36" << std::endl;

  // Filter on the crab's PSR B name (no-op).
  database.filterName("PsR b0531+21");
  num_eph = database.getNumEph();
  if (36 != num_eph)
    err() << "there are " << num_eph << " ephemerides for the crab, not 36" << std::endl;

  // Write this output to form basis for comparing future tests.
  std::string outfile("crab_db.fits");
  remove(outfile.c_str());
  database.save(outfile, m_creator, m_author);

  // Check the result against its reference file in data/outref/ directory.
  tester.checkOutputFits(outfile, prependOutrefPath(outfile));
}

void PulsarDbTestApp::testTime() {
  setMethod("testTime");
  PulsarDbAppTester tester(*this);

  // Create pulsar database object.
  PulsarDb database(m_tpl_file);
  database.load(m_in_file);

  // Give a time filter.
  database.filterInterval(53400., 53800.);

  // Make sure the test data has expected size.
  int num_eph = database.getNumEph();
  if (406 != num_eph)
    err() << "after filterInterval(53400., 53800.) there are " << num_eph << " ephemerides, not 406" << std::endl;

  // Save the result for basis of comparing future test output.
  std::string outfile("time_db.fits");
  remove(outfile.c_str());
  database.save(outfile, m_creator, m_author);

  // Check the result against its reference file in data/outref/ directory.
  tester.checkOutputFits(outfile, prependOutrefPath(outfile));
}

void PulsarDbTestApp::testBadInterval() {
  setMethod("testBadInterval");

  // Create pulsar database object.
  PulsarDb database(m_tpl_file);
  database.load(m_in_file);

  try {
    // Invalid interval, with start time later than stop time.
    database.filterInterval(54500., 54499.);
    err() << "filterInterval(" << 54500. << ", " <<
      54499. << ") did not throw an exception" << std::endl;
  } catch (const std::exception &) {
    // This is fine.
  }
}

void PulsarDbTestApp::testSolarEph() {
  setMethod("testSolarEph");
  PulsarDbAppTester tester(*this);

  // Create pulsar database object.
  PulsarDb database(m_tpl_file);
  database.load(m_in_file);

  // Give a time filter.
  database.filterSolarEph("JPL DE405");

  // Make sure the test data has expected size (SPIN_PARAMETERS extension).
  int num_eph = database.getNumEph();
  if (5 != num_eph)
    err() << "after filterSolarEph(\"JPL DE405\") there are " << num_eph << " spin ephemerides, not 5" << std::endl;

  // Make sure the test data has expected size (ORBITAL_PARAMETERS extension).
  num_eph = database.getNumEph(false);
  if (4 != num_eph)
    err() << "after filterSolarEph(\"JPL DE405\") there are " << num_eph << " orbital ephemerides, not 4" << std::endl;

  // Save the result for basis of comparing future test output.
  std::string outfile("solar_db.fits");
  remove(outfile.c_str());
  database.save(outfile, m_creator, m_author);

  // Check the result against its reference file in data/outref/ directory.
  tester.checkOutputFits(outfile, prependOutrefPath(outfile));
}

void PulsarDbTestApp::testExpression() {
  setMethod("testExpression");
  PulsarDbAppTester tester(*this);

  // Create pulsar database object.
  PulsarDb database(m_tpl_file);
  database.load(m_in_file);

  // Filter using more complex criteria.
  database.filterExpression("F2 != 0.");

  // Confirm that the correct number of ephemerides were found.
  int num_eph = database.getNumEph();
  if (837 != num_eph)
    err() << "found " << num_eph << " ephemerides with F2 != 0., not 837" << std::endl;

  // Test saving this for basis of comparing future test output.
  std::string outfile("f2not0_db.fits");
  remove(outfile.c_str());
  database.save(outfile, m_creator, m_author);

  // Check the result against its reference file in data/outref/ directory.
  tester.checkOutputFits(outfile, prependOutrefPath(outfile));
}

void PulsarDbTestApp::testAppend() {
  setMethod("testAppend");
  PulsarDbAppTester tester(*this);

  // Load the same data base twice.
  PulsarDb database(m_tpl_file);
  database.load(m_in_file);
  database.load(m_in_file);

  // Save the result for basis of comparing future test output.
  std::string outfile("twice_db.fits");
  remove(outfile.c_str());
  database.save(outfile, m_creator, m_author);

  // Check the result against its reference file in data/outref/ directory.
  tester.checkOutputFits(outfile, prependOutrefPath(outfile));


  // Load a FITS database with TELESCOP=FERMI.
  std::string filename(prependDataPath("groD4-dc2v6.fits"));
  try {
    database.load(filename);
  } catch (const std::exception & x) {
    err() << "PulsarDb::load method threw an exception for FITS file \"" << filename << "\": " << std::endl << x.what() << std::endl;
  }

  // Load a TEXT database with TELESCOP=FERMI.
  filename = prependDataPath("psrdb_spin_fermi.txt");
  try {
    database.load(filename);
  } catch (const std::exception & x) {
    err() << "PulsarDb::load method threw an exception for TEXT file \"" << filename << "\": " << std::endl << x.what() << std::endl;
  }
}

void PulsarDbTestApp::testTextPulsarDb() {
  setMethod("testTextPulsarDb");
  PulsarDbAppTester tester(*this);

  // Ingest one of each type of table.
  PulsarDb database(m_tpl_file);
  database.load(prependDataPath("psrdb_spin.txt"));
  database.load(prependDataPath("psrdb_binary.txt"));
  database.load(prependDataPath("psrdb_remark.txt"));
  database.load(prependDataPath("psrdb_obs.txt"));
  database.load(prependDataPath("psrdb_name.txt"));

  // Save all tables into one FITS file.
  std::string filename1("psrdb_all.fits");
  remove(filename1.c_str());
  database.save(filename1, m_creator, m_author);

  // Check the result against its reference file in data/outref/ directory.
  tester.checkOutputFits(filename1, prependOutrefPath(filename1));

  // Copy an existing FITS database.
  PulsarDb fits_psrdb(m_tpl_file);
  fits_psrdb.load("crab_db.fits");

  // Append all tables into the FITS database.
  fits_psrdb.load(filename1);

  // Save all tables into another FITS file.
  std::string filename2 = "psrdb_append.fits";
  remove(filename2.c_str());
  fits_psrdb.save(filename2, m_creator, m_author);

  // Check the result against its reference file in data/outref/ directory.
  tester.checkOutputFits(filename2, prependOutrefPath(filename2));

  // Test the getter of history records.
  std::list<std::string> expected_command;
  expected_command.push_back("Load FITSDB AUTHOR='Anonymous Tester' DATE=");
  expected_command.push_back("Load FITSDB AUTHOR='Anonymous Tester' DATE=");
  std::list<std::string> expected_ancestry;
  expected_ancestry.push_back("PULSARDB AUTHOR='Anonymous Tester' DATE=");
  expected_ancestry.push_back("* Load FITSDB AUTHOR='' DATE=");
  expected_ancestry.push_back("* Filter by pulsar name 'CRab'");
  expected_ancestry.push_back("* Filter by pulsar name 'PsR b0531+21'");
  expected_ancestry.push_back("PULSARDB AUTHOR='Anonymous Tester' DATE=");
  expected_ancestry.push_back("* Load TEXTDB SPIN_PARAMETERS(FREQ) FILENAME='psrdb_spin.txt'");
  expected_ancestry.push_back("* Load TEXTDB ORBITAL_PARAMETERS(DD) FILENAME='psrdb_binary.txt'");
  expected_ancestry.push_back("* Load TEXTDB REMARKS FILENAME='psrdb_remark.txt'");
  expected_ancestry.push_back("* Load TEXTDB OBSERVERS FILENAME='psrdb_obs.txt'");
  expected_ancestry.push_back("* Load TEXTDB ALTERNATIVE_NAMES FILENAME='psrdb_name.txt'");
  std::list<std::string> result_command;
  std::list<std::string> result_ancestry;
  fits_psrdb.getHistory(result_command, result_ancestry);

  // Compare the command history.
  if (result_command.size() != expected_command.size()) {
    err() << "PulsarDb::getHistory returned " << result_command.size() << " commands, not " <<
      expected_command.size() << " as expected." << std::endl;
  } else {
    std::list<std::string>::const_iterator exp_itor = expected_command.begin();
    std::list<std::string>::const_iterator res_itor = result_command.begin();
    int line_index = 0;
    for (; exp_itor != expected_command.end() && res_itor != result_command.end(); ++exp_itor, ++res_itor, ++line_index) {
      std::string res_string = res_itor->substr(0, exp_itor->size());
      if (res_string != *exp_itor) {
        err() << "PulsarDb::getHistory returned '" << *res_itor << "', not '" << *exp_itor <<
          "' as expected for command No. " << line_index + 1 << "." << std::endl;
      }
    }
  }

  // Compare the ancestry records.
  if (result_ancestry.size() != expected_ancestry.size()) {
    err() << "PulsarDb::getHistory returned " << result_ancestry.size() << " ancestry records, not " <<
      expected_ancestry.size() << " as expected." << std::endl;
  } else {
    std::list<std::string>::const_iterator exp_itor = expected_ancestry.begin();
    std::list<std::string>::const_iterator res_itor = result_ancestry.begin();
    int line_index = 0;
    for (; exp_itor != expected_ancestry.end() && res_itor != result_ancestry.end(); ++exp_itor, ++res_itor, ++line_index) {
      std::string res_string = res_itor->substr(0, exp_itor->size());
      if (res_string != *exp_itor) {
        err() << "PulsarDb::getHistory returned '" << *res_itor << "', not '" << *exp_itor <<
          "' as expected for ancestry record No. " << line_index + 1 << "." << std::endl;
      }
    }
  }

  // Test detection of badly-formatted tables.
  std::list<std::string> filename_list;
  filename_list.push_back("no_such_file.txt");
  filename_list.push_back("baddb_noextname.txt");
  filename_list.push_back("baddb_nosuchextname.txt");
  filename_list.push_back("baddb_nosuchephstyle.txt");
  filename_list.push_back("baddb_noheader.txt");
  filename_list.push_back("baddb_nosuchfield.txt");
  filename_list.push_back("baddb_toolittlefield.txt");
  filename_list.push_back("baddb_toomanyfield.txt");
  filename_list.push_back("baddb_unbalancedquote.txt");
  for (std::list<std::string>::const_iterator itor = filename_list.begin(); itor != filename_list.end(); ++itor) {
    const std::string & filename(*itor);
    try {
      database.load(prependDataPath(filename));
      err() << "PulsarDb::load method did not throw an exception for TEXT file \"" << filename << "\"" << std::endl;
    } catch (const std::exception &) {
      // This is fine.
    }
  }

  // Test no detection of poorly-formatted, but legal tables.
  filename_list.clear();
  filename_list.push_back("okdb_extraspace.txt");
  filename_list.push_back("okdb_extraquote.txt");
  for (std::list<std::string>::const_iterator itor = filename_list.begin(); itor != filename_list.end(); ++itor) {
    const std::string & filename(*itor);
    try {
      database.load(prependDataPath(filename));
    } catch (const std::exception & x) {
      err() << "PulsarDb::load method threw an exception for TEXT file \"" << filename << "\": " << x.what() << std::endl;
    }
  }
}

class TextEph: public FormattedEph {
  public:
    TextEph() {}
    virtual ~TextEph() {}

    void append(const std::string & par_name, const std::string & par_value) {
      m_par_list.push_back(std::make_pair(par_name, par_value));
    }

    template <typename DataType>
    void append(const std::string & par_name, const DataType & par_value) {
      std::ostringstream oss;
      oss.precision(std::numeric_limits<double>::digits10);
      oss << par_value;
      append(par_name, oss.str());
    }

    void clear() {
      m_par_list.clear();
    }

    virtual st_stream::OStream & write(st_stream::OStream & os) const {
      for (plist_type::const_iterator itor = m_par_list.begin(); itor != m_par_list.end(); ++itor) {
        // Place a carriage return at the end of line, except for the last line.
        if (itor != m_par_list.begin()) os << std::endl;

        // Write parameter value.
        const std::string & par_name(itor->first);
        const std::string & par_value(itor->second);
        if (("Valid Since" == par_name) || ("Valid Until" == par_name)) {
          os << format(par_name, par_value, " : ");
        } else {
          os << format(par_name, par_value);
        }
      }
      return os;
    }

  private:
    typedef std::list<std::pair<std::string, std::string> > plist_type;
    plist_type m_par_list;
};

void PulsarDbTestApp::testFrequencyEph() {
  setMethod("testFrequencyEph");

  // Prepare absolute times and tolerance for testing computations.
  AbsoluteTime since("TT", 51910, 0.);
  AbsoluteTime until("TT", 51910, 1.);
  AbsoluteTime epoch("TT", 51910, 123.456789);
  double epsilon = 1.e-8;

  // Create a frequency ephemeris.
  FrequencyEph eph1("TDB", since, until, epoch, 22., 45., 0.11, 1.125e-2, -2.25e-4, 6.75e-6);

  // Create a time to pick an ephemeris.
  AbsoluteTime pick_time("TT", 51910, 223.456789);

  // Test pulse phase computations.
  double phase = eph1.calcPulsePhase(pick_time);
  if (fabs(phase/.235 - 1.) > epsilon)
    err() << "FrequencyEph::calcPulsePhase produced phase == " << phase << " not .235" << std::endl;
 
  // Test pulse phase computations, with a non-zero global phase offset.
  phase = eph1.calcPulsePhase(pick_time, 0.1234);
  if (fabs(phase/.3584 - 1.) > epsilon)
    err() << "FrequencyEph::calcPulsePhase produced phase == " << phase << " not .3584" << std::endl;
 
  // Change ephemeris to produce a noticeable effect.
  epoch = AbsoluteTime("TT", 51910, 123.4567891234567);
  FrequencyEph eph2("TDB", since, until, epoch, 22., 45., .11, 1.125e-2, -2.25e-4, 13.5e-6);
  AbsoluteTime ev_time("TT", 51910, 223.4567891234567);

  // Test coordinates.
  std::pair<double, double> computed_ra_dec = eph2.calcSkyPosition(ev_time);
  double computed_ra = computed_ra_dec.first;
  double correct_ra = 22.;
  if (fabs(computed_ra - correct_ra) > epsilon) {
    err() << "FrequencyEph::calcSkyPosition produced ra == " << computed_ra << " not " << correct_ra << std::endl;
  }
  
  double computed_dec = computed_ra_dec.second;
  double correct_dec = 45.;
  if (fabs(computed_dec - correct_dec) > epsilon) {
    err() << "FrequencyEph::calcSkyPosition produced dec == " << computed_dec << " not " << correct_dec << std::endl;
  }
  
  // Test phase computation.
  double computed_phi0 = eph2.calcPulsePhase(ev_time);
  if (fabs(computed_phi0/.36 - 1.) > epsilon)
    err() << "FrequencyEph::calcPulsePhase produced phi0 == " << computed_phi0 << " not .36" << std::endl;
 
  // Test frequency computation.
  double computed_f0 = eph2.calcFrequency(ev_time, 0);
  double correct_f0 = 5.625e-2;
  if (fabs(computed_f0 - correct_f0) > epsilon) {
    err() << "FrequencyEph::calcFrequency produced f0 == " << computed_f0 << " not " << correct_f0 << std::endl;
  }
  
  double computed_f1 = eph2.calcFrequency(ev_time, 1);
  double correct_f1 = 11.25e-4;
  if (fabs(computed_f1 - correct_f1) > epsilon) {
    err() << "FrequencyEph::calcFrequency produced f1 == " << computed_f1 << " not " << correct_f1 << std::endl;
  }
  
  double computed_f2 = eph2.calcFrequency(ev_time, 2);
  double correct_f2 = 13.5e-6;
  if (fabs(computed_f2 - correct_f2) > epsilon) {
    err() << "FrequencyEph::calcFrequency produced f2 == " << computed_f2 << " not " << correct_f2 << std::endl;
  }

  // Test the constructor that takes numerical arguments.
  FrequencyEph eph3("TDB", AbsoluteTime("TDB", 12345, 51840.), AbsoluteTime("TDB", 23456, 60480.), AbsoluteTime("TDB", 34567, 69120.),
    11., 22., 33., 44., 55., 66.);
  TextEph t_eph;
  t_eph.append("Valid Since", "12345.6 MJD (TDB)");
  t_eph.append("Valid Until", "23456.7 MJD (TDB)");
  t_eph.append("Epoch",       "34567.8 MJD (TDB)");
  t_eph.append("RA",   "11");
  t_eph.append("Dec",  "22");
  t_eph.append("Phi0", "33");
  t_eph.append("F0",   "44");
  t_eph.append("F1",   "55");
  t_eph.append("F2",   "66");
  checkEphParameter(getMethod() + "_numeric", eph3, t_eph);

  // Test the constructor that takes a FITS record.
  tip::TipFile tip_file = tip::IFileSvc::instance().createMemFile(getMethod() + ".fits", prependDataPath("test_FrequencyEph.tpl"));
  tip::Table * table = tip_file.editTable("1");
  table->setNumRecords(1);
  tip::Header & header(table->getHeader());
  tip::Table::Record & record(*(table->begin()));
  record["VALID_SINCE"].set(54321);
  record["VALID_UNTIL"].set(65432);
  record["EPOCH_INT"].set(76543);
  record["EPOCH_FRAC"].set(.5);
  record["TOABARY_INT"].set(76543);
  record["TOABARY_FRAC"].set(.50001); // 0.864 seconds off.
  record["RA"].set(1.1);
  record["DEC"].set(2.2);
  record["F0"].set(3.3);
  record["F1"].set(4.4);
  record["F2"].set(5.5);
  FrequencyEph eph4(record, header);
  t_eph.clear();
  t_eph.append("Valid Since", "54321 MJD (TDB)");
  t_eph.append("Valid Until", "65433 MJD (TDB)"); // The end of 65432 == the begining of 65433.
  t_eph.append("Epoch",       "76543.5 MJD (TDB)");
  t_eph.append("RA",   "1.1");
  t_eph.append("Dec",  "2.2");
  double dt = .864;
  double phi0 = -(3.3*dt + 4.4/2.*dt*dt + 5.5/6.*dt*dt*dt);
  while (phi0 < 0.) ++phi0;
  while (phi0 >= 1.) --phi0;
  t_eph.append("Phi0", phi0);
  t_eph.append("F0",   "3.3");
  t_eph.append("F1",   "4.4");
  t_eph.append("F2",   "5.5");
  checkEphParameter(getMethod() + "_fits", eph4, t_eph);
}

class SimplePeriodEph {
  public:
    SimplePeriodEph(double phi0, double p0, double p1, double p2): m_phi0(phi0), m_p0(p0), m_p1(p1), m_p2(p2) {}

    double calcPulsePhase(double dt, double step, double global_offset) const {
      // Check and normalize step size.
      if (step == 0.) throw std::runtime_error("Bad test parameter given: step size for differenciation is zero");
      step = std::fabs(step);

      // Set integration boundary.
      double t0 = 0.;
      double t1 = 0.;
      if (dt < 0.) {
        t0 = dt;
        t1 = 0.;
      } else if (dt > 0.) {
        t0 = 0;
        t1 = dt;
      } else {
        t0 = 0.;
        t1 = 0.;
      }
      std::list<double> tx;
      for (int ii=0; ii < static_cast<int>(std::ceil(dt/step)) - 1; ++ii) tx.push_back(step*ii);
      tx.push_back(t1);

      // Perform integration.
      double t_left = t0;
      double v_left = 1. / calcPeriod(t0);
      double t_right = 0.;
      double v_right = 0.;
      double phase = m_phi0 + global_offset;
      for (std::list<double>::const_iterator itor = tx.begin(); itor != tx.end(); ++itor) {
        t_right = *itor;
        v_right = 1. / calcPeriod(t_right);
        phase += (v_left + v_right) * (t_right - t_left) / 2.;

        if (phase < 0. || phase >= 1.) {
          double int_part; // ignored, needed for modf.
          phase = std::modf(phase, &int_part);
          if (phase < 0.) ++phase;
        }

        t_left = t_right;
        v_left = v_right;
      }
      return phase;
    }

    double calcFrequency(double dt, double step, int order) const {
      // Check and normalize step size.
      if (step == 0.) throw std::runtime_error("Bad test parameter given: step size for differenciation is zero");
      step = std::fabs(step);

      // Perform differenciation recursively.
      if (order < 0) throw std::runtime_error("Bad test parameter given: a negative derivative order");
      if (order == 0) return 1. / calcPeriod(dt);
      else return (calcFrequency(dt+step/2., step, order-1) - calcFrequency(dt-step/2., step, order-1)) / step;
    }

  private:
    double calcPeriod(double dt) const { return m_p0 + m_p1*dt + m_p2/2.*dt*dt; }
    double m_phi0, m_p0, m_p1, m_p2;
};

void PulsarDbTestApp::testPeriodEph() {
  setMethod("testPeriodEph");

  AbsoluteTime since("TDB", 51910, 0.);
  AbsoluteTime until("TDB", 51910, 1.);
  AbsoluteTime epoch("TDB", 51910, 123.456789);

  // Create a frequency ephemeris.
  FrequencyEph f_eph("TDB", since, until, epoch, 22., 45., 0.875, 1.125e-2, -2.25e-4, 6.75e-6);

  // Create a period ephemeris.
  // This is a set of values known to be the inverses of the frequency coefficients above.
  PeriodEph p_eph("TDB", since, until, epoch, 22., 45., 0.875, 88.8888888888888888888889,
    1.777777777777777777777778, 0.0177777777777777777777778);

  // Compare frequency & period.
  double epsilon = 1.e-8;
  const double nano_sec = 1.e-9;
  ElapsedTime tolerance("TDB", Duration(nano_sec, "Sec"));

  if (!f_eph.getEpoch().equivalentTo(p_eph.getEpoch(), tolerance))
    err() << "FrequencyEph and PeriodEph give different values for epoch" << std::endl;

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
        err() << "FrequencyEph and PeriodEph give absolutely different values for " << field[ii] <<
          " (" << value1[ii] << " and " << value2[ii] << ")" << std::endl;
    } else if (fabs(value1[ii] / value2[ii] - 1.) > epsilon) {
      err() << "FrequencyEph and PeriodEph give fractionally different values for " << field[ii] <<
        " (" << value1[ii] << " and " << value2[ii] << ")" << std::endl;
    }
  }

  // Test frequency computation away from the reference epoch.
  double time_since_epoch = 1000.;
  double step_size = 100.;
  AbsoluteTime abs_time = epoch + ElapsedTime("TDB", Duration(time_since_epoch, "Sec"));
  double ra = 22.;
  double dec = 45.;
  double phi0 = 0.875;
  double p0 = 1.23456789;
  double p1 = 9.87654321e-5;
  double p2 = 1.357902468e-10;
  PeriodEph eph1("TDB", since, until, epoch, ra, dec, phi0, p0, p1, p2);
  SimplePeriodEph s_eph1(phi0, p0, p1, p2);

  epsilon = 1.e-3; // Note: Need a loose tolerance because SimplePeriodEph::calcFrequency is not that precise.
  double result = 0.;
  double expected = 0.;
  bool test_failed = true;
  for (int order = 0; order < 5; ++order) {
    result = eph1.calcFrequency(abs_time, order);
    expected = s_eph1.calcFrequency(time_since_epoch, step_size, order);
    if (0. == result || 0. == expected) test_failed = fabs(result + expected) > std::numeric_limits<double>::epsilon();
    else test_failed = std::fabs(result / expected - 1.) > epsilon;
    if (test_failed) {
      err() << "PeriodEph::calcFrequency(abs_time, " << order << ") returned " << result <<
        ", not " << expected << " as expected" << std::endl;
    }
  }

  // Test pulse phase computation away from the reference epoch, with and without a non-zero global offset.
  step_size = time_since_epoch / 1000.;
  double global_phase_offset = 0.34567;

  double phase_tolerance = 1.e-5;
  result = eph1.calcPulsePhase(abs_time, global_phase_offset);
  expected = s_eph1.calcPulsePhase(time_since_epoch, step_size, global_phase_offset);
  test_failed = (fabs(result - expected) > phase_tolerance && fabs(result - expected + 1) > phase_tolerance
    && fabs(result - expected - 1) > phase_tolerance);
  if (test_failed) {
    err() << "PeriodEph::calcPulsePhase(abs_time, " << global_phase_offset << ") returned " << result <<
      ", not close enough to " << expected << " as expected" << std::endl;
  }

  global_phase_offset = 0.;
  result = eph1.calcPulsePhase(abs_time, global_phase_offset);
  expected = s_eph1.calcPulsePhase(time_since_epoch, step_size, global_phase_offset);
  test_failed = (fabs(result - expected) > phase_tolerance && fabs(result - expected + 1) > phase_tolerance
    && fabs(result - expected - 1) > phase_tolerance);
  if (test_failed) {
    err() << "PeriodEph::calcPulsePhase(abs_time, " << global_phase_offset << ") returned " << result <<
      ", not close enough to " << expected << " as expected" << std::endl;
  }

  // Test pulse phase computation for artificial parameters to cover all cases of frequency integration.
  phase_tolerance = 1.e-3;
  double good_period_par[][3] = {
    {0.123456789, 0., 0.}, {123.456789, -1.23456789e-2, 0.},
    {2., 1.e-5, 1.e-10}, {2., 1.9e-5, 1.e-10}, {2., 2.e-5, 1.e-10}
    // Note: Branch for 2 * p0p2 == p1^2 is almost impossible to test; difficult to achieve the exact equality.
    //       Parameters {2., 2.e-5, 1.e-10} and {1., 1.e-5, 2.e-10} were tried on Linux, but either of them
    //       did not achieve the equality as desired.
  };
  for (size_t ii = 0; ii != sizeof (good_period_par)/sizeof(double[3]); ++ii) {
    p0 = good_period_par[ii][0];
    p1 = good_period_par[ii][1];
    p2 = good_period_par[ii][2];
    PeriodEph eph2("TDB", since, until, epoch, ra, dec, phi0, p0, p1, p2);
    SimplePeriodEph s_eph2(phi0, p0, p1, p2);
    result = eph2.calcPulsePhase(abs_time);
    expected = s_eph2.calcPulsePhase(time_since_epoch, step_size, 0.);
    test_failed = (fabs(result - expected) > phase_tolerance && fabs(result - expected + 1) > phase_tolerance
      && fabs(result - expected - 1) > phase_tolerance);
    if (test_failed) {
      err() << "PeriodEph::calcPulsePhase(abs_time, 0.) returned " << result <<
        ", not close enough to " << expected << " as expected" << std::endl;
    }
  }

  // Test pulse phase computation for problematic parameters.
  double bad_period_par[][3] = {
    {0., 0., 0.}, {123.456789, -1.23456789, 0.}, {2., -2.e-2, 1.e-8}
    // Note: Branch for 2 * p0p2 == p1^2 is almost impossible to test; difficult to achieve the exact equality.
    //       Parameters {2., -2.e-2, 1.e-4} and {2., 2.e-2, -1.e-4} were tried on Linux, but either of them
    //       did not achieve the equality as desired.
  };
  for (size_t ii = 0; ii != sizeof (bad_period_par)/sizeof(double[3]); ++ii) {
    p0 = bad_period_par[ii][0];
    p1 = bad_period_par[ii][1];
    p2 = bad_period_par[ii][2];
    PeriodEph eph3("TDB", since, until, epoch, ra, dec, phi0, p0, p1, p2);
    try {
      eph3.calcPulsePhase(abs_time);
      err() << "PeriodEph::calcPulsePhase(abs_time) did not throw an exception for p0=" << p0 <<
        ", p1=" << p1 << ", p2=" << p2 << std::endl;
    } catch (const std::exception &) {
      // This is fine.
    }
  }

  // Test the constructor that takes numerical arguments.
  PeriodEph eph4("TDB", AbsoluteTime("TDB", 12345, 51840.), AbsoluteTime("TDB", 23456, 60480.), AbsoluteTime("TDB", 34567, 69120.),
    11., 22., 33., 44., 55., 66.);
  TextEph t_eph;
  t_eph.append("Valid Since", "12345.6 MJD (TDB)");
  t_eph.append("Valid Until", "23456.7 MJD (TDB)");
  t_eph.append("Epoch",       "34567.8 MJD (TDB)");
  t_eph.append("RA",   "11");
  t_eph.append("Dec",  "22");
  t_eph.append("Phi0", "33");
  t_eph.append("P0",   "44");
  t_eph.append("P1",   "55");
  t_eph.append("P2",   "66");
  checkEphParameter(getMethod() + "_numeric", eph4, t_eph);

  // Test the constructor that takes a FITS record.
  tip::TipFile tip_file = tip::IFileSvc::instance().createMemFile(getMethod() + ".fits", prependDataPath("test_PeriodEph.tpl"));
  tip::Table * table = tip_file.editTable("1");
  table->setNumRecords(1);
  tip::Header & header(table->getHeader());
  tip::Table::Record & record(*(table->begin()));
  record["VALID_SINCE"].set(54321);
  record["VALID_UNTIL"].set(65432);
  record["EPOCH_INT"].set(76543);
  record["EPOCH_FRAC"].set(.5);
  record["TOABARY_INT"].set(76543);
  record["TOABARY_FRAC"].set(.50001); // 0.864 seconds off.
  record["RA"].set(1.1);
  record["DEC"].set(2.2);
  record["P0"].set(3.3);
  record["P1"].set(4.4);
  record["P2"].set(5.5);
  PeriodEph eph5(record, header);
  t_eph.clear();
  t_eph.append("Valid Since", "54321 MJD (TDB)");
  t_eph.append("Valid Until", "65433 MJD (TDB)"); // The end of 65432 == the begining of 65433.
  t_eph.append("Epoch",       "76543.5 MJD (TDB)");
  t_eph.append("RA",   "1.1");
  t_eph.append("Dec",  "2.2");
  double dt = .864;
  double sqrt_term = std::sqrt(2.*3.3*5.5 - 4.4*4.4);
  phi0 = -(2. / sqrt_term * std::atan(sqrt_term*dt/(2.*3.3+4.4*dt)));
  while (phi0 < 0.) ++phi0;
  while (phi0 >= 1.) --phi0;
  t_eph.append("Phi0", phi0);
  t_eph.append("P0",   "3.3");
  t_eph.append("P1",   "4.4");
  t_eph.append("P2",   "5.5");
  checkEphParameter(getMethod() + "_fits", eph5, t_eph);
}

void PulsarDbTestApp::testSimpleDdEph() {
  setMethod("testSimpleDdEph");

  AbsoluteTime t0("TDB", 51910, 123.456789);
  SimpleDdEph eph1("TDB", 1000., .2, 0., 0., 0., 0., 0., 0., t0, 0., 0., 0.);
  AbsoluteTime ev_time("TDB", 51910, 223.456789);

  double epsilon = 1.e-8;

  // Test time system getter.
  std::string sys_name = eph1.getSystem().getName();
  if ("TDB" != sys_name) {
    err() << "SimpleDdEph::getSystem returned \"" << sys_name << "\", not \"TDB\"" << std::endl;
  }

  // Test orbital phase computations.
  double phase = eph1.calcOrbitalPhase(ev_time);
  if (fabs(phase/.099 - 1.) > epsilon)
    err() << "SimpleDdEph::calcOrbitalPhase produced phase == " << phase << " not .099" << std::endl;

  // Test orbital phase computations, with a non-zero global phase offset.
  phase = eph1.calcOrbitalPhase(ev_time, 0.1234);
  if (fabs(phase/.2224 - 1.) > epsilon)
    err() << "SimpleDdEph::calcOrbitalPhase produced phase == " << phase << " not .2224" << std::endl;

  // Test binary modulation and demodulation.
  // Binary parameters: (PB, PBDOT, A1, XDOT, ECC, ECCDOT, OM, OMDOT, T0, GAMMA, SHAPIRO_R, SHAPIRO_S)
  // = (27906.980897, -2.43e-12, 2.3417598, 0.0, 0.61713101, 0.0, 220.142729, 4.22662, 45888.1172487, 0.004295, 0.0, 0.0)
  // Note: OMDOT has been changed to 4.22662*365.25/365.0 to adopt the bug fix (by M. Hirayama on March 16th, 2010).
  AbsoluteTime abs_t0("TDB", Mjd(45888, .1172487));
  SimpleDdEph eph2("TDB", 27906.980897, -2.43e-12, 2.3417598, 0.0, 0.61713101, 0.0, 220.142729, 4.22662 * 365.25 / 365.0,
                   abs_t0, 0.004295, 0.0, 0.0);

  // MJD's: { {original-MJD, modulated-MJD}, ... }
  Mjd mjd_test_values[][2] = {
    { Mjd(45988, .00001172486599971307e+04),
      Mjd(45988, .00001172823346967320e+04) },
    { Mjd(45988, .00001519708822161192e+04),
      Mjd(45988, .00001520057480382526e+04) },
    { Mjd(45988, .00001866931044423836e+04),
      Mjd(45988, .00001867233097334768e+04) },
    { Mjd(45988, .00002214153266613721e+04),
      Mjd(45988, .00002214303721880313e+04) },
    { Mjd(45988, .00002561375488876365e+04),
      Mjd(45988, .00002561254377882811e+04) },
    { Mjd(45988, .00002908597711066250e+04),
      Mjd(45988, .00002908532530402717e+04) },
    { Mjd(45988, .00003255819933328894e+04),
      Mjd(45988, .00003255877721779576e+04) },
    { Mjd(45988, .00003603042155518779e+04),
      Mjd(45988, .00003603213319299883e+04) },
    { Mjd(45988, .00003950264377781423e+04),
      Mjd(45988, .00003950526724878252e+04) },
    { Mjd(45988, .00004297486599971307e+04),
      Mjd(45988, .00004297811300504257e+04) },
    { Mjd(45988, .00004644708822161192e+04),
      Mjd(45988, .00004645058955188297e+04) }
  };

  // Permitted difference is 100 ns.
  double delta = 100. * 1.e-9;
  ElapsedTime tolerance("TDB", Duration(delta, "Sec"));

  std::cerr.precision(24);
  for (size_t ii = 0; ii != sizeof(mjd_test_values)/sizeof(Mjd[2]); ++ii) {
    AbsoluteTime tdb_mjd("TDB", mjd_test_values[ii][0]);
    AbsoluteTime expected_tdb_mjd("TDB", mjd_test_values[ii][1]);
    AbsoluteTime original_tdb_mjd(tdb_mjd);
    eph2.modulateBinary(tdb_mjd);
    if (!tdb_mjd.equivalentTo(expected_tdb_mjd, tolerance)) {
      err() << "Binary modulation of " << original_tdb_mjd.represent("TDB", MjdFmt) << " was computed to be " <<
        tdb_mjd.represent("TDB", MjdFmt) << ", not " << expected_tdb_mjd.represent("TDB", MjdFmt) << ", as expected." << std::endl;
    }
  }

  for (size_t ii = 0; ii != sizeof(mjd_test_values)/sizeof(Mjd[2]); ++ii) {
    AbsoluteTime tdb_mjd("TDB", mjd_test_values[ii][1]);
    AbsoluteTime expected_tdb_mjd("TDB", mjd_test_values[ii][0]);
    AbsoluteTime original_tdb_mjd(tdb_mjd);
    eph2.demodulateBinary(tdb_mjd);
    if (!tdb_mjd.equivalentTo(expected_tdb_mjd, tolerance)) {
      err() << "Binary demodulation of " << original_tdb_mjd.represent("TDB", MjdFmt) << " was computed to be " <<
        tdb_mjd.represent("TDB", MjdFmt) << ", not " << expected_tdb_mjd.represent("TDB", MjdFmt) << ", as expected." << std::endl;
    }
  }

  // Test the constructor that takes numerical arguments.
  SimpleDdEph eph3("TDB", 12345.6789, 1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7, AbsoluteTime("TDB", 12345, 51840.), 8.8, 9.9, 11.1);
  TextEph t_eph;
  t_eph.append("PB",        "12345.6789");
  t_eph.append("PBDOT",     "1.1");
  t_eph.append("A1",        "2.2");
  t_eph.append("XDOT",      "3.3");
  t_eph.append("ECC",       "4.4");
  t_eph.append("ECCDOT",    "5.5");
  t_eph.append("OM",        6.6 * SimpleDdEph::s_rad_per_deg);
  t_eph.append("OMDOT",     7.7 * SimpleDdEph::s_rad_year_per_deg_sec);
  t_eph.append("T0",        "12345.6");
  t_eph.append("GAMMA",     "8.8");
  t_eph.append("SHAPIRO_R", "9.9e-06");
  t_eph.append("SHAPIRO_S", "11.1");
  checkEphParameter(getMethod() + "_numeric", eph3, t_eph);

  // Test the constructor that takes a FITS record.
  tip::TipFile tip_file = tip::IFileSvc::instance().createMemFile(getMethod() + ".fits", prependDataPath("test_SimpleDdEph.tpl"));
  tip::Table * table = tip_file.editTable("1");
  table->setNumRecords(1);
  tip::Header & header(table->getHeader());
  tip::Table::Record & record(*(table->begin()));
  record["PB"].set(12345.6789);
  record["PBDOT"].set(1.1);
  record["A1"].set(2.2);
  record["XDOT"].set(3.3);
  record["ECC"].set(4.4);
  record["ECCDOT"].set(5.5);
  record["OM"].set(6.6);
  record["OMDOT"].set(7.7);
  record["T0"].set(12345.6);
  record["GAMMA"].set(8.8);
  record["SHAPIRO_R"].set(9.9);
  record["SHAPIRO_S"].set(11.1);
  SimpleDdEph eph4(record, header);
  checkEphParameter(getMethod() + "_fits", eph4, t_eph);
} 

void PulsarDbTestApp::testPdotCanceler() {
  setMethod("testPdotCanceler");

  // Set test parameters.
  AbsoluteTime origin("TT", 51910, 123.4567891234567);
  AbsoluteTime ev_time1("TT", 51910, 223.4567891234567);
  AbsoluteTime ev_time2("TT", 51910, 223.4567891234567);
  double f0 = 1.125e-2;
  double f1 = -2.25e-4;
  double f2 = 13.5e-6;
  AbsoluteTime correct_time("TT", 51910, 323.4567891234567);
  ElapsedTime tolerance("TT", Duration(1.e-6, "Sec")); // 1 microsecond.

  // Test PdotCanceler created from literal numbers.
  std::vector<double> fdot_ratio(2);
  fdot_ratio[0] = f1 / f0;
  fdot_ratio[1] = f2 / f0;
  PdotCanceler canceler1("TDB", origin, fdot_ratio);

  canceler1.cancelPdot(ev_time1);
  if (!correct_time.equivalentTo(ev_time1, tolerance)) {
    err() << "After constructed from literal numbers, PdotCanceler::cancelPdot produced pdot-corrected time == "
      << ev_time1 << " not " << correct_time << std::endl;
  }

  // Test PdotCanceler created from a PulsarEph object.
  AbsoluteTime since("TT", 51910, 0.);
  AbsoluteTime until("TT", 51910, 1.);
  FrequencyEph f_eph("TDB", since, until, origin, 22., 45., .11, f0, f1, f2);
  PdotCanceler canceler2(origin, f_eph, 2);

  canceler2.cancelPdot(ev_time2);
  if (!correct_time.equivalentTo(ev_time2, tolerance)) {
    err() << "After constructed from PulsarEph, PdotCanceler::cancelPdot produced pdot-corrected time == "
      << ev_time2 << " not " << correct_time << std::endl;
  }

  // Test time system getter.
  PdotCanceler canceler3("tdB", origin, fdot_ratio);
  std::string result = canceler3.getSystem().getName();
  std::string expected = "TDB";
  if (result != expected) {
    err() << "After constructed from literal numbers, PdotCanceler::getSystem returned \"" << result << "\", not \""
      << expected << "\"" << std::endl;
  }
}

void PulsarDbTestApp::testChooser() {
  setMethod("testChooser");
  PulsarDbAppTester tester(*this);

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
    err() << "there are " << num_eph << " ephemerides for " << pulsar_name << ", not 8" << std::endl;

  // Write this output to form basis for comparing future tests.
  std::string outfile("chooser_db.fits");
  remove(outfile.c_str());
  database.save(outfile, m_creator, m_author);

  // Check the result against its reference file in data/outref/ directory.
  tester.checkOutputFits(outfile, prependOutrefPath(outfile));

  AbsoluteTime pick_time("TDB", Mjd1(54012.5));
  AbsoluteTime expected_epoch("TDB", 0, 0.);

  PulsarEphCont eph_cont;

  database.getEph(eph_cont);

  StrictEphChooser chooser;

  // Test one with no tiebreaking needed.
  const PulsarEph * chosen = &chooser.choose(eph_cont, pick_time);
  expected_epoch.set("TDB", Mjd1(54262.));
  ElapsedTime tolerance("TDB", Duration(1.e-9, "Sec")); // 1 nanosecond.
  if (!expected_epoch.equivalentTo(chosen->getEpoch(), tolerance))
    err() << "for time " << pick_time << ", chooser chose ephemeris with EPOCH == " << chosen->getEpoch() <<
      ", not " << expected_epoch << " as expected." << std::endl;

  // Test one with tiebreaking.
  pick_time.set("TDB", Mjd1(53545.5));
  chosen = &chooser.choose(eph_cont, pick_time);
  expected_epoch.set("TDB", Mjd1(53891.));
  if (!expected_epoch.equivalentTo(chosen->getEpoch(), tolerance))
    err() << "for time " << pick_time << ", chooser chose ephemeris with EPOCH == " << chosen->getEpoch() <<
      ", not " << expected_epoch << " as expected." << std::endl;

  // Test one which is too early.
  pick_time.set("TDB", Mjd1(53544.5));
  try {
    chosen = &chooser.choose(eph_cont, pick_time);
    err() << "for time " << pick_time << ", chooser chose ephemeris with EPOCH == " << chosen->getEpoch() << std::endl;
  } catch (const std::runtime_error &) {
    // This is to be expected.
  }

  // Test one which is too late.
  pick_time.set("TDB", Mjd1(55579.5));
  try {
    chosen = &chooser.choose(eph_cont, pick_time);
    err() << "for time " << pick_time << ", chooser chose ephemeris with EPOCH == " << chosen->getEpoch() << std::endl;
  } catch (const std::runtime_error &) {
    // This is to be expected.
  }

  // Try one which is too late, but without being strict about validity.
  pick_time.set("TDB", Mjd1(55579.5));
  try {
    chosen = &(SloppyEphChooser().choose(eph_cont, pick_time));
  } catch (const std::runtime_error &) {
    err() << "for time " << pick_time << ", chooser did not find an ephemeris even with strict_validity == false" <<
      std::endl;
  }

  // Make a selection which will result in an empty container of ephemerides.
  database.filterName("Aunt Gertrude");
  if (0 != database.getNumEph()) {
    err() << "What? Aunt Gertrude is a PULSAR???" << std::endl;
  }
  
  database.getEph(eph_cont);

  // Try to choose an ephemeris from the empty set.
  pick_time.set("TDB", Mjd1(55579.5));
  try {
    chosen = &chooser.choose(eph_cont, pick_time);
    err() << "chooser chose ephemeris from an empty set of candidates." << std::endl;
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
  pick_time.set("TDB", Mjd1(52500.));
  try {
    const OrbitalEph & orbital_eph = chooser.choose(orbital_cont, pick_time);
    AbsoluteTime expected_t0("TDB", Mjd1(52060.84100795));
    ElapsedTime tolerance("TT", Duration(1.e-9, "Sec")); // 1 nanosecond.
    if (!orbital_eph.t0().equivalentTo(expected_t0, tolerance)) {
      err() << "chooser chose orbital ephemeris with time " << orbital_eph.t0() << ", not " <<
        expected_t0 << std::endl;
    }
  } catch (const std::runtime_error & x) {
    err() << "for time " << pick_time << ", chooser had trouble choosing orbital eph: " <<
    x.what() << std::endl;
  }

  // Clean up.
  for (PulsarEphCont::reverse_iterator itor = eph_cont.rbegin(); itor != eph_cont.rend(); ++itor) delete *itor;
  eph_cont.clear();

  // Test tiebreaking with different tolerances.
  long origin = 51910;
  AbsoluteTime valid_since("TDB", origin, 101.);
  AbsoluteTime valid_until("TDB", origin, 200.);
  AbsoluteTime epoch("TDB", origin, 150.);
  eph_cont.push_back(new FrequencyEph("TT", valid_since, valid_until, epoch, 22., 45., 0., 1., 0., 0.));
  valid_since = AbsoluteTime("TDB", origin, 100.);
  eph_cont.push_back(new FrequencyEph("TDB", valid_since, valid_until, epoch, 22., 45., 0., 1., 0., 0.));

  StrictEphChooser strict_chooser(ElapsedTime("TDB", Duration(.99, "Sec")));
  pick_time = AbsoluteTime("TDB", origin, 120.);
  chosen = &strict_chooser.choose(eph_cont, pick_time);
  if ("TT" != chosen->getSystem().getName())
    err() << "for time " << pick_time << " with tolerance .99, chooser chose eph with " <<
      chosen->getSystem().getName() << ", not TT as expected" << std::endl;

  strict_chooser = StrictEphChooser(ElapsedTime("TDB", Duration(1.01, "Sec")));
  chosen = &strict_chooser.choose(eph_cont, pick_time);
  if ("TDB" != chosen->getSystem().getName())
    err() << "for time " << pick_time << " with tolerance 1.01, chooser chose eph with " <<
      chosen->getSystem().getName() << ", not TDB as expected" << std::endl;

  // Clean up.
  for (PulsarEphCont::reverse_iterator itor = eph_cont.rbegin(); itor != eph_cont.rend(); ++itor) delete *itor;
  eph_cont.clear();

  // Test sloppy chooser around two disjoint ephemerides.
  epoch.set("TDB", Mjd1(51910.));
  eph_cont.push_back(new FrequencyEph("TT", AbsoluteTime("TDB", Mjd1(51910.)), AbsoluteTime("TDB", Mjd1(51920.)),
    epoch, 22., 45., 0., 1., 0., 0.));
  eph_cont.push_back(new FrequencyEph("TDB", AbsoluteTime("TDB", Mjd1(51930.)), AbsoluteTime("TDB", Mjd1(51940.)),
    epoch, 22., 45., 0., 1., 0., 0.));

  SloppyEphChooser sloppy_chooser;
  pick_time.set("TDB", Mjd1(51905.));
  chosen = &sloppy_chooser.choose(eph_cont, pick_time);
  if ("TT" != chosen->getSystem().getName())
    err() << "for time before either ephemeris: " << pick_time << ", chooser chose eph with " <<
      chosen->getSystem().getName() << ", not TT as expected" << std::endl;

  pick_time.set("TDB", Mjd1(51915.));
  chosen = &sloppy_chooser.choose(eph_cont, pick_time);
  if ("TT" != chosen->getSystem().getName())
    err() << "for time during first ephemeris: " << pick_time << ", chooser chose eph with " <<
      chosen->getSystem().getName() << ", not TT as expected" << std::endl;

  pick_time.set("TDB", Mjd1(51921.));
  chosen = &sloppy_chooser.choose(eph_cont, pick_time);
  if ("TT" != chosen->getSystem().getName())
    err() << "for time shortly after first ephemeris: " << pick_time << ", chooser chose eph with " <<
      chosen->getSystem().getName() << ", not TT as expected" << std::endl;

  pick_time.set("TDB", Mjd1(51929.));
  chosen = &sloppy_chooser.choose(eph_cont, pick_time);
  if ("TDB" != chosen->getSystem().getName())
    err() << "for time shortly before second ephemeris: " << pick_time << ", chooser chose eph with " <<
      chosen->getSystem().getName() << ", not TDB as expected" << std::endl;

  pick_time.set("TDB", Mjd1(51935.));
  chosen = &sloppy_chooser.choose(eph_cont, pick_time);
  if ("TDB" != chosen->getSystem().getName())
    err() << "for time during second ephemeris: " << pick_time << ", chooser chose eph with " <<
      chosen->getSystem().getName() << ", not TDB as expected" << std::endl;

  pick_time.set("TDB", Mjd1(51945.));
  chosen = &sloppy_chooser.choose(eph_cont, pick_time);
  if ("TDB" != chosen->getSystem().getName())
    err() << "for time after second ephemeris: " << pick_time << ", chooser chose eph with " <<
      chosen->getSystem().getName() << ", not TDB as expected" << std::endl;

  // Test choice from prehistory, say 100000 years before origin of MJD.
  pick_time.set("TDB", Mjd1(-36525000.));
  try {
    chosen = &sloppy_chooser.choose(eph_cont, pick_time);
    if ("TT" != chosen->getSystem().getName())
      err() << "for time a long time before the first ephemeris: " << pick_time <<
        ", chooser chose eph with " << chosen->getSystem().getName() << ", not TT as expected" << std::endl;
  } catch (const std::exception & x) {
    err() << "for time a long time before the first ephemeris: " << pick_time <<
      ", chooser threw exception: " << std::endl << x.what() << std::endl;
  }

  // Test choice from far future, say 100000 years after origin of MJD.
  pick_time.set("TDB", Mjd1(36525000.));
  try {
    chosen = &sloppy_chooser.choose(eph_cont, pick_time);
    if ("TDB" != chosen->getSystem().getName())
      err() << "for time a long time after the second ephemeris: " << pick_time <<
        ", chooser chose eph with " << chosen->getSystem().getName() << ", not TDB as expected" << std::endl;
  } catch (const std::exception & x) {
    err() << "for time a long time after the second ephemeris: " << pick_time <<
      ", chooser threw exception: " << std::endl << x.what() << std::endl;
  }

  // Clean up.
  for (PulsarEphCont::reverse_iterator itor = eph_cont.rbegin(); itor != eph_cont.rend(); ++itor) delete *itor;
  eph_cont.clear();
  for (OrbitalEphCont::reverse_iterator itor = orbital_cont.rbegin(); itor != orbital_cont.rend(); ++itor) delete *itor;
  orbital_cont.clear();
}

void PulsarDbTestApp::testEphComputer() {
  setMethod("testEphComputer");

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

  AbsoluteTime expected_gtdb("TDB", 54101, 100.);
  ElapsedTime tolerance("TDB", Duration(1.e-9, "Sec")); // 1 nanosecond.
  const PulsarEph & eph(chooser.choose(eph_cont, expected_gtdb));
  PdotCanceler canceler(eph.getEpoch(), eph, 2);
  canceler.cancelPdot(expected_gtdb);

  // Repeat computations using the EphComputer class, and compare results.
  // Create the computer.
  EphComputer computer;

  // Load the data from the database.
  computer.load(database);

  // Test cancelPdot, and compare result to previous result.
  AbsoluteTime gtdb("TDB", 54101, 100.);
  std::vector<double> fdot_ratio(2, 0.);
  double f0 = eph.calcFrequency(eph.getEpoch(), 0);
  fdot_ratio[0] = eph.calcFrequency(eph.getEpoch(), 1) / f0;
  fdot_ratio[1] = eph.calcFrequency(eph.getEpoch(), 2) / f0;
  computer.setPdotCancelParameter(eph.getSystem().getName(), eph.getEpoch(), fdot_ratio);
  computer.cancelPdot(gtdb);
  if (!expected_gtdb.equivalentTo(gtdb, tolerance))
    err() << "Given time system name, time origin, and fdot ratios, " <<
      "EphComputer::cancelPdot returned absolute time " << gtdb << ", not " << expected_gtdb << ", as expected." << std::endl;

  // Test cancelPdot after setting pdot parameters in a different way, and compare result to previous result.
  gtdb = AbsoluteTime("TDB", 54101, 100.);
  computer.setPdotCancelParameter(eph.getEpoch(), eph, 2);
  computer.cancelPdot(gtdb);
  if (!expected_gtdb.equivalentTo(gtdb, tolerance))
    err() << "Given time origin, pulsar ephemeris, and maximum derivative order, " <<
      "EphComputer::cancelPdot returned absolute time " << gtdb << ", not " << expected_gtdb << ", as expected." << std::endl;

  // Test cancelPdot after setting pdot parameters in yet another way, and compare result to previous result.
  gtdb = AbsoluteTime("TDB", 54101, 100.);
  computer.setPdotCancelParameter(eph.getEpoch(), 2);
  computer.cancelPdot(gtdb);
  if (!expected_gtdb.equivalentTo(gtdb, tolerance))
    err() << "Given time origin and maximum derivative order, " <<
      "EphComputer::cancelPdot returned absolute time " << gtdb << ", not " << expected_gtdb << ", as expected." << std::endl;

  // Test calcPulsePhase, by comparing it with PulsarEph::calcPulsePhase.
  double expected_pulse_phase = eph.calcPulsePhase(expected_gtdb);
  double pulse_phase = computer.calcPulsePhase(expected_gtdb);
  if (expected_pulse_phase != pulse_phase)
    err() << "EphComputer::calcPulsePhase returned phase " << pulse_phase << ", not " <<
      expected_pulse_phase << ", as expected." << std::endl;

  // Test calcPulsePhase, by comparing it with PulsarEph::calcPulsePhase, with a non-zero global phase offset.
  expected_pulse_phase = eph.calcPulsePhase(expected_gtdb, 0.1234);
  pulse_phase = computer.calcPulsePhase(expected_gtdb, 0.1234);
  if (expected_pulse_phase != pulse_phase)
    err() << "EphComputer::calcPulsePhase returned phase " << pulse_phase << ", not " <<
      expected_pulse_phase << ", as expected." << std::endl;

  // Test calcSkyPosition, by comparing it with FrequencyEph::calcSkyPosition.
  std::pair<double, double> expected_ra_dec = eph.calcSkyPosition(expected_gtdb);
  std::pair<double, double> ra_dec = computer.calcSkyPosition(expected_gtdb);
  if (expected_ra_dec != ra_dec)
    err() << "EphComputer::calcSkyPosition returned (RA, Dec) = (" << ra_dec.first << ", " << ra_dec.second <<
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
    err() << "EphComputer::calcOrbitalPhase returned phase " << orbital_phase << ", not " <<
      expected_orbital_phase << ", as expected." << std::endl;

  // Test calcOrbitalPhase, by comparing it with OrbitalEph::calcOrbitalPhase, with a non-zero global phase offset.
  expected_orbital_phase = orbital_eph.calcOrbitalPhase(expected_gtdb, 0.1234);
  orbital_phase = computer.calcOrbitalPhase(expected_gtdb, 0.1234);
  if (expected_orbital_phase != orbital_phase)
    err() << "EphComputer::calcOrbitalPhase returned phase " << orbital_phase << ", not " <<
      expected_orbital_phase << ", as expected." << std::endl;

  // Test binary modulation/demodulation.
  // First perform computations without the computer.
  orbital_eph.modulateBinary(expected_gtdb);

  // Then perform computations with the computer.
  computer.modulateBinary(gtdb);
  if (!expected_gtdb.equivalentTo(gtdb, tolerance))
    err() << "After EphComputer::modulateBinary, absolute time was " << gtdb << ", not " <<
      expected_gtdb << ", as expected." << std::endl;

  // First perform computations without the computer.
  expected_gtdb = gtdb;
  orbital_eph.demodulateBinary(expected_gtdb);
  
  // Then perform computations with the computer.
  computer.demodulateBinary(gtdb);
  if (!expected_gtdb.equivalentTo(gtdb, tolerance))
    err() << "After EphComputer::demodulateBinary, absolute time was " << gtdb << ", not " <<
      expected_gtdb << ", as expected." << std::endl;

  // Prepare for tests of loading methods.
  EphComputer computer2;
  PulsarDb database3(m_tpl_file);
  database3.load(m_in_file);

  // Test emptyness of the pulsar ephemeris container.
  if (0 != computer2.getNumPulsarEph())
    err() << "After creating a new ephemeris computer, there were " << computer.getNumPulsarEph() <<
      " pulsar ephemeri(de)s, not 0 as expected." << std::endl;

  // Test emptyness of the orbital ephemeris container.
  if (0 != computer2.getNumOrbitalEph())
    err() << "After creating a new ephemeris computer, there were " << computer.getNumOrbitalEph() <<
      " orbital ephemeri(de)s, not 0 as expected." << std::endl;

  // Test emptyness of the ephemeris remark container.
  if (0 != computer2.getNumEphRemark())
    err() << "After creating a new ephemeris computer, there were " << computer.getNumEphRemark() <<
      " ephemeris remark(s), not 0 as expected." << std::endl;

  // Register PulsarEph and OrbitalEph subclasses for various ephemeris models.
  database3.registerPulsarEph<FrequencyEph>("FREQ");
  database3.registerOrbitalEph<SimpleDdEph>("DD");

  // Load everything in this database at a time.
  computer2.load(database3);

  if (1141 != computer2.getNumPulsarEph())
    err() << "After loading all data from pulsar database, there were " << computer2.getNumPulsarEph() <<
      " spin pulsar ephemeri(de)s, not 1141 as expected." << std::endl;

  if (19 != computer2.getNumOrbitalEph())
    err() << "After loading all data from pulsar database, there were " << computer2.getNumOrbitalEph() <<
      " orbital ephemeri(de)s, not 19 as expected." << std::endl;

  if (12 != computer2.getNumEphRemark())
    err() << "After loading all data from pulsar database, there were " << computer2.getNumEphRemark() <<
      " ephemeris remark(s), not 12 as expected." << std::endl;

  // Create a new ephemeris computer for tests of type-specific loaders.
  EphComputer computer3;

  // Load just the spin parameters from this database.
  computer3.loadPulsarEph(database3);

  if (1141 != computer3.getNumPulsarEph())
    err() << "After loading spin pulsar ephemerides, there were " << computer3.getNumPulsarEph() <<
      " spin pulsar ephemeri(de)s, not 1141 as expected." << std::endl;

  if (0 != computer3.getNumOrbitalEph())
    err() << "After loading spin pulsar ephemerides, there were " << computer3.getNumOrbitalEph() <<
      " orbital ephemeri(de)s, not 0 as expected." << std::endl;

  if (0 != computer3.getNumEphRemark())
    err() << "After loading spin pulsar ephemerides, there were " << computer3.getNumEphRemark() <<
      " ephemeris remark(s), not 0 as expected." << std::endl;

  // Load just the orbital parameters from this database.
  computer3.loadOrbitalEph(database3);

  if (1141 != computer3.getNumPulsarEph())
    err() << "After loading orbital ephemerides, there were " << computer3.getNumPulsarEph() <<
      " spin pulsar ephemeri(de)s, not 1141 as expected." << std::endl;

  if (19 != computer3.getNumOrbitalEph())
    err() << "After loading orbital ephemerides, there were " << computer3.getNumOrbitalEph() <<
      " orbital ephemeri(de)s, not 19 as expected." << std::endl;

  if (0 != computer3.getNumEphRemark())
    err() << "After loading orbital ephemerides, there were " << computer3.getNumEphRemark() <<
      " ephemeris remark(s), not 0 as expected." << std::endl;

  // Load just the ephemeris remarks from this database.
  computer3.loadEphRemark(database3);

  if (1141 != computer3.getNumPulsarEph())
    err() << "After loading ephemeris remarks, there were " << computer3.getNumPulsarEph() <<
      " spin pulsar ephemeri(de)s, not 1141 as expected." << std::endl;

  if (19 != computer3.getNumOrbitalEph())
    err() << "After loading ephemeris remarks, there were " << computer3.getNumOrbitalEph() <<
      " orbital ephemeri(de)s, not 19 as expected." << std::endl;

  if (12 != computer3.getNumEphRemark())
    err() << "After loading ephemeris remarks, there were " << computer3.getNumEphRemark() <<
      " ephemeris remark(s), not 12 as expected." << std::endl;

  // Prepare variables for loading single ephemeris entry.
  AbsoluteTime since("TDB", 51910, 100.);
  AbsoluteTime until("TDB", 51910, 300.);
  AbsoluteTime epoch("TDB", 51910, 200.);

  // Load just one set of spin parameters.
  computer3.loadPulsarEph(FrequencyEph("TDB", since, until, epoch, 0., 0., 0., 1., 0., 0.));

  if (1142 != computer3.getNumPulsarEph())
    err() << "After loading one spin pulsar ephemeris, there were " << computer3.getNumPulsarEph() <<
      " spin pulsar ephemeri(de)s, not 1142 as expected." << std::endl;

  if (19 != computer3.getNumOrbitalEph())
    err() << "After loading one spin pulsar ephemeris, there were " << computer3.getNumOrbitalEph() <<
      " orbital ephemeri(de)s, not 19 as expected." << std::endl;

  if (12 != computer3.getNumEphRemark())
    err() << "After loading one spin pulsar ephemeris, there were " << computer3.getNumEphRemark() <<
      " ephemeris remark(s), not 12 as expected." << std::endl;

  // Load just one set of orbital parameters.
  computer3.loadOrbitalEph(SimpleDdEph("TDB", 1000., .2, 0., 0., 0., 0., 0., 0., epoch, 0., 0., 0.));

  if (1142 != computer3.getNumPulsarEph())
    err() << "After loading one orbital ephemeris, there were " << computer3.getNumPulsarEph() <<
      " spin pulsar ephemeri(de)s, not 1142 as expected." << std::endl;

  if (20 != computer3.getNumOrbitalEph())
    err() << "After loading one orbital ephemeris, there were " << computer3.getNumOrbitalEph() <<
      " orbital ephemeri(de)s, not 20 as expected." << std::endl;

  if (12 != computer3.getNumEphRemark())
    err() << "After loading one orbital ephemeris, there were " << computer3.getNumEphRemark() <<
      " ephemeris remark(s), not 12 as expected." << std::endl;

  // Load just one ephemeris remark.
  computer3.loadEphRemark(EphStatus(since, until, Remarked, "This is a remark."));

  if (1142 != computer3.getNumPulsarEph())
    err() << "After loading one ephemeris remark, there were " << computer3.getNumPulsarEph() <<
      " spin pulsar ephemeri(de)s, not 1142 as expected." << std::endl;

  if (20 != computer3.getNumOrbitalEph())
    err() << "After loading one ephemeris remark, there were " << computer3.getNumOrbitalEph() <<
      " orbital ephemeri(de)s, not 20 as expected." << std::endl;

  if (13 != computer3.getNumEphRemark())
    err() << "After loading one ephemeris remark, there were " << computer3.getNumEphRemark() <<
      " ephemeris remark(s), not 13 as expected." << std::endl;

  // Load non-remark ephemeris status, which should not be loaded.
  computer3.loadEphRemark(EphStatus(since, until, Unavailable, "No data"));
  computer3.loadEphRemark(EphStatus(since, until, Extrapolated, "Ephemeris gap"));

  if (1142 != computer3.getNumPulsarEph())
    err() << "After loading non-remark ephemeris status, there were " << computer3.getNumPulsarEph() <<
      " spin pulsar ephemeri(de)s, not 1142 as expected." << std::endl;

  if (20 != computer3.getNumOrbitalEph())
    err() << "After loading non-remark ephemeris status, there were " << computer3.getNumOrbitalEph() <<
      " orbital ephemeri(de)s, not 20 as expected." << std::endl;

  if (13 != computer3.getNumEphRemark())
    err() << "After loading non-remark ephemeris status, there were " << computer3.getNumEphRemark() <<
      " ephemeris remark(s), not 13 as expected." << std::endl;

  // Prepare variables for tests of examinePulsarEph method.
  EphComputer computer4;
  EphStatusCont eph_status_cont;
  AbsoluteTime abs_time_01("TDB", 51910, 100.);
  // AbsoluteTime abs_time_02("TDB", 51910, 200.);
  AbsoluteTime abs_time_03("TDB", 51910, 300.);
  AbsoluteTime abs_time_04("TDB", 51910, 400.);
  AbsoluteTime abs_time_05("TDB", 51910, 500.);
  AbsoluteTime abs_time_06("TDB", 51910, 600.);

  // Load two ephemerides with a gap between them.
  computer4.loadPulsarEph(FrequencyEph("TDB", abs_time_01, abs_time_03, epoch, 0., 0., 0., 1., 0., 0.));
  computer4.loadPulsarEph(FrequencyEph("TDB", abs_time_04, abs_time_06, epoch, 0., 0., 0., 1., 0., 0.));

  // Test of delegation of ephemeris gap detection.
  computer4.examinePulsarEph(abs_time_03, abs_time_05, eph_status_cont);
  if (1 != eph_status_cont.size()) {
    err() << "EphComputer::examinePulsarEph method returned " << eph_status_cont.size() <<
      " ephemeris status, not 1, when one ephemeris gap exists" << std::endl;
  } else {
    const EphStatus & eph_status = *(eph_status_cont.begin());
    const EphStatusCodeType & result_code = eph_status.getStatusCode();
    EphStatusCodeType expected_code = Unavailable;
    if (result_code != expected_code) {
      err() << "EphComputer::examinePulsarEph method returned an EphStatus object with status code " << result_code <<
        ", not " << expected_code << ", when one ephemeris gap exists" << std::endl;
    }

    const AbsoluteTime & result_since = eph_status.getEffectiveSince();
    const AbsoluteTime & expected_since = abs_time_03;
    if (!result_since.equivalentTo(expected_since, tolerance)) {
      err() << "EphComputer::examinePulsarEph method returned an EphStatus object effective since " << result_since <<
        ", not " << expected_since << " as expected" << std::endl;
    }

    const AbsoluteTime & result_until = eph_status.getEffectiveUntil();
    const AbsoluteTime & expected_until = abs_time_04;
    if (!result_until.equivalentTo(expected_until, tolerance)) {
      err() << "EphComputer::examinePulsarEph method returned an EphStatus object effective until " << result_until <<
        ", not " << expected_until << " as expected" << std::endl;
    }
  }

  // Prepare for tests of summarizeTimeSystem method.
  EphComputer computer5;

  // Test time system summary with empty computer.
  std::string result_spin;
  std::string result_orbital;
  std::string result_pdot;
  std::string expected_spin("None");
  std::string expected_orbital("None");
  std::string expected_pdot("None");
  computer5.summarizeTimeSystem(result_spin, result_orbital, result_pdot);
  if (result_spin != expected_spin) {
    err () << "EphComputer::summarizeTimeSystem returned \"" << result_spin << "\" for spin ephemerides, not \"" << expected_spin <<
      "\", as expected" << std::endl;
  }
  if (result_orbital != expected_orbital) {
    err () << "EphComputer::summarizeTimeSystem returned \"" << result_orbital << "\" for orbital ephemerides, not \"" <<
      expected_orbital << "\", as expected" << std::endl;
  }
  if (result_pdot != expected_pdot) {
    err () << "EphComputer::summarizeTimeSystem returned \"" << result_pdot << "\" for pdot cancellation, not \"" << expected_pdot <<
      "\", as expected" << std::endl;
  }

  // Load lots of ephemeris data with various time systems.
  computer5.load(database3);
  computer5.loadPulsarEph(FrequencyEph("TT", abs_time_01, abs_time_03, epoch, 0., 0., 0., 1., 0., 0.));
  computer5.loadPulsarEph(PeriodEph("TDB", abs_time_04, abs_time_06, epoch, 0., 0., 0., 1., 0., 0.));
  computer5.loadOrbitalEph(SimpleDdEph("TDB", 1000., .2, 0., 0., 0., 0., 0., 0., epoch, 0., 0., 0.));
  computer5.loadOrbitalEph(SimpleDdEph("TAI", 1000., .2, 0., 0., 0., 0., 0., 0., epoch, 0., 0., 0.));
  computer5.loadOrbitalEph(SimpleDdEph("TT", 1000., .2, 0., 0., 0., 0., 0., 0., epoch, 0., 0., 0.));
  computer5.setPdotCancelParameter("TAI", abs_time_01, fdot_ratio);

  // Test time system summary with filled ephemeris computer.
  expected_spin = "TDB(1142) TT(1)";
  expected_orbital = "TAI(1) TDB(20) TT(1)";
  expected_pdot = "TAI";
  computer5.summarizeTimeSystem(result_spin, result_orbital, result_pdot);
  if (result_spin != expected_spin) {
    err () << "EphComputer::summarizeTimeSystem returned \"" << result_spin << "\" for spin ephemerides, not \"" << expected_spin <<
      "\", as expected" << std::endl;
  }
  if (result_orbital != expected_orbital) {
    err () << "EphComputer::summarizeTimeSystem returned \"" << result_orbital << "\" for orbital ephemerides, not \"" <<
      expected_orbital << "\", as expected" << std::endl;
  }
  if (result_pdot != expected_pdot) {
    err () << "EphComputer::summarizeTimeSystem returned \"" << result_pdot << "\" for pdot cancellation, not \"" << expected_pdot <<
      "\", as expected" << std::endl;
  }

  // Clean up.
  for (PulsarEphCont::reverse_iterator itor = eph_cont.rbegin(); itor != eph_cont.rend(); ++itor) delete *itor;
  eph_cont.clear();
  for (OrbitalEphCont::reverse_iterator itor = orbital_eph_cont.rbegin(); itor != orbital_eph_cont.rend(); ++itor) delete *itor;
  orbital_eph_cont.clear();
}

void PulsarDbTestApp::testEphGetter() {
  setMethod("testEphGetter");

  // Get access to database.
  PulsarDb database(m_tpl_file);
  database.load(m_in_file);

  // Register PulsarEph and OrbitalEph subclasses for various ephemeris models.
  database.registerPulsarEph<FrequencyEph>("FREQ");
  database.registerOrbitalEph<SimpleDdEph>("DD");

  PulsarEphCont pulsar_eph_cont;
  database.getEph(pulsar_eph_cont);
  if (pulsar_eph_cont.size() != size_t(database.getNumEph()))
    err() << "PulsarDb::getEph(PulsarEphCont &) got " << pulsar_eph_cont.size() << " ephemerides, not " <<
      database.getNumEph() << ", as expected." << std::endl;

  OrbitalEphCont orbital_eph_cont;
  database.getEph(orbital_eph_cont);
  size_t expected_orbital = 19;
  if (orbital_eph_cont.size() != expected_orbital) 
    err() << "PulsarDb::getEph(OrbitalEphCont &) got " << orbital_eph_cont.size() << " ephemerides, not " <<
      expected_orbital << ", as expected." << std::endl;

  EphStatusCont eph_status_cont;
  database.getRemark(eph_status_cont);
  size_t expected_remark = 12;
  if (eph_status_cont.size() != expected_remark)
    err() << "PulsarDb::getRemark(EphStatusCont &) got " << eph_status_cont.size() << " remarks, not " <<
      expected_remark << ", as expected." << std::endl;

  // Clean up.
  for (PulsarEphCont::reverse_iterator itor = pulsar_eph_cont.rbegin(); itor != pulsar_eph_cont.rend(); ++itor) delete *itor;
  pulsar_eph_cont.clear();
  for (OrbitalEphCont::reverse_iterator itor = orbital_eph_cont.rbegin(); itor != orbital_eph_cont.rend(); ++itor) delete *itor;
  orbital_eph_cont.clear();
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
      static const AbsoluteTime s_bogus_time("TDB", 0, 0.);
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
    BogusOrbitalEphBase(): OrbitalEph(timeSystem::ElapsedTime("TDB", Duration(10.e-9, "Sec")), 100) {}

    virtual const timeSystem::TimeSystem & getSystem() const { return TimeSystem::getSystem("TDB"); }

    virtual const EphRoutingInfo & getRoutingInfo() const {
      static const EphRoutingInfo s_bogus_info;
      return s_bogus_info;
    }

    virtual const AbsoluteTime & t0() const {
      static const AbsoluteTime s_bogus_absolute_time("TDB", 0, 0.);
      return s_bogus_absolute_time;
    }
    virtual double calcOrbitalPhase(const AbsoluteTime & /* ev_time */, double /* phase_offset */ = 0.) const {
      return 0.;
    }
    virtual ElapsedTime calcOrbitalDelay(const timeSystem::AbsoluteTime & /* ev_time */) const {
      static const ElapsedTime s_bogus_elapsed_time("TDB", Duration::zero());
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

void PulsarDbTestApp::testMultipleEphModel() {
  setMethod("testMultipleEphModel");

  std::auto_ptr<PulsarDb> database(0);

  // Test rejection of a template file w/o EPHSTYLE in SPIN_PARAMETERS extension.
  std::string tpl_file = prependDataPath("test_pulsarDb_badspin.tpl");
  try {
    database.reset(new PulsarDb(tpl_file));
    err() << "PulsarDb::PulsarDb(\"" << tpl_file << "\") did not throw an exception" << std::endl;
  } catch (const std::exception &) {
    // This is fine.
  }

  // Test rejection of a template file w/o EPHSTYLE in ORBITAL_PARAMETERS extension.
  tpl_file = prependDataPath("test_pulsarDb_badorbital.tpl");
  try {
    database.reset(new PulsarDb(tpl_file));
    err() << "PulsarDb::PulsarDb(\"" << tpl_file << "\") did not throw an exception" << std::endl;
  } catch (const std::exception &) {
    // This is fine.
  }

  // Test successful creation of a PulsarDb object with a normal template.
  tpl_file = prependDataPath("test_pulsarDb.tpl");
  try {
    database.reset(new PulsarDb(tpl_file));
  } catch (const std::exception & x) {
    err() << "PulsarDb::PulsarDb(\"" << tpl_file << "\") threw exception: " << std::endl <<
      x.what() << std::endl;
  }

  // Test successful creation of a PulsarDb object with an unusual template.
  // --- First generation SPIN_PARAMETERS tables duplicated.
  tpl_file = prependDataPath("test_pulsarDb_g0g1g2g1.tpl");
  try {
    database.reset(new PulsarDb(tpl_file));
  } catch (const std::exception & x) {
    err() << "PulsarDb::PulsarDb(\"" << tpl_file << "\") threw exception: " << std::endl <<
      x.what() << std::endl;
  }

  // Test successful creation of a PulsarDb object with an unusual template.
  // --- First generation ORBITAL_PARAMETERS tables duplicated.
  tpl_file = prependDataPath("test_pulsarDb_g2g1g0g1.tpl");
  try {
    database.reset(new PulsarDb(tpl_file));
  } catch (const std::exception & x) {
    err() << "PulsarDb::PulsarDb(\"" << tpl_file << "\") threw exception: " << std::endl <<
      x.what() << std::endl;
  }

  // Test successful creation of a PulsarDb object with an unusual template.
  // --- No duplication of tables, but tables are not orderd by generation.
  tpl_file = prependDataPath("test_pulsarDb_g2g1g2g0.tpl");
  try {
    database.reset(new PulsarDb(tpl_file));
  } catch (const std::exception & x) {
    err() << "PulsarDb::PulsarDb(\"" << tpl_file << "\") threw exception: " << std::endl <<
      x.what() << std::endl;
  }

  // Test loading ephemerides from FITS database files in the current format.
  std::string test_subject("current FITS");
  database.reset(new PulsarDb(tpl_file));
  tpl_file = prependDataPath("test_pulsarDb.tpl");
  bool load_original = false;
  bool expected_to_fail = false;
  testLoadingFits(test_subject, *database, tpl_file, load_original, expected_to_fail);

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
  checkEphRouting(test_subject, *database, expected_route_dict);

  // Test loading ephemerides from FITS database files in the original format.
  test_subject = "original FITS";
  database.reset(new PulsarDb(tpl_file));
  load_original = true;
  expected_to_fail = false;
  testLoadingFits(test_subject, *database, tpl_file, load_original, expected_to_fail);

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
  checkEphRouting(test_subject, *database, expected_route_dict);

  // Test loading ephemerides from text database files in the current format.
  test_subject = "current TEXT";
  database.reset(new PulsarDb(tpl_file));
  load_original = false;
  expected_to_fail = false;
  testLoadingText(test_subject, *database, "SPIN_PARAMETERS", "MODEL1", "VALUE1", load_original, expected_to_fail);
  testLoadingText(test_subject, *database, "SPIN_PARAMETERS", "MODEL2", "VALUE2", load_original, expected_to_fail);
  testLoadingText(test_subject, *database, "ORBITAL_PARAMETERS", "MODEL1", "VALUE3", load_original, expected_to_fail);
  testLoadingText(test_subject, *database, "ORBITAL_PARAMETERS", "MODEL2", "VALUE4", load_original, expected_to_fail);

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
  checkEphRouting(test_subject, *database, expected_route_dict);

  // Test loading ephemerides from FITS database files in the original format.
  test_subject = "original TEXT";
  database.reset(new PulsarDb(tpl_file));
  load_original = true;
  expected_to_fail = false;
  testLoadingText(test_subject, *database, "SPIN_PARAMETERS", "MODEL1", "VALUE1", load_original, expected_to_fail);
  testLoadingText(test_subject, *database, "SPIN_PARAMETERS", "MODEL2", "VALUE2", load_original, expected_to_fail);
  testLoadingText(test_subject, *database, "ORBITAL_PARAMETERS", "MODEL1", "VALUE3", load_original, expected_to_fail);
  testLoadingText(test_subject, *database, "ORBITAL_PARAMETERS", "MODEL2", "VALUE4", load_original, expected_to_fail);

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
  checkEphRouting(test_subject, *database, expected_route_dict);

  // Test loading ephemerides from FITS database files in the current format into a badly formatted database.
  test_subject = "current FITS into two SPIN";
  database.reset(new PulsarDb(prependDataPath("test_pulsarDb_g0g1g2g1.tpl")));
  load_original = false;
  expected_to_fail = true;
  testLoadingFits(test_subject, *database, prependDataPath("test_pulsarDb.tpl"), load_original, expected_to_fail);

  // Test loading ephemerides from FITS database files in the current format into a badly formatted database.
  test_subject = "current FITS into two ORBITAL";
  database.reset(new PulsarDb(prependDataPath("test_pulsarDb_g2g1g0g1.tpl")));
  load_original = false;
  expected_to_fail = true;
  testLoadingFits(test_subject, *database, prependDataPath("test_pulsarDb.tpl"), load_original, expected_to_fail);

  // Test loading ephemerides from FITS database files in the current format into an unusually formatted database.
  test_subject = "current FITS into reverse-ordered";
  database.reset(new PulsarDb(prependDataPath("test_pulsarDb_g2g1g2g0.tpl")));
  load_original = false;
  expected_to_fail = false;
  testLoadingFits(test_subject, *database, prependDataPath("test_pulsarDb.tpl"), load_original,
    expected_to_fail);

  // Test getting ephemerides that were loaded from FITS database files in the current format.
  database->registerPulsarEph<BogusPulsarEph<1> >("MODEL1");
  database->registerPulsarEph<BogusPulsarEph<2> >("MODEL2");
  database->registerOrbitalEph<BogusOrbitalEph<1> >("MODEL1");
  database->registerOrbitalEph<BogusOrbitalEph<2> >("MODEL2");
  expected_route_dict.clear();
  expected_route_dict["VALUE1"] = EphRoutingInfo("VALUE1", "SPIN_PARAMETERS", "MODEL1", "BogusPulsarEph1");
  expected_route_dict["VALUE2"] = EphRoutingInfo("VALUE2", "SPIN_PARAMETERS", "MODEL2", "BogusPulsarEph2");
  expected_route_dict["VALUE3"] = EphRoutingInfo("VALUE3", "ORBITAL_PARAMETERS", "MODEL1", "BogusOrbitalEph1");
  expected_route_dict["VALUE4"] = EphRoutingInfo("VALUE4", "ORBITAL_PARAMETERS", "MODEL2", "BogusOrbitalEph2");
  checkEphRouting(test_subject, *database, expected_route_dict);

  // Test loading ephemerides from FITS database files in the original format into a badly formatted database.
  test_subject = "original FITS into two SPIN";
  database.reset(new PulsarDb(prependDataPath("test_pulsarDb_g0g1g2g1.tpl")));
  load_original = true;
  expected_to_fail = true;
  testLoadingFits(test_subject, *database, prependDataPath("test_pulsarDb.tpl"), load_original, expected_to_fail);

  // Test loading ephemerides from FITS database files in the original format into a badly formatted database.
  test_subject = "original FITS into two ORBITAL";
  database.reset(new PulsarDb(prependDataPath("test_pulsarDb_g2g1g0g1.tpl")));
  load_original = true;
  expected_to_fail = true;
  testLoadingFits(test_subject, *database, prependDataPath("test_pulsarDb.tpl"), load_original, expected_to_fail);

  // Test loading ephemerides from FITS database files in the original format into an unusually formatted database.
  test_subject = "original FITS into reverse-ordered";
  database.reset(new PulsarDb(prependDataPath("test_pulsarDb_g2g1g2g0.tpl")));
  load_original = true;
  expected_to_fail = false;
  testLoadingFits(test_subject, *database, prependDataPath("test_pulsarDb.tpl"), load_original,
    expected_to_fail);

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
  checkEphRouting(test_subject, *database, expected_route_dict);

  // Test loading ephemerides from TEXT database files in the current format into a badly formatted database.
  test_subject = "current TEXT into two SPIN";
  database.reset(new PulsarDb(prependDataPath("test_pulsarDb_g0g1g2g1.tpl")));
  load_original = false;
  expected_to_fail = true;
  testLoadingText(test_subject, *database, "SPIN_PARAMETERS", "MODEL1", "VALUE1", load_original, expected_to_fail);
  testLoadingText(test_subject, *database, "SPIN_PARAMETERS", "MODEL2", "VALUE2", load_original, expected_to_fail);
  expected_to_fail = false;
  testLoadingText(test_subject, *database, "ORBITAL_PARAMETERS", "MODEL1", "VALUE3", load_original, expected_to_fail);
  testLoadingText(test_subject, *database, "ORBITAL_PARAMETERS", "MODEL2", "VALUE4", load_original, expected_to_fail);

  // Test loading ephemerides from TEXT database files in the current format into a badly formatted database.
  database.reset(new PulsarDb(prependDataPath("test_pulsarDb_g2g1g0g1.tpl")));
  load_original = false;
  expected_to_fail = false;
  testLoadingText(test_subject, *database, "SPIN_PARAMETERS", "MODEL1", "VALUE1", load_original, expected_to_fail);
  testLoadingText(test_subject, *database, "SPIN_PARAMETERS", "MODEL2", "VALUE2", load_original, expected_to_fail);
  expected_to_fail = true;
  testLoadingText(test_subject, *database, "ORBITAL_PARAMETERS", "MODEL1", "VALUE3", load_original, expected_to_fail);
  testLoadingText(test_subject, *database, "ORBITAL_PARAMETERS", "MODEL2", "VALUE4", load_original, expected_to_fail);

  // Test loading ephemerides from TEXT database files in the current format into an unusually formatted database.
  test_subject = "current TEXT into reverse-ordered";
  database.reset(new PulsarDb(prependDataPath("test_pulsarDb_g2g1g2g0.tpl")));
  load_original = false;
  expected_to_fail = false;
  testLoadingText(test_subject, *database, "SPIN_PARAMETERS", "MODEL1", "VALUE1", load_original, expected_to_fail);
  testLoadingText(test_subject, *database, "SPIN_PARAMETERS", "MODEL2", "VALUE2", load_original, expected_to_fail);
  testLoadingText(test_subject, *database, "ORBITAL_PARAMETERS", "MODEL1", "VALUE3", load_original, expected_to_fail);
  testLoadingText(test_subject, *database, "ORBITAL_PARAMETERS", "MODEL2", "VALUE4", load_original, expected_to_fail);

  // Test getting ephemerides that were loaded from TEXT database files in the current format.
  database->registerPulsarEph<BogusPulsarEph<1> >("MODEL1");
  database->registerPulsarEph<BogusPulsarEph<2> >("MODEL2");
  database->registerOrbitalEph<BogusOrbitalEph<1> >("MODEL1");
  database->registerOrbitalEph<BogusOrbitalEph<2> >("MODEL2");
  expected_route_dict.clear();
  expected_route_dict["VALUE1"] = EphRoutingInfo("VALUE1", "SPIN_PARAMETERS", "MODEL1", "BogusPulsarEph1");
  expected_route_dict["VALUE2"] = EphRoutingInfo("VALUE2", "SPIN_PARAMETERS", "MODEL2", "BogusPulsarEph2");
  expected_route_dict["VALUE3"] = EphRoutingInfo("VALUE3", "ORBITAL_PARAMETERS", "MODEL1", "BogusOrbitalEph1");
  expected_route_dict["VALUE4"] = EphRoutingInfo("VALUE4", "ORBITAL_PARAMETERS", "MODEL2", "BogusOrbitalEph2");
  checkEphRouting(test_subject, *database, expected_route_dict);

  // Test loading ephemerides from TEXT database files in the original format into a badly formatted database.
  test_subject = "original TEXT into two SPIN";
  database.reset(new PulsarDb(prependDataPath("test_pulsarDb_g0g1g2g1.tpl")));
  load_original = true;
  expected_to_fail = true;
  testLoadingText(test_subject, *database, "SPIN_PARAMETERS", "MODEL1", "VALUE1", load_original, expected_to_fail);
  testLoadingText(test_subject, *database, "SPIN_PARAMETERS", "MODEL2", "VALUE2", load_original, expected_to_fail);
  expected_to_fail = false;
  testLoadingText(test_subject, *database, "ORBITAL_PARAMETERS", "MODEL1", "VALUE3", load_original, expected_to_fail);
  testLoadingText(test_subject, *database, "ORBITAL_PARAMETERS", "MODEL2", "VALUE4", load_original, expected_to_fail);

  // Test loading ephemerides from TEXT database files in the original format into a badly formatted database.
  test_subject = "original TEXT into two ORBITAL";
  database.reset(new PulsarDb(prependDataPath("test_pulsarDb_g2g1g0g1.tpl")));
  load_original = true;
  expected_to_fail = false;
  testLoadingText(test_subject, *database, "SPIN_PARAMETERS", "MODEL1", "VALUE1", load_original, expected_to_fail);
  testLoadingText(test_subject, *database, "SPIN_PARAMETERS", "MODEL2", "VALUE2", load_original, expected_to_fail);
  expected_to_fail = true;
  testLoadingText(test_subject, *database, "ORBITAL_PARAMETERS", "MODEL1", "VALUE3", load_original, expected_to_fail);
  testLoadingText(test_subject, *database, "ORBITAL_PARAMETERS", "MODEL2", "VALUE4", load_original, expected_to_fail);

  // Test loading ephemerides from TEXT database files in the original format into an unusually formatted database.
  test_subject = "original TEXT into reverse-ordered";
  database.reset(new PulsarDb(prependDataPath("test_pulsarDb_g2g1g2g0.tpl")));
  load_original = true;
  expected_to_fail = false;
  testLoadingText(test_subject, *database, "SPIN_PARAMETERS", "MODEL1", "VALUE1", load_original, expected_to_fail);
  testLoadingText(test_subject, *database, "SPIN_PARAMETERS", "MODEL2", "VALUE2", load_original, expected_to_fail);
  testLoadingText(test_subject, *database, "ORBITAL_PARAMETERS", "MODEL1", "VALUE3", load_original, expected_to_fail);
  testLoadingText(test_subject, *database, "ORBITAL_PARAMETERS", "MODEL2", "VALUE4", load_original, expected_to_fail);

  // Test getting ephemerides that were loaded from TEXT database files in the original format.
  database->registerPulsarEph<BogusPulsarEph<1> >("MODEL1");
  database->registerPulsarEph<BogusPulsarEph<2> >("MODEL2");
  database->registerOrbitalEph<BogusOrbitalEph<1> >("MODEL1");
  database->registerOrbitalEph<BogusOrbitalEph<2> >("MODEL2");
  expected_route_dict.clear();
  expected_route_dict["VALUE1"] = EphRoutingInfo("VALUE1", "SPIN_PARAMETERS", "MODEL1", "BogusPulsarEph1");
  expected_route_dict["VALUE2"] = EphRoutingInfo("VALUE2", "SPIN_PARAMETERS", "MODEL1", "BogusPulsarEph1");
  expected_route_dict["VALUE3"] = EphRoutingInfo("VALUE3", "ORBITAL_PARAMETERS", "MODEL1", "BogusOrbitalEph1");
  expected_route_dict["VALUE4"] = EphRoutingInfo("VALUE4", "ORBITAL_PARAMETERS", "MODEL1", "BogusOrbitalEph1");
  checkEphRouting(test_subject, *database, expected_route_dict);
}

void PulsarDbTestApp::testEphStatus() {
  setMethod("testEphStatus");

  // Prepare variables for tests of getters.
  AbsoluteTime expected_since("TDB", 51910, 0.);
  AbsoluteTime expected_until("TDB", 51911, 0.);
  EphStatusCodeType expected_code(Remarked);
  std::string expected_description("Bogus description for testing purpose");

  // Create an EphStatus object.
  EphStatus eph_status(expected_since, expected_until, expected_code, expected_description);

  // Test getter for "effective since" time.
  ElapsedTime tolerance("TDB", Duration(1.e-9, "Sec")); // 1 nanosecond.
  AbsoluteTime result_since = eph_status.getEffectiveSince();
  if (!result_since.equivalentTo(expected_since, tolerance)) {
    err() << "EphStatus::getEffectiveSince() returned " << result_since << ", not " << expected_since <<
      " as expected." << std::endl;
  }

  // Test getter for "effective until" time.
  AbsoluteTime result_until = eph_status.getEffectiveUntil();
  if (!result_until.equivalentTo(expected_until, tolerance)) {
    err() << "EphStatus::getEffectiveUntil() returned " << result_until << ", not " << expected_until <<
      " as expected." << std::endl;
  }

  // Test getter for status code.
  EphStatusCodeType result_code = eph_status.getStatusCode();
  if (result_code != expected_code) {
    err() << "EphStatus::getStatusCode() returned " << result_code << ", not " << expected_code <<
      " as expected." << std::endl;
  }

  // Test getter for description.
  std::string result_description = eph_status.getDescription();
  if (result_description != expected_description) {
    err() << "EphStatus::getDescription() returned \"" << result_description << "\", not \"" << expected_description <<
      "\" as expected." << std::endl;
  }

  // Prepare variables for tests of effectiveness checker.
  AbsoluteTime abs_time_earlier("TDB", 51909, 86399.);
  AbsoluteTime abs_time_during("TDB", 51910, 43200.);
  AbsoluteTime abs_time_later("TDB", 51911, 1.);
  AbsoluteTime abs_time_latest("TDB", 51911, 1.);

  // Test determination of status effectiveness.
  bool result_effective = eph_status.effectiveBetween(abs_time_earlier, abs_time_during);
  bool expected_effective = true;
  if (result_effective != expected_effective) {
    err() << "EphStatus::effectiveBetween(" << abs_time_earlier << ", " << abs_time_during << ") returned " <<
      result_effective << ", not " << expected_effective << " as expected." << std::endl;
  }
  result_effective = eph_status.effectiveBetween(abs_time_during, abs_time_later);
  expected_effective = true;
  if (result_effective != expected_effective) {
    err() << "EphStatus::effectiveBetween(" << abs_time_during << ", " << abs_time_later << ") returned " <<
      result_effective << ", not " << expected_effective << " as expected." << std::endl;
  }
  result_effective = eph_status.effectiveBetween(abs_time_later, abs_time_latest);
  expected_effective = false;
  if (result_effective != expected_effective) {
    err() << "EphStatus::effectiveBetween(" << abs_time_later << ", " << abs_time_latest << ") returned " <<
      result_effective << ", not " << expected_effective << " as expected." << std::endl;
  }
  result_effective = eph_status.effectiveBetween(abs_time_earlier, abs_time_later);
  expected_effective = true;
  if (result_effective != expected_effective) {
    err() << "EphStatus::effectiveBetween(" << abs_time_earlier << ", " << abs_time_later << ") returned " <<
      result_effective << ", not " << expected_effective << " as expected." << std::endl;
  }

  // Test handling of never-effective status.
  EphStatus eph_status_never(expected_until, expected_since, expected_code, expected_description);
  result_effective = eph_status_never.effectiveBetween(abs_time_earlier, abs_time_later);
  expected_effective = false;
  if (result_effective != expected_effective) {
    err() << "EphStatus::effectiveBetween(" << abs_time_earlier << ", " << abs_time_later << ") returned " <<
      result_effective << " for never-effective ephemeris status, not " << expected_effective << " as expected." << std::endl;
  }

  // Test EphStatus::report method.
  std::string result_report = EphStatus(expected_since, expected_until, Unavailable, "").report("TDB", MjdFmt);
  std::string expected_report = "No ephemeris available since 51910 MJD (TDB) until 51911 MJD (TDB)";
  if (result_report != expected_report) {
    err() << "EphStatus::report(\"TDB\", MjdFmt) returned \"" << result_report << "\", not \"" <<
      expected_report << "\", expected." << std::endl;
  }
  result_report = EphStatus(expected_since, expected_until, Unavailable, "No data").report("TDB", MjdFmt);
  expected_report = "No ephemeris available (No data) since 51910 MJD (TDB) until 51911 MJD (TDB)";
  if (result_report != expected_report) {
    err() << "EphStatus::report(\"TDB\", MjdFmt) returned \"" << result_report << "\", not \"" <<
      expected_report << "\", expected." << std::endl;
  }
  result_report = EphStatus(expected_since, expected_until, Extrapolated, "").report("TDB", MjdFmt);
  expected_report = "Ephemeris to be extrapolated since 51910 MJD (TDB) until 51911 MJD (TDB)";
  if (result_report != expected_report) {
    err() << "EphStatus::report(\"TDB\", MjdFmt) returned \"" << result_report << "\", not \"" <<
      expected_report << "\", expected." << std::endl;
  }
  result_report = EphStatus(expected_since, expected_until, Extrapolated, "Ephemeris gap").report("TDB", MjdFmt);
  expected_report = "Ephemeris to be extrapolated (Ephemeris gap) since 51910 MJD (TDB) until 51911 MJD (TDB)";
  if (result_report != expected_report) {
    err() << "EphStatus::report(\"TDB\", MjdFmt) returned \"" << result_report << "\", not \"" <<
      expected_report << "\", expected." << std::endl;
  }
  result_report = EphStatus(expected_since, expected_until, Remarked, "Test remark entry").report("TDB", MjdFmt);
  expected_report = "Remarked \"Test remark entry\" since 51910 MJD (TDB) until 51911 MJD (TDB)";
  if (result_report != expected_report) {
    err() << "EphStatus::report(\"TDB\", MjdFmt) returned \"" << result_report << "\", not \"" <<
      expected_report << "\", expected." << std::endl;
  }

  // Prepare variables for tests of ephemeris status computations.
  StrictEphChooser strict_chooser;
  SloppyEphChooser sloppy_chooser;
  PulsarEphCont pulsar_eph_cont;
  AbsoluteTime abs_time_01("TDB", 51910, 53100.);
  AbsoluteTime abs_time_02("TDB", 51910, 53200.);
  AbsoluteTime abs_time_03("TDB", 51910, 53300.);
  AbsoluteTime abs_time_04("TDB", 51910, 53400.);
  AbsoluteTime abs_time_05("TDB", 51910, 53500.);
  AbsoluteTime abs_time_06("TDB", 51910, 53600.);
  AbsoluteTime abs_time_07("TDB", 51910, 53700.);
  AbsoluteTime abs_time_08("TDB", 51910, 53800.);
  AbsoluteTime abs_time_09("TDB", 51910, 53900.);
  AbsoluteTime abs_time_10("TDB", 51910, 54000.);
  EphStatusCont eph_status_cont;

  // Test for detection of ephemeris unavailability by a strict chooser.
  strict_chooser.examine(pulsar_eph_cont, abs_time_04, abs_time_08, eph_status_cont);
  if (1 != eph_status_cont.size()) {
    err() << "StrictEphChooser::examine method returned " << eph_status_cont.size() <<
      " ephemeris status, not 1, when no pulsar ephemeris is passed" << std::endl;
  } else {
    const EphStatusCodeType & result_code = eph_status_cont.begin()->getStatusCode();
    EphStatusCodeType expected_code = Unavailable;
    if (result_code != expected_code) {
      err() << "StrictEphChooser::examine method returned an EphStatus object with status code " << result_code <<
        ", not " << expected_code << ", when no pulsar ephemeris is passed" << std::endl;
    }
  }

  // Test for detection of ephemeris unavailability by a sloppy chooser.
  sloppy_chooser.examine(pulsar_eph_cont, abs_time_04, abs_time_08, eph_status_cont);
  if (1 != eph_status_cont.size()) {
    err() << "SloppyEphChooser::examine method returned " << eph_status_cont.size() <<
      " ephemeris status, not 1, when no pulsar ephemeris is passed" << std::endl;
  } else {
    const EphStatusCodeType & result_code = eph_status_cont.begin()->getStatusCode();
    EphStatusCodeType expected_code = Unavailable;
    if (result_code != expected_code) {
      err() << "SloppyEphChooser::examine method returned an EphStatus object with status code " << result_code <<
        ", not " << expected_code << ", when no pulsar ephemeris is passed" << std::endl;
    }
  }

  // Prepare variables for tests of ephemeris gap detection.
  pulsar_eph_cont.push_back(new FrequencyEph("TDB", abs_time_01, abs_time_02, abs_time_02, 0., 0., 0., 1., 0., 0.));
  pulsar_eph_cont.push_back(new FrequencyEph("TDB", abs_time_01, abs_time_03, abs_time_02, 0., 0., 0., 1., 0., 0.));
  pulsar_eph_cont.push_back(new FrequencyEph("TDB", abs_time_02, abs_time_04, abs_time_02, 0., 0., 0., 1., 0., 0.));
  pulsar_eph_cont.push_back(new FrequencyEph("TDB", abs_time_03, abs_time_06, abs_time_02, 0., 0., 0., 1., 0., 0.));
  pulsar_eph_cont.push_back(new FrequencyEph("TDB", abs_time_07, abs_time_10, abs_time_02, 0., 0., 0., 1., 0., 0.));
  pulsar_eph_cont.push_back(new FrequencyEph("TDB", abs_time_08, abs_time_09, abs_time_02, 0., 0., 0., 1., 0., 0.));

  // Test detection of ephemeris gaps by a strict chooser.
  strict_chooser.examine(pulsar_eph_cont, abs_time_04, abs_time_08, eph_status_cont);
  if (1 != eph_status_cont.size()) {
    err() << "StrictEphChooser::examine method returned " << eph_status_cont.size() <<
      " ephemeris status, not 1, when one pulsar ephemeris gap exists" << std::endl;
  } else {
    const EphStatus & eph_status = *(eph_status_cont.begin());
    const EphStatusCodeType & result_code = eph_status.getStatusCode();
    EphStatusCodeType expected_code = Unavailable;
    if (result_code != expected_code) {
      err() << "StrictEphChooser::examine method returned an EphStatus object with status code " << result_code <<
        ", not " << expected_code << ", when one pulsar ephemeris gap exists" << std::endl;
    }

    const AbsoluteTime & result_since = eph_status.getEffectiveSince();
    const AbsoluteTime & expected_since = abs_time_06;
    if (!result_since.equivalentTo(expected_since, tolerance)) {
      err() << "StrictEphChooser::examine method returned an EphStatus object effective since " << result_since <<
        ", not " << expected_since << " as expected" << std::endl;
    }

    const AbsoluteTime & result_until = eph_status.getEffectiveUntil();
    const AbsoluteTime & expected_until = abs_time_07;
    if (!result_until.equivalentTo(expected_until, tolerance)) {
      err() << "StrictEphChooser::examine method returned an EphStatus object effective until " << result_until <<
        ", not " << expected_until << " as expected" << std::endl;
    }
  }

  // Test detection of ephemeris gaps by a sloppy chooser.
  sloppy_chooser.examine(pulsar_eph_cont, abs_time_04, abs_time_08, eph_status_cont);
  if (1 != eph_status_cont.size()) {
    err() << "SloppyEphChooser::examine method returned " << eph_status_cont.size() <<
      " ephemeris status, not 1, when one pulsar ephemeris gap exists" << std::endl;
  } else {
    const EphStatus & eph_status = *(eph_status_cont.begin());
    const EphStatusCodeType & result_code = eph_status.getStatusCode();
    EphStatusCodeType expected_code = Extrapolated;
    if (result_code != expected_code) {
      err() << "SloppyEphChooser::examine method returned an EphStatus object with status code " << result_code <<
        ", not " << expected_code << ", when one pulsar ephemeris gap exists" << std::endl;
    }

    const AbsoluteTime & result_since = eph_status.getEffectiveSince();
    const AbsoluteTime & expected_since = abs_time_06;
    if (!result_since.equivalentTo(expected_since, tolerance)) {
      err() << "SloppyEphChooser::examine method returned an EphStatus object effective since " << result_since <<
        ", not " << expected_since << " as expected" << std::endl;
    }

    const AbsoluteTime & result_until = eph_status.getEffectiveUntil();
    const AbsoluteTime & expected_until = abs_time_07;
    if (!result_until.equivalentTo(expected_until, tolerance)) {
      err() << "SloppyEphChooser::examine method returned an EphStatus object effective until " << result_until <<
        ", not " << expected_until << " as expected" << std::endl;
    }
  }

  // Clean up.
  for (PulsarEphCont::reverse_iterator itor = pulsar_eph_cont.rbegin(); itor != pulsar_eph_cont.rend(); ++itor) delete *itor;
  pulsar_eph_cont.clear();

  // Prepare variables for tests of ephemeris gap detection.
  pulsar_eph_cont.push_back(new FrequencyEph("TDB", abs_time_07, abs_time_10, abs_time_02, 0., 0., 0., 1., 0., 0.));

  // Test detection of ephemeris gaps at the beginning of the time interval of interest, by a strict chooser.
  strict_chooser.examine(pulsar_eph_cont, abs_time_04, abs_time_08, eph_status_cont);
  if (1 != eph_status_cont.size()) {
    err() << "StrictEphChooser::examine method returned " << eph_status_cont.size() <<
      " ephemeris status, not 1, when no pulsar ephemeris is available at the beginning of the given time interval" << std::endl;
  } else {
    const EphStatus & eph_status = *(eph_status_cont.begin());
    const EphStatusCodeType & result_code = eph_status.getStatusCode();
    EphStatusCodeType expected_code = Unavailable;
    if (result_code != expected_code) {
      err() << "StrictEphChooser::examine method returned an EphStatus object with status code " << result_code <<
        ", not " << expected_code << ", when no pulsar ephemeris is available at the beginning of the given time interval" << std::endl;
    }

    const AbsoluteTime & result_since = eph_status.getEffectiveSince();
    const AbsoluteTime & expected_since = abs_time_04;
    if (!result_since.equivalentTo(expected_since, tolerance)) {
      err() << "StrictEphChooser::examine method returned an EphStatus object effective since " << result_since <<
        ", not " << expected_since << " as expected" << std::endl;
    }

    const AbsoluteTime & result_until = eph_status.getEffectiveUntil();
    const AbsoluteTime & expected_until = abs_time_07;
    if (!result_until.equivalentTo(expected_until, tolerance)) {
      err() << "StrictEphChooser::examine method returned an EphStatus object effective until " << result_until <<
        ", not " << expected_until << " as expected" << std::endl;
    }
  }

  // Test detection of ephemeris gaps at the beginning of the time interval of interest, by a sloppy chooser.
  sloppy_chooser.examine(pulsar_eph_cont, abs_time_04, abs_time_08, eph_status_cont);
  if (1 != eph_status_cont.size()) {
    err() << "SloppyEphChooser::examine method returned " << eph_status_cont.size() <<
      " ephemeris status, not 1, when no pulsar ephemeris is available at the beginning of the given time interval" << std::endl;
  } else {
    const EphStatus & eph_status = *(eph_status_cont.begin());
    const EphStatusCodeType & result_code = eph_status.getStatusCode();
    EphStatusCodeType expected_code = Extrapolated;
    if (result_code != expected_code) {
      err() << "SloppyEphChooser::examine method returned an EphStatus object with status code " << result_code <<
        ", not " << expected_code << ", when no pulsar ephemeris is available at the beginning of the given time interval" << std::endl;
    }

    const AbsoluteTime & result_since = eph_status.getEffectiveSince();
    const AbsoluteTime & expected_since = abs_time_04;
    if (!result_since.equivalentTo(expected_since, tolerance)) {
      err() << "SloppyEphChooser::examine method returned an EphStatus object effective since " << result_since <<
        ", not " << expected_since << " as expected" << std::endl;
    }

    const AbsoluteTime & result_until = eph_status.getEffectiveUntil();
    const AbsoluteTime & expected_until = abs_time_07;
    if (!result_until.equivalentTo(expected_until, tolerance)) {
      err() << "SloppyEphChooser::examine method returned an EphStatus object effective until " << result_until <<
        ", not " << expected_until << " as expected" << std::endl;
    }
  }

  // Clean up.
  for (PulsarEphCont::reverse_iterator itor = pulsar_eph_cont.rbegin(); itor != pulsar_eph_cont.rend(); ++itor) delete *itor;
  pulsar_eph_cont.clear();

  // Prepare variables for tests of ephemeris gap detection.
  pulsar_eph_cont.push_back(new FrequencyEph("TDB", abs_time_03, abs_time_06, abs_time_02, 0., 0., 0., 1., 0., 0.));

  // Test detection of ephemeris gaps at the end of the time interval of interest, by a strict chooser.
  strict_chooser.examine(pulsar_eph_cont, abs_time_04, abs_time_08, eph_status_cont);
  if (1 != eph_status_cont.size()) {
    err() << "StrictEphChooser::examine method returned " << eph_status_cont.size() <<
      " ephemeris status, not 1, when no pulsar ephemeris is available at the end of the given time interval" << std::endl;
  } else {
    const EphStatus & eph_status = *(eph_status_cont.begin());
    const EphStatusCodeType & result_code = eph_status.getStatusCode();
    EphStatusCodeType expected_code = Unavailable;
    if (result_code != expected_code) {
      err() << "StrictEphChooser::examine method returned an EphStatus object with status code " << result_code <<
        ", not " << expected_code << ", when no pulsar ephemeris is available at the end of the given time interval" << std::endl;
    }

    const AbsoluteTime & result_since = eph_status.getEffectiveSince();
    const AbsoluteTime & expected_since = abs_time_06;
    if (!result_since.equivalentTo(expected_since, tolerance)) {
      err() << "StrictEphChooser::examine method returned an EphStatus object effective since " << result_since <<
        ", not " << expected_since << " as expected" << std::endl;
    }

    const AbsoluteTime & result_until = eph_status.getEffectiveUntil();
    const AbsoluteTime & expected_until = abs_time_08;
    if (!result_until.equivalentTo(expected_until, tolerance)) {
      err() << "StrictEphChooser::examine method returned an EphStatus object effective until " << result_until <<
        ", not " << expected_until << " as expected" << std::endl;
    }
  }

  // Test detection of ephemeris gaps at the end of the time interval of interest, by a sloppy chooser.
  sloppy_chooser.examine(pulsar_eph_cont, abs_time_04, abs_time_08, eph_status_cont);
  if (1 != eph_status_cont.size()) {
    err() << "SloppyEphChooser::examine method returned " << eph_status_cont.size() <<
      " ephemeris status, not 1, when no pulsar ephemeris is available at the end of the given time interval" << std::endl;
  } else {
    const EphStatus & eph_status = *(eph_status_cont.begin());
    const EphStatusCodeType & result_code = eph_status.getStatusCode();
    EphStatusCodeType expected_code = Extrapolated;
    if (result_code != expected_code) {
      err() << "SloppyEphChooser::examine method returned an EphStatus object with status code " << result_code <<
        ", not " << expected_code << ", when no pulsar ephemeris is available at the end of the given time interval" << std::endl;
    }

    const AbsoluteTime & result_since = eph_status.getEffectiveSince();
    const AbsoluteTime & expected_since = abs_time_06;
    if (!result_since.equivalentTo(expected_since, tolerance)) {
      err() << "SloppyEphChooser::examine method returned an EphStatus object effective since " << result_since <<
        ", not " << expected_since << " as expected" << std::endl;
    }

    const AbsoluteTime & result_until = eph_status.getEffectiveUntil();
    const AbsoluteTime & expected_until = abs_time_08;
    if (!result_until.equivalentTo(expected_until, tolerance)) {
      err() << "SloppyEphChooser::examine method returned an EphStatus object effective until " << result_until <<
        ", not " << expected_until << " as expected" << std::endl;
    }
  }

  // Clean up.
  for (PulsarEphCont::reverse_iterator itor = pulsar_eph_cont.rbegin(); itor != pulsar_eph_cont.rend(); ++itor) delete *itor;
  pulsar_eph_cont.clear();

  // Prepare variables for tests of ephemeris remark selection.
  EphComputer eph_computer;
  eph_computer.loadEphRemark(EphStatus(abs_time_03, abs_time_05, Remarked, "Remark No. 1"));
  eph_computer.loadEphRemark(EphStatus(abs_time_05, abs_time_07, Remarked, "Remark No. 2"));
  eph_computer.loadEphRemark(EphStatus(abs_time_07, abs_time_09, Remarked, "Remark No. 3"));
  eph_computer.loadEphRemark(EphStatus(abs_time_09, abs_time_10, Remarked, "Remark No. 4"));

  // Test selection of ephemeris remarks by an ephemeris computer.
  eph_computer.getEphRemark(abs_time_04, abs_time_08, eph_status_cont);
  if (3 != eph_status_cont.size()) {
    err() << "EphComputer::getEphRemark method returned " << eph_status_cont.size() <<
      " ephemeris remark(s), not 3 as expected" << std::endl;
  } else {
    EphStatusCont::iterator status_itor = eph_status_cont.begin();
    const EphStatus & eph_status_1 = *status_itor;
    ++status_itor;
    const EphStatus & eph_status_2 = *status_itor;
    ++status_itor;
    const EphStatus & eph_status_3 = *status_itor;

    const EphStatusCodeType & result_code_1 = eph_status_1.getStatusCode();
    const EphStatusCodeType & result_code_2 = eph_status_2.getStatusCode();
    const EphStatusCodeType & result_code_3 = eph_status_3.getStatusCode();
    EphStatusCodeType expected_code = Remarked;
    if (result_code_1 != expected_code || result_code_2 != expected_code || result_code_3 != expected_code) {
      err() << "EphComputer::getEphRemark method returned three EphStatus objects with status codes of " <<
        result_code_1 << ", " << result_code_2 << ", and " << result_code_3 << ", respectively, not all " << expected_code <<
        " as expected" << std::endl;
    }

    const AbsoluteTime & result_since_1 = eph_status_1.getEffectiveSince();
    const AbsoluteTime & result_since_2 = eph_status_2.getEffectiveSince();
    const AbsoluteTime & result_since_3 = eph_status_3.getEffectiveSince();
    const AbsoluteTime & expected_since_1 = abs_time_03;
    const AbsoluteTime & expected_since_2 = abs_time_05;
    const AbsoluteTime & expected_since_3 = abs_time_07;
    if (!result_since_1.equivalentTo(expected_since_1, tolerance) || !result_since_2.equivalentTo(expected_since_2, tolerance) ||
        !result_since_3.equivalentTo(expected_since_3, tolerance)) {
      err() << "EphComputer::getEphRemark method returned three EphStatus objects effective since " <<
        result_since_1 << ", " << result_since_2 << ", and " << result_since_3 << ", respectively, not " << expected_since_1 <<
        ", " << expected_since_2 << ", and " << expected_since_3 << " as expected" << std::endl;
    }

    const AbsoluteTime & result_until_1 = eph_status_1.getEffectiveUntil();
    const AbsoluteTime & result_until_2 = eph_status_2.getEffectiveUntil();
    const AbsoluteTime & result_until_3 = eph_status_3.getEffectiveUntil();
    const AbsoluteTime & expected_until_1 = abs_time_05;
    const AbsoluteTime & expected_until_2 = abs_time_07;
    const AbsoluteTime & expected_until_3 = abs_time_09;
    if (!result_until_1.equivalentTo(expected_until_1, tolerance) || !result_until_2.equivalentTo(expected_until_2, tolerance) ||
        !result_until_3.equivalentTo(expected_until_3, tolerance)) {
      err() << "EphComputer::getEphRemark method returned three EphStatus objects effective until " <<
        result_until_1 << ", " << result_until_2 << ", and " << result_until_3 << ", respectively, not " << expected_until_1 <<
        ", " << expected_until_2 << ", and " << expected_until_3 << " as expected" << std::endl;
    }
  }
}

void PulsarDbTestApp::testBinary() {
  setMethod("testBinary");

  // Create pulsar database object.
  PulsarDb database(m_tpl_file);
  database.load(m_in_file);

  // Test a non-binary pulsar case.
  if (database.isBinary("cRAb"))
    err() << "PulsarDb::isBinary method determined pulsar \"cRAb\" is a binary pulsar." << std::endl;

  // Test a binary pulsar case.
  if (!database.isBinary("PSR J1834-0010"))
    err() << "PulsarDb::isBinary method determined pulsar \"PSR J1834-0010\" is not a binary pulsar." << std::endl;

  // Test no ephemeris case.
  if (database.isBinary("No_such_pulsar"))
    err() << "PulsarDb::isBinary method determined non-existing pulsar \"No_such_pulsar\" is a binary pulsar." << std::endl;

  // Load the wrong spin ephemeris for the Crab.
  database.load(prependDataPath("wrongdb_spin.txt"));
  try {
    database.isBinary("cRAb");
    err() << "PulsarDb::isBinary method did not throw an exception for the Crab pulsar with a wrong binary flag" << std::endl;
  } catch (const std::exception &) {
    // This is fine.
  }

  // Refresh the database contents.
  database.filterName("No_such_pulsar");
  database.load(m_in_file);

  // Load the wrong orbital ephemeris for the Crab.
  database.load(prependDataPath("wrongdb_binary.txt"));
  try {
    database.isBinary("cRAb");
    err() << "PulsarDb::isBinary method did not throw an exception for the Crab pulsar with a wrong orbital ephemeris" << std::endl;
  } catch (const std::exception &) {
    // This is fine.
  }
}

void PulsarDbTestApp::testPulsarDbApp() {
  setMethod("testPulsarDbApp");

  // Create an application tester object.
  PulsarDbAppTester app_tester(*this);

  // Prepare variables to create application objects.
  std::list<std::string> test_name_cont;
  test_name_cont.push_back("par1");
  test_name_cont.push_back("par2");
  test_name_cont.push_back("par3");
  test_name_cont.push_back("par4");
  test_name_cont.push_back("par5");

  // Loop over parameter sets.
  for (std::list<std::string>::const_iterator test_itor = test_name_cont.begin(); test_itor != test_name_cont.end(); ++test_itor) {
    const std::string & test_name = *test_itor;
    std::string out_file(getMethod() + "_" + test_name + ".fits");
    std::string out_file_ref(prependOutrefPath(out_file));

    // Set default parameters.
    st_app::AppParGroup pars(app_tester.getName());
    pars["psrdbfile"] = "";
    pars["outfile"] = "";
    pars["filter"] = "NONE";
    pars["psrname"] = "ANY";
    pars["tstart"] = 0.;
    pars["tstop"] = 1.e5;
    pars["solareph"] = "JPL DE405";
    pars["author"] = "Anonymous User";
    pars["leapsecfile"] = "DEFAULT";
    pars["chatter"] = 2;
    pars["clobber"] = "yes";
    pars["debug"] = "no";
    pars["gui"] = "no";
    pars["mode"] = "ql";

    // Set test-specific parameters.
    if ("par1" == test_name) {
      // Test filtering by official pulsar name.
      pars["psrdbfile"] = m_in_file;
      pars["outfile"] = out_file;
      pars["filter"] = "NAME";
      pars["psrname"] = "PSR J0323+3944";
      pars["author"] = "Anonymous Tester";

    } else if ("par2" == test_name) {
      // Test filtering by pulsar's alternative name.
      pars["psrdbfile"] = m_in_file;
      pars["outfile"] = out_file;
      pars["filter"] = "NAME";
      pars["psrname"] = "Crab";
      pars["author"] = "Anonymous Tester";

    } else if ("par3" == test_name) {
      // Test filtering by time interval.
      pars["psrdbfile"] = m_in_file;
      pars["outfile"] = out_file;
      pars["filter"] = "TIME";
      pars["tstart"] = 53400.;
      pars["tstop"] = 53800.;
      pars["author"] = "Anonymous Tester";

    } else if ("par4" == test_name) {
      // Test filtering by time interval.
      pars["psrdbfile"] = m_in_file;
      pars["outfile"] = out_file;
      pars["filter"] = "SOLAREPH";
      pars["solareph"] = "JPL DE405";
      pars["author"] = "Anonymous Tester";

    } else if ("par5" == test_name) {
      // Test creation of pulsar database file.
      std::string summary_file("psrdb_all.txt");
      remove(summary_file.c_str());
      std::ofstream ofs(summary_file.c_str());
      ofs << prependDataPath("psrdb_spin.txt") << std::endl;
      ofs << prependDataPath("psrdb_binary.txt") << std::endl;
      ofs << prependDataPath("psrdb_remark.txt") << std::endl;
      ofs << prependDataPath("psrdb_obs.txt") << std::endl;
      ofs << prependDataPath("psrdb_name.txt") << std::endl;
      ofs.close();
      pars["psrdbfile"] = "@" + summary_file;
      pars["outfile"] = out_file;
      pars["filter"] = "NONE";
      pars["author"] = "Anonymous Tester";

    } else {
      // Skip this iteration.
      continue;
    }

    // Remove output FITS file.
    remove(out_file.c_str());

    // Test the application.
    app_tester.test(pars, "", "", out_file, out_file_ref);
  }
}

void PulsarDbTestApp::testEphComputerApp() {
  setMethod("testEphComputerApp");

  // Create an application tester object.
  EphComputerAppTester app_tester(*this);

  // Prepare variables to create application objects.
  std::list<std::string> test_name_cont;
  test_name_cont.push_back("par1");
  test_name_cont.push_back("par2");
  test_name_cont.push_back("par3");
  test_name_cont.push_back("par4");
  test_name_cont.push_back("par5");
  test_name_cont.push_back("par6");
  test_name_cont.push_back("par7");
  test_name_cont.push_back("par8");
  test_name_cont.push_back("par9");
  test_name_cont.push_back("par10");
  test_name_cont.push_back("par11");
  test_name_cont.push_back("par12");
  test_name_cont.push_back("par13");
  test_name_cont.push_back("par14");
  test_name_cont.push_back("par15");

  // Prepare files to be used in the tests.
  std::string master_pulsardb = prependDataPath("master_pulsardb_v2.fits");
  std::string leap_file = prependDataPath("gtephem_leapsec.fits");
  std::string remark_file = prependDataPath("gtephem_remark.txt");
  std::string summary_file = getMethod() + "_summary.txt";
  std::ofstream ofs_summary(summary_file.c_str());
  ofs_summary << master_pulsardb << std::endl;
  ofs_summary << remark_file << std::endl;
  ofs_summary << prependDataPath("psrdb_user.fits") << std::endl;
  ofs_summary.close();

  // Loop over parameter sets.
  for (std::list<std::string>::const_iterator test_itor = test_name_cont.begin(); test_itor != test_name_cont.end(); ++test_itor) {
    const std::string & test_name = *test_itor;
    std::string log_file(getMethod() + "_" + test_name + ".log");
    std::string log_file_ref(getMethod() + "_" + test_name + ".ref");
    bool ignore_exception(false);

    // Set default parameters.
    st_app::AppParGroup pars(app_tester.getName());
    pars["psrdbfile"] = "";
    pars["psrname"] = "ANY";
    pars["reftime"] = 0.;
    pars["timeformat"] = "MJD";
    pars["timesys"] = "TDB";
    pars["strict"] = "no";
    pars["solareph"] = "JPL DE405";
    pars["matchsolareph"] = "NONE";
    pars["leapsecfile"] = "DEFAULT";
    pars["reportephstatus"] = "yes";
    pars["chatter"] = 2;
    pars["clobber"] = "yes";
    pars["debug"] = "no";
    pars["gui"] = "no";
    pars["mode"] = "ql";

    // Set test-specific parameters.
    if ("par1" == test_name) {
      // Test standard computation.
      pars["psrdbfile"] = master_pulsardb;
      pars["psrname"] = "PSR B0540-69";
      pars["reftime"] = 212380785.922;
      pars["timeformat"] = "FERMI";
      pars["timesys"] = "TDB";
      log_file_ref = prependOutrefPath(log_file);

    } else if ("par2" == test_name) {
      // Test strict=yes.
      pars["psrdbfile"] = master_pulsardb;
      pars["psrname"] = "PSR B0540-69";
      pars["reftime"] = 212502400.0;
      pars["timeformat"] = "FERMI";
      pars["timesys"] = "TDB";
      pars["strict"] = "yes";
      log_file_ref = prependOutrefPath(log_file);

    } else if ("par3" == test_name) {
      // Test chatter=3.
      pars["psrname"] = "PSR B0540-69";
      pars["reftime"] = 212502400.0;
      pars["psrdbfile"] = master_pulsardb;
      pars["chatter"] = 3;
      pars["timeformat"] = "FERMI";
      pars["timesys"] = "TDB";
      log_file_ref = prependOutrefPath(log_file);

    } else if ("par4" == test_name) {
      // Test an orbital ephemeris.
      pars["psrname"] = "PSR J1834-0010";
      pars["reftime"] = 212502400.0;
      pars["psrdbfile"] = master_pulsardb;
      pars["timeformat"] = "FERMI";
      pars["timesys"] = "TDB";
      log_file_ref = prependOutrefPath(log_file);

    } else if ("par5" == test_name) {
      // Test chatter=3 with an orbital ephemeris.
      pars["psrname"] = "PSR J1834-0010";
      pars["reftime"] = 212502400.0;
      pars["psrdbfile"] = master_pulsardb;
      pars["chatter"] = 3;
      pars["timeformat"] = "FERMI";
      pars["timesys"] = "TDB";
      log_file_ref = prependOutrefPath(log_file);

    } else if ("par6" == test_name) {
      // Test a reference time during a leap second.
      pars["psrname"] = "PSR B0540-69";
      pars["reftime"] = 157766336.0;
      pars["psrdbfile"] = master_pulsardb;
      pars["timeformat"] = "FERMI";
      pars["timesys"] = "UTC";
      pars["leapsecfile"] = leap_file;
      log_file_ref = prependOutrefPath(log_file);

    } else if ("par7" == test_name) {
      // Create a reference output for testing ISO 8601 format.
      pars["psrname"] = "PSR B0540-69";
      pars["reftime"] = "54368.11169104166666666666"; // Need to be a string to preserve precision.
      pars["psrdbfile"] = master_pulsardb;
      pars["strict"] = "yes";
      pars["timeformat"] = "MJD";
      pars["timesys"] = "TDB";
      log_file_ref = prependOutrefPath(log_file);

    } else if ("par8" == test_name) {
      // Test ISO 8601 calendar format.
      pars["psrname"] = "PSR B0540-69";
      pars["reftime"] = "2007-09-25T02:40:50.106";
      pars["psrdbfile"] = master_pulsardb;
      pars["strict"] = "yes";
      pars["timeformat"] = "ISO";
      pars["timesys"] = "TDB";
      log_file_ref = prependOutrefPath(log_file);

    } else if ("par9" == test_name) {
      // Test ISO 8601 week-day format.
      pars["psrname"] = "PSR B0540-69";
      pars["reftime"] = "2007-W39-2T02:40:50.106";
      pars["psrdbfile"] = master_pulsardb;
      pars["strict"] = "yes";
      pars["timeformat"] = "ISO";
      pars["timesys"] = "TDB";
      log_file_ref = prependOutrefPath(log_file);

    } else if ("par10" == test_name) {
      // Test ISO 8601 ordinal-day format.
      pars["psrname"] = "PSR B0540-69";
      pars["reftime"] = "2007-268T02:40:50.106";
      pars["psrdbfile"] = master_pulsardb;
      pars["strict"] = "yes";
      pars["timeformat"] = "ISO";
      pars["timesys"] = "TDB";
      log_file_ref = prependOutrefPath(log_file);

    } else if ("par11" == test_name) {
      // Test detection of empty ephemeris database.
      pars["psrname"] = "PSR B0540-69";
      pars["reftime"] = 212380785.922;
      pars["psrdbfile"] = "NONE";
      pars["timeformat"] = "FERMI";
      pars["timesys"] = "TDB";

      remove(log_file_ref.c_str());
      std::ofstream ofs(log_file_ref.c_str());
      std::runtime_error error("No spin ephemeris is in the database");
      app_tester.writeException(ofs, error);
      ofs.close();

      ignore_exception = true;

    } else if ("par12" == test_name) {
      // Test detection of unknown pulsar name.
      pars["psrname"] = "No Such Pulsar";
      pars["reftime"] = 212380785.922;
      pars["psrdbfile"] = master_pulsardb;
      pars["timeformat"] = "FERMI";
      pars["timesys"] = "TDB";

      remove(log_file_ref.c_str());
      std::ofstream ofs(log_file_ref.c_str());
      std::runtime_error error("No spin ephemeris is available for pulsar \"No Such Pulsar\" in the database");
      app_tester.writeException(ofs, error);
      ofs.close();

      ignore_exception = true;

    } else if ("par13" == test_name) {
      // Test reporting no spin ephemeris available for a given solar system ephemeris.
      pars["psrname"] = "PSR B0540-69";
      pars["reftime"] = 212380785.922;
      pars["psrdbfile"] = master_pulsardb;
      pars["timeformat"] = "FERMI";
      pars["timesys"] = "TDB";
      pars["solareph"] = "JPL DE405";
      pars["matchsolareph"] = "ALL";

      remove(log_file_ref.c_str());
      std::ofstream ofs(log_file_ref.c_str());
      std::runtime_error error("No spin ephemeris is available for solar system ephemeris \"JPL DE405\" for pulsar \"PSR B0540-69\" in the database");
      app_tester.writeException(ofs, error);
      ofs.close();

      ignore_exception = true;

    } else if ("par14" == test_name) {
      // Test reporting no orbital ephemeris available for a given solar system ephemeris.
      // Note: This test requires a spin ephemeris with DE200 and an orbital ephemeris with DE405.
      pars["psrname"] = "PSR J0700+6418";
      pars["reftime"] = 212502400.0;
      pars["psrdbfile"] = m_in_file;
      pars["timeformat"] = "FERMI";
      pars["timesys"] = "TDB";
      pars["solareph"] = "JPL DE200";
      pars["matchsolareph"] = "ALL";

      remove(log_file_ref.c_str());
      std::ofstream ofs(log_file_ref.c_str());
      std::runtime_error error("No orbital ephemeris is available for solar system ephemeris \"JPL DE200\" for pulsar \"PSR J0700+6418\" in the database");
      app_tester.writeException(ofs, error);
      ofs.close();

      ignore_exception = true;

    } else if ("par15" == test_name) {
      // Test reporting of ephemeris remarks.
      pars["psrname"] = "PSR B0540-69";
      pars["reftime"] = 212380785.922;
      pars["psrdbfile"] = "@" + summary_file;
      pars["timeformat"] = "FERMI";
      pars["timesys"] = "TDB";
      log_file_ref = prependOutrefPath(log_file);

    } else if ("par16" == test_name) {
      // Test no reporting of ephemeris remarks with reportephstatus=no.
      pars["psrname"] = "PSR B0540-69";
      pars["reftime"] = 212380785.922;
      pars["psrdbfile"] = "@" + summary_file;
      pars["timeformat"] = "FERMI";
      pars["timesys"] = "TDB";
      pars["reportephstatus"] = "no";
      log_file_ref = prependOutrefPath(log_file);

    } else if ("par17" == test_name) {
      // Test reporting of ephemeris database creation history.
      pars["psrname"] = "PSR B0540-69";
      pars["reftime"] = 212380785.922;
      pars["chatter"] = 4;
      pars["psrdbfile"] = "@" + summary_file;
      pars["timeformat"] = "FERMI";
      pars["timesys"] = "TDB";
      pars["reportephstatus"] = "no";
      log_file_ref = prependOutrefPath(log_file);

    } else {
      // Skip this iteration.
      continue;
    }

    // Test the application.
    app_tester.test(pars, log_file, log_file_ref, "", "", ignore_exception);
  }
}

void PulsarDbTestApp::testLoadingFits(const std::string & test_subject, PulsarDb & database, const std::string & tpl_file,
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
    if (load_original) {
      table->getHeader().erase("EPHSTYLE");
      table->getHeader().erase("PDBTGEN");
    }
  }

  // Test loading ephemerides.
  try {
    // Load the text database.
    database.load(filename);
    if (expected_to_fail) {
      err() << "Loading " << test_subject << ": PulsarDb::load method did not throw exception for FITS file \"" << filename <<
        "\"" << std::endl;
    } else {
      // This is fine.
    }
  } catch (const std::exception & x) {
    if (expected_to_fail) {
      // This is fine.
    } else {
      err() << "Loading " << test_subject << ": PulsarDb::load method threw exception for FITS file \"" << filename <<
        "\": " << std::endl << x.what() << std::endl;
    }
  }
}

void PulsarDbTestApp::testLoadingText(const std::string & test_subject, PulsarDb & database, const std::string & ext_name,
      const std::string & eph_style, const std::string & string_value, bool load_original, bool expected_to_fail) {
  std::string filename("testdb.txt");

  // Create a text database file.
  remove(filename.c_str());
  std::ofstream ofs(filename.c_str());
  ofs << ext_name << std::endl;
  if (!load_original) ofs << "EPHSTYLE = " << eph_style << std::endl;
  ofs << "STRING_VALUE" << std::endl;
  ofs << string_value << std::endl;
  ofs.close();

  // Load the text database.
  try {
    database.load(filename);
    if (expected_to_fail) {
      err() << "Loading " << test_subject << ": PulsarDb::load method did not throw exception for text file \"" << filename <<
        "\" with EXTNAME=" << ext_name << ", EPHSTYLE=" << eph_style << ", STRING_VALUE=" << string_value <<
        ": " << std::endl;
    } else {
      // This is fine.
    }
  } catch (const std::exception & x) {
    if (expected_to_fail) {
      // This is fine.
    } else {
      err() << "Loading " << test_subject << ": PulsarDb::load method threw exception for text file \"" << filename <<
        "\" with EXTNAME=" << ext_name << ", EPHSTYLE=" << eph_style << ", STRING_VALUE=" << string_value <<
        ": " << std::endl << x.what() << std::endl;
    }
  }
}

void PulsarDbTestApp::checkEphRouting(const std::string & test_subject, const PulsarDb & database,
  const std::map<std::string, EphRoutingInfo> & expected_route_dict) {
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
      err() << "Checking " << test_subject << ": PulsarDb::getEph(PulsarEphCont &) returned an object of an unregistered class" <<
        std::endl;
      continue;
    } else {
      returned_route_list.push_back(eph->getRoutingInfo());
    }
  }
  for (OrbitalEphCont::const_iterator itor = orbital_eph_cont.begin(); itor != orbital_eph_cont.end(); ++itor) {
    BogusOrbitalEphBase * eph(dynamic_cast<BogusOrbitalEphBase *>(*itor));
    if (eph == 0) {
      err() << "Checking " << test_subject << ": PulsarDb::getEph(OrbitalEphCont &) returned an object of an unregistered class" <<
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
      err() << "Checking " << test_subject << ": Ephemeris with value \"" << string_value <<
        "\" was not returned by PulsarDb::getEph method" << std::endl;
    } else if (num_eph_found > 1) {
      err() << "Checking " << test_subject << ": Ephemeris with value \"" << string_value <<
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
      err() << "Checking " << test_subject << ": PulsarDb::getEph method returned an unexpected ephemeris data: " <<
        string_value << std::endl;

    } else {
      // Check an extension value of this ephemeris entry.
      std::string expected_ext_name = expected_itor->second.getExtensionName();
      if (ext_name != expected_ext_name) {
      err() << "Checking " << test_subject << ": Ephemeris with value \"" << string_value <<
        "\" was loaded into an extension with EXTNAME=" << ext_name << ", not " << expected_ext_name <<
        " as expected" << std::endl;
      }

      // Check EPHSTYLE keyword value of an extension that this ephemeris entry was coming through.
      std::string expected_model_name = expected_itor->second.getModelName();
      if (model_name != expected_model_name) {
        err() << "Checking " << test_subject << ": Ephemeris with value \"" << string_value <<
          "\" was loaded into an extension with EPHSTYLE=" << model_name << ", not " << expected_model_name <<
          " as expected" << std::endl;
      }

      // Check a class name that this ephemeris entry was passed to.
      std::string expected_class_name = expected_itor->second.getClassName();
      if (class_name != expected_class_name) {
        err() << "Checking " << test_subject << ": Ephemeris with value \"" << string_value <<
          "\" was passed to " << class_name << " class, not " << expected_class_name <<
          " as expected" << std::endl;
      }
    }
  }

  // Clean up.
  for (PulsarEphCont::reverse_iterator itor = pulsar_eph_cont.rbegin(); itor != pulsar_eph_cont.rend(); ++itor) delete *itor;
  pulsar_eph_cont.clear();
  for (OrbitalEphCont::reverse_iterator itor = orbital_eph_cont.rbegin(); itor != orbital_eph_cont.rend(); ++itor) delete *itor;
  orbital_eph_cont.clear();
}

void PulsarDbTestApp::checkEphParameter(const std::string & base_name, const FormattedEph & eph, const TextEph & text_eph) {
  // Prepare files for text outputs.
  std::string outfile(base_name + ".out");
  std::string reffile(base_name + ".ref");
  remove(outfile.c_str());
  remove(reffile.c_str());

  // Prepare an output stream to capture a text output from ephemeris objects.
  st_stream::OStream st_os(false);

  // Write out an ephemeris being tested.
  std::ofstream ofs_out(outfile.c_str());
  st_os.connect(ofs_out);
  eph.write(st_os);
  st_os.disconnect(ofs_out);
  ofs_out.close();

  // Write out a reference ephemeris.
  std::ofstream ofs_ref(reffile.c_str());
  st_os.connect(ofs_ref);
  text_eph.write(st_os);
  st_os.disconnect(ofs_ref);
  ofs_ref.close();

  // Compare the text outputs.
  EphComputerAppTester tester(*this);
  tester.checkOutputText(outfile, reffile);
}

st_app::StAppFactory<PulsarDbTestApp> g_factory("test_pulsarDb");
