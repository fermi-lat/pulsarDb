/** \file test_pulsarDb.cxx
    \brief Test code for pulsarDb package.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/

/// $Header$

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

#include "pulsarDb/BtModelEph.h"
#include "pulsarDb/Ell1ModelEph.h"
#include "pulsarDb/EphChooser.h"
#include "pulsarDb/EphComputer.h"
#include "pulsarDb/EphComputerApp.h"
#include "pulsarDb/EphStatus.h"
#include "pulsarDb/FrequencyEph.h"
#include "pulsarDb/HighPrecisionEph.h"
#include "pulsarDb/MssModelEph.h"
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
#include "timeSystem/SourcePosition.h"
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
    "POS_EPOCH_INT" == column_name || "FREQ_EPOCH_INT" == column_name || "TOAGEO_INT" == column_name || "TOABARY_INT" == column_name ||
    "OBSERVER_CODE" == column_name || "BINARY_FLAG" == column_name || "SOLAR_SYSTEM_EPHEMERIS" == column_name ||
    "DESCRIPTION" == column_name || "OBSERVATORY" == column_name || "CONTACT_PERSON" == column_name || "REFERENCE" == column_name ||
    "ALTNAME" == column_name) {
    // Require an exact match as character strings.
    std::string out_value;
    std::string ref_value;
    out_cell.get(out_value);
    ref_cell.get(ref_value);
    verified = (out_value == ref_value);
    if (!verified) error_stream << "Character string \"" << out_value << "\" not identical to \"" << ref_value << "\"";

  } else if ("FREQ_PARAMETERS" == column_name || "WAVE_SINE" == column_name || "WAVE_COSINE" == column_name ||
    "GLITCH_PARAMETERS" == column_name || "GLITCH_DIMENSIONS" == column_name || "3J_COLUMN" == column_name ||
    "3D_COLUMN" == column_name || "PJ_COLUMN" == column_name || "PD_COLUMN" == column_name) {
    // Extract cell values as arrays of floating-point numbers.
    std::vector<double> out_array;
    std::vector<double> ref_array;
    out_cell.get(out_array);
    ref_cell.get(ref_array);
    bool difference_found = false;

    // Require an exact match in the number of elements.
    if (ref_array.size() != out_array.size()) {
      error_stream << "Number of elements is " << out_array.size() << ", not " << ref_array.size();
      difference_found = true;

    } else {
      // Compare as floating-point numbers.
      double rel_tol = std::numeric_limits<double>::epsilon() * 1000.;
      for (std::vector<double>::size_type idx = 0; idx < ref_array.size(); ++idx) {
        // NOTE: This requires out_value be exactly zero when ref_value is zero.
        if (std::fabs(out_array[idx] - ref_array[idx]) > std::fabs(ref_array[idx]) * rel_tol) {
          error_stream << std::endl << "Element " << idx + 1 << " of " << ref_array.size() << " is " << out_array[idx] <<
            ", not " << ref_array[idx];
          difference_found = true;
        }
      }
    }
    verified = !difference_found;

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

    } else if ("EPOCH_FRAC" == column_name || "POS_EPOCH_FRAC" == column_name || "FREQ_EPOCH_FRAC" == column_name ||
      "TOAGEO_FRAC" == column_name || "TOABARY_FRAC" == column_name) {
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

    /// \brief Test FormattedEph class.
    virtual void testFormattedEph();

    /// \brief Test FrequencyEph class (a subclass of PulsarEph class).
    virtual void testFrequencyEph();

    /// \brief Test PeriodEph class (a subclass of PulsarEph class).
    virtual void testPeriodEph();

    /// \brief Test HighPrecisionEph class (a subclass of PulsarEph class).
    virtual void testHighPrecisionEph();

    /// \brief Test SimpleDdEph class (a subclass of OrbitalEph class).
    virtual void testSimpleDdEph();

    /// \brief Test BtModelEph class (a subclass of OrbitalEph class).
    virtual void testBtModelEph();

    /// \brief Test Ell1ModelEph class (a subclass of OrbitalEph class).
    virtual void testEll1ModelEph();

    /// \brief Test MssModelEph class (a subclass of OrbitalEph class).
    virtual void testMssModelEph();

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
  m_in_file = prependDataPath("testpsrdb_dbmanip.fits");

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
  testFormattedEph();
  testFrequencyEph();
  testPeriodEph();
  testHighPrecisionEph();
  testSimpleDdEph();
  testBtModelEph();
  testEll1ModelEph();
  testMssModelEph();
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
  if (1155 != num_eph)
    err() << "there are initially " << num_eph << " ephemerides, not 1155" << std::endl;

  // Perform filtering which doesn't exclude any ephemerides.
  database.filterInterval(0., 1.e6);

  // Make sure the test data still has expected size.
  num_eph = database.getNumEph();
  if (1155 != num_eph)
    err() << "After no-op filterInterval there are " << num_eph << " ephemerides, not 1155" << std::endl;

  // Another no-op is to use "any".
  database.filterName("aNy");

  // Make sure the test data still has expected size.
  num_eph = database.getNumEph();
  if (1155 != num_eph)
    err() << "After no-op filterName there are " << num_eph << " ephemerides, not 1155" << std::endl;

  // Another no-op is to use a silly expression
  database.filterExpression("#row<1142");

  // Make sure the test data still has expected size.
  num_eph = database.getNumEph();
  if (1155 != num_eph)
    err() << "After no-op filter there are " << num_eph << " ephemerides, not 1155" << std::endl;

  // Another no-op is to use a blank string.
  database.filterExpression("\t  \t");

  // Make sure the test data still has expected size.
  num_eph = database.getNumEph();
  if (1155 != num_eph)
    err() << "After no-op filter there are " << num_eph << " ephemerides, not 1155" << std::endl;

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
  if (4 != num_eph)
    err() << "there are " << num_eph << " ephemerides for PSR J0323+3944, not 4" << std::endl;

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
  if (38 != num_eph)
    err() << "there are " << num_eph << " ephemerides for the crab, not 38" << std::endl;

  // Filter on the crab's PSR B name (no-op).
  database.filterName("PsR b0531+21");
  num_eph = database.getNumEph();
  if (38 != num_eph)
    err() << "there are " << num_eph << " ephemerides for the crab, not 38" << std::endl;

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
  if (408 != num_eph)
    err() << "after filterInterval(53400., 53800.) there are " << num_eph << " ephemerides, not 408" << std::endl;

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
  if (7 != num_eph)
    err() << "after filterSolarEph(\"JPL DE405\") there are " << num_eph << " spin ephemerides, not 7" << std::endl;

  // Make sure the test data has expected size (ORBITAL_PARAMETERS extension).
  num_eph = database.getNumEph(false);
  if (7 != num_eph)
    err() << "after filterSolarEph(\"JPL DE405\") there are " << num_eph << " orbital ephemerides, not 7" << std::endl;

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
  database.filterExpression("OBSERVER_CODE == \"J\"");

  // Confirm that the correct number of ephemerides were found.
  int num_eph = database.getNumEph();
  if (235 != num_eph)
    err() << "found " << num_eph << " ephemerides with OBSERVER_CODE == \"J\", not 235" << std::endl;

  // Test saving this for basis of comparing future test output.
  std::string outfile("obscodeJ_db.fits");
  remove(outfile.c_str());
  database.save(outfile, m_creator, m_author);

  // Check the result against its reference file in data/outref/ directory.
  tester.checkOutputFits(outfile, prependOutrefPath(outfile));
}

void PulsarDbTestApp::testAppend() {
  setMethod("testAppend");
  PulsarDbAppTester tester(*this);

  std::auto_ptr<PulsarDb> database(0);

  // Open an existing FITS database, created by filtering by name.
  database.reset(new PulsarDb(m_tpl_file));
  database->load(prependDataPath("testpsrdb_crab.fits"));

  // Append all tables in another database, originating from text files.
  database->load(prependDataPath("testpsrdb_text.fits"));

  // Save all tables into a new FITS file.
  std::string outfile1 = "psrdb_append.fits";
  remove(outfile1.c_str());
  database->save(outfile1, m_creator, m_author);

  // Check the result against its reference file in data/outref/ directory.
  tester.checkOutputFits(outfile1, prependOutrefPath(outfile1));

  // Test the getter of history records.
  std::list<std::string> expected_command;
  expected_command.push_back("Load FITSDB AUTHOR='Anonymous Tester' DATE=");
  expected_command.push_back("Load FITSDB AUTHOR='Anonymous Tester' DATE=");
  std::list<std::string> expected_ancestry;
  expected_ancestry.push_back("PULSARDB AUTHOR='Anonymous Tester' DATE=");
  expected_ancestry.push_back("* Load FITSDB AUTHOR='Masaharu Hirayama' DATE=");
  expected_ancestry.push_back("* Filter by pulsar name 'CRab'");
  expected_ancestry.push_back("* Filter by pulsar name 'PsR b0531+21'");
  expected_ancestry.push_back("PULSARDB AUTHOR='Masaharu Hirayama' DATE=");
  expected_ancestry.push_back("* Load FITSDB AUTHOR='Masaharu Hirayama' DATE=");
  expected_ancestry.push_back("* Load TEXTDB ORBITAL_PARAMETERS(ELL1) FILENAME='testpsrdb_binary_ell1.t");
  expected_ancestry.push_back("xt'");
  expected_ancestry.push_back("* Load TEXTDB ORBITAL_PARAMETERS(MSS) FILENAME='testpsrdb_binary_mss.txt");
  expected_ancestry.push_back("'");
  expected_ancestry.push_back("PULSARDB AUTHOR='Masaharu Hirayama' DATE=");
  expected_ancestry.push_back("* Load FITSDB AUTHOR='' DATE=");
  expected_ancestry.push_back("* Load TEXTDB SPIN_PARAMETERS(PER) FILENAME='testpsrdb_spin_per.txt'");
  expected_ancestry.push_back("* Load TEXTDB SPIN_PARAMETERS(HP) FILENAME='testpsrdb_spin_hp.txt'");
  expected_ancestry.push_back("* Load TEXTDB ORBITAL_PARAMETERS(BT) FILENAME='testpsrdb_binary_bt.txt'");
  expected_ancestry.push_back("PULSARDB AUTHOR='Anonymous Tester' DATE=");
  expected_ancestry.push_back("* Load TEXTDB SPIN_PARAMETERS(FREQ) FILENAME='psrdb_spin.txt'");
  expected_ancestry.push_back("* Load TEXTDB SPIN_PARAMETERS(PER) FILENAME='psrdb_spin_per.txt'");
  expected_ancestry.push_back("* Load TEXTDB SPIN_PARAMETERS(HP) FILENAME='psrdb_spin_hp.txt'");
  expected_ancestry.push_back("* Load TEXTDB ORBITAL_PARAMETERS(DD) FILENAME='psrdb_binary.txt'");
  expected_ancestry.push_back("* Load TEXTDB ORBITAL_PARAMETERS(BT) FILENAME='psrdb_binary_bt.txt'");
  expected_ancestry.push_back("* Load TEXTDB REMARKS FILENAME='psrdb_remark.txt'");
  expected_ancestry.push_back("* Load TEXTDB OBSERVERS FILENAME='psrdb_obs.txt'");
  expected_ancestry.push_back("* Load TEXTDB ALTERNATIVE_NAMES FILENAME='psrdb_name.txt'");
  std::list<std::string> result_command;
  std::list<std::string> result_ancestry;
  database->getHistory(result_command, result_ancestry);

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

  // Test loading the same data base twice.
  database.reset(new PulsarDb(m_tpl_file));
  database->load(m_in_file);
  database->load(m_in_file);

  // Save the result for basis of comparing future test output.
  std::string outfile2("twice_db.fits");
  remove(outfile2.c_str());
  database->save(outfile2, m_creator, m_author);

  // Check the result against its reference file in data/outref/ directory.
  tester.checkOutputFits(outfile2, prependOutrefPath(outfile2));

  // Test loading databases with TELESCOP=FERMI.
  database.reset(new PulsarDb(m_tpl_file));

  // Test loading a FITS database with TELESCOP=FERMI.
  std::string filename(prependDataPath("testpsrdb_dbmanip_fermi.fits"));
  try {
    database->load(filename);
  } catch (const std::exception & x) {
    err() << "PulsarDb::load method threw an exception for FITS file \"" << filename << "\": " << std::endl << x.what() << std::endl;
  }

  // Test loading a TEXT database with TELESCOP=FERMI.
  filename = prependDataPath("psrdb_spin_fermi.txt");
  try {
    database->load(filename);
  } catch (const std::exception & x) {
    err() << "PulsarDb::load method threw an exception for TEXT file \"" << filename << "\": " << std::endl << x.what() << std::endl;
  }

  // Save the result for basis of comparing future test output.
  std::string outfile3("psrdb_append_fermi.fits");
  remove(outfile3.c_str());
  database->save(outfile3, m_creator, m_author);

  // Check the result against its reference file in data/outref/ directory.
  tester.checkOutputFits(outfile3, prependOutrefPath(outfile3));
}

void PulsarDbTestApp::testTextPulsarDb() {
  setMethod("testTextPulsarDb");
  PulsarDbAppTester tester(*this);

  // Ingest one of each type of table.
  PulsarDb database(m_tpl_file);
  database.load(prependDataPath("psrdb_spin.txt"));
  database.load(prependDataPath("psrdb_spin_per.txt"));
  database.load(prependDataPath("psrdb_spin_hp.txt"));
  database.load(prependDataPath("psrdb_binary.txt"));
  database.load(prependDataPath("psrdb_binary_bt.txt"));
  database.load(prependDataPath("psrdb_remark.txt"));
  database.load(prependDataPath("psrdb_obs.txt"));
  database.load(prependDataPath("psrdb_name.txt"));

  // Save all tables into one FITS file.
  std::string filename1("psrdb_all.fits");
  remove(filename1.c_str());
  database.save(filename1, m_creator, m_author);

  // Check the result against its reference file in data/outref/ directory.
  tester.checkOutputFits(filename1, prependOutrefPath(filename1));

  // Test detection of error(s) in loading badly-formatted tables.
  PulsarDb bad_database(prependDataPath("test_TextPulsarDb.tpl"));
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
  filename_list.push_back("baddb_escapedquote.txt");
  filename_list.push_back("baddb_nestedvector.txt");
  filename_list.push_back("baddb_quotedvector.txt");
  filename_list.push_back("baddb_vectoredscalar.txt");
  filename_list.push_back("baddb_badparenthesis1.txt");
  filename_list.push_back("baddb_badparenthesis2.txt");
  filename_list.push_back("baddb_badparenthesis3.txt");
  filename_list.push_back("baddb_badparenthesis4.txt");
  filename_list.push_back("baddb_badparenthesis5.txt");
  filename_list.push_back("baddb_badparenthesis6.txt");
  filename_list.push_back("baddb_badparenthesis7.txt");
  filename_list.push_back("baddb_badparenthesis8.txt");
  for (std::list<std::string>::const_iterator itor = filename_list.begin(); itor != filename_list.end(); ++itor) {
    const std::string & filename(*itor);
    try {
      bad_database.load(prependDataPath(filename));
      err() << "PulsarDb::load method did not throw an exception for TEXT file \"" << filename << "\"" << std::endl;
    } catch (const std::exception & x) {
      // This is fine.
    }
  }

  // Test no detection of errors in loading poorly-formatted, but legal tables.
  std::list<std::string> basename_list;
  basename_list.push_back("okdb_extraspace");
  basename_list.push_back("okdb_extraquote");
  basename_list.push_back("okdb_extraescape");
  basename_list.push_back("okdb_notallcolumn1");
  basename_list.push_back("okdb_notallcolumn2");
  basename_list.push_back("okdb_emptyfield");
  for (std::list<std::string>::const_iterator itor = basename_list.begin(); itor != basename_list.end(); ++itor) {
    const std::string & basename(*itor);
    PulsarDb ok_database(prependDataPath("test_TextPulsarDb.tpl"));
    std::string filename(basename + ".txt");
    bool successfully_loaded = false;
    try {
      ok_database.load(prependDataPath(filename));
      successfully_loaded = true;
    } catch (const std::exception & x) {
      err() << "PulsarDb::load method threw an exception for TEXT file \"" << filename << "\": " << std::endl << x.what() << std::endl;
    }

    // Write out to a FITS file and compare it with a reference.
    if (successfully_loaded) {
      std::string outfile(basename + ".fits");
      remove(outfile.c_str());
      ok_database.save(outfile, m_creator, m_author);
      tester.checkOutputFits(outfile, prependOutrefPath(outfile));
    }
  }
}

class FormattedEphTester: public FormattedEph {
  public:
    FormattedEphTester(PulsarTestApp & test_app): m_test_app(&test_app) {}
    virtual ~FormattedEphTester() {}

    virtual st_stream::OStream & write(st_stream::OStream & os) const { return os; }

    template <typename DataType>
    void testFormat(const std::string & param_name, const DataType & param_obj, const std::string & param_unit) const {
      st_stream::OStream st_os(false);

      // Test format method without a separator specified.
      std::ostringstream oss1;
      st_os.connect(oss1);
      st_os << format(param_name, param_obj, param_unit);
      st_os.disconnect(oss1);
      std::string result = oss1.str();

      std::ostringstream oss2;
      oss2.width(16); oss2 << param_name << " = " << param_obj;
      if (!param_unit.empty()) oss2 << " " << param_unit;
      std::string expected = oss2.str();

      if (result != expected) {
        m_test_app->err() << "FormattedEph::format(\"" << param_name << "\", " << param_obj << "\", \"" << param_unit <<
          "\") returned a formatter that writes out \"" << result << "\", not \"" << expected << "\"" << std::endl;
      }

      // Test format method with a separator specified.
      std::ostringstream oss3;
      st_os.connect(oss3);
      st_os << format(param_name, param_obj, param_unit, " : ");
      st_os.disconnect(oss3);
      result = oss3.str();

      std::ostringstream oss4;
      oss4.width(16); oss4 << param_name << " : " << param_obj;
      if (!param_unit.empty()) oss4 << " " << param_unit;
      expected = oss4.str();

      if (result != expected) {
        m_test_app->err() << "FormattedEph::format(\"" << param_name << "\", " << param_obj << "\", \"" << param_unit <<
          "\", \":\") returned a formatter that writes out \"" << result << "\", not \"" << expected << "\"" << std::endl;
      }
    }

    template <typename DataType>
    void testReadScalarColumn(const tip::Table::ConstRecord & record, const std::string & field_name, const DataType & default_value,
      const DataType & expected_value, bool expected_to_fail) const {
      DataType data_value;

      // Test read method without a default value.
      bool exception_thrown = false;
      bool return_value = true;
      try {
        return_value = read(record, field_name, data_value);
        if (expected_to_fail) {
          m_test_app->err() << "FormattedEph::read(record, \"" << field_name << "\", data_value) did not throw an exception" <<
            std::endl;
        }
      } catch (const CellReadError & x) {
        exception_thrown = true;
        if (!expected_to_fail) {
          m_test_app->err() << "FormattedEph::read(record, \"" << field_name <<
            "\", data_value) unexpectedly threw an exception: " << std::endl << x.what() << std::endl;
        }
      } catch (const std::exception & x) {
        exception_thrown = true;
        m_test_app->err() << "FormattedEph::read(record, \"" << field_name <<
          "\", data_value) threw a standard exception: " << std::endl << x.what() << std::endl;
      }
      if (!expected_to_fail && !exception_thrown) {
        if (data_value != expected_value) {
          m_test_app->err() << "FormattedEph::read(record, \"" << field_name << "\", data_value) returned with data_value = " <<
            represent(data_value) << ", not " << represent(expected_value) << std::endl;
        }
        if (false != return_value) {
          m_test_app->err() << "FormattedEph::read(record, \"" << field_name <<
            "\", data_value) returned a logical TRUE, not FALSE" << std::endl;
        }
      }

      // Test read method with a default value.
      exception_thrown = false;
      return_value = false;
      try {
        return_value = read(record, field_name, data_value, default_value);
      } catch (const std::exception & x) {
        exception_thrown = true;
        m_test_app->err() << "FormattedEph::read(record, \"" << field_name << "\", data_value, " << represent(default_value) <<
          ") unexpectedly threw an exception: " << std::endl << x.what() << std::endl;
      }
      if (!exception_thrown) {
        if(data_value != expected_value) {
          m_test_app->err() << "FormattedEph::read(record, \"" << field_name << "\", data_value, " << represent(default_value) <<
            ") returned with data_value = " << represent(data_value) << ", not " << represent(expected_value) << std::endl;
        }
        if (return_value != expected_to_fail) {
          m_test_app->err() << "FormattedEph::read(record, \"" << field_name << "\", data_value, " << represent(default_value) <<
            ") returned a logical " << (return_value ? "TRUE" : "FALSE") << ", not " <<
            (expected_to_fail ? "TRUE" : "FALSE") << std::endl;
        }
      }
    }

    template <typename DataType>
    void testReadVectorColumn(const tip::Table::ConstRecord & record, const std::string & field_name,
      const std::vector<DataType> & default_array, const DataType & default_value, const std::vector<DataType> & expected_array,
      const std::vector<bool> & expected_flag, bool expected_to_fail) const {
      std::vector<DataType> data_array;
      std::vector<bool> flag_array;

      // Test read method for a scalar column.
      const std::vector<DataType> & expected_array_scalar = (expected_to_fail ? default_array : expected_array);
      testReadScalarColumn(record, field_name, default_array, expected_array_scalar, expected_to_fail);

      // Test read method with vector arguments.
      bool exception_thrown = false;
      bool return_value = false;
      try {
        return_value = read(record, field_name, data_array, default_value, flag_array);
      } catch (const std::exception & x) {
        exception_thrown = true;
        m_test_app->err() << "FormattedEph::read(record, \"" << field_name << "\", data_array, " << default_value <<
          ", flag_array) unexpectedly threw an exception: " << std::endl << x.what() << std::endl;
      }
      if (!exception_thrown) {
        if (data_array != expected_array) {
          m_test_app->err() << "FormattedEph::read(record, \"" << field_name << "\", data_array, " << default_value <<
            ", flag_array) returned with data_array = " << represent(data_array) << ", not " <<
            represent(expected_array) << std::endl;
        }
        if (flag_array != expected_flag) {
          m_test_app->err() << "FormattedEph::read(record, \"" << field_name << "\", data_array, " << default_value <<
            ", flag_array) returned with flag_array = " << represent(flag_array) << ", not " <<
            represent(expected_flag) << std::endl;
        }
        if (return_value != expected_to_fail) {
          m_test_app->err() << "FormattedEph::read(record, \"" << field_name << "\", data_array, " << default_value <<
            ", flag_array) returned a logical " << (return_value ? "TRUE" : "FALSE") << ", not " <<
            (expected_to_fail ? "TRUE" : "FALSE") << std::endl;
        }
      }
    }

    void testTrimPhaseValue(double phase_value, double phase_offset, double expected_with_offset, double expected_no_offset,
      double tolerance) const {
      // Test trimPhaseValue method with offset argument.
      double phase_result = trimPhaseValue(phase_value, phase_offset);
      if (std::fabs(expected_with_offset - phase_result) > tolerance) {
        m_test_app->err() << "FormattedEph::trimPhaseValue(" << phase_value << ", " << phase_offset << ") returned " <<
          phase_result << ", not " << expected_with_offset << std::endl;
      }

      // Test trimPhaseValue method without offset argument.
      phase_result = trimPhaseValue(phase_value);
      if (std::fabs(expected_no_offset - phase_result) > tolerance) {
        m_test_app->err() << "FormattedEph::trimPhaseValue(" << phase_value << ", " << phase_offset << ") returned " <<
          phase_result << ", not " << expected_no_offset << std::endl;
      }
    }

  private:
    PulsarTestApp * m_test_app;

    template <typename DataType>
    std::string represent(const DataType & data_value) const {
      std::ostringstream oss;
      oss << data_value;
      return oss.str();
    }

    template <typename DataType>
    std::string represent(const std::vector<DataType> & data_array) const {
      std::ostringstream oss;
      oss << "(";
      for (typename std::vector<DataType>::const_iterator itor = data_array.begin(); itor != data_array.end(); ++itor) {
        if (data_array.begin() != itor) oss << ", ";
        oss << *itor;
      }
      oss << ")";
      return oss.str();
    }

    std::string represent(const std::vector<bool> & flag_array) const {
      std::ostringstream oss;
      oss << "(";
      for (std::vector<bool>::const_iterator itor = flag_array.begin(); itor != flag_array.end(); ++itor) {
        if (flag_array.begin() != itor) oss << ", ";
        oss << (*itor ? "TRUE" : "FALSE");
      }
      oss << ")";
      return oss.str();
    }
};

void PulsarDbTestApp::testFormattedEph() {
  setMethod("testFormattedEph");

  // Create a tester object.
  FormattedEphTester tester(*this);

  // Test format method for basic data types.
  float  float_value  = 1.1; tester.testFormat("float",  float_value,  "unit for float");
  double double_value = 2.2; tester.testFormat("double", double_value, "unit for double");
  char   char_value   = 'c'; tester.testFormat("char",   char_value,   "unit for char");
  signed char  signed_char_value  = 's'; tester.testFormat("signed char",  signed_char_value,  "unit for signed char");
  signed short signed_short_value = -11; tester.testFormat("signed short", signed_short_value, "unit for signed short");
  signed int   signed_int_value   = -22; tester.testFormat("signed int",   signed_int_value,   "unit for signed int");
  signed long  signed_long_value  = -33; tester.testFormat("signed long",  signed_long_value,  "unit for signed long");
  unsigned char  unsigned_char_value  = 'u'; tester.testFormat("unsigned char",  unsigned_char_value,  "unit for unsigned char");
  unsigned short unsigned_short_value = 11;  tester.testFormat("unsigned short", unsigned_short_value, "unit for unsigned short");
  unsigned int   unsigned_int_value   = 22;  tester.testFormat("unsigned int",   unsigned_int_value,   "unit for unsigned int");
  unsigned long  unsigned_long_value  = 33;  tester.testFormat("unsigned long",  unsigned_long_value,  "unit for unsigned long");

  // Test format method for compound data types.
  std::string string_value = "string value"; tester.testFormat("string", string_value, "unit for std::string");
  AbsoluteTime abs_time("TDB", 51910, 12345.6789); tester.testFormat("absolute time", abs_time, "unit for unsigned AbsoluteTime");

  // Open a FITS extension to test reading column values.
  std::string filename(prependDataPath("column_samples.fits"));
  std::auto_ptr<const tip::Table> table(tip::IFileSvc::instance().readTable(filename, "SCALAR"));

  // Set test constants.
  bool content_1L = true;
  bool default_1L = false;
  char content_1B = 12;
  char default_1B = 34;
  int content_1I = 123;
  int default_1I = 456;
  long content_1J = 1234;
  long default_1J = 5678;
  float content_1E = 12.345;
  float default_1E = 543.21;
  double content_1D = 1234.56789;
  double default_1D = 98765.4321;

  // Test read method for various data types.
  tip::Table::ConstIterator record_itor = table->begin();
  bool expected_to_fail = false;
  tester.testReadScalarColumn(*record_itor, "1L_COLUMN", default_1L, content_1L, expected_to_fail);
  tester.testReadScalarColumn(*record_itor, "1B_COLUMN", default_1B, content_1B, expected_to_fail);
  tester.testReadScalarColumn(*record_itor, "1I_COLUMN", default_1I, content_1I, expected_to_fail);
  tester.testReadScalarColumn(*record_itor, "1J_COLUMN", default_1J, content_1J, expected_to_fail);
  tester.testReadScalarColumn(*record_itor, "1E_COLUMN", default_1E, content_1E, expected_to_fail);
  tester.testReadScalarColumn(*record_itor, "1D_COLUMN", default_1D, content_1D, expected_to_fail);

  // Test read method for non-existing column.
  expected_to_fail = true;
  tester.testReadScalarColumn(*record_itor, "NO_SUCH_COLUMN", default_1L, default_1L, expected_to_fail);
  tester.testReadScalarColumn(*record_itor, "NO_SUCH_COLUMN", default_1B, default_1B, expected_to_fail);
  tester.testReadScalarColumn(*record_itor, "NO_SUCH_COLUMN", default_1I, default_1I, expected_to_fail);
  tester.testReadScalarColumn(*record_itor, "NO_SUCH_COLUMN", default_1J, default_1J, expected_to_fail);
  tester.testReadScalarColumn(*record_itor, "NO_SUCH_COLUMN", default_1E, default_1E, expected_to_fail);
  tester.testReadScalarColumn(*record_itor, "NO_SUCH_COLUMN", default_1D, default_1D, expected_to_fail);

  // Test read method for NULL detection.
  ++record_itor;
  expected_to_fail = true;
  tester.testReadScalarColumn(*record_itor, "1L_COLUMN", default_1L, default_1L, expected_to_fail);
  tester.testReadScalarColumn(*record_itor, "1B_COLUMN", default_1B, default_1B, expected_to_fail);
  tester.testReadScalarColumn(*record_itor, "1I_COLUMN", default_1I, default_1I, expected_to_fail);
  tester.testReadScalarColumn(*record_itor, "1J_COLUMN", default_1J, default_1J, expected_to_fail);
  tester.testReadScalarColumn(*record_itor, "1E_COLUMN", default_1E, default_1E, expected_to_fail);
  tester.testReadScalarColumn(*record_itor, "1D_COLUMN", default_1D, default_1D, expected_to_fail);

  // Open a FITS extension to test reading vector-column values.
  table.reset(tip::IFileSvc::instance().readTable(filename, "VECTOR"));
  record_itor = table->begin();

  // Set test constants.
  std::vector<bool> null_array(5);
  null_array[0] = false;
  null_array[1] = true;
  null_array[2] = true;
  null_array[3] = false;
  null_array[4] = false;
  std::vector<bool> allok_array(5, false);
  std::vector<bool> empty_array(0);

  std::vector<bool> content_5L(5);
  content_5L[0] = true;
  content_5L[1] = false;
  content_5L[2] = true;
  content_5L[3] = true;
  content_5L[4] = false;
  std::vector<bool> defarray_5L(3);
  defarray_5L[0] = true;
  defarray_5L[1] = true;
  defarray_5L[2] = false;
  bool defvalue_5L = false;
  std::vector<bool> empty_5L(0);

  std::vector<char> content_5B(5);
  content_5B[0] = 12;
  content_5B[1] = 23;
  content_5B[2] = 34;
  content_5B[3] = 45;
  content_5B[4] = 56;
  std::vector<char> defarray_5B(3);
  defarray_5B[0] = 11;
  defarray_5B[1] = 22;
  defarray_5B[2] = 33;
  char defvalue_5B = 67;
  std::vector<char> empty_5B(0);

  std::vector<int> content_5I(5);
  content_5I[0] = 123;
  content_5I[1] = 234;
  content_5I[2] = 345;
  content_5I[3] = 456;
  content_5I[4] = 567;
  std::vector<int> defarray_5I(3);
  defarray_5I[0] = 111;
  defarray_5I[1] = 222;
  defarray_5I[2] = 333;
  int defvalue_5I = 678;
  std::vector<int> empty_5I(0);

  std::vector<long> content_5J(5);
  content_5J[0] = 1234;
  content_5J[1] = 2345;
  content_5J[2] = 3456;
  content_5J[3] = 4567;
  content_5J[4] = 5678;
  std::vector<long> defarray_5J(3);
  defarray_5J[0] = 1111;
  defarray_5J[1] = 2222;
  defarray_5J[2] = 3333;
  long defvalue_5J = 6789;
  std::vector<long> empty_5J(0);

  std::vector<float> content_5E(5);
  content_5E[0] = 12.345;
  content_5E[1] = 23.456;
  content_5E[2] = 34.567;
  content_5E[3] = 45.678;
  content_5E[4] = 56.789;
  std::vector<float> defarray_5E(3);
  defarray_5E[0] = 11.111;
  defarray_5E[1] = 22.222;
  defarray_5E[2] = 33.333;
  float defvalue_5E = 67.891;
  std::vector<float> empty_5E(0);

  std::vector<double> content_5D(5);
  content_5D[0] = 1234.56789;
  content_5D[1] = 2345.67891;
  content_5D[2] = 3456.78912;
  content_5D[3] = 4567.89123;
  content_5D[4] = 5678.91234;
  std::vector<double> defarray_5D(3);
  defarray_5D[0] = 1111.11111;
  defarray_5D[1] = 2222.22222;
  defarray_5D[2] = 3333.33333;
  double defvalue_5D = 6789.12345;
  std::vector<double> empty_5D(0);

  // Test read method for vector columns.
  expected_to_fail = false;
  tester.testReadVectorColumn(*record_itor, "5L_COLUMN", defarray_5L, defvalue_5L, content_5L, allok_array, expected_to_fail);
  tester.testReadVectorColumn(*record_itor, "5B_COLUMN", defarray_5B, defvalue_5B, content_5B, allok_array, expected_to_fail);
  tester.testReadVectorColumn(*record_itor, "5I_COLUMN", defarray_5I, defvalue_5I, content_5I, allok_array, expected_to_fail);
  tester.testReadVectorColumn(*record_itor, "5J_COLUMN", defarray_5J, defvalue_5J, content_5J, allok_array, expected_to_fail);
  tester.testReadVectorColumn(*record_itor, "5E_COLUMN", defarray_5E, defvalue_5E, content_5E, allok_array, expected_to_fail);
  tester.testReadVectorColumn(*record_itor, "5D_COLUMN", defarray_5D, defvalue_5D, content_5D, allok_array, expected_to_fail);

  // Test read method for non-existing vector columns.
  expected_to_fail = true;
  tester.testReadVectorColumn(*record_itor, "NO_SUCH_COLUMN", defarray_5L, defvalue_5L, empty_5L, empty_array, expected_to_fail);
  tester.testReadVectorColumn(*record_itor, "NO_SUCH_COLUMN", defarray_5B, defvalue_5B, empty_5B, empty_array, expected_to_fail);
  tester.testReadVectorColumn(*record_itor, "NO_SUCH_COLUMN", defarray_5I, defvalue_5I, empty_5I, empty_array, expected_to_fail);
  tester.testReadVectorColumn(*record_itor, "NO_SUCH_COLUMN", defarray_5J, defvalue_5J, empty_5J, empty_array, expected_to_fail);
  tester.testReadVectorColumn(*record_itor, "NO_SUCH_COLUMN", defarray_5E, defvalue_5E, empty_5E, empty_array, expected_to_fail);
  tester.testReadVectorColumn(*record_itor, "NO_SUCH_COLUMN", defarray_5D, defvalue_5D, empty_5D, empty_array, expected_to_fail);

  // Test read method for NULL detection for vector columns.
  ++record_itor;
  expected_to_fail = true;

  content_5L[0] = false;
  content_5L[1] = defvalue_5L;
  content_5L[2] = defvalue_5L;
  content_5L[3] = false;
  content_5L[4] = true;

  content_5B[0] = 21;
  content_5B[1] = defvalue_5B;
  content_5B[2] = defvalue_5B;
  content_5B[3] = 54;
  content_5B[4] = 65;

  content_5I[0] = 321;
  content_5I[1] = defvalue_5I;
  content_5I[2] = defvalue_5I;
  content_5I[3] = 654;
  content_5I[4] = 765;

  content_5J[0] = 4321;
  content_5J[1] = defvalue_5J;
  content_5J[2] = defvalue_5J;
  content_5J[3] = 7654;
  content_5J[4] = 8765;

  content_5E[0] = 54.321;
  content_5E[1] = defvalue_5E;
  content_5E[2] = defvalue_5E;
  content_5E[3] = 87.654;
  content_5E[4] = 98.765;

  content_5D[0] = 9876.54321;
  content_5D[1] = defvalue_5D;
  content_5D[2] = defvalue_5D;
  content_5D[3] = 3219.87654;
  content_5D[4] = 4321.98765;

  tester.testReadVectorColumn(*record_itor, "5L_COLUMN", defarray_5L, defvalue_5L, content_5L, null_array, expected_to_fail);
  tester.testReadVectorColumn(*record_itor, "5B_COLUMN", defarray_5B, defvalue_5B, content_5B, null_array, expected_to_fail);
  tester.testReadVectorColumn(*record_itor, "5I_COLUMN", defarray_5I, defvalue_5I, content_5I, null_array, expected_to_fail);
  tester.testReadVectorColumn(*record_itor, "5J_COLUMN", defarray_5J, defvalue_5J, content_5J, null_array, expected_to_fail);
  tester.testReadVectorColumn(*record_itor, "5E_COLUMN", defarray_5E, defvalue_5E, content_5E, null_array, expected_to_fail);
  tester.testReadVectorColumn(*record_itor, "5D_COLUMN", defarray_5D, defvalue_5D, content_5D, null_array, expected_to_fail);

  // Open a FITS extension to test reading variable-length-column values.
  table.reset(tip::IFileSvc::instance().readTable(filename, "VARRAY"));
  record_itor = table->begin();

  // Set test constants.
  null_array.resize(7);
  null_array[0] = false;
  null_array[1] = false;
  null_array[2] = true;
  null_array[3] = true;
  null_array[4] = false;
  null_array[5] = false;
  null_array[6] = false;
  allok_array.resize(3);

  std::vector<bool> content_PL(3);
  content_PL[0] = true;
  content_PL[1] = false;
  content_PL[2] = true;

  std::vector<char> content_PB(3);
  content_PB[0] = 12;
  content_PB[1] = 23;
  content_PB[2] = 34;

  std::vector<int> content_PI(3);
  content_PI[0] = 123;
  content_PI[1] = 234;
  content_PI[2] = 345;

  std::vector<long> content_PJ(3);
  content_PJ[0] = 1234;
  content_PJ[1] = 2345;
  content_PJ[2] = 3456;

  std::vector<float> content_PE(3);
  content_PE[0] = 12.345;
  content_PE[1] = 23.456;
  content_PE[2] = 34.567;

  std::vector<double> content_PD(3);
  content_PD[0] = 1234.56789;
  content_PD[1] = 2345.67891;
  content_PD[2] = 3456.78912;

  // Test read method for variable-length columns.
  expected_to_fail = false;
  tester.testReadVectorColumn(*record_itor, "PL_COLUMN", content_5L, defvalue_5L, content_PL, allok_array, expected_to_fail);
  tester.testReadVectorColumn(*record_itor, "PB_COLUMN", content_5B, defvalue_5B, content_PB, allok_array, expected_to_fail);
  tester.testReadVectorColumn(*record_itor, "PI_COLUMN", content_5I, defvalue_5I, content_PI, allok_array, expected_to_fail);
  tester.testReadVectorColumn(*record_itor, "PJ_COLUMN", content_5J, defvalue_5J, content_PJ, allok_array, expected_to_fail);
  tester.testReadVectorColumn(*record_itor, "PE_COLUMN", content_5E, defvalue_5E, content_PE, allok_array, expected_to_fail);
  tester.testReadVectorColumn(*record_itor, "PD_COLUMN", content_5D, defvalue_5D, content_PD, allok_array, expected_to_fail);

  // Test read method for non-existing variable-length columns.
  expected_to_fail = true;
  tester.testReadVectorColumn(*record_itor, "NO_SUCH_COLUMN", content_5L, defvalue_5L, empty_5L, empty_array, expected_to_fail);
  tester.testReadVectorColumn(*record_itor, "NO_SUCH_COLUMN", content_5B, defvalue_5B, empty_5B, empty_array, expected_to_fail);
  tester.testReadVectorColumn(*record_itor, "NO_SUCH_COLUMN", content_5I, defvalue_5I, empty_5I, empty_array, expected_to_fail);
  tester.testReadVectorColumn(*record_itor, "NO_SUCH_COLUMN", content_5J, defvalue_5J, empty_5J, empty_array, expected_to_fail);
  tester.testReadVectorColumn(*record_itor, "NO_SUCH_COLUMN", content_5E, defvalue_5E, empty_5E, empty_array, expected_to_fail);
  tester.testReadVectorColumn(*record_itor, "NO_SUCH_COLUMN", content_5D, defvalue_5D, empty_5D, empty_array, expected_to_fail);

  // Test read method for NULL detection for vector columns.
  ++record_itor;
  expected_to_fail = true;

  content_PL.resize(7);
  content_PL[0] = false;
  content_PL[1] = true;
  content_PL[2] = defvalue_5L;
  content_PL[3] = defvalue_5L;
  content_PL[4] = false;
  content_PL[5] = false;
  content_PL[6] = true;

  content_PB.resize(7);
  content_PB[0] = 21;
  content_PB[1] = 32;
  content_PB[2] = defvalue_5B;
  content_PB[3] = defvalue_5B;
  content_PB[4] = 43;
  content_PB[5] = 54;
  content_PB[6] = 65;

  content_PI.resize(7);
  content_PI[0] = 321;
  content_PI[1] = 432;
  content_PI[2] = defvalue_5I;
  content_PI[3] = defvalue_5I;
  content_PI[4] = 543;
  content_PI[5] = 654;
  content_PI[6] = 765;

  content_PJ.resize(7);
  content_PJ[0] = 4321;
  content_PJ[1] = 5432;
  content_PJ[2] = defvalue_5J;
  content_PJ[3] = defvalue_5J;
  content_PJ[4] = 6543;
  content_PJ[5] = 7654;
  content_PJ[6] = 8765;

  content_PE.resize(7);
  content_PE[0] = 54.321;
  content_PE[1] = 65.432;
  content_PE[2] = defvalue_5E;
  content_PE[3] = defvalue_5E;
  content_PE[4] = 76.543;
  content_PE[5] = 87.654;
  content_PE[6] = 98.765;

  content_PD.resize(7);
  content_PD[0] = 9876.54321;
  content_PD[1] = 1987.65432;
  content_PD[2] = defvalue_5D;
  content_PD[3] = defvalue_5D;
  content_PD[4] = 2198.76543;
  content_PD[5] = 3219.87654;
  content_PD[6] = 4321.98765;

  tester.testReadVectorColumn(*record_itor, "PL_COLUMN", content_5L, defvalue_5L, content_PL, null_array, expected_to_fail);
  tester.testReadVectorColumn(*record_itor, "PB_COLUMN", content_5B, defvalue_5B, content_PB, null_array, expected_to_fail);
  tester.testReadVectorColumn(*record_itor, "PI_COLUMN", content_5I, defvalue_5I, content_PI, null_array, expected_to_fail);
  tester.testReadVectorColumn(*record_itor, "PJ_COLUMN", content_5J, defvalue_5J, content_PJ, null_array, expected_to_fail);
  tester.testReadVectorColumn(*record_itor, "PE_COLUMN", content_5E, defvalue_5E, content_PE, null_array, expected_to_fail);
  tester.testReadVectorColumn(*record_itor, "PD_COLUMN", content_5D, defvalue_5D, content_PD, null_array, expected_to_fail);

  // Test trimPhaseValue method.
  double tolerance = 1.e-6;
  tester.testTrimPhaseValue(0.00345,         0.12,         0.12345, 0.00345, tolerance);
  tester.testTrimPhaseValue(0.00345 +     1, 0.12,         0.12345, 0.00345, tolerance);
  tester.testTrimPhaseValue(0.00345,         0.12 +     1, 0.12345, 0.00345, tolerance);
  tester.testTrimPhaseValue(0.00345 +     1, 0.12 +     1, 0.12345, 0.00345, tolerance);
  tester.testTrimPhaseValue(0.00345 +    21, 0.12,         0.12345, 0.00345, tolerance);
  tester.testTrimPhaseValue(0.00345,         0.12 +    21, 0.12345, 0.00345, tolerance);
  tester.testTrimPhaseValue(0.00345 +    21, 0.12 +    21, 0.12345, 0.00345, tolerance);
  tester.testTrimPhaseValue(0.00345 +   321, 0.12,         0.12345, 0.00345, tolerance);
  tester.testTrimPhaseValue(0.00345,         0.12 +   321, 0.12345, 0.00345, tolerance);
  tester.testTrimPhaseValue(0.00345 +   321, 0.12 +   321, 0.12345, 0.00345, tolerance);
  tester.testTrimPhaseValue(0.00345 +  4321, 0.12,         0.12345, 0.00345, tolerance);
  tester.testTrimPhaseValue(0.00345,         0.12 +  4321, 0.12345, 0.00345, tolerance);
  tester.testTrimPhaseValue(0.00345 +  4321, 0.12 +  4321, 0.12345, 0.00345, tolerance);
  tester.testTrimPhaseValue(0.00345 + 54321, 0.12,         0.12345, 0.00345, tolerance);
  tester.testTrimPhaseValue(0.00345,         0.12 + 54321, 0.12345, 0.00345, tolerance);
  tester.testTrimPhaseValue(0.00345 + 54321, 0.12 + 54321, 0.12345, 0.00345, tolerance);
  tester.testTrimPhaseValue(0.00345 -     1, 0.12,         0.12345, 0.00345, tolerance);
  tester.testTrimPhaseValue(0.00345,         0.12 -     1, 0.12345, 0.00345, tolerance);
  tester.testTrimPhaseValue(0.00345 -     1, 0.12 -     1, 0.12345, 0.00345, tolerance);
  tester.testTrimPhaseValue(0.00345 -    21, 0.12,         0.12345, 0.00345, tolerance);
  tester.testTrimPhaseValue(0.00345,         0.12 -    21, 0.12345, 0.00345, tolerance);
  tester.testTrimPhaseValue(0.00345 -    21, 0.12 -    21, 0.12345, 0.00345, tolerance);
  tester.testTrimPhaseValue(0.00345 -   321, 0.12,         0.12345, 0.00345, tolerance);
  tester.testTrimPhaseValue(0.00345,         0.12 -   321, 0.12345, 0.00345, tolerance);
  tester.testTrimPhaseValue(0.00345 -   321, 0.12 -   321, 0.12345, 0.00345, tolerance);
  tester.testTrimPhaseValue(0.00345 -  4321, 0.12,         0.12345, 0.00345, tolerance);
  tester.testTrimPhaseValue(0.00345,         0.12 -  4321, 0.12345, 0.00345, tolerance);
  tester.testTrimPhaseValue(0.00345 -  4321, 0.12 -  4321, 0.12345, 0.00345, tolerance);
  tester.testTrimPhaseValue(0.00345 - 54321, 0.12,         0.12345, 0.00345, tolerance);
  tester.testTrimPhaseValue(0.00345,         0.12 - 54321, 0.12345, 0.00345, tolerance);
  tester.testTrimPhaseValue(0.00345 - 54321, 0.12 - 54321, 0.12345, 0.00345, tolerance);
}

class TextEph: public FormattedEph {
  public:
    TextEph() {}
    virtual ~TextEph() {}

    void append(const std::string & par_name, const std::string & par_value, const std::string & par_unit) {
      m_par_list.push_back(std::make_pair(par_name, std::make_pair(par_value, par_unit)));
    }

    template <typename DataType>
    void append(const std::string & par_name, const DataType & par_value, const std::string & par_unit) {
      std::ostringstream oss;
      oss.precision(std::numeric_limits<double>::digits10);
      oss << par_value;
      append(par_name, oss.str(), par_unit);
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
        const std::string & par_value(itor->second.first);
        const std::string & par_unit(itor->second.second);
        if (("Valid Since" == par_name) || ("Valid Until" == par_name)) {
          os << format(par_name, par_value, par_unit, " : ");
        } else {
          os << format(par_name, par_value, par_unit);
        }
      }
      return os;
    }

  private:
    typedef std::list<std::pair<std::string, std::pair<std::string, std::string> > > plist_type;
    plist_type m_par_list;
};

void PulsarDbTestApp::testFrequencyEph() {
  setMethod("testFrequencyEph");

  // Prepare absolute times and tolerance for testing computations.
  AbsoluteTime since("TT", 51910, 0.);
  AbsoluteTime until("TT", 51910, 1.);
  AbsoluteTime epoch("TT", 51910, 123.456789);
  AbsoluteTime pick_time("TT", 51910, 223.456789);
  ElapsedTime tolerance("TDB", Duration(1.e-9, "Sec")); // 1 nanosecond.
  double epsilon = 1.e-8;

  // Create a frequency ephemeris.
  std::auto_ptr<FrequencyEph> eph(new FrequencyEph("TDB", since, until, epoch, 22., 45., 0.11, 1.125e-2, -2.25e-4, 6.75e-6));

  // Test time system getter.
  std::string sys_name = eph->getSystem().getName();
  if ("TDB" != sys_name) {
    err() << "FrequencyEph::getSystem returned \"" << sys_name << "\", not \"TDB\"" << std::endl;
  }

  // Test valid-since getter.
  const AbsoluteTime & since_returned = eph->getValidSince();
  if (!since_returned.equivalentTo(since, tolerance)) {
    err() << "FrequencyEph::getValidSince returned AbsoluteTime(" << since_returned << "), not AbsoluteTime(" << since <<
      ")" << std::endl;
  }

  // Test valid-until getter.
  const AbsoluteTime & until_returned = eph->getValidUntil();
  if (!until_returned.equivalentTo(until, tolerance)) {
    err() << "FrequencyEph::getValidUntil returned AbsoluteTime(" << until_returned << "), not AbsoluteTime(" << until <<
      ")" << std::endl;
  }

  // Test epoch getter.
  const AbsoluteTime & epoch_returned = eph->getEpoch();
  if (!epoch_returned.equivalentTo(epoch, tolerance)) {
    err() << "FrequencyEph::getEpoch returned AbsoluteTime(" << epoch_returned << "), not AbsoluteTime(" << epoch <<
      ")" << std::endl;
  }

  // Test remark getter.
  if (eph->getRemark().size()) {
    err() << "FrequencyEph::getRemark returned a non-empty container of remarks." << std::endl;
  }

  // Test pulse phase computations.
  double phase = eph->calcPulsePhase(pick_time);
  if (std::fabs(phase/.235 - 1.) > epsilon)
    err() << "FrequencyEph::calcPulsePhase produced phase == " << phase << ", not .235" << std::endl;
 
  // Test pulse phase computations, with a non-zero global phase offset.
  phase = eph->calcPulsePhase(pick_time, 0.1234);
  if (std::fabs(phase/.3584 - 1.) > epsilon)
    err() << "FrequencyEph::calcPulsePhase produced phase == " << phase << ", not .3584" << std::endl;
 
  // Change ephemeris to produce a noticeable effect.
  epoch = AbsoluteTime("TT", 51910, 123.4567891234567);
  eph.reset(new FrequencyEph("TDB", since, until, epoch, 22., 45., .11, 1.125e-2, -2.25e-4, 13.5e-6));
  AbsoluteTime ev_time("TT", 51910, 223.4567891234567);

  // Test the pulsar position.
  SourcePosition computed_pos = eph->calcPosition(ev_time);
  SourcePosition correct_pos = SourcePosition(22., 45.);
  std::string coord_name("XYZ");
  for (std::size_t ii = 0; ii < 3; ++ii) {
    double computed_coord = computed_pos.getDirection()[ii];
    double correct_coord = correct_pos.getDirection()[ii];
    if (std::fabs(computed_coord - correct_coord) > epsilon) {
      err() << "FrequencyEph::calcPosition returned " << coord_name[ii] << "=" << computed_coord << ", not " <<
        correct_coord << "." << std::endl;
    }
  }
  if (computed_pos.hasDistance() != false) {
    err() << "FrequencyEph::calcPosition returned a source position with its distance known." << std::endl;
  }

  // Test phase computation.
  double computed_phi0 = eph->calcPulsePhase(ev_time);
  if (std::fabs(computed_phi0/.36 - 1.) > epsilon)
    err() << "FrequencyEph::calcPulsePhase produced phi0 == " << computed_phi0 << ", not .36" << std::endl;
 
  // Test frequency computation.
  double computed_f0 = eph->calcFrequency(ev_time, 0);
  double correct_f0 = 5.625e-2;
  if (std::fabs(computed_f0 - correct_f0) > epsilon) {
    err() << "FrequencyEph::calcFrequency produced f0 == " << computed_f0 << ", not " << correct_f0 << std::endl;
  }
  
  double computed_f1 = eph->calcFrequency(ev_time, 1);
  double correct_f1 = 11.25e-4;
  if (std::fabs(computed_f1 - correct_f1) > epsilon) {
    err() << "FrequencyEph::calcFrequency produced f1 == " << computed_f1 << ", not " << correct_f1 << std::endl;
  }
  
  double computed_f2 = eph->calcFrequency(ev_time, 2);
  double correct_f2 = 13.5e-6;
  if (std::fabs(computed_f2 - correct_f2) > epsilon) {
    err() << "FrequencyEph::calcFrequency produced f2 == " << computed_f2 << ", not " << correct_f2 << std::endl;
  }

  // Test the constructor that takes numerical arguments.
  eph.reset(new FrequencyEph("TDB", AbsoluteTime("TDB", 12345, 51840.), AbsoluteTime("TDB", 23456, 60480.),
    AbsoluteTime("TDB", 34567, 69120.), 11., 22., 33., 44., 55., 66.));
  TextEph t_eph;
  t_eph.append("Valid Since", "12345.6 MJD (TDB)", "");
  t_eph.append("Valid Until", "23456.7 MJD (TDB)", "");
  t_eph.append("Epoch",       "34567.8 MJD (TDB)", "");
  t_eph.append("RA",   "11", "degrees");
  t_eph.append("Dec",  "22", "degrees");
  t_eph.append("Phi0", "33", "");
  t_eph.append("F0",   "44", "s**(-1)");
  t_eph.append("F1",   "55", "s**(-2)");
  t_eph.append("F2",   "66", "s**(-3)");
  checkEphParameter(getMethod() + "_numeric", *eph, t_eph);

  // Test the constructor that takes a FITS record.
  std::string test_tpl("test_FrequencyEph.tpl");
  remove(test_tpl.c_str());
  std::ofstream ofs(test_tpl.c_str());
  ofs << "\\include " << prependDataPath("PulsarDb_primary.tpl") << std::endl;
  ofs << "\\include " << prependDataPath("PulsarDb_spin_freq.tpl") << std::endl;
  ofs.close();

  tip::IFileSvc::instance().createFile(getMethod() + ".fits", test_tpl);
  tip::TipFile tip_file = tip::IFileSvc::instance().openFile(getMethod() + ".fits");

  std::auto_ptr<tip::Table> table(tip_file.editTable("1"));
  tip::Header & header(table->getHeader());
  table->setNumRecords(1);
  tip::Table::Record record(table.get(), 0);
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
  eph.reset(new FrequencyEph(record, header));
  t_eph.clear();
  t_eph.append("Valid Since", "54321 MJD (TDB)",   "");
  t_eph.append("Valid Until", "65433 MJD (TDB)",   ""); // The end of 65432 == the begining of 65433.
  t_eph.append("Epoch",       "76543.5 MJD (TDB)", "");
  t_eph.append("RA",   "1.1", "degrees");
  t_eph.append("Dec",  "2.2", "degrees");
  double dt = .864;
  double phi0 = -(3.3*dt + 4.4/2.*dt*dt + 5.5/6.*dt*dt*dt);
  while (phi0 < 0.) ++phi0;
  while (phi0 >= 1.) --phi0;
  t_eph.append("Phi0", phi0, "");
  t_eph.append("F0",   "3.3", "s**(-1)");
  t_eph.append("F1",   "4.4", "s**(-2)");
  t_eph.append("F2",   "5.5", "s**(-3)");
  checkEphParameter(getMethod() + "_fits", *eph, t_eph);
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

  // Prepare absolute times and tolerance for testing computations.
  AbsoluteTime since("TDB", 51910, 0.);
  AbsoluteTime until("TDB", 51910, 1.);
  AbsoluteTime epoch("TDB", 51910, 123.456789);
  ElapsedTime tolerance("TDB", Duration(1.e-9, "Sec")); // 1 nanosecond.
  double epsilon = 1.e-8;

  // Create a frequency ephemeris.
  FrequencyEph f_eph("TDB", since, until, epoch, 22., 45., 0.875, 1.125e-2, -2.25e-4, 6.75e-6);

  // Create a period ephemeris.
  // This is a set of values known to be the inverses of the frequency coefficients above.
  PeriodEph p_eph("TDB", since, until, epoch, 22., 45., 0.875, 88.8888888888888888888889,
    1.777777777777777777777778, 0.0177777777777777777777778);

  // Compare frequency & period.
  if (!f_eph.getEpoch().equivalentTo(p_eph.getEpoch(), tolerance))
    err() << "FrequencyEph and PeriodEph give different values for epoch" << std::endl;

  char * field[] = { "phi0", "f0", "f1" , "f2" };
  double value1[] = { f_eph.calcPulsePhase(epoch), f_eph.calcFrequency(epoch, 0), f_eph.calcFrequency(epoch, 1),
                      f_eph.calcFrequency(epoch, 2) };
  double value2[] = { p_eph.calcPulsePhase(epoch), p_eph.calcFrequency(epoch, 0), p_eph.calcFrequency(epoch, 1),
                      p_eph.calcFrequency(epoch, 2) };
  for (std::size_t ii = 0; ii != sizeof(value1) / sizeof(double); ++ii) {
    if (0. == value1[ii] || 0. == value2[ii]) {
      if (std::fabs(value1[ii] + value2[ii]) > std::numeric_limits<double>::epsilon())
        err() << "FrequencyEph and PeriodEph give absolutely different values for " << field[ii] <<
          " (" << value1[ii] << " and " << value2[ii] << ")" << std::endl;
    } else if (std::fabs(value1[ii] / value2[ii] - 1.) > epsilon) {
      err() << "FrequencyEph and PeriodEph give fractionally different values for " << field[ii] <<
        " (" << value1[ii] << " and " << value2[ii] << ")" << std::endl;
    }
  }

  // Create period ephemeris for regular testings.
  double ra = 22.;
  double dec = 45.;
  double phi0 = 0.875;
  double p0 = 1.23456789;
  double p1 = 9.87654321e-5;
  double p2 = 1.357902468e-10;
  std::auto_ptr<PeriodEph> eph(new PeriodEph("TDB", since, until, epoch, ra, dec, phi0, p0, p1, p2));

  // Test time system getter.
  std::string sys_name = eph->getSystem().getName();
  if ("TDB" != sys_name) {
    err() << "PeriodEph::getSystem returned \"" << sys_name << "\", not \"TDB\"" << std::endl;
  }

  // Test valid-since getter.
  const AbsoluteTime & since_returned = eph->getValidSince();
  if (!since_returned.equivalentTo(since, tolerance)) {
    err() << "PeriodEph::getValidSince returned AbsoluteTime(" << since_returned << "), not AbsoluteTime(" << since <<
      ")" << std::endl;
  }

  // Test valid-until getter.
  const AbsoluteTime & until_returned = eph->getValidUntil();
  if (!until_returned.equivalentTo(until, tolerance)) {
    err() << "PeriodEph::getValidUntil returned AbsoluteTime(" << until_returned << "), not AbsoluteTime(" << until <<
      ")" << std::endl;
  }

  // Test epoch getter.
  const AbsoluteTime & epoch_returned = eph->getEpoch();
  if (!epoch_returned.equivalentTo(epoch, tolerance)) {
    err() << "PeriodEph::getEpoch returned AbsoluteTime(" << epoch_returned << "), not AbsoluteTime(" << epoch <<
      ")" << std::endl;
  }

  // Test remark getter.
  if (eph->getRemark().size()) {
    err() << "PeriodEph::getRemark returned a non-empty container of remarks." << std::endl;
  }

  // Test the pulsar position.
  SourcePosition computed_pos = eph->calcPosition(epoch);
  SourcePosition correct_pos = SourcePosition(22., 45.);
  std::string coord_name("XYZ");
  for (std::size_t ii = 0; ii < 3; ++ii) {
    double computed_coord = computed_pos.getDirection()[ii];
    double correct_coord = correct_pos.getDirection()[ii];
    if (std::fabs(computed_coord - correct_coord) > epsilon) {
      err() << "PeriodEph::calcPosition returned " << coord_name[ii] << "=" << computed_coord << ", not " <<
        correct_coord << "." << std::endl;
    }
  }
  if (computed_pos.hasDistance() != false) {
    err() << "PeriodEph::calcPosition returned a source position with its distance known." << std::endl;
  }

  // Test frequency computation away from the reference epoch.
  double time_since_epoch = 1000.;
  double step_size = 100.;
  AbsoluteTime abs_time = epoch + ElapsedTime("TDB", Duration(time_since_epoch, "Sec"));
  std::auto_ptr<SimplePeriodEph> s_eph(new SimplePeriodEph(phi0, p0, p1, p2));
  epsilon = 1.e-3; // Note: Need a loose tolerance because SimplePeriodEph::calcFrequency is not that precise.
  double result = 0.;
  double expected = 0.;
  bool test_failed = true;
  for (int order = 0; order < 5; ++order) {
    result = eph->calcFrequency(abs_time, order);
    expected = s_eph->calcFrequency(time_since_epoch, step_size, order);
    if (0. == result || 0. == expected) test_failed = std::fabs(result + expected) > std::numeric_limits<double>::epsilon();
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
  result = eph->calcPulsePhase(abs_time, global_phase_offset);
  expected = s_eph->calcPulsePhase(time_since_epoch, step_size, global_phase_offset);
  test_failed = (std::fabs(result - expected) > phase_tolerance && std::fabs(result - expected + 1) > phase_tolerance
    && std::fabs(result - expected - 1) > phase_tolerance);
  if (test_failed) {
    err() << "PeriodEph::calcPulsePhase(abs_time, " << global_phase_offset << ") returned " << result <<
      ", not close enough to " << expected << " as expected" << std::endl;
  }

  global_phase_offset = 0.;
  result = eph->calcPulsePhase(abs_time, global_phase_offset);
  expected = s_eph->calcPulsePhase(time_since_epoch, step_size, global_phase_offset);
  test_failed = (std::fabs(result - expected) > phase_tolerance && std::fabs(result - expected + 1) > phase_tolerance
    && std::fabs(result - expected - 1) > phase_tolerance);
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
  for (std::size_t ii = 0; ii != sizeof (good_period_par)/sizeof(double[3]); ++ii) {
    p0 = good_period_par[ii][0];
    p1 = good_period_par[ii][1];
    p2 = good_period_par[ii][2];
    eph.reset(new PeriodEph("TDB", since, until, epoch, ra, dec, phi0, p0, p1, p2));
    s_eph.reset(new SimplePeriodEph(phi0, p0, p1, p2));
    result = eph->calcPulsePhase(abs_time);
    expected = s_eph->calcPulsePhase(time_since_epoch, step_size, 0.);
    test_failed = (std::fabs(result - expected) > phase_tolerance && std::fabs(result - expected + 1) > phase_tolerance
      && std::fabs(result - expected - 1) > phase_tolerance);
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
  for (std::size_t ii = 0; ii != sizeof (bad_period_par)/sizeof(double[3]); ++ii) {
    p0 = bad_period_par[ii][0];
    p1 = bad_period_par[ii][1];
    p2 = bad_period_par[ii][2];
    eph.reset(new PeriodEph("TDB", since, until, epoch, ra, dec, phi0, p0, p1, p2));
    try {
      eph->calcPulsePhase(abs_time);
      err() << "PeriodEph::calcPulsePhase(abs_time) did not throw an exception for p0=" << p0 <<
        ", p1=" << p1 << ", p2=" << p2 << std::endl;
    } catch (const std::exception &) {
      // This is fine.
    }
  }

  // Test the constructor that takes numerical arguments.
  eph.reset(new PeriodEph("TDB", AbsoluteTime("TDB", 12345, 51840.), AbsoluteTime("TDB", 23456, 60480.),
    AbsoluteTime("TDB", 34567, 69120.), 11., 22., 33., 44., 55., 66.));
  TextEph t_eph;
  t_eph.append("Valid Since", "12345.6 MJD (TDB)", "");
  t_eph.append("Valid Until", "23456.7 MJD (TDB)", "");
  t_eph.append("Epoch",       "34567.8 MJD (TDB)", "");
  t_eph.append("RA",   "11", "degrees");
  t_eph.append("Dec",  "22", "degrees");
  t_eph.append("Phi0", "33", "");
  t_eph.append("P0",   "44", "s");
  t_eph.append("P1",   "55", "");
  t_eph.append("P2",   "66", "s**(-1)");
  checkEphParameter(getMethod() + "_numeric", *eph, t_eph);

  // Test the constructor that takes a FITS record.
  std::string test_tpl("test_PeriodEph.tpl");
  remove(test_tpl.c_str());
  std::ofstream ofs(test_tpl.c_str());
  ofs << "\\include " << prependDataPath("PulsarDb_primary.tpl") << std::endl;
  ofs << "\\include " << prependDataPath("PulsarDb_spin_per.tpl") << std::endl;
  ofs.close();

  tip::IFileSvc::instance().createFile(getMethod() + ".fits", test_tpl);
  tip::TipFile tip_file = tip::IFileSvc::instance().openFile(getMethod() + ".fits");

  std::auto_ptr<tip::Table> table(tip_file.editTable("1"));
  tip::Header & header(table->getHeader());
  table->setNumRecords(1);
  tip::Table::Record record(table.get(), 0);
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
  eph.reset(new PeriodEph(record, header));
  t_eph.clear();
  t_eph.append("Valid Since", "54321 MJD (TDB)",   "");
  t_eph.append("Valid Until", "65433 MJD (TDB)",   ""); // The end of 65432 == the begining of 65433.
  t_eph.append("Epoch",       "76543.5 MJD (TDB)", "");
  t_eph.append("RA",   "1.1", "degrees");
  t_eph.append("Dec",  "2.2", "degrees");
  double dt = .864;
  double sqrt_term = std::sqrt(2.*3.3*5.5 - 4.4*4.4);
  phi0 = -(2. / sqrt_term * std::atan(sqrt_term*dt/(2.*3.3+4.4*dt)));
  while (phi0 < 0.) ++phi0;
  while (phi0 >= 1.) --phi0;
  t_eph.append("Phi0", phi0,  "");
  t_eph.append("P0",   "3.3", "s");
  t_eph.append("P1",   "4.4", "");
  t_eph.append("P2",   "5.5", "s**(-1)");
  checkEphParameter(getMethod() + "_fits", *eph, t_eph);
}

void PulsarDbTestApp::testHighPrecisionEph() {
  setMethod("testHighPrecisionEph");
  EphComputerAppTester tester(*this);

  // Prepare variables for tests.
  AbsoluteTime since("TDB", 51910, 0.);
  AbsoluteTime until("TDB", 51910, 1.);
  AbsoluteTime epoch("TDB", 51910, 123.456789);
  AbsoluteTime ev_time("TDB", 51910, 0.);
  ElapsedTime tolerance("TDB", Duration(1.e-9, "Sec")); // 1 nanosecond.
  const HighPrecisionEph::freq_type freq_empty(0);
  const HighPrecisionEph::wave_type wave_empty(0);
  const HighPrecisionEph::glitch_type glitch_empty(0);
  double elapsed = 12.345 * SecPerDay(); // 12.345 days in seconds.
  ev_time = epoch + ElapsedTime("TDB", Duration(elapsed, "Sec"));

  // Create a high-precision ephemeris.
  std::auto_ptr<HighPrecisionEph> eph(new HighPrecisionEph("TDB", since, until, epoch, 0., 0., 0., 0., 0., -1.,
    epoch, freq_empty, 1., wave_empty, wave_empty, glitch_empty));

  // Test time system getter.
  std::string sys_name = eph->getSystem().getName();
  if ("TDB" != sys_name) {
    err() << "HighPrecisionEph::getSystem returned \"" << sys_name << "\", not \"TDB\"" << std::endl;
  }

  // Test valid-since getter.
  const AbsoluteTime & since_returned = eph->getValidSince();
  if (!since_returned.equivalentTo(since, tolerance)) {
    err() << "HighPrecisionEph::getValidSince returned AbsoluteTime(" << since_returned << "), not AbsoluteTime(" << since <<
      ")" << std::endl;
  }

  // Test valid-until getter.
  const AbsoluteTime & until_returned = eph->getValidUntil();
  if (!until_returned.equivalentTo(until, tolerance)) {
    err() << "HighPrecisionEph::getValidUntil returned AbsoluteTime(" << until_returned << "), not AbsoluteTime(" << until <<
      ")" << std::endl;
  }

  // Test epoch getter.
  const AbsoluteTime & epoch_returned = eph->getEpoch();
  if (!epoch_returned.equivalentTo(epoch, tolerance)) {
    err() << "HighPrecisionEph::getEpoch returned AbsoluteTime(" << epoch_returned << "), not AbsoluteTime(" << epoch <<
      ")" << std::endl;
  }

  // Test remark getter, for an empty result.
  // Note: Tests with glitches are included below, following a test of each type of constructor.
  if (eph->getRemark().size()) {
    err() << "HighPrecisionEph::getRemark returned a non-empty container of remarks." << std::endl;
  }

  // Test pulse phase computation, only with frequency parameters.
  HighPrecisionEph::freq_type freq_pars(20);
  double factor = 1.;
  int nth_term = 1;
  for (HighPrecisionEph::freq_type::iterator itor = freq_pars.begin(); itor != freq_pars.end(); ++itor, ++nth_term) {
    *itor = 1.001 * nth_term * factor;
    factor *= nth_term / elapsed;
  }

  eph.reset(new HighPrecisionEph("TDB", since, until, epoch, 0., 0., 0., 0., 0., -1.,
    epoch, freq_pars, 1., wave_empty, wave_empty, glitch_empty));
  double phase = eph->calcPulsePhase(ev_time);
  double correct_phase_freq = freq_pars.size() * (freq_pars.size() + 1) / 2000.;
  double correct_phase = correct_phase_freq;
  double epsilon = 1.e-6;
  if (std::fabs(phase - correct_phase) > epsilon) {
    err() << "HighPrecisionEph::calcPulsePhase returned " << phase << ", not " << correct_phase << ", for " <<
      freq_pars.size() << "-term frequency parameters." << std::endl;
  }

  // Test pulse phase computation with the sine component of wave parameters, as well as frequency parameters.
  double wave_angle = 5.4321;
  double wave_omega = wave_angle / (elapsed / SecPerDay());
  HighPrecisionEph::wave_type wave_sine(20);
  nth_term = 1;
  for (HighPrecisionEph::wave_type::iterator itor = wave_sine.begin(); itor != wave_sine.end(); ++itor, ++nth_term) {
    *itor = 1.0001 * nth_term / std::sin(wave_angle * nth_term) / freq_pars[1];
  }

  eph.reset(new HighPrecisionEph("TDB", since, until, epoch, 0., 0., 0., 0., 0., -1.,
    epoch, freq_pars, wave_omega, wave_sine, wave_empty, glitch_empty));
  phase = eph->calcPulsePhase(ev_time);
  double correct_phase_sine = wave_sine.size() * (wave_sine.size() + 1) / 20000.;
  correct_phase = correct_phase_freq + correct_phase_sine;
  epsilon = 1.e-7;
  if (std::fabs(phase - correct_phase) > epsilon) {
    err() << "HighPrecisionEph::calcPulsePhase returned " << phase << ", not " << correct_phase << ", for " <<
      wave_sine.size() << "-term wave parameters (sine components only) on top of frequency parameters." << std::endl;
  }

  // Test pulse phase computation, only with the cosine component of wave parameters.
  HighPrecisionEph::wave_type wave_cosine(20);
  nth_term = 1;
  for (HighPrecisionEph::wave_type::iterator itor = wave_cosine.begin(); itor != wave_cosine.end(); ++itor, ++nth_term) {
    *itor = 1.00001 * nth_term / std::cos(wave_angle * nth_term) / freq_pars[1];
  }

  eph.reset(new HighPrecisionEph("TDB", since, until, epoch, 0., 0., 0., 0., 0., -1.,
    epoch, freq_pars, wave_omega, wave_empty, wave_cosine, glitch_empty));
  phase = eph->calcPulsePhase(ev_time);
  double correct_phase_cosine = wave_cosine.size() * (wave_cosine.size() + 1) / 200000.;
  correct_phase = correct_phase_freq + correct_phase_cosine;
  epsilon = 1.e-8;
  if (std::fabs(phase - correct_phase) > epsilon) {
    err() << "HighPrecisionEph::calcPulsePhase returned " << phase << ", not " << correct_phase << ", for " <<
      wave_cosine.size() << "-term wave parameters (cosine components only) on top of frequency parameters." << std::endl;
  }

  // Test pulse phase computation, only with glitch parameters.
  int num_past_glitch = 5;
  HighPrecisionEph::glitch_type glitch_list(2 * num_past_glitch);
  int nth_glitch = 1;
  int num_all_term = 0;
  for (HighPrecisionEph::glitch_type::iterator itor = glitch_list.begin(); itor != glitch_list.end(); ++itor, ++nth_glitch) {
    int sign = (nth_glitch % 2 ? +1 : -1);
    int num_jump_term = (nth_glitch + 1) / 2;
    int num_decay_term = (nth_glitch + 1) / 2;
    double dt_glitch = sign * nth_glitch * 12345.6789;
    itor->m_epoch = ev_time - ElapsedTime("TDB", Duration(dt_glitch, "Sec"));
    itor->m_perm_jump.resize(num_jump_term);
    factor = 1.;
    nth_term = num_all_term + 1;
    for (HighPrecisionEph::jump_type::iterator jump_itor = itor->m_perm_jump.begin(); jump_itor != itor->m_perm_jump.end();
      ++jump_itor, ++nth_term) {
      *jump_itor = sign * 1.000001 * nth_term * factor;
      factor *= (nth_term - num_all_term) / dt_glitch;
    }
    itor->m_decay_comp.resize(num_decay_term);
    for (HighPrecisionEph::decay_type::iterator decay_itor = itor->m_decay_comp.begin(); decay_itor != itor->m_decay_comp.end();
      ++decay_itor, ++nth_term) {
      double decay_time = nth_term * 1.2345; // in days.
      decay_itor->first = sign * 1.000001 * nth_term / (decay_time * SecPerDay())
        / (1. - std::exp(-dt_glitch / SecPerDay() / decay_time));
      decay_itor->second = decay_time;
    }
    if (nth_glitch % 2 == 0) num_all_term += num_jump_term + num_decay_term;
  }

  eph.reset(new HighPrecisionEph("TDB", since, until, epoch, 0., 0., 0., 0., 0., -1.,
    epoch, freq_empty, 1., wave_empty, wave_empty, glitch_list));
  phase = eph->calcPulsePhase(ev_time);
  double correct_phase_glitch = num_all_term * (num_all_term + 1) / 2000000.;
  correct_phase = correct_phase_glitch;
  epsilon = 1.e-9;
  if (std::fabs(phase - correct_phase) > epsilon) {
    err() << "HighPrecisionEph::calcPulsePhase returned " << phase << ", not " << correct_phase << ", for " <<
      glitch_list.size() << " glitches." << std::endl;
  }

  // Test pulse phase computation with all components.
  eph.reset(new HighPrecisionEph("TDB", since, until, epoch, 0., 0., 0., 0., 0., -1.,
    epoch, freq_pars, wave_omega, wave_sine, wave_cosine, glitch_list));
  phase = eph->calcPulsePhase(ev_time);
  correct_phase = correct_phase_freq + correct_phase_sine + correct_phase_cosine + correct_phase_glitch;
  epsilon = 1.e-9;
  if (std::fabs(phase - correct_phase) > epsilon) {
    err() << "HighPrecisionEph::calcPulsePhase returned " << phase << ", not " << correct_phase <<
      ", for all components combined." << std::endl;
  }

  // Test frequency computations for various orders of time-derivatives.
  epsilon = std::numeric_limits<double>::epsilon() * 100.;
  for (int nth_derivative = 0; nth_derivative < 10; ++nth_derivative) {
    // Test frequency computation, only with frequency parameters.
    factor = 1.;
    nth_term = 1;
    for (HighPrecisionEph::freq_type::iterator itor = freq_pars.begin(); itor != freq_pars.end(); ++itor, ++nth_term) {
      *itor = 1.1 * nth_term * factor;
      if (nth_term > nth_derivative + 1) factor *= (nth_term - nth_derivative - 1) / elapsed;
    }

    eph.reset(new HighPrecisionEph("TDB", since, until, epoch, 0., 0., 0., 0., 0., -1.,
      epoch, freq_pars, 1., wave_empty, wave_empty, glitch_empty));
    double frequency = eph->calcFrequency(ev_time, nth_derivative);
    double correct_frequency_freq = 1.1 * (freq_pars.size() - nth_derivative - 1) * (freq_pars.size() + nth_derivative + 2) / 2.;
    double correct_frequency = correct_frequency_freq;
    if (std::fabs(frequency/correct_frequency - 1.) > epsilon) {
      err() << "HighPrecisionEph::calcFrequency(ev_time, " << nth_derivative << ") returned " << frequency << ", not " <<
        correct_frequency << ", for " << freq_pars.size() << "-term frequency parameters." << std::endl;
    }

    // Test frequency computation with the sine component of wave parameters, as well as frequency parameters.
    nth_term = 1;
    for (HighPrecisionEph::wave_type::iterator itor = wave_sine.begin(); itor != wave_sine.end(); ++itor, ++nth_term) {
      double dsin_dt = 1.;
      for (int ii = 0; ii < nth_derivative + 1; ++ii) dsin_dt *= wave_omega / SecPerDay() * nth_term;
      switch (nth_derivative % 4){
      case 0: dsin_dt *= +std::cos(wave_angle * nth_term); break;
      case 1: dsin_dt *= -std::sin(wave_angle * nth_term); break;
      case 2: dsin_dt *= -std::cos(wave_angle * nth_term); break;
      case 3: dsin_dt *= +std::sin(wave_angle * nth_term); break;
      default: break;
      }
      *itor = 1.11 * nth_term / dsin_dt / freq_pars[1];
    }

    eph.reset(new HighPrecisionEph("TDB", since, until, epoch, 0., 0., 0., 0., 0., -1.,
      epoch, freq_pars, wave_omega, wave_sine, wave_empty, glitch_empty));
    frequency = eph->calcFrequency(ev_time, nth_derivative);
    double correct_frequency_sine = 1.11 * wave_sine.size() * (wave_sine.size() + 1) / 2.;
    correct_frequency = correct_frequency_freq + correct_frequency_sine;
    if (std::fabs(frequency/correct_frequency - 1.) > epsilon) {
      err() << "HighPrecisionEph::calcFrequency(ev_time, " << nth_derivative << ") returned " << frequency << ", not " <<
        correct_frequency << ", for " << wave_sine.size() <<
        "-term wave parameters (sine components only) on top of frequency parameters." << std::endl;
    }

    // Test frequency computation with the cosine component of wave parameters, as well as frequency parameters.
    nth_term = 1;
    for (HighPrecisionEph::wave_type::iterator itor = wave_cosine.begin(); itor != wave_cosine.end(); ++itor, ++nth_term) {
      double dcos_dt = 1.;
      for (int ii = 0; ii < nth_derivative + 1; ++ii) dcos_dt *= wave_omega / SecPerDay() * nth_term;
      switch (nth_derivative % 4){
      case 0: dcos_dt *= -std::sin(wave_angle * nth_term); break;
      case 1: dcos_dt *= -std::cos(wave_angle * nth_term); break;
      case 2: dcos_dt *= +std::sin(wave_angle * nth_term); break;
      case 3: dcos_dt *= +std::cos(wave_angle * nth_term); break;
      default: break;
      }
      *itor = 1.111 * nth_term / dcos_dt / freq_pars[1];
    }

    eph.reset(new HighPrecisionEph("TDB", since, until, epoch, 0., 0., 0., 0., 0., -1.,
      epoch, freq_pars, wave_omega, wave_empty, wave_cosine, glitch_empty));
    frequency = eph->calcFrequency(ev_time, nth_derivative);
    double correct_frequency_cosine = 1.111 * wave_cosine.size() * (wave_cosine.size() + 1) / 2.;
    correct_frequency = correct_frequency_freq + correct_frequency_cosine;
    if (std::fabs(frequency/correct_frequency - 1.) > epsilon) {
      err() << "HighPrecisionEph::calcFrequency(ev_time, " << nth_derivative << ") returned " << frequency << ", not " <<
        correct_frequency << ", for " << wave_cosine.size() <<
        "-term wave parameters (cosine components only) on top of frequency parameters." << std::endl;
    }

    // Test frequency computation, only with glitch parameters.
    int nth_glitch = 1;
    int num_all_term = 0;
    int rejected_term = 0;
    for (HighPrecisionEph::glitch_type::iterator itor = glitch_list.begin(); itor != glitch_list.end(); ++itor, ++nth_glitch) {
      int sign = (nth_glitch % 2 ? +1 : -1);
      int num_jump_term = (nth_glitch + 1) / 2;
      int num_decay_term = (nth_glitch + 1) / 2;
      double dt_glitch = sign * nth_glitch * 12345.6789;
      itor->m_epoch = ev_time - ElapsedTime("TDB", Duration(dt_glitch, "Sec"));
      itor->m_perm_jump.resize(num_jump_term);
      factor = 1.;
      nth_term = num_all_term + 1;
      for (HighPrecisionEph::jump_type::iterator jump_itor = itor->m_perm_jump.begin(); jump_itor != itor->m_perm_jump.end();
           ++jump_itor, ++nth_term) {
        *jump_itor = sign * 1.1111 * nth_term * factor;
        if (nth_term - num_all_term > nth_derivative + 1) factor *= (nth_term - num_all_term - nth_derivative - 1) / dt_glitch;
        else rejected_term += nth_term;
      }
      itor->m_decay_comp.resize(num_decay_term);
      for (HighPrecisionEph::decay_type::iterator decay_itor = itor->m_decay_comp.begin(); decay_itor != itor->m_decay_comp.end();
           ++decay_itor, ++nth_term) {
        double decay_time = nth_term * 1.2345; // in days.
        decay_itor->first = sign * 1.1111 * nth_term / std::exp(-dt_glitch / SecPerDay() / decay_time);
        for (int ii = 0; ii < nth_derivative; ++ii) decay_itor->first *= -decay_time * SecPerDay();
        decay_itor->second = decay_time;
      }
      if (nth_glitch % 2 == 0) num_all_term += num_jump_term + num_decay_term;
    }

    eph.reset(new HighPrecisionEph("TDB", since, until, epoch, 0., 0., 0., 0., 0., -1.,
      epoch, freq_empty, 1., wave_empty, wave_empty, glitch_list));
    frequency = eph->calcFrequency(ev_time, nth_derivative);
    double correct_frequency_glitch = 1.1111 * (num_all_term * (num_all_term + 1) - rejected_term) / 2.;
    correct_frequency = correct_frequency_glitch;
    if (std::fabs(frequency/correct_frequency - 1.) > epsilon) {
      err() << "HighPrecisionEph::calcFrequency(ev_time, " << nth_derivative << ") returned " << frequency << ", not " <<
        correct_frequency << ", for " << glitch_list.size() << " glitches." << std::endl;
    }

    // Test frequency computation with all components.
    eph.reset(new HighPrecisionEph("TDB", since, until, epoch, 0., 0., 0., 0., 0., -1.,
      epoch, freq_pars, wave_omega, wave_sine, wave_cosine, glitch_list));
    frequency = eph->calcFrequency(ev_time, nth_derivative);
    correct_frequency = correct_frequency_freq + correct_frequency_sine + correct_frequency_cosine + correct_frequency_glitch;
    if (std::fabs(frequency/correct_frequency - 1.) > epsilon) {
      err() << "HighPrecisionEph::calcFrequency(ev_time, " << nth_derivative << ") returned " << frequency << ", not " <<
        correct_frequency << ", for all components combined." << std::endl;
    }
  }

  // Test the pulsar position.
  double parsec = 3.26 * 365.25 * 86400.; // Approximately 1 parsec in light-seconds.
  const double sqrt14 = std::sqrt(14.);
  double par_list[][8] = {
    // Contains physical parameters of the source position and the proper motion, which are:
    //   RA & Dec (degree), distance (light-secod), dX, dY, & dZ (radian), displacement due to proper motion (light-second),
    //   and an elapsed time (second).
    { 60., +30., 4.*parsec, 1./sqrt14, 3./sqrt14, 2./sqrt14, 3.e-5*parsec, 3.e+7},
    { 60., -30., 5.*parsec, 3./sqrt14, 1./sqrt14, 2./sqrt14, 2.e-5*parsec, 3.e+7},
    {150., +30., 6.*parsec, 2./sqrt14, 3./sqrt14, 1./sqrt14, 1.e-5*parsec, 3.e+7},
    {150., -30., 7.*parsec, 3./sqrt14, 2./sqrt14, 1./sqrt14, 4.e-5*parsec, 3.e+7},
    {240., +30., 8.*parsec, 1./sqrt14, 2./sqrt14, 3./sqrt14, 5.e-5*parsec, 3.e+7},
    {240., -30., 3.*parsec, 2./sqrt14, 1./sqrt14, 3./sqrt14, 6.e-5*parsec, 3.e+7},
    {330., +30., 2.*parsec, 1./sqrt14, 3./sqrt14, 2./sqrt14, 7.e-5*parsec, 3.e+7},
    {330., -30., 1.*parsec, 3./sqrt14, 1./sqrt14, 2./sqrt14, 8.e-5*parsec, 3.e+7},
    // Repeat the parameters for a different period of time.
    { 60., +30., 4.*parsec, 1./sqrt14, 3./sqrt14, 2./sqrt14, 3.e-5*parsec, 7.e+7},
    { 60., -30., 5.*parsec, 3./sqrt14, 1./sqrt14, 2./sqrt14, 2.e-5*parsec, 7.e+7},
    {150., +30., 6.*parsec, 2./sqrt14, 3./sqrt14, 1./sqrt14, 1.e-5*parsec, 7.e+7},
    {150., -30., 7.*parsec, 3./sqrt14, 2./sqrt14, 1./sqrt14, 4.e-5*parsec, 7.e+7},
    {240., +30., 8.*parsec, 1./sqrt14, 2./sqrt14, 3./sqrt14, 5.e-5*parsec, 7.e+7},
    {240., -30., 3.*parsec, 2./sqrt14, 1./sqrt14, 3./sqrt14, 6.e-5*parsec, 7.e+7},
    {330., +30., 2.*parsec, 1./sqrt14, 3./sqrt14, 2./sqrt14, 7.e-5*parsec, 7.e+7},
    {330., -30., 1.*parsec, 3./sqrt14, 1./sqrt14, 2./sqrt14, 8.e-5*parsec, 7.e+7}
  };
  std::string coord_name("XYZ");
  epsilon = std::numeric_limits<double>::epsilon() * 100.;
  for (std::size_t ii = 0; ii != sizeof(par_list)/sizeof(par_list[0]); ++ii) {
    // Copy parameters of physical model.
    double * par = par_list[ii];
    double ra = par[0];
    double dec = par[1];
    double initial_distance = par[2];
    std::vector<double> pm_direction(par+3, par+6);
    double pm_distance = par[6];
    double elapsed_second = par[7];
    ev_time = epoch + ElapsedTime("TDB", Duration(elapsed_second, "Sec"));

    // Compute common parameters of the physical model.
    const double deg_per_rad = 180. / 3.14159265358979323846;
    const double mas_year_per_rad_sec = deg_per_rad * 3600. * 1000. * 365.25 * 86400.;
    std::vector<double> initial_direction(3);
    initial_direction[0] = std::cos(dec/deg_per_rad) * std::cos(ra/deg_per_rad);
    initial_direction[1] = std::cos(dec/deg_per_rad) * std::sin(ra/deg_per_rad);
    initial_direction[2] = std::sin(dec/deg_per_rad);
    std::vector<double> final_position(3);
    double final_distance = 0.;
    for (std::size_t jj = 0; jj < 3; ++jj) {
      final_position[jj] = initial_direction[jj] * initial_distance + pm_direction[jj] * pm_distance;
      final_distance += final_position[jj] * final_position[jj];
    }
    final_distance = std::sqrt(final_distance);
    std::vector<double> final_direction(3);
    for (std::size_t jj = 0; jj < 3; ++jj) final_direction[jj] = final_position[jj] / final_distance;
    SourcePosition ra_vel_srcpos(ra + 90., 0.);
    const std::vector<double> & ra_direction = ra_vel_srcpos.getDirection();
    SourcePosition dec_vel_srcpos(ra, dec + 90.);
    const std::vector<double> & dec_direction = dec_vel_srcpos.getDirection();
    
    // Case 1: stationary pulsar at an unknown distance.
    eph.reset(new HighPrecisionEph("TDB", since, until, epoch, ra, dec, 0., 0., 0., -1.,
      epoch, freq_empty, 1., wave_empty, wave_empty, glitch_empty));
    SourcePosition computed_pos = eph->calcPosition(ev_time);
    for (std::size_t jj = 0; jj < 3; ++jj) {
      double computed_coord = computed_pos.getDirection()[jj];
      double correct_coord = initial_direction[jj];
      if (std::fabs(computed_coord - correct_coord) > epsilon) {
        err() << "HighPrecisionEph::calcPosition returned " << coord_name[jj] << "=" << computed_coord << ", not " <<
          correct_coord << " (RA: " << ra << ", Dec: " << dec << ", RA velocity: " << 0. << ", Dec velocity: " << 0. <<
          ", radial velocity: " << 0. << ", parallax: unknown)." << std::endl;
      }
    }
    if (computed_pos.hasDistance() != false) {
      err() << "HighPrecisionEph::calcPosition returned a source position with its distance known (RA: " << ra <<
        ", Dec: " << dec << ", RA velocity: " << 0. << ", Dec velocity: " << 0. << ", radial velocity: " << 0. <<
        ", parallax: unknown)." << std::endl;
    }

    // Case 2: stationary pulsar at a known distance.
    double parallax = 1.49597870691e+11 / 2.99792458e+8 / initial_distance; // in radians.
    parallax *= deg_per_rad * 3600. * 1000.; // in milliarcseconds.

    eph.reset(new HighPrecisionEph("TDB", since, until, epoch, ra, dec, 0., 0., 0., parallax,
      epoch, freq_empty, 1., wave_empty, wave_empty, glitch_empty));
    computed_pos = eph->calcPosition(ev_time);
    for (std::size_t jj = 0; jj < 3; ++jj) {
      double computed_coord = computed_pos.getDirection()[jj];
      double correct_coord = initial_direction[jj];
      if (std::fabs(computed_coord - correct_coord) > epsilon) {
        err() << "HighPrecisionEph::calcPosition returned " << coord_name[jj] << "=" << computed_coord << ", not " <<
          correct_coord << " (RA: " << ra << ", Dec: " << dec << ", RA velocity: " << 0. << ", Dec velocity: " << 0. <<
          ", radial velocity: " << 0. << ", parallax: " << parallax << ")." << std::endl;
      }
    }
    if (computed_pos.hasDistance() != true) {
      err() << "HighPrecisionEph::calcPosition returned a source position with its distance unknown (RA: " << ra <<
        ", Dec: " << dec << ", RA velocity: " << 0. << ", Dec velocity: " << 0. << ", radial velocity: " << 0. <<
        ", parallax: " << parallax << ")." << std::endl;
    } else {
      double computed_dis = computed_pos.getDistance();
      double correct_dis = initial_distance;
      if (std::fabs(computed_dis/correct_dis - 1.) > epsilon) {
        err() << "HighPrecisionEph::calcPosition returned a source position with the distance of " << computed_dis <<
          ", not " << correct_dis << " (RA: " << ra << ", Dec: " << dec << ", RA velocity: " << 0. << ", Dec velocity: " <<
          0. << ", radial velocity: " << 0. << ", parallax: " << parallax << ")." << std::endl;
      }
    }

    // Case 3: pulsar moving along the line of signt at an unknown distance.
    double distance_change_rate = (final_distance - initial_distance) * 2.99792458e+5 / elapsed_second; // in km/s.

    eph.reset(new HighPrecisionEph("TDB", since, until, epoch, ra, dec, 0., 0., distance_change_rate, -1.,
      epoch, freq_empty, 1., wave_empty, wave_empty, glitch_empty));
    computed_pos = eph->calcPosition(ev_time);
    for (std::size_t jj = 0; jj < 3; ++jj) {
      double computed_coord = computed_pos.getDirection()[jj];
      double correct_coord = initial_direction[jj];
      if (std::fabs(computed_coord - correct_coord) > epsilon) {
        err() << "HighPrecisionEph::calcPosition returned " << coord_name[jj] << "=" << computed_coord << ", not " <<
          correct_coord << " (RA: " << ra << ", Dec: " << dec << ", RA velocity: " << 0. << ", Dec velocity: " << 0. <<
          ", radial velocity: " << distance_change_rate << ", parallax: unknown)." << std::endl;
      }
    }
    if (computed_pos.hasDistance() != false) {
      err() << "HighPrecisionEph::calcPosition returned a source position with its distance known (RA: " << ra <<
        ", Dec: " << dec << ", RA velocity: " << 0. << ", Dec velocity: " << 0. << ", radial velocity: " <<
        distance_change_rate << ", parallax: unknown)." << std::endl;
    }

    // Case 4: pulsar moving along the line of signt at a known distance.
    eph.reset(new HighPrecisionEph("TDB", since, until, epoch, ra, dec, 0., 0., distance_change_rate, parallax,
      epoch, freq_empty, 1., wave_empty, wave_empty, glitch_empty));
    computed_pos = eph->calcPosition(ev_time);
    for (std::size_t jj = 0; jj < 3; ++jj) {
      double computed_coord = computed_pos.getDirection()[jj];
      double correct_coord = initial_direction[jj];
      if (std::fabs(computed_coord - correct_coord) > epsilon) {
        err() << "HighPrecisionEph::calcPosition returned " << coord_name[jj] << "=" << computed_coord << ", not " <<
          correct_coord << " (RA: " << ra << ", Dec: " << dec << ", RA velocity: " << 0. << ", Dec velocity: " << 0. <<
          ", radial velocity: " << distance_change_rate << ", parallax: " << parallax << ")." << std::endl;
      }
    }
    if (computed_pos.hasDistance() != true) {
      err() << "HighPrecisionEph::calcPosition returned a source position with its distance unknown (RA: " << ra <<
        ", Dec: " << dec << ", RA velocity: " << 0. << ", Dec velocity: " << 0. << ", radial velocity: " <<
        distance_change_rate << ", parallax: " << parallax << ")." << std::endl;
    } else {
      double computed_dis = computed_pos.getDistance();
      double correct_dis = final_distance;
      if (std::fabs(computed_dis/correct_dis - 1.) > epsilon) {
        err() << "HighPrecisionEph::calcPosition returned a source position with the distance of " << computed_dis <<
          ", not " << correct_dis << " (RA: " << ra << ", Dec: " << dec << ", RA velocity: " << 0. << ", Dec velocity: " <<
          0. << ", radial velocity: " << distance_change_rate << ", parallax: " << parallax << ")." << std::endl;
      }
    }

    // Case 5: general moving pulsar at an unknown distance.
    double stretch_factor = 0.;
    for (std::size_t jj = 0; jj < 3; ++jj) stretch_factor += initial_direction[jj] * final_direction[jj];
    std::vector<double> pm_vector_trans(3);
    for (std::size_t jj = 0; jj < 3; ++jj) pm_vector_trans[jj] = final_direction[jj] / stretch_factor - initial_direction[jj];
    double ra_vel = 0.;
    for (std::size_t jj = 0; jj < 3; ++jj) ra_vel += pm_vector_trans[jj] * ra_direction[jj] / std::cos(dec/deg_per_rad) / elapsed_second;
    ra_vel *= mas_year_per_rad_sec; // in milliarcseconds per Julian year (365.25 days).
    double dec_vel = 0.;
    for (std::size_t jj = 0; jj < 3; ++jj) dec_vel += pm_vector_trans[jj] * dec_direction[jj] / elapsed_second;
    dec_vel *= mas_year_per_rad_sec; // in milliarcseconds per Julian year (365.25 days).
    double radial_vel = distance_change_rate;

    eph.reset(new HighPrecisionEph("TDB", since, until, epoch, ra, dec, ra_vel, dec_vel, radial_vel, -1.,
      epoch, freq_empty, 1., wave_empty, wave_empty, glitch_empty));
    computed_pos = eph->calcPosition(ev_time);
    for (std::size_t jj = 0; jj < 3; ++jj) {
      double computed_coord = computed_pos.getDirection()[jj];
      double correct_coord = final_direction[jj];
      if (std::fabs(computed_coord - correct_coord) > epsilon) {
        err() << "HighPrecisionEph::calcPosition returned " << coord_name[jj] << "=" << computed_coord << ", not " <<
          correct_coord << " (RA: " << ra << ", Dec: " << dec << ", RA velocity: " << ra_vel << ", Dec velocity: " <<
          dec_vel << ", radial velocity: " << radial_vel << ", parallax: unknown)." << std::endl;
      }
    }
    if (computed_pos.hasDistance() != false) {
      err() << "HighPrecisionEph::calcPosition returned a source position with its distance known (RA: " << ra <<
        ", Dec: " << dec << ", RA velocity: " << ra_vel << ", Dec velocity: " << dec_vel << ", radial velocity: " <<
        radial_vel << ", parallax: unknown)." << std::endl;
    }

    // Case 6: general moving pulsar at a known distance.
    double pm_speed_radial = pm_distance * 2.99792458e+5 / elapsed_second; // in km/s.
    double pm_speed_trans = pm_distance / initial_distance / elapsed_second; // in radians per second.
    pm_speed_trans *= mas_year_per_rad_sec; // in milliarcseconds per Julian year (365.25 days).
    ra_vel = 0.;
    for (std::size_t jj = 0; jj < 3; ++jj) ra_vel += pm_direction[jj] * pm_speed_trans * ra_direction[jj] / std::cos(dec/deg_per_rad);
    dec_vel = 0.;
    for (std::size_t jj = 0; jj < 3; ++jj) dec_vel += pm_direction[jj] * pm_speed_trans * dec_direction[jj];
    radial_vel = 0.;
    for (std::size_t jj = 0; jj < 3; ++jj) radial_vel += pm_direction[jj] * pm_speed_radial * initial_direction[jj];

    eph.reset(new HighPrecisionEph("TDB", since, until, epoch, ra, dec, ra_vel, dec_vel, radial_vel, parallax,
      epoch, freq_empty, 1., wave_empty, wave_empty, glitch_empty));
    computed_pos = eph->calcPosition(ev_time);
    for (std::size_t jj = 0; jj < 3; ++jj) {
      double computed_coord = computed_pos.getDirection()[jj];
      double correct_coord = final_direction[jj];
      if (std::fabs(computed_coord - correct_coord) > epsilon) {
        err() << "HighPrecisionEph::calcPosition returned " << coord_name[jj] << "=" << computed_coord << ", not " <<
          correct_coord << " (RA: " << ra << ", Dec: " << dec << ", RA velocity: " << ra_vel << ", Dec velocity: " <<
          dec_vel << ", radial velocity: " << radial_vel << ", parallax: " << parallax << ")." << std::endl;
      }
    }
    if (computed_pos.hasDistance() != true) {
      err() << "HighPrecisionEph::calcPosition returned a source position with its distance unknown (RA: " << ra <<
        ", Dec: " << dec << ", RA velocity: " << ra_vel << ", Dec velocity: " << dec_vel << ", radial velocity: " <<
        radial_vel << ", parallax: " << parallax << ")." << std::endl;
    } else {
      double computed_dis = computed_pos.getDistance();
      double correct_dis = final_distance;
      if (std::fabs(computed_dis/correct_dis - 1.) > epsilon) {
        err() << "HighPrecisionEph::calcPosition returned a source position with the distance of " << computed_dis <<
          ", not " << correct_dis << " (RA: " << ra << ", Dec: " << dec << ", RA velocity: " << ra_vel << ", Dec velocity: " <<
          dec_vel << ", radial velocity: " << radial_vel << ", parallax: " << parallax << ")." << std::endl;
      }
    }
  }

  // Test the constructor that takes numerical arguments.
  freq_pars.clear();
  freq_pars.push_back(.123456789);
  freq_pars.push_back(1.23456789);
  freq_pars.push_back(12.3456789);
  freq_pars.push_back(123.456789);
  freq_pars.push_back(1234.56789);
  freq_pars.push_back(12345.6789);
  wave_sine.clear();
  wave_sine.push_back(.1122334455);
  wave_sine.push_back(1.122334455);
  wave_sine.push_back(11.22334455);
  wave_cosine.clear();
  wave_cosine.push_back(.9988776655);
  wave_cosine.push_back(9.988776655);
  wave_cosine.push_back(99.88776655);
  wave_cosine.push_back(998.8776655);
  wave_cosine.push_back(9988.776655);
  glitch_list.clear();
  glitch_list.resize(2);
  HighPrecisionEph::glitch_type::iterator glitch_itor = glitch_list.begin();
  glitch_itor->m_epoch = AbsoluteTime("TDB", 54321,  8640.);
  glitch_itor->m_perm_jump.push_back(.987654321);
  glitch_itor->m_perm_jump.push_back(9.87654321);
  glitch_itor->m_perm_jump.push_back(98.7654321);
  ++glitch_itor;
  glitch_itor->m_epoch = AbsoluteTime("TDB", 65432, 17280.);
  glitch_itor->m_decay_comp.push_back(std::make_pair(0.00054321, 123.45));
  glitch_itor->m_decay_comp.push_back(std::make_pair(0.0056789, 54.321));
  eph.reset(new HighPrecisionEph("TDB", AbsoluteTime("TDB", 12345, 51840.), AbsoluteTime("TDB", 23456, 60480.),
    AbsoluteTime("TDB", 34567, 69120.), 11., 22., 33., 44., 55., 66., AbsoluteTime("TDB", 45678, 77760.),
    freq_pars, 77., wave_sine, wave_cosine, glitch_list));
  TextEph t_eph;
  t_eph.append("Valid Since",     "12345.6 MJD (TDB)", "");
  t_eph.append("Valid Until",     "23456.7 MJD (TDB)", "");
  t_eph.append("Position Epoch",  "34567.8 MJD (TDB)", "");
  t_eph.append("RA",              "11", "degrees");
  t_eph.append("Dec",             "22", "degrees");
  t_eph.append("RA Velocity",     "33", "milliarcseconds/yr");
  t_eph.append("Dec Velocity",    "44", "milliarcseconds/yr");
  t_eph.append("Radial Velocity", "55", "km/s");
  t_eph.append("Annual Parallax", "66", "milliarcseconds");
  t_eph.append("Frequency Epoch", "45678.9 MJD (TDB)", "");
  t_eph.append("Phi0",            "0.123456789", "");
  t_eph.append("F0",              "1.23456789",  "s**(-1)");
  t_eph.append("F1",              "12.3456789",  "s**(-2)");
  t_eph.append("F2",              "123.456789",  "s**(-3)");
  t_eph.append("F3",              "1234.56789",  "s**(-4)");
  t_eph.append("F4",              "12345.6789",  "s**(-5)");
  t_eph.append("Wave Frequency",  "77", "radians/day");
  t_eph.append("Sin1",  "0.1122334455", "s");
  t_eph.append("Cos1",  "0.9988776655", "s");
  t_eph.append("Sin2",  "1.122334455",  "s");
  t_eph.append("Cos2",  "9.988776655",  "s");
  t_eph.append("Sin3",  "11.22334455",  "s");
  t_eph.append("Cos3",  "99.88776655",  "s");
  t_eph.append("Cos4",  "998.8776655",  "s");
  t_eph.append("Cos5",  "9988.776655",  "s");
  t_eph.append("Glitch Epoch", "54321.1 MJD (TDB)", "");
  t_eph.append("dPhi0", "0.987654321", "");
  t_eph.append("dF0",   "9.87654321",  "s**(-1)");
  t_eph.append("dF1",   "98.7654321",  "s**(-2)");
  t_eph.append("Glitch Epoch", "65432.2 MJD (TDB)", "");
  t_eph.append("Amplitude",  "0.00054321", "s**(-1)");
  t_eph.append("Decay Time", "123.45",     "days");
  t_eph.append("Amplitude",  "0.0056789",  "s**(-1)");
  t_eph.append("Decay Time", "54.321",     "days");
  checkEphParameter(getMethod() + "_numeric", *eph, t_eph);

  // Test glitch reporting capability, with the constructor with numerical arguments.
  EphStatusCont remark_list = eph->getRemark();
  if (2 != remark_list.size()) {
    err() << "HighPrecisionEph object constructed with numerical arguments reported " << remark_list.size() <<
      " ephemeris remark(s), not 2." << std::endl;
  } else {
    std::list<std::string> report_list;
    report_list.push_back("Remarked \"Glitch observed at 54321.1 MJD (TDB)\" since 54321.1 MJD (TDB) until 23456.7 MJD (TDB)");
    report_list.push_back("Remarked \"Glitch observed at 65432.2 MJD (TDB)\" since 65432.2 MJD (TDB) until 23456.7 MJD (TDB)");
    std::list<std::string>::const_iterator str_itor = report_list.begin();
    int remark_number = 1;
    for (EphStatusCont::const_iterator rem_itor = remark_list.begin(); rem_itor != remark_list.end() && str_itor != report_list.end();
      ++rem_itor, ++str_itor, ++remark_number) {
      std::ostringstream oss;
      oss << rem_itor->report("TDB", MjdFmt);
      std::ostringstream dummy_oss;
      if (!tester.verify(oss.str(), *str_itor, dummy_oss)) {
        err() << "HighPrecisionEph object constructed with numerical arguments reported '" << oss.str() <<
          "' as ephemeris remark No. " << remark_number << ", not '" << *str_itor << "'." << std::endl;
      }
    }
  }

  // Test the constructor that takes a FITS record.
  std::string test_tpl("test_HighPrecisionEph.tpl");
  remove(test_tpl.c_str());
  std::ofstream ofs(test_tpl.c_str());
  ofs << "\\include " << prependDataPath("PulsarDb_primary.tpl") << std::endl;
  ofs << "\\include " << prependDataPath("PulsarDb_spin_hp.tpl") << std::endl;
  ofs.close();

  tip::IFileSvc::instance().createFile(getMethod() + ".fits", test_tpl);
  tip::TipFile tip_file = tip::IFileSvc::instance().openFile(getMethod() + ".fits");

  std::auto_ptr<tip::Table> table(tip_file.editTable("1"));
  tip::Header & header(table->getHeader());
  table->setNumRecords(1);
  tip::Table::Record record(table.get(), 0);
  record["VALID_SINCE"].set(12345);
  record["VALID_UNTIL"].set(23456);
  record["POS_EPOCH_INT"].set(34567);
  record["POS_EPOCH_FRAC"].set(.8);
  record["RA"].set(1.1);
  record["DEC"].set(2.2);
  record["RA_VELOCITY"].set(3.3);
  record["DEC_VELOCITY"].set(4.4);
  record["RADIAL_VELOCITY"].set(5.5);
  record["PARALLAX"].set(6.6);
  record["FREQ_EPOCH_INT"].set(45678);
  record["FREQ_EPOCH_FRAC"].set(.9);
  record["TOABARY_INT"].set(45678);
  record["TOABARY_FRAC"].set(.9001); // 8.64 seconds off.
  std::vector<double> pars;
  pars.push_back(1.23456789);
  pars.push_back(12.3456789);
  pars.push_back(123.456789);
  pars.push_back(1234.56789);
  pars.push_back(12345.6789);
  record["FREQ_PARAMETERS"].set(pars);
  record["WAVE_OMEGA"].set(7.7);
  pars.clear();
  pars.push_back(.1122334455);
  pars.push_back(1.122334455);
  pars.push_back(11.22334455);
  record["WAVE_SINE"].set(pars);
  pars.clear();
  pars.push_back(.9988776655);
  pars.push_back(9.988776655);
  pars.push_back(99.88776655);
  pars.push_back(998.8776655);
  pars.push_back(9988.776655);
  record["WAVE_COSINE"].set(pars);
  pars.clear();
  pars.push_back(54321.1);
  pars.push_back(0.987654321);
  pars.push_back(9.87654321);
  pars.push_back(98.7654321);
  pars.push_back(65432.2);
  pars.push_back(0.00054321);
  pars.push_back(123.45);
  pars.push_back(0.0056789);
  pars.push_back(54.321);
  record["GLITCH_PARAMETERS"].set(pars);
  pars.clear();
  pars.push_back(1);
  pars.push_back(3);
  pars.push_back(0);
  pars.push_back(1);
  pars.push_back(0);
  pars.push_back(4);
  record["GLITCH_DIMENSIONS"].set(pars);
  eph.reset(new HighPrecisionEph(record, header));
  t_eph.clear();
  t_eph.append("Valid Since", "12345 MJD (TDB)", "");
  t_eph.append("Valid Until", "23457 MJD (TDB)", ""); // The end of 23456 == the begining of 23457.
  t_eph.append("Position Epoch",  "34567.8 MJD (TDB)", "");
  t_eph.append("RA",              "1.1", "degrees");
  t_eph.append("Dec",             "2.2", "degrees");
  t_eph.append("RA Velocity",     "3.3", "milliarcseconds/yr");
  t_eph.append("Dec Velocity",    "4.4", "milliarcseconds/yr");
  t_eph.append("Radial Velocity", "5.5", "km/s");
  t_eph.append("Annual Parallax", "6.6", "milliarcseconds");
  t_eph.append("Frequency Epoch", "45678.9 MJD (TDB)", "");
  double dt = 8.64;
  double int_part = 0.;
  double phi_at_toa = 0.;
  phi_at_toa += std::modf(1.23456789*dt, &int_part);
  phi_at_toa += std::modf(12.3456789/2.*dt*dt, &int_part);
  phi_at_toa += std::modf(123.456789/6.*dt*dt*dt, &int_part);
  phi_at_toa += std::modf(1234.56789/24.*dt*dt*dt*dt, &int_part);
  phi_at_toa += std::modf(12345.6789/120.*dt*dt*dt*dt*dt, &int_part);
  double f0 = 1.23456789;
  double omega_at_toa = 7.7 * dt / SecPerDay();
  phi_at_toa += std::modf(.1122334455 * f0 * std::sin(1. * omega_at_toa), &int_part);
  phi_at_toa += std::modf(1.122334455 * f0 * std::sin(2. * omega_at_toa), &int_part);
  phi_at_toa += std::modf(11.22334455 * f0 * std::sin(3. * omega_at_toa), &int_part);
  phi_at_toa += std::modf(.9988776655 * f0 * std::cos(1. * omega_at_toa), &int_part);
  phi_at_toa += std::modf(9.988776655 * f0 * std::cos(2. * omega_at_toa), &int_part);
  phi_at_toa += std::modf(99.88776655 * f0 * std::cos(3. * omega_at_toa), &int_part);
  phi_at_toa += std::modf(998.8776655 * f0 * std::cos(4. * omega_at_toa), &int_part);
  phi_at_toa += std::modf(9988.776655 * f0 * std::cos(5. * omega_at_toa), &int_part);
  double phi0 = -phi_at_toa;
  while (phi0 < 0.) ++phi0;
  while (phi0 >= 1.) --phi0;
  t_eph.append("Phi0", phi0, "");
  t_eph.append("F0",   "1.23456789", "s**(-1)");
  t_eph.append("F1",   "12.3456789", "s**(-2)");
  t_eph.append("F2",   "123.456789", "s**(-3)");
  t_eph.append("F3",   "1234.56789", "s**(-4)");
  t_eph.append("F4",   "12345.6789", "s**(-5)");
  t_eph.append("Wave Frequency", "7.7", "radians/day");
  t_eph.append("Sin1",  "0.1122334455", "s");
  t_eph.append("Cos1",  "0.9988776655", "s");
  t_eph.append("Sin2",  "1.122334455",  "s");
  t_eph.append("Cos2",  "9.988776655",  "s");
  t_eph.append("Sin3",  "11.22334455",  "s");
  t_eph.append("Cos3",  "99.88776655",  "s");
  t_eph.append("Cos4",  "998.8776655",  "s");
  t_eph.append("Cos5",  "9988.776655",  "s");
  t_eph.append("Glitch Epoch", "54321.1 MJD (TDB)", "");
  t_eph.append("dPhi0", "0.987654321", "");
  t_eph.append("dF0",   "9.87654321",  "s**(-1)");
  t_eph.append("dF1",   "98.7654321",  "s**(-2)");
  t_eph.append("Glitch Epoch", "65432.2 MJD (TDB)", "");
  t_eph.append("Amplitude",  "0.00054321", "s**(-1)");
  t_eph.append("Decay Time", "123.45",     "days");
  t_eph.append("Amplitude",  "0.0056789",  "s**(-1)");
  t_eph.append("Decay Time", "54.321",     "days");
  checkEphParameter(getMethod() + "_fits", *eph, t_eph);

  // Test glitch reporting capability, with the constructor with a FITS record.
  remark_list = eph->getRemark();
  if (2 != remark_list.size()) {
    err() << "HighPrecisionEph object constructed with a FITS record reported " << remark_list.size() <<
      " ephemeris remark(s), not 2." << std::endl;
  } else {
    std::list<std::string> report_list;
    report_list.push_back("Remarked \"Glitch observed at 54321.1 MJD (TDB)\" since 54321.1 MJD (TDB) until 23457 MJD (TDB)");
    report_list.push_back("Remarked \"Glitch observed at 65432.2 MJD (TDB)\" since 65432.2 MJD (TDB) until 23457 MJD (TDB)");
    std::list<std::string>::const_iterator str_itor = report_list.begin();
    int remark_number = 1;
    for (EphStatusCont::const_iterator rem_itor = remark_list.begin(); rem_itor != remark_list.end() && str_itor != report_list.end();
      ++rem_itor, ++str_itor, ++remark_number) {
      std::ostringstream oss;
      oss << rem_itor->report("TDB", MjdFmt);
      std::ostringstream dummy_oss;
      if (!tester.verify(oss.str(), *str_itor, dummy_oss)) {
        err() << "HighPrecisionEph object constructed with a FITS record reported '" << oss.str() <<
          "' as ephemeris remark No. " << remark_number << ", not '" << *str_itor << "'." << std::endl;
      }
    }
  }

  // Test detection of glitch parameter errors.
  // TODO: Use tip::TableCell::setNull method once it is implemented in tip.
  const double null_value = std::numeric_limits<double>::min();
  pars.clear();
  pars.push_back(54321.1);
  pars.push_back(0.987654321);
  pars.push_back(9.87654321);
  pars.push_back(98.7654321);
  pars.push_back(null_value);
  pars.push_back(0.00054321);
  pars.push_back(123.45);
  pars.push_back(0.0056789);
  pars.push_back(54.321);
  record["GLITCH_PARAMETERS"].set(pars);
  try {
    HighPrecisionEph eph_obj(record, header);
    err() << "HighPrecisionEph constructor did not throw an exception for an undefined glitch epoch" << std::endl;
  } catch (const std::exception &) {
    // This is fine.
  }

  pars.clear();
  pars.push_back(54321.1);
  pars.push_back(0.987654321);
  pars.push_back(9.87654321);
  pars.push_back(98.7654321);
  pars.push_back(65432.2);
  pars.push_back(0.00054321);
  pars.push_back(null_value);
  pars.push_back(0.0056789);
  pars.push_back(54.321);
  record["GLITCH_PARAMETERS"].set(pars);
  try {
    HighPrecisionEph eph_obj(record, header);
    err() << "HighPrecisionEph constructor did not throw an exception for an undefined glitch decay time" << std::endl;
  } catch (const std::exception &) {
    // This is fine.
  }

  pars.clear();
  pars.push_back(1);
  pars.push_back(3);
  pars.push_back(0);
  pars.push_back(1);
  pars.push_back(0);
  pars.push_back(5);
  record["GLITCH_DIMENSIONS"].set(pars);
  try {
    HighPrecisionEph eph_obj(record, header);
    err() << "HighPrecisionEph constructor did not throw an exception for unmatched GLITCH_DIMENSIONS (1, 3, 0, 1, 0, 5)" <<
      std::endl;
  } catch (const std::exception &) {
    // This is fine.
  }

  pars.clear();
  pars.push_back(2);
  pars.push_back(2);
  pars.push_back(0);
  pars.push_back(1);
  pars.push_back(0);
  pars.push_back(4);
  record["GLITCH_DIMENSIONS"].set(pars);
  try {
    HighPrecisionEph eph_obj(record, header);
    err() << "HighPrecisionEph constructor did not throw an exception for having two fields for a glitch epoch" <<
      std::endl;
  } catch (const std::exception &) {
    // This is fine.
  }
}

void PulsarDbTestApp::testSimpleDdEph() {
  setMethod("testSimpleDdEph");

  // Prepare variables for testing ephemeris computations.
  AbsoluteTime t0("TDB", 51910, 123.456789);
  std::auto_ptr<SimpleDdEph> eph(new SimpleDdEph("TDB", 1000., .2, 0., 0., 0., 0., 0., 0., t0, 0., 0., 0.));
  AbsoluteTime ev_time("TDB", 51910, 223.456789);
  double epsilon = 1.e-8;

  // Test time system getter.
  std::string sys_name = eph->getSystem().getName();
  if ("TDB" != sys_name) {
    err() << "SimpleDdEph::getSystem returned \"" << sys_name << "\", not \"TDB\"" << std::endl;
  }

  // Test orbital phase computations.
  double phase = eph->calcOrbitalPhase(ev_time);
  double phase_expected = .099;
  if (std::fabs(phase/phase_expected - 1.) > epsilon) {
    err() << "SimpleDdEph::calcOrbitalPhase produced phase == " << phase << ", not " << phase_expected << std::endl;
  }

  // Test orbital phase computations, with a non-zero global phase offset.
  phase = eph->calcOrbitalPhase(ev_time, 0.1234);
  phase_expected += 0.1234;
  if (std::fabs(phase/phase_expected - 1.) > epsilon) {
    err() << "SimpleDdEph::calcOrbitalPhase produced phase == " << phase << ", not " << phase_expected << std::endl;
  }

  // Test orbital delay computations.
  // Note: Below the equations in Taylor et al. (ApJ, 345, 434-450, 1989) are computed in reverse.
  int turn_array[] = {0, 1, -1, 2, -2, 5, -5, 20, -20, 100, -100};
  std::list<int> turn_list(turn_array, turn_array + sizeof(turn_array)/sizeof(int));
  for (std::list<int>::const_iterator turn_itor = turn_list.begin(); turn_itor != turn_list.end(); ++turn_itor) {
    // 1) Set constants.
    double par_PB = 10. * 86400.; // 10 days in seconds.
    double max_elapsed = par_PB * (std::fabs(*turn_itor) + 1.);
    double par_XDOT = -0.44 / max_elapsed;
    double par_ECCDOT = -0.00055 / max_elapsed;
    double par_OMDOT = 6.6 / max_elapsed;
    double par_GAMMA = 0.0011;
    double par_SHAPIRO_R = 0.0022; // in micro-seconds.
    double par_SHAPIRO_S = 0.0033;

    // 2) Set values to variables that appear in equations 8, 9, and 10.
    double var_e = 0.6; // 1-e^2 = 0.64, sqrt(1-e^2) = 0.8.
    double var_w = 30. / 180. * M_PI; // 30 degrees in radian (sin w = 1/2, cos w = sqrt(3)/2).
    double var_u = 60. / 180. * M_PI; // 60 degrees in radian (sin u = sqrt(3)/2, cos u = 1/2, tan u/2 = sqrt(3)/3).
    double var_x = 10.; // 10 light-seconds.
    double var_r = par_SHAPIRO_R * 1.e-6; // micro-seconds to seconds.

    // 3) Compute orbital delay and true anomaly.
    double delay = var_x * 1./2. * (1./2. - 0.6) + var_x * 0.8 * 3./4.; // Equation 8.
    delay += par_GAMMA * std::sqrt(3.)/2.; // Equation 9.
    delay += -2. * var_r * std::log(1. - 0.6*1./2. - par_SHAPIRO_S*(1./2.*(1./2. - 0.6) + 0.8*3./4.)); // Equation 10.
    double true_anomaly = 2. * std::atan(2. * std::sqrt(3.)/3.); // Equation 13.

    // 4) Modify variables for additional turns of the binary system.
    var_u += 2. * M_PI * (*turn_itor);
    true_anomaly += 2. * M_PI * (*turn_itor);

    // 5) Compute back orbital parameters so as to reproduce values set in step 2.
    double elapsed = var_u / 2. / M_PI * par_PB; // Equation 12.
    double par_PBDOT = 0.6 * std::sqrt(3.)/2. / M_PI * (par_PB/elapsed) * (par_PB/elapsed); // Equation 12.
    double par_A1 = var_x - par_XDOT * elapsed;
    double par_ECC = var_e - par_ECCDOT * elapsed;
    double var_wdot = par_OMDOT / 180.0 * M_PI / 365.25 / 86400.; // degrees-per-year to radian-per-second.
    double var_k = var_wdot / 2. / M_PI * par_PB;
    double par_OM = var_w - var_k * true_anomaly; // Equation 14.
    par_OM *= 180. / M_PI; // radians to degrees.
    ev_time = t0 + ElapsedTime("TDB", Duration(elapsed, "Sec"));

    // 6) Create an ephemeris object and test its calcOrbitalDelay method.
    eph.reset(new SimpleDdEph("TDB", par_PB, par_PBDOT, par_A1, par_XDOT, par_ECC, par_ECCDOT, par_OM, par_OMDOT, t0,
      par_GAMMA, par_SHAPIRO_R, par_SHAPIRO_S));
    ElapsedTime delay_result = eph->calcOrbitalDelay(ev_time);
    ElapsedTime delay_expected("TDB", Duration(delay, "Sec"));
    Duration delay_tolerance(std::numeric_limits<double>::epsilon() * 86400. * 100., "Sec");
    if (!delay_result.getDuration().equivalentTo(delay_expected.getDuration(), delay_tolerance))
      err() << "SimpleDdEph::calcOrbitalDelay(" << ev_time << ") returned " << delay_result << " after " << *turn_itor <<
        " turn(s), not equivalent to " << delay_expected << " as expected." << std::endl;
  }

  // Test binary modulation and demodulation.
  // Binary parameters: (PB, PBDOT, A1, XDOT, ECC, ECCDOT, OM, OMDOT, T0, GAMMA, SHAPIRO_R, SHAPIRO_S)
  // = (27906.980897, -2.43e-12, 2.3417598, 0.0, 0.61713101, 0.0, 220.142729, 4.22662, 45888.1172487, 0.004295, 0.0, 0.0)
  // Note: OMDOT has been changed to 4.22662*365.25/365.0 to adopt the bug fix (by M. Hirayama on March 16th, 2010).
  AbsoluteTime abs_t0("TDB", Mjd(45888, .1172487));
  eph.reset(new SimpleDdEph("TDB", 27906.980897, -2.43e-12, 2.3417598, 0.0, 0.61713101, 0.0, 220.142729, 4.22662 * 365.25 / 365.0,
    abs_t0, 0.004295, 0.0, 0.0));

  // MJD's: { {original-MJD, modulated-MJD}, ... }
  Mjd mjd_test_values[][2] = {
    { Mjd(45988, .1172486599971307), Mjd(45988, .1172823346967320) },
    { Mjd(45988, .1519708822161192), Mjd(45988, .1520057480382526) },
    { Mjd(45988, .1866931044423836), Mjd(45988, .1867233097334768) },
    { Mjd(45988, .2214153266613721), Mjd(45988, .2214303721880313) },
    { Mjd(45988, .2561375488876365), Mjd(45988, .2561254377882811) },
    { Mjd(45988, .2908597711066250), Mjd(45988, .2908532530402717) },
    { Mjd(45988, .3255819933328894), Mjd(45988, .3255877721779576) },
    { Mjd(45988, .3603042155518779), Mjd(45988, .3603213319299883) },
    { Mjd(45988, .3950264377781423), Mjd(45988, .3950526724878252) },
    { Mjd(45988, .4297486599971307), Mjd(45988, .4297811300504257) },
    { Mjd(45988, .4644708822161192), Mjd(45988, .4645058955188297) }
  };

  // Permitted difference is 100 ns.
  double delta = 100. * 1.e-9;
  ElapsedTime tolerance("TDB", Duration(delta, "Sec"));

  setPrecision(24);
  for (std::size_t ii = 0; ii != sizeof(mjd_test_values)/sizeof(Mjd[2]); ++ii) {
    AbsoluteTime tdb_mjd("TDB", mjd_test_values[ii][0]);
    AbsoluteTime expected_tdb_mjd("TDB", mjd_test_values[ii][1]);
    AbsoluteTime original_tdb_mjd(tdb_mjd);
    eph->modulateBinary(tdb_mjd);
    if (!tdb_mjd.equivalentTo(expected_tdb_mjd, tolerance)) {
      err() << "Binary modulation of " << original_tdb_mjd.represent("TDB", MjdFmt) << " by SimpleDdEph was computed to be " <<
        tdb_mjd.represent("TDB", MjdFmt) << ", not " << expected_tdb_mjd.represent("TDB", MjdFmt) << ", as expected." << std::endl;
    }
  }

  for (std::size_t ii = 0; ii != sizeof(mjd_test_values)/sizeof(Mjd[2]); ++ii) {
    AbsoluteTime tdb_mjd("TDB", mjd_test_values[ii][1]);
    AbsoluteTime expected_tdb_mjd("TDB", mjd_test_values[ii][0]);
    AbsoluteTime original_tdb_mjd(tdb_mjd);
    eph->demodulateBinary(tdb_mjd);
    if (!tdb_mjd.equivalentTo(expected_tdb_mjd, tolerance)) {
      err() << "Binary demodulation of " << original_tdb_mjd.represent("TDB", MjdFmt) << " by SimpleDdEph was computed to be " <<
        tdb_mjd.represent("TDB", MjdFmt) << ", not " << expected_tdb_mjd.represent("TDB", MjdFmt) << ", as expected." << std::endl;
    }
  }

  // Test the constructor that takes numerical arguments.
  const double rad_per_deg = M_PI / 180.;
  const double rad_year_per_deg_sec = rad_per_deg / (365.25 * 86400.);
  eph.reset(new SimpleDdEph("TDB", 12345.6789, 1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7, AbsoluteTime("TDB", 12345, 51840.), 8.8, 9.9, 11.1));
  TextEph t_eph;
  t_eph.append("PB",        "12345.6789",               "s");
  t_eph.append("PBDOT",     "1.1",                      "");
  t_eph.append("A1",        "2.2",                      "lt-s");
  t_eph.append("XDOT",      "3.3",                      "lt-s/s");
  t_eph.append("ECC",       "4.4",                      "");
  t_eph.append("ECCDOT",    "5.5",                      "s**(-1)");
  t_eph.append("OM",        6.6 * rad_per_deg,          "radians");
  t_eph.append("OMDOT",     7.7 * rad_year_per_deg_sec, "radians/s");
  t_eph.append("T0",        "12345.6 MJD (TDB)",        "");
  t_eph.append("GAMMA",     "8.8",                      "s");
  t_eph.append("SHAPIRO_R", "9.9e-06",                  "s");
  t_eph.append("SHAPIRO_S", "11.1",                     "");
  checkEphParameter(getMethod() + "_numeric", *eph, t_eph);

  // Test the constructor that takes a FITS record.
  std::string test_tpl("test_SimpleDdEph.tpl");
  remove(test_tpl.c_str());
  std::ofstream ofs(test_tpl.c_str());
  ofs << "\\include " << prependDataPath("PulsarDb_primary.tpl") << std::endl;
  ofs << "\\include " << prependDataPath("PulsarDb_orbital_dd.tpl") << std::endl;
  ofs.close();

  tip::IFileSvc::instance().createFile(getMethod() + ".fits", test_tpl);
  tip::TipFile tip_file = tip::IFileSvc::instance().openFile(getMethod() + ".fits");

  std::auto_ptr<tip::Table> table(tip_file.editTable("1"));
  tip::Header & header(table->getHeader());
  table->setNumRecords(1);
  tip::Table::Record record(table.get(), 0);
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
  eph.reset(new SimpleDdEph(record, header));
  checkEphParameter(getMethod() + "_fits", *eph, t_eph);
} 

void PulsarDbTestApp::testBtModelEph() {
  setMethod("testBtModelDdEph");

  // Prepare variables for testing ephemeris computations.
  AbsoluteTime t0("TDB", 51910, 123.456789);
  std::auto_ptr<BtModelEph> eph(new BtModelEph("TDB", 1000., .2, 0., 0., 0., 0., 0., 0., t0, 0.));
  AbsoluteTime ev_time("TDB", 51910, 223.456789);
  double epsilon = 1.e-8;

  // Test time system getter.
  std::string sys_name = eph->getSystem().getName();
  if ("TDB" != sys_name) {
    err() << "BtModelEph::getSystem returned \"" << sys_name << "\", not \"TDB\"" << std::endl;
  }

  // Test orbital phase computations.
  double phase = eph->calcOrbitalPhase(ev_time);
  double phase_expected = 100./1010.;
  if (std::fabs(phase/phase_expected - 1.) > epsilon) {
    err() << "BtModelEph::calcOrbitalPhase produced phase == " << phase << ", not " << phase_expected << std::endl;
  }

  // Test orbital phase computations, with a non-zero global phase offset.
  phase = eph->calcOrbitalPhase(ev_time, 0.1234);
  phase_expected += 0.1234;
  if (std::fabs(phase/phase_expected - 1.) > epsilon) {
    err() << "BtModelEph::calcOrbitalPhase produced phase == " << phase << ", not " << phase_expected << std::endl;
  }

  // Test orbital delay computations.
  // Note: Below the equations in Blandford et al. (ApJ 205, 580-591, 1976) are computed in reverse.
  int turn_array[] = {0, 1, -1, 2, -2, 5, -5, 20, -20, 100, -100};
  std::list<int> turn_list(turn_array, turn_array + sizeof(turn_array)/sizeof(int));
  for (std::list<int>::const_iterator turn_itor = turn_list.begin(); turn_itor != turn_list.end(); ++turn_itor) {
    // 1) Set constants.
    double var_PB = 10. * 86400.; // 10 days in seconds.
    double max_elapsed = var_PB * (std::fabs(*turn_itor) + 1.);
    double par_PBDOT = -0.33 * 86400. / max_elapsed;
    double par_XDOT = -0.44 / max_elapsed;
    double par_ECCDOT = -0.00055 / max_elapsed;
    double par_OMDOT = 6.6 / max_elapsed;
    double par_GAMMA = 0.0011;

    // 2) Set values to variables that appear in equations 2.27, 2.30 and 2.31.
    double var_e = 0.6; // 1-e^2 = 0.64, sqrt(1-e^2) = 0.8.
    double var_w = 30. / 180. * M_PI; // 30 degrees in radian (sin w = 1/2, cos w = sqrt(3)/2).
    double var_large_e = 60. / 180. * M_PI; // 60 degrees in radian (sin E = sqrt(3)/2, cos E = 1/2).
    double var_x = 10.; // 10 light-seconds.

    // 3) Compute orbital delay and true anomaly.
    double alpha = var_x * 1./2.; // Equation 2.31, the first definition.
    double beta = 0.8 * var_x * std::sqrt(3.)/2.; // Equation 2.31, the second definition.
    double delay = alpha * (1./2. - 0.6) + (beta + par_GAMMA) * std::sqrt(3.)/2.; // Equation 2.30.

    // 4) Modify variables for additional turns of the binary system.
    var_large_e += 2. * M_PI * (*turn_itor);

    // 5) Compute back orbital parameters so as to reproduce values set in step 2.
    double elapsed = (var_large_e - var_e * std::sqrt(3.)/2.) / 2. / M_PI * var_PB; // Equation 2.27, where sigma = -T0*2pi/var_PB.
    double par_PB = var_PB - par_PBDOT * elapsed / 2.; // Equation 2.38.
    double par_A1 = var_x - par_XDOT * elapsed; // Equation 2.38.
    double par_ECC = var_e - par_ECCDOT * elapsed; // Equation 2.38.
    double var_wdot = par_OMDOT / 180.0 * M_PI / 365.25 / 86400.; // degrees-per-year to radian-per-second.
    double par_OM = var_w - var_wdot * elapsed; // Equation 2.38.
    par_OM *= 180. / M_PI; // radians to degrees.
    ev_time = t0 + ElapsedTime("TDB", Duration(elapsed, "Sec"));

    // 6) Create an ephemeris object and test its calcOrbitalDelay method.
    eph.reset(new BtModelEph("TDB", par_PB, par_PBDOT, par_A1, par_XDOT, par_ECC, par_ECCDOT, par_OM, par_OMDOT, t0, par_GAMMA));
    ElapsedTime delay_result = eph->calcOrbitalDelay(ev_time);
    ElapsedTime delay_expected("TDB", Duration(delay, "Sec"));
    Duration delay_tolerance(std::numeric_limits<double>::epsilon() * 86400. * 100., "Sec");
    if (!delay_result.getDuration().equivalentTo(delay_expected.getDuration(), delay_tolerance))
      err() << "BtModelEph::calcOrbitalDelay(" << ev_time << ") returned " << delay_result << " after " << *turn_itor <<
        " turn(s), not equivalent to " << delay_expected << " as expected." << std::endl;
  }

  // Test binary modulation and demodulation.
  // Binary parameters: (PB, PBDOT, A1, XDOT, ECC, ECCDOT, OM, OMDOT, T0, GAMMA)
  // = (27906.980897, -2.43e-12, 2.3417598, 0.0, 0.61713101, 0.0, 220.142729, 4.22662, 45888.1172487, 0.004295)
  AbsoluteTime abs_t0("TDB", Mjd(45888, .1172487));
  eph.reset(new BtModelEph("TDB", 27906.980897, -2.43e-12, 2.3417598, 0.0, 0.61713101, 0.0, 220.142729, 4.22662, abs_t0, 0.004295));

  // MJD's: { {original-MJD, modulated-MJD}, ... }
  Mjd mjd_test_values[][2] = {
    { Mjd(45988, .1172486599971307), Mjd(45988, .117282334454636) },
    { Mjd(45988, .1519708822161192), Mjd(45988, .152005747963997) },
    { Mjd(45988, .1866931044423836), Mjd(45988, .186723309730255) },
    { Mjd(45988, .2214153266613721), Mjd(45988, .221430372200973) },
    { Mjd(45988, .2561375488876365), Mjd(45988, .256125437773890) },
    { Mjd(45988, .2908597711066250), Mjd(45988, .290853252365181) },
    { Mjd(45988, .3255819933328894), Mjd(45988, .325587771318563) },
    { Mjd(45988, .3603042155518779), Mjd(45988, .360321331166131) },
    { Mjd(45988, .3950264377781423), Mjd(45988, .395052671940559) },
    { Mjd(45988, .4297486599971307), Mjd(45988, .429781129741675) },
    { Mjd(45988, .4644708822161192), Mjd(45988, .464505895402304) }
  };

  // Permitted difference is 100 ns.
  double delta = 100. * 1.e-9;
  ElapsedTime tolerance("TDB", Duration(delta, "Sec"));

  setPrecision(24);
  for (std::size_t ii = 0; ii != sizeof(mjd_test_values)/sizeof(Mjd[2]); ++ii) {
    AbsoluteTime tdb_mjd("TDB", mjd_test_values[ii][0]);
    AbsoluteTime expected_tdb_mjd("TDB", mjd_test_values[ii][1]);
    AbsoluteTime original_tdb_mjd(tdb_mjd);
    eph->modulateBinary(tdb_mjd);
    if (!tdb_mjd.equivalentTo(expected_tdb_mjd, tolerance)) {
      err() << "Binary modulation of " << original_tdb_mjd.represent("TDB", MjdFmt) << " by BtModelEph was computed to be " <<
        tdb_mjd.represent("TDB", MjdFmt) << ", not " << expected_tdb_mjd.represent("TDB", MjdFmt) << ", as expected." << std::endl;
    }
  }

  for (std::size_t ii = 0; ii != sizeof(mjd_test_values)/sizeof(Mjd[2]); ++ii) {
    AbsoluteTime tdb_mjd("TDB", mjd_test_values[ii][1]);
    AbsoluteTime expected_tdb_mjd("TDB", mjd_test_values[ii][0]);
    AbsoluteTime original_tdb_mjd(tdb_mjd);
    eph->demodulateBinary(tdb_mjd);
    if (!tdb_mjd.equivalentTo(expected_tdb_mjd, tolerance)) {
      err() << "Binary demodulation of " << original_tdb_mjd.represent("TDB", MjdFmt) << " by BtModelEph was computed to be " <<
        tdb_mjd.represent("TDB", MjdFmt) << ", not " << expected_tdb_mjd.represent("TDB", MjdFmt) << ", as expected." << std::endl;
    }
  }

  // Test the constructor that takes numerical arguments.
  const double rad_per_deg = M_PI / 180.;
  const double rad_year_per_deg_sec = rad_per_deg / (365.25 * 86400.);
  eph.reset(new BtModelEph("TDB", 12345.6789, 1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7, AbsoluteTime("TDB", 12345, 51840.), 8.8));
  TextEph t_eph;
  t_eph.append("PB",        "12345.6789",               "s");
  t_eph.append("PBDOT",     "1.1",                      "");
  t_eph.append("A1",        "2.2",                      "lt-s");
  t_eph.append("XDOT",      "3.3",                      "lt-s/s");
  t_eph.append("ECC",       "4.4",                      "");
  t_eph.append("ECCDOT",    "5.5",                      "s**(-1)");
  t_eph.append("OM",        6.6 * rad_per_deg,          "radians");
  t_eph.append("OMDOT",     7.7 * rad_year_per_deg_sec, "radians/s");
  t_eph.append("T0",        "12345.6 MJD (TDB)",        "");
  t_eph.append("GAMMA",     "8.8",                      "s");
  checkEphParameter(getMethod() + "_numeric", *eph, t_eph);

  // Test the constructor that takes a FITS record.
  std::string test_tpl("test_BtModelEph.tpl");
  remove(test_tpl.c_str());
  std::ofstream ofs(test_tpl.c_str());
  ofs << "\\include " << prependDataPath("PulsarDb_primary.tpl") << std::endl;
  ofs << "\\include " << prependDataPath("PulsarDb_orbital_bt.tpl") << std::endl;
  ofs.close();

  tip::IFileSvc::instance().createFile(getMethod() + ".fits", test_tpl);
  tip::TipFile tip_file = tip::IFileSvc::instance().openFile(getMethod() + ".fits");

  std::auto_ptr<tip::Table> table(tip_file.editTable("1"));
  tip::Header & header(table->getHeader());
  table->setNumRecords(1);
  tip::Table::Record record(table.get(), 0);
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
  eph.reset(new BtModelEph(record, header));
  checkEphParameter(getMethod() + "_fits", *eph, t_eph);
}

void PulsarDbTestApp::testEll1ModelEph() {
  setMethod("testEll1ModelDdEph");

  // Prepare variables for testing ephemeris computations.
  AbsoluteTime tasc("TDB", 51910, 123.456789);
  std::auto_ptr<Ell1ModelEph> eph(new Ell1ModelEph("TDB", 1000., .2, 0., 0., 0., 0., 0., 0., tasc, 0., 0.));
  AbsoluteTime ev_time("TDB", 51910, 223.456789);
  double epsilon = 1.e-8;

  // Test time system getter.
  std::string sys_name = eph->getSystem().getName();
  if ("TDB" != sys_name) {
    err() << "Ell1ModelEph::getSystem returned \"" << sys_name << "\", not \"TDB\"" << std::endl;
  }

  // Test orbital phase computations.
  double pb = 1000.;
  double pbdot = .2;
  double ecc = 3.45678912e-3;
  double eccdot = 4.56789123e-6;
  double omega = 123.456789 / 180. * M_PI;
  double omdot = 2.34567891e-3 / 180. * M_PI;

  // Try three cases for zero and non-zero eccentricity parameters.
  // Note: For consistency in parametrization, omega or omdot must be zero when ecc or eccdot is zero (0).
  //       See the comments in the implementation of Ell1ModelEph class for more explanations.
  for (int ii=0; ii<3; ii++) {
    if (ii > 0) ecc = omdot = 0.;
    if (ii > 1) eccdot = omega = 0.;

    double nb = 2. * M_PI / pb; // Eq. A4 in Lange, et al. (MNRAS 326, 274)
    double nbdot = - pbdot / pb * nb;
    double nbbar = nb + omdot - nbdot * omega / (nb + omdot); // Eq. A11 in Lange, et al. (MNRAS 326, 274)
    double elapsed_second = 100.;
    double phase_expected = (nbbar + nbdot / 2.0 * elapsed_second) * elapsed_second / (2. * M_PI);

    double eps1 = ecc * std::sin(omega); // Eq. A8 in Lange, et al. (MNRAS 326, 274)
    double eps1dot = eccdot * std::sin(omega) + ecc * std::cos(omega) * omdot; // Eq. A14 in Lange, et al. (MNRAS 326, 274)
    double eps2 = ecc * std::cos(omega); // Eq. A8 in Lange, et al. (MNRAS 326, 274)
    double eps2dot = eccdot * std::cos(omega) - ecc * std::sin(omega) * omdot; // Eq. A15 in Lange, et al. (MNRAS 326, 274)
    eph.reset(new Ell1ModelEph("TDB", pb, pbdot, 0., 0., eps1, eps1dot, eps2, eps2dot, tasc, 0., 0.));
    double phase = eph->calcOrbitalPhase(ev_time);
    if (std::fabs(phase/phase_expected - 1.) > epsilon) {
      err() << "Ell1ModelEph::calcOrbitalPhase produced phase == " << phase << ", not " << phase_expected << " (ECC=" <<
        ecc << ", ECCDOT=" << eccdot << ")" << std::endl;
    }

    // Test orbital phase computations, with a non-zero global phase offset.
    phase = eph->calcOrbitalPhase(ev_time, 0.1234);
    phase_expected += 0.1234;
    if (std::fabs(phase/phase_expected - 1.) > epsilon) {
      err() << "Ell1ModelEph::calcOrbitalPhase produced phase == " << phase << ", not " << phase_expected << " (ECC=" <<
        ecc << ", ECCDOT=" << eccdot << ")" << std::endl;
    }
  }

  // Test orbital delay computations.
  // Note: Below the equations in Lange, et al. (MNRAS 326, 274, 2001) are computed in reverse.
  int turn_array[] = {0, 1, -1, 2, -2, 5, -5, 20, -20, 100, -100};
  std::list<int> turn_list(turn_array, turn_array + sizeof(turn_array)/sizeof(int));
  for (std::list<int>::const_iterator turn_itor = turn_list.begin(); turn_itor != turn_list.end(); ++turn_itor) {
    // 1) Set constants.
    double par_PB = 10. * 86400.; // 10 days in seconds.
    double max_elapsed = par_PB * (std::fabs(*turn_itor) + 1.);
    double par_XDOT = -0.44 / max_elapsed;
    double par_EPS1DOT = 5.5e-6 / max_elapsed;
    double par_EPS2DOT = 6.6e-6 / max_elapsed;
    double par_SHAPIRO_R = 0.0022; // in micro-seconds.
    double par_SHAPIRO_S = 0.0033;

    // 2) Set values to variables that appear in equations A6 and A16.
    double var_eta = 7.7e-3;
    double var_kappa = 8.8e-3;
    double var_large_phi = 60. / 180. * M_PI; // 60 degrees in radian (sin phi = sin 2*phi = sqrt(3)/2, cos 2*phi = -1/2).
    double var_x = 10.; // 10 light-seconds.
    double var_r = par_SHAPIRO_R * 1.e-6; // micro-seconds to seconds.

    // 3) Compute orbital delay.
    double delay = var_x * (std::sqrt(3.)/2. + var_kappa / 2. * std::sqrt(3.)/2. + var_eta / 4.); // Equation A6.
    delay += -2. * var_r * std::log(1. - par_SHAPIRO_S*std::sqrt(3.)/2.); // Equation A16.

    // 4) Modify variables for additional turns of the binary system.
    var_large_phi += 2. * M_PI * (*turn_itor);

    // 5) Compute back orbital parameters so as to reproduce values set in step 2.
    double var_nb = 2. * M_PI / par_PB; // Equation A4.
    double elapsed = 1.001 * var_large_phi / var_nb; // First term of Equation A9.
    double par_EPS1 = var_eta - par_EPS1DOT * elapsed; // Equation A10.
    double par_EPS2 = var_kappa - par_EPS2DOT * elapsed; // Equation A10.
    double var_e0_squared = par_EPS1 * par_EPS1 + par_EPS2 * par_EPS2; // Derived from equation A8.
    double var_w0 = std::atan2(par_EPS1, par_EPS2); // Derived from equation A8.
    if (var_w0 < 0.) var_w0 += 2. * M_PI;
    double var_wdot = (par_EPS1DOT * par_EPS2 - par_EPS2DOT * par_EPS1) / var_e0_squared; // Derived from equation A8.
    double var_nbdot = (var_large_phi - (var_nb + var_wdot) * elapsed)
      / (elapsed * (elapsed / 2. - var_w0 / (var_nb + var_wdot))); // Derived from equation A9.
    double par_PBDOT = - var_nbdot / var_nb * par_PB; // Derived from equation A4.
    double par_A1 = var_x - par_XDOT * elapsed; // Equation A10.
    ev_time = tasc + ElapsedTime("TDB", Duration(elapsed, "Sec"));

    // 6) Create an ephemeris object and test its calcOrbitalDelay method.
    eph.reset(new Ell1ModelEph("TDB", par_PB, par_PBDOT, par_A1, par_XDOT, par_EPS1, par_EPS1DOT, par_EPS2, par_EPS2DOT, tasc,
      par_SHAPIRO_R, par_SHAPIRO_S));
    ElapsedTime delay_result = eph->calcOrbitalDelay(ev_time);
    ElapsedTime delay_expected("TDB", Duration(delay, "Sec"));
    Duration delay_tolerance(std::numeric_limits<double>::epsilon() * 86400. * 100., "Sec");
    if (!delay_result.getDuration().equivalentTo(delay_expected.getDuration(), delay_tolerance))
      err() << "Ell1ModelEph::calcOrbitalDelay(" << ev_time << ") returned " << delay_result << " after " << *turn_itor <<
        " turn(s), not equivalent to " << delay_expected << " as expected." << std::endl;
  }

  // Test binary modulation and demodulation.
  // Binary parameters: (PB, PBDOT, A1, XDOT, ECC, ECCDOT, OM, OMDOT, T0, GAMMA)
  // = (27906.980897, -2.43e-12, 2.3417598, 0.0, 0.61713101e-3, 0.0, 220.142729, 4.22662, 45888.1172487, 0.004295)
  // and all other parameters set to 0.0 (zero).
  // Note: The following conversions were made by hand.
  //   EPS1 = ECC * sin(OM) = -3.9786059746415216e-4
  //   EPS1DOT = ECCDOT * sin(OM) + ECC * cos(OM) * OMDOT = -1.10277737e-12
  //   EPS2 = ECC * cos(OM) = -4.7176013872421223e-4
  //   EPS2DOT = ECCDOT * cos(OM) - ECC * sin(OM) * OMDOT = -9.3003123e-13
  //   TASC = T0 - OM / (2*pi/PB + OMDOT) = T0 - 17065.15319093980494836272 seconds
  AbsoluteTime abs_t0("TDB", Mjd(45888, .1172487));
  AbsoluteTime abs_tasc = abs_t0 - ElapsedTime("TDB", Duration(17065.15319093980494836272, "Sec"));
  eph.reset(new Ell1ModelEph("TDB", 27906.980897, -2.43e-12, 2.3417598, 0.0, -3.9786059746415216e-4, -1.10277737e-12,
    -4.7176013872421223e-4, -9.3003123e-13, abs_tasc, 0.0, 0.0));

  // MJD's: { {original-MJD, modulated-MJD}, ... }
  Mjd mjd_test_values[][2] = {
    { Mjd(45988, .1172486599971307), Mjd(45988, .117274986866046) },
    { Mjd(45988, .1519708822161192), Mjd(45988, .151995444036652) },
    { Mjd(45988, .1866931044423836), Mjd(45988, .186705113506475) },
    { Mjd(45988, .2214153266613721), Mjd(45988, .221409499721483) },
    { Mjd(45988, .2561375488876365), Mjd(45988, .256116442226323) },
    { Mjd(45988, .2908597711066250), Mjd(45988, .290832661977126) },
    { Mjd(45988, .3255819933328894), Mjd(45988, .325560792708099) },
    { Mjd(45988, .3603042155518779), Mjd(45988, .360298227751942) },
    { Mjd(45988, .3950264377781423), Mjd(45988, .395038283532021) },
    { Mjd(45988, .4297486599971307), Mjd(45988, .429773139267359) },
    { Mjd(45988, .4644708822161192), Mjd(45988, .464497254754357) }
  };

  // Permitted difference is 100 ns.
  double delta = 100. * 1.e-9;
  ElapsedTime tolerance("TDB", Duration(delta, "Sec"));

  setPrecision(24);
  for (std::size_t ii = 0; ii != sizeof(mjd_test_values)/sizeof(Mjd[2]); ++ii) {
    AbsoluteTime tdb_mjd("TDB", mjd_test_values[ii][0]);
    AbsoluteTime expected_tdb_mjd("TDB", mjd_test_values[ii][1]);
    AbsoluteTime original_tdb_mjd(tdb_mjd);
    eph->modulateBinary(tdb_mjd);
    if (!tdb_mjd.equivalentTo(expected_tdb_mjd, tolerance)) {
      err() << "Binary modulation of " << original_tdb_mjd.represent("TDB", MjdFmt) << " by Ell1ModelEph was computed to be " <<
        tdb_mjd.represent("TDB", MjdFmt) << ", not " << expected_tdb_mjd.represent("TDB", MjdFmt) << ", as expected." << std::endl;
    }
  }

  for (std::size_t ii = 0; ii != sizeof(mjd_test_values)/sizeof(Mjd[2]); ++ii) {
    AbsoluteTime tdb_mjd("TDB", mjd_test_values[ii][1]);
    AbsoluteTime expected_tdb_mjd("TDB", mjd_test_values[ii][0]);
    AbsoluteTime original_tdb_mjd(tdb_mjd);
    eph->demodulateBinary(tdb_mjd);
    if (!tdb_mjd.equivalentTo(expected_tdb_mjd, tolerance)) {
      err() << "Binary demodulation of " << original_tdb_mjd.represent("TDB", MjdFmt) << " by Ell1ModelEph was computed to be " <<
        tdb_mjd.represent("TDB", MjdFmt) << ", not " << expected_tdb_mjd.represent("TDB", MjdFmt) << ", as expected." << std::endl;
    }
  }

  // Test the constructor that takes numerical arguments.
  eph.reset(new Ell1ModelEph("TDB", 12345.6789, 1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7, AbsoluteTime("TDB", 12345, 51840.), 8.8, 9.9));
  TextEph t_eph;
  t_eph.append("PB",        "12345.6789",        "s");
  t_eph.append("PBDOT",     "1.1",               "");
  t_eph.append("A1",        "2.2",               "lt-s");
  t_eph.append("XDOT",      "3.3",               "lt-s/s");
  t_eph.append("EPS1",       "4.4",              "");
  t_eph.append("EPS1DOT",   "5.5",               "s**(-1)");
  t_eph.append("EPS2",      "6.6",               "");
  t_eph.append("EPS2DOT",   "7.7",               "s**(-1)");
  t_eph.append("TASC",      "12345.6 MJD (TDB)", "");
  t_eph.append("SHAPIRO_R", "8.8e-06",           "s");
  t_eph.append("SHAPIRO_S", "9.9",               "");
  checkEphParameter(getMethod() + "_numeric", *eph, t_eph);

  // Test the constructor that takes a FITS record.
  std::string test_tpl("test_Ell1ModelEph.tpl");
  remove(test_tpl.c_str());
  std::ofstream ofs(test_tpl.c_str());
  ofs << "\\include " << prependDataPath("PulsarDb_primary.tpl") << std::endl;
  ofs << "\\include " << prependDataPath("PulsarDb_orbital_ell1.tpl") << std::endl;
  ofs.close();

  tip::IFileSvc::instance().createFile(getMethod() + ".fits", test_tpl);
  tip::TipFile tip_file = tip::IFileSvc::instance().openFile(getMethod() + ".fits");

  std::auto_ptr<tip::Table> table(tip_file.editTable("1"));
  tip::Header & header(table->getHeader());
  table->setNumRecords(1);
  tip::Table::Record record(table.get(), 0);
  record["PB"].set(12345.6789);
  record["PBDOT"].set(1.1);
  record["A1"].set(2.2);
  record["XDOT"].set(3.3);
  record["EPS1"].set(4.4);
  record["EPS1DOT"].set(5.5);
  record["EPS2"].set(6.6);
  record["EPS2DOT"].set(7.7);
  record["TASC"].set(12345.6);
  record["SHAPIRO_R"].set(8.8);
  record["SHAPIRO_S"].set(9.9);
  eph.reset(new Ell1ModelEph(record, header));
  checkEphParameter(getMethod() + "_fits", *eph, t_eph);
}

void PulsarDbTestApp::testMssModelEph() {
  setMethod("testMssModelDdEph");

  // Prepare variables for testing ephemeris computations.
  AbsoluteTime t0("TDB", 51910, 123.456789);
  std::auto_ptr<MssModelEph> eph(new MssModelEph("TDB", 1000., .2, 0., 0., 0., 0., 0., 0., 0., 0., t0, 0., 0., 0., 0., 0., 0., 0.));
  AbsoluteTime ev_time("TDB", 51910, 223.456789);
  double epsilon = 1.e-8;

  // Test time system getter.
  std::string sys_name = eph->getSystem().getName();
  if ("TDB" != sys_name) {
    err() << "MssModelEph::getSystem returned \"" << sys_name << "\", not \"TDB\"" << std::endl;
  }

  // Test orbital phase computations.
  double phase = eph->calcOrbitalPhase(ev_time);
  double phase_expected = .099;
  if (std::fabs(phase/phase_expected - 1.) > epsilon) {
    err() << "MssModelEph::calcOrbitalPhase produced phase == " << phase << ", not " << phase_expected << std::endl;
  }

  // Test orbital phase computations, with a non-zero global phase offset.
  phase = eph->calcOrbitalPhase(ev_time, 0.1234);
  phase_expected += 0.1234;
  if (std::fabs(phase/phase_expected - 1.) > epsilon) {
    err() << "MssModelEph::calcOrbitalPhase produced phase == " << phase << ", not " << phase_expected << std::endl;
  }

  // Test orbital delay computations.
  // Note: Below the equations in Wex (MNRAS 298, 67-77, 1998) and Taylor et al. (ApJ, 345, 434-450, 1989) are computed in reverse.
  int turn_array[] = {0, 1, -1, 2, -2, 5, -5, 20, -20, 100, -100};
  std::list<int> turn_list(turn_array, turn_array + sizeof(turn_array)/sizeof(int));
  for (std::list<int>::const_iterator turn_itor = turn_list.begin(); turn_itor != turn_list.end(); ++turn_itor) {
    // 1) Set constants.
    double par_PB = 10. * 86400.; // 10 days in seconds.
    double max_elapsed = par_PB * (std::fabs(*turn_itor) + 1.);
    double par_XDOT = -0.44 / max_elapsed;
    double par_X2DOT = 0.3355 / max_elapsed / max_elapsed;
    double par_ECCDOT = -0.00055 / max_elapsed;
    double par_OMDOT = 6.6 / max_elapsed;
    double par_OM2DOT = 5.577 / max_elapsed / max_elapsed;
    double par_DELTA_R = 0.05;
    double par_DELTA_THETA = 1./3.;
    double par_GAMMA = 0.0011;
    double par_SHAPIRO_R = 0.0022; // in micro-seconds.
    double par_SHAPIRO_S = 0.0033;
    double par_ABERRATION_A = 0.00000432198765;
    double par_ABERRATION_B = 0.00000321987654;

    // 2) Set values to variables that appear in equations 8, 9, and 10 in Taylor et al. (1989).
    double var_e = 0.6; // 1-e^2 = 0.64, sqrt(1-e^2) = 0.8, e*(1+delta_r) = 0.63, sqrt(1-e^2*(1+delta_theta)^2) = 0.6.
    double var_w = 30. / 180. * M_PI; // 30 degrees in radian (sin w = 1/2, cos w = sqrt(3)/2).
    double var_u = 60. / 180. * M_PI; // 60 degrees in radian (sin u = sqrt(3)/2, cos u = 1/2, tan u/2 = sqrt(3)/3).
    double var_x = 10.; // 10 light-seconds.
    double var_r = par_SHAPIRO_R * 1.e-6; // micro-seconds to seconds.

    // 3) Compute orbital delay and true anomaly.
    double delay = var_x * 1./2. * (1./2. - 0.63) + var_x * 0.6 * 3./4.; // Equation 8 in Taylor et al. (1989).
    delay += par_GAMMA * std::sqrt(3.)/2.; // Equation 9 in Taylor et al. (1989).
    delay += -2. * var_r * std::log(1. - 0.6*1./2. - par_SHAPIRO_S*(1./2.*(1./2. - 0.6) + 0.8*3./4.)); // Equation 10 in Taylor et al. (1989).
    double true_anomaly = 2. * std::atan(2. * std::sqrt(3.)/3.); // Equation 13 in Taylor et al. (1989).
    delay += par_ABERRATION_A * (std::sin(var_w + true_anomaly) + 0.6 * 1./2.); // Equation 11 in Taylor et al. (1989).
    delay += par_ABERRATION_B * (std::cos(var_w + true_anomaly) + 0.6 * std::sqrt(3.)/2.); // Equation 11 in Taylor et al. (1989).

    // 4) Modify variables for additional turns of the binary system.
    var_u += 2. * M_PI * (*turn_itor);
    true_anomaly += 2. * M_PI * (*turn_itor);

    // 5) Compute back orbital parameters so as to reproduce values set in step 2.
    double elapsed = var_u / 2. / M_PI * par_PB; // Equation 12.
    double par_PBDOT = 0.6 * std::sqrt(3.)/2. / M_PI * (par_PB/elapsed) * (par_PB/elapsed); // Equation 12.
    double var_xi = par_XDOT / 2. / M_PI * par_PB; // Equation 64 in Wex (1998).
    double par_A1 = var_x - var_xi * true_anomaly - par_X2DOT / 2. * elapsed * elapsed; // Equations 64 and 73 in Wex (1998).
    double par_ECC = var_e - par_ECCDOT * elapsed;
    double var_wdot = par_OMDOT / 180.0 * M_PI / 365.25 / 86400.; // degrees-per-year to radian-per-second.
    double var_k = var_wdot / 2. / M_PI * par_PB; // Equation 65 in Wex (1998).
    double var_w2dot = par_OM2DOT / 180.0 * M_PI / 365.25 / 86400. / 365.25 / 86400.; // degrees-per-year2 to radian-per-second2.
    double par_OM = var_w - var_k * true_anomaly - var_w2dot / 2. * elapsed * elapsed; // Equations 65 and 74 in Wex (1998).
    par_OM *= 180. / M_PI; // radians to degrees.
    ev_time = t0 + ElapsedTime("TDB", Duration(elapsed, "Sec"));

    // 6) Create an ephemeris object and test its calcOrbitalDelay method.
    eph.reset(new MssModelEph("TDB", par_PB, par_PBDOT, par_A1, par_XDOT, par_X2DOT, par_ECC, par_ECCDOT, par_OM, par_OMDOT,
      par_OM2DOT, t0, par_DELTA_R, par_DELTA_THETA, par_GAMMA, par_SHAPIRO_R, par_SHAPIRO_S, par_ABERRATION_A, par_ABERRATION_B));
    ElapsedTime delay_result = eph->calcOrbitalDelay(ev_time);
    ElapsedTime delay_expected("TDB", Duration(delay, "Sec"));
    Duration delay_tolerance(std::numeric_limits<double>::epsilon() * 86400. * 100., "Sec");
    if (!delay_result.getDuration().equivalentTo(delay_expected.getDuration(), delay_tolerance))
      err() << "MssModelEph::calcOrbitalDelay(" << ev_time << ") returned " << delay_result << " after " << *turn_itor <<
        " turn(s), not equivalent to " << delay_expected << " as expected." << std::endl;
  }

  // Test binary modulation and demodulation.
  // Binary parameters: (PB, PBDOT, A1, XDOT, ECC, ECCDOT, OM, OMDOT, T0, GAMMA)
  // = (27906.980897, -2.43e-12, 2.3417598, 0.0, 0.61713101e-3, 0.0, 220.142729, 4.22662, 45888.1172487, 0.004295)
  // and all other parameters set to 0.0 (zero).
  AbsoluteTime abs_t0("TDB", Mjd(45888, .1172487));
  eph.reset(new MssModelEph("TDB", 27906.980897, -2.43e-12, 2.3417598, 0.0, 0.0, 0.61713101, 0.0, 220.142729, 4.22662, 0.0,
    abs_t0, 0.0, 0.0, 0.004295, 0.0, 0.0, 0.0, 0.0));

  // MJD's: { {original-MJD, modulated-MJD}, ... }
  Mjd mjd_test_values[][2] = {
    { Mjd(45988, .1172486599971307), Mjd(45988, .117282334337111) },
    { Mjd(45988, .1519708822161192), Mjd(45988, .152005747822377) },
    { Mjd(45988, .1866931044423836), Mjd(45988, .186723309695442) },
    { Mjd(45988, .2214153266613721), Mjd(45988, .221430372328996) },
    { Mjd(45988, .2561375488876365), Mjd(45988, .25612543777892 ) },
    { Mjd(45988, .2908597711066250), Mjd(45988, .2908532526978  ) },
    { Mjd(45988, .3255819933328894), Mjd(45988, .32558777169579 ) },
    { Mjd(45988, .3603042155518779), Mjd(45988, .360321331412386) },
    { Mjd(45988, .3950264377781423), Mjd(45988, .39505267200415 ) },
    { Mjd(45988, .4297486599971307), Mjd(45988, .429781129654144) },
    { Mjd(45988, .4644708822161192), Mjd(45988, .464505895254603) }
  };

  // Permitted difference is 100 ns.
  double delta = 100. * 1.e-9;
  ElapsedTime tolerance("TDB", Duration(delta, "Sec"));

  setPrecision(24);
  for (std::size_t ii = 0; ii != sizeof(mjd_test_values)/sizeof(Mjd[2]); ++ii) {
    AbsoluteTime tdb_mjd("TDB", mjd_test_values[ii][0]);
    AbsoluteTime expected_tdb_mjd("TDB", mjd_test_values[ii][1]);
    AbsoluteTime original_tdb_mjd(tdb_mjd);
    eph->modulateBinary(tdb_mjd);
    if (!tdb_mjd.equivalentTo(expected_tdb_mjd, tolerance)) {
      err() << "Binary modulation of " << original_tdb_mjd.represent("TDB", MjdFmt) << " by MssModelEph was computed to be " <<
        tdb_mjd.represent("TDB", MjdFmt) << ", not " << expected_tdb_mjd.represent("TDB", MjdFmt) << ", as expected." << std::endl;
    }
  }

  for (std::size_t ii = 0; ii != sizeof(mjd_test_values)/sizeof(Mjd[2]); ++ii) {
    AbsoluteTime tdb_mjd("TDB", mjd_test_values[ii][1]);
    AbsoluteTime expected_tdb_mjd("TDB", mjd_test_values[ii][0]);
    AbsoluteTime original_tdb_mjd(tdb_mjd);
    eph->demodulateBinary(tdb_mjd);
    if (!tdb_mjd.equivalentTo(expected_tdb_mjd, tolerance)) {
      err() << "Binary demodulation of " << original_tdb_mjd.represent("TDB", MjdFmt) << " by MssModelEph was computed to be " <<
        tdb_mjd.represent("TDB", MjdFmt) << ", not " << expected_tdb_mjd.represent("TDB", MjdFmt) << ", as expected." << std::endl;
    }
  }

  // Test the constructor that takes numerical arguments.
  const double rad_per_deg = M_PI / 180.;
  const double rad_year_per_deg_sec = rad_per_deg / (365.25 * 86400.);
  const double rad_year2_per_deg_sec2 = rad_year_per_deg_sec / (365.25 * 86400.);
  eph.reset(new MssModelEph("TDB", 12345.6789, 1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7, 8.8, 9.9, AbsoluteTime("TDB", 12345, 51840.),
    11.1, 22.2, 33.3, 44.4, 55.5, 66.6, 77.7));
  TextEph t_eph;
  t_eph.append("PB",           "12345.6789",                 "s");
  t_eph.append("PBDOT",        "1.1",                        "");
  t_eph.append("A1",           "2.2",                        "lt-s");
  t_eph.append("XDOT",         "3.3",                        "lt-s/s");
  t_eph.append("X2DOT",        "4.4",                        "lt-s/s**2");
  t_eph.append("ECC",          "5.5",                        "");
  t_eph.append("ECCDOT",       "6.6",                        "s**(-1)");
  t_eph.append("OM",           7.7 * rad_per_deg,            "radians");
  t_eph.append("OMDOT",        8.8 * rad_year_per_deg_sec,   "radians/s");
  t_eph.append("OM2DOT",       9.9 * rad_year2_per_deg_sec2, "radians/s**2");
  t_eph.append("T0",           "12345.6 MJD (TDB)",          "");
  t_eph.append("DELTA_R",      "11.1",                       "");
  t_eph.append("DELTA_THETA",  "22.2",                       "");
  t_eph.append("GAMMA",        "33.3",                       "s");
  t_eph.append("SHAPIRO_R",    "44.4e-06",                   "s");
  t_eph.append("SHAPIRO_S",    "55.5",                       "");
  t_eph.append("ABERRATION_A", "66.6",                       "s");
  t_eph.append("ABERRATION_B", "77.7",                       "s");
  checkEphParameter(getMethod() + "_numeric", *eph, t_eph);

  // Test the constructor that takes a FITS record.
  std::string test_tpl("test_MssModelEph.tpl");
  remove(test_tpl.c_str());
  std::ofstream ofs(test_tpl.c_str());
  ofs << "\\include " << prependDataPath("PulsarDb_primary.tpl") << std::endl;
  ofs << "\\include " << prependDataPath("PulsarDb_orbital_mss.tpl") << std::endl;
  ofs.close();

  tip::IFileSvc::instance().createFile(getMethod() + ".fits", test_tpl);
  tip::TipFile tip_file = tip::IFileSvc::instance().openFile(getMethod() + ".fits");

  std::auto_ptr<tip::Table> table(tip_file.editTable("1"));
  tip::Header & header(table->getHeader());
  table->setNumRecords(1);
  tip::Table::Record record(table.get(), 0);
  record["PB"].set(12345.6789);
  record["PBDOT"].set(1.1);
  record["A1"].set(2.2);
  record["XDOT"].set(3.3);
  record["X2DOT"].set(4.4);
  record["ECC"].set(5.5);
  record["ECCDOT"].set(6.6);
  record["OM"].set(7.7);
  record["OMDOT"].set(8.8);
  record["OM2DOT"].set(9.9);
  record["T0"].set(12345.6);
  record["DELTA_R"].set(11.1);
  record["DELTA_THETA"].set(22.2);
  record["GAMMA"].set(33.3);
  record["SHAPIRO_R"].set(44.4);
  record["SHAPIRO_S"].set(55.5);
  record["ABERRATION_A"].set(66.6);
  record["ABERRATION_B"].set(77.7);
  eph.reset(new MssModelEph(record, header));
  checkEphParameter(getMethod() + "_fits", *eph, t_eph);
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
      << ev_time1 << ", not " << correct_time << std::endl;
  }

  // Test PdotCanceler created from a PulsarEph object.
  AbsoluteTime since("TT", 51910, 0.);
  AbsoluteTime until("TT", 51910, 1.);
  FrequencyEph f_eph("TDB", since, until, origin, 22., 45., .11, f0, f1, f2);
  PdotCanceler canceler2(origin, f_eph, 2);

  canceler2.cancelPdot(ev_time2);
  if (!correct_time.equivalentTo(ev_time2, tolerance)) {
    err() << "After constructed from PulsarEph, PdotCanceler::cancelPdot produced pdot-corrected time == "
      << ev_time2 << ", not " << correct_time << std::endl;
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
  database.registerPulsarEph<PeriodEph>("PER");
  database.registerPulsarEph<HighPrecisionEph>("HP");
  database.registerOrbitalEph<SimpleDdEph>("DD");
  database.registerOrbitalEph<BtModelEph>("BT");
  database.registerOrbitalEph<Ell1ModelEph>("ELL1");
  database.registerOrbitalEph<MssModelEph>("MSS");

  std::string pulsar_name = "PSR J0139+5814";

  // Select pulsar we want to use.
  database.filterName(pulsar_name);

  // Confirm that the correct number of ephemerides were found.
  int num_eph = database.getNumEph();
  if (10 != num_eph)
    err() << "there are " << num_eph << " ephemerides for " << pulsar_name << ", not 10" << std::endl;

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
  database2.registerPulsarEph<PeriodEph>("PER");
  database2.registerPulsarEph<HighPrecisionEph>("HP");
  database2.registerOrbitalEph<SimpleDdEph>("DD");
  database2.registerOrbitalEph<BtModelEph>("BT");
  database2.registerOrbitalEph<Ell1ModelEph>("ELL1");
  database2.registerOrbitalEph<MssModelEph>("MSS");

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
    err() << "for time " << pick_time << ", chooser had trouble choosing orbital eph: " << std::endl <<
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
  database.registerPulsarEph<PeriodEph>("PER");
  database.registerPulsarEph<HighPrecisionEph>("HP");
  database.registerOrbitalEph<SimpleDdEph>("DD");
  database.registerOrbitalEph<BtModelEph>("BT");
  database.registerOrbitalEph<Ell1ModelEph>("ELL1");
  database.registerOrbitalEph<MssModelEph>("MSS");

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

  // Test calcPosition, by comparing it with FrequencyEph::calcPosition.
  SourcePosition expected_src_pos = eph.calcPosition(expected_gtdb);
  SourcePosition src_pos = computer.calcPosition(expected_gtdb);
  std::string coord_name("XYZ");
  for (std::size_t ii = 0; ii < 3; ++ii) {
    double expected_coord = expected_src_pos.getDirection()[ii];
    double coord = src_pos.getDirection()[ii];
    if (expected_coord != coord) {
      err() << "EphComputer::calcPosition returned " << coord_name[ii] << "=" << coord << ", not " <<
        expected_coord << "." << std::endl;
    }
  }
  if (src_pos.hasDistance() != false) {
    err() << "PeriodEph::calcPosition returned a source position with its distance known" << std::endl;
  }

  // Test binary modulation/demodulation.
  // Get new independent access to database, to keep independent from the tests above.
  PulsarDb database2(m_tpl_file);
  database2.load(m_in_file);

  // Register PulsarEph and OrbitalEph subclasses for various ephemeris models.
  database2.registerPulsarEph<FrequencyEph>("FREQ");
  database2.registerPulsarEph<PeriodEph>("PER");
  database2.registerPulsarEph<HighPrecisionEph>("HP");
  database2.registerOrbitalEph<SimpleDdEph>("DD");
  database2.registerOrbitalEph<BtModelEph>("BT");
  database2.registerOrbitalEph<Ell1ModelEph>("ELL1");
  database2.registerOrbitalEph<MssModelEph>("MSS");

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
  database3.registerPulsarEph<PeriodEph>("PER");
  database3.registerPulsarEph<HighPrecisionEph>("HP");
  database3.registerOrbitalEph<SimpleDdEph>("DD");
  database3.registerOrbitalEph<BtModelEph>("BT");
  database3.registerOrbitalEph<Ell1ModelEph>("ELL1");
  database3.registerOrbitalEph<MssModelEph>("MSS");

  // Load everything in this database at a time.
  computer2.load(database3);

  if (1155 != computer2.getNumPulsarEph())
    err() << "After loading all data from pulsar database, there were " << computer2.getNumPulsarEph() <<
      " spin pulsar ephemeri(de)s, not 1155 as expected." << std::endl;

  if (25 != computer2.getNumOrbitalEph())
    err() << "After loading all data from pulsar database, there were " << computer2.getNumOrbitalEph() <<
      " orbital ephemeri(de)s, not 25 as expected." << std::endl;

  if (26 != computer2.getNumEphRemark())
    err() << "After loading all data from pulsar database, there were " << computer2.getNumEphRemark() <<
      " ephemeris remark(s), not 26 as expected." << std::endl;

  // Create a new ephemeris computer for tests of type-specific loaders.
  EphComputer computer3;

  // Load just the spin parameters from this database.
  computer3.loadPulsarEph(database3);

  if (1155 != computer3.getNumPulsarEph())
    err() << "After loading spin pulsar ephemerides, there were " << computer3.getNumPulsarEph() <<
      " spin pulsar ephemeri(de)s, not 1155 as expected." << std::endl;

  if (0 != computer3.getNumOrbitalEph())
    err() << "After loading spin pulsar ephemerides, there were " << computer3.getNumOrbitalEph() <<
      " orbital ephemeri(de)s, not 0 as expected." << std::endl;

  if (14 != computer3.getNumEphRemark())
    err() << "After loading spin pulsar ephemerides, there were " << computer3.getNumEphRemark() <<
      " ephemeris remark(s), not 14 as expected." << std::endl;

  // Load just the orbital parameters from this database.
  computer3.loadOrbitalEph(database3);

  if (1155 != computer3.getNumPulsarEph())
    err() << "After loading orbital ephemerides, there were " << computer3.getNumPulsarEph() <<
      " spin pulsar ephemeri(de)s, not 1155 as expected." << std::endl;

  if (25 != computer3.getNumOrbitalEph())
    err() << "After loading orbital ephemerides, there were " << computer3.getNumOrbitalEph() <<
      " orbital ephemeri(de)s, not 25 as expected." << std::endl;

  if (14 != computer3.getNumEphRemark())
    err() << "After loading orbital ephemerides, there were " << computer3.getNumEphRemark() <<
      " ephemeris remark(s), not 14 as expected." << std::endl;

  // Load just the ephemeris remarks from this database.
  computer3.loadEphRemark(database3);

  if (1155 != computer3.getNumPulsarEph())
    err() << "After loading ephemeris remarks, there were " << computer3.getNumPulsarEph() <<
      " spin pulsar ephemeri(de)s, not 1155 as expected." << std::endl;

  if (25 != computer3.getNumOrbitalEph())
    err() << "After loading ephemeris remarks, there were " << computer3.getNumOrbitalEph() <<
      " orbital ephemeri(de)s, not 25 as expected." << std::endl;

  if (26 != computer3.getNumEphRemark())
    err() << "After loading ephemeris remarks, there were " << computer3.getNumEphRemark() <<
      " ephemeris remark(s), not 26 as expected." << std::endl;

  // Prepare variables for loading single ephemeris entry.
  AbsoluteTime since("TDB", 51910, 100.);
  AbsoluteTime until("TDB", 51910, 300.);
  AbsoluteTime epoch("TDB", 51910, 200.);

  // Load just one set of spin parameters.
  computer3.loadPulsarEph(FrequencyEph("TDB", since, until, epoch, 0., 0., 0., 1., 0., 0.));

  if (1156 != computer3.getNumPulsarEph())
    err() << "After loading one spin pulsar ephemeris, there were " << computer3.getNumPulsarEph() <<
      " spin pulsar ephemeri(de)s, not 1156 as expected." << std::endl;

  if (25 != computer3.getNumOrbitalEph())
    err() << "After loading one spin pulsar ephemeris, there were " << computer3.getNumOrbitalEph() <<
      " orbital ephemeri(de)s, not 25 as expected." << std::endl;

  if (26 != computer3.getNumEphRemark())
    err() << "After loading one spin pulsar ephemeris, there were " << computer3.getNumEphRemark() <<
      " ephemeris remark(s), not 26 as expected." << std::endl;

  // Load just one set of orbital parameters.
  computer3.loadOrbitalEph(SimpleDdEph("TDB", 1000., .2, 0., 0., 0., 0., 0., 0., epoch, 0., 0., 0.));

  if (1156 != computer3.getNumPulsarEph())
    err() << "After loading one orbital ephemeris, there were " << computer3.getNumPulsarEph() <<
      " spin pulsar ephemeri(de)s, not 1156 as expected." << std::endl;

  if (26 != computer3.getNumOrbitalEph())
    err() << "After loading one orbital ephemeris, there were " << computer3.getNumOrbitalEph() <<
      " orbital ephemeri(de)s, not 26 as expected." << std::endl;

  if (26 != computer3.getNumEphRemark())
    err() << "After loading one orbital ephemeris, there were " << computer3.getNumEphRemark() <<
      " ephemeris remark(s), not 26 as expected." << std::endl;

  // Load just one ephemeris remark.
  computer3.loadEphRemark(EphStatus(since, until, Remarked, "This is a remark."));

  if (1156 != computer3.getNumPulsarEph())
    err() << "After loading one ephemeris remark, there were " << computer3.getNumPulsarEph() <<
      " spin pulsar ephemeri(de)s, not 1156 as expected." << std::endl;

  if (26 != computer3.getNumOrbitalEph())
    err() << "After loading one ephemeris remark, there were " << computer3.getNumOrbitalEph() <<
      " orbital ephemeri(de)s, not 26 as expected." << std::endl;

  if (27 != computer3.getNumEphRemark())
    err() << "After loading one ephemeris remark, there were " << computer3.getNumEphRemark() <<
      " ephemeris remark(s), not 27 as expected." << std::endl;

  // Load non-remark ephemeris status, which should not be loaded.
  computer3.loadEphRemark(EphStatus(since, until, Unavailable, "No data"));
  computer3.loadEphRemark(EphStatus(since, until, Extrapolated, "Ephemeris gap"));

  if (1156 != computer3.getNumPulsarEph())
    err() << "After loading non-remark ephemeris status, there were " << computer3.getNumPulsarEph() <<
      " spin pulsar ephemeri(de)s, not 1156 as expected." << std::endl;

  if (26 != computer3.getNumOrbitalEph())
    err() << "After loading non-remark ephemeris status, there were " << computer3.getNumOrbitalEph() <<
      " orbital ephemeri(de)s, not 26 as expected." << std::endl;

  if (27 != computer3.getNumEphRemark())
    err() << "After loading non-remark ephemeris status, there were " << computer3.getNumEphRemark() <<
      " ephemeris remark(s), not 27 as expected." << std::endl;

  // Load just one set of spin parameters, with multiple remarks included.
  const HighPrecisionEph::freq_type freq_empty(0);
  const HighPrecisionEph::wave_type wave_empty(0);
  HighPrecisionEph::glitch_type glitch_list(2);
  HighPrecisionEph::glitch_type::iterator glitch_itor = glitch_list.begin();
  glitch_itor->m_epoch = AbsoluteTime("TDB", 54321,  8640.);
  ++glitch_itor;
  glitch_itor->m_epoch = AbsoluteTime("TDB", 65432, 17280.);
  computer3.loadPulsarEph(HighPrecisionEph("TDB", since, until, epoch, 0., 0., 0., 0., 0., 0., epoch,
    freq_empty, 1., wave_empty, wave_empty, glitch_list));
  if (1157 != computer3.getNumPulsarEph())
    err() << "After loading one spin pulsar ephemeris, there were " << computer3.getNumPulsarEph() <<
      " spin pulsar ephemeri(de)s, not 1157 as expected." << std::endl;

  if (26 != computer3.getNumOrbitalEph())
    err() << "After loading one spin pulsar ephemeris, there were " << computer3.getNumOrbitalEph() <<
      " orbital ephemeri(de)s, not 26 as expected." << std::endl;

  if (29 != computer3.getNumEphRemark())
    err() << "After loading one spin pulsar ephemeris, there were " << computer3.getNumEphRemark() <<
      " ephemeris remark(s), not 29 as expected." << std::endl;

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
  expected_spin = "TDB(1156) TT(1)";
  expected_orbital = "TAI(1) TDB(26) TT(1)";
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
  database.registerPulsarEph<PeriodEph>("PER");
  database.registerPulsarEph<HighPrecisionEph>("HP");
  database.registerOrbitalEph<SimpleDdEph>("DD");
  database.registerOrbitalEph<BtModelEph>("BT");
  database.registerOrbitalEph<Ell1ModelEph>("ELL1");
  database.registerOrbitalEph<MssModelEph>("MSS");

  PulsarEphCont pulsar_eph_cont;
  database.getEph(pulsar_eph_cont);
  std::size_t expected_spin = 1155;
  if (pulsar_eph_cont.size() != expected_spin)
    err() << "PulsarDb::getEph(PulsarEphCont &) got " << pulsar_eph_cont.size() << " ephemerides, not " <<
      expected_spin << ", as expected." << std::endl;

  OrbitalEphCont orbital_eph_cont;
  database.getEph(orbital_eph_cont);
  std::size_t expected_orbital = 25;
  if (orbital_eph_cont.size() != expected_orbital) 
    err() << "PulsarDb::getEph(OrbitalEphCont &) got " << orbital_eph_cont.size() << " ephemerides, not " <<
      expected_orbital << ", as expected." << std::endl;

  EphStatusCont eph_status_cont;
  database.getRemark(eph_status_cont);
  std::size_t expected_remark = 12;
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
    virtual const EphStatusCont & getRemark() const { static const EphStatusCont s_remark_cont; return s_remark_cont; }
    virtual PulsarEph * clone() const { return new BogusPulsarEphBase(*this); }
    virtual double calcPulsePhase(const AbsoluteTime & /* ev_time */, double /* phase_offset */ = 0.) const {
      return 0.;
    }
    virtual double calcFrequency(const AbsoluteTime & /* ev_time */, int /* derivative_order */ = 0) const { return 0.; }
    virtual SourcePosition calcPosition(const AbsoluteTime & /* ev_time */) const {
      return SourcePosition(0., 0.); 
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
    BogusOrbitalEphBase(): OrbitalEph(ElapsedTime("TDB", Duration(10.e-9, "Sec")), 100) {}

    virtual const TimeSystem & getSystem() const { return TimeSystem::getSystem("TDB"); }

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
    virtual ElapsedTime calcOrbitalDelay(const AbsoluteTime & /* ev_time */) const {
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
  AbsoluteTime abs_time_earliest("TDB", 51909, 0.);
  AbsoluteTime abs_time_earlier("TDB", 51909, 86399.);
  AbsoluteTime abs_time_during1("TDB", 51910, 21600.);
  AbsoluteTime abs_time_during2("TDB", 51910, 64800.);
  AbsoluteTime abs_time_later("TDB", 51911, 1.);
  AbsoluteTime abs_time_latest("TDB", 51912, 0.);

  // Test determination of status effectiveness.
  bool result_effective = eph_status.effectiveBetween(abs_time_earliest, abs_time_earlier);
  bool expected_effective = false;
  if (result_effective != expected_effective) {
    err() << "EphStatus::effectiveBetween(" << abs_time_earliest << ", " << abs_time_earlier << ") returned " <<
      result_effective << ", not " << expected_effective << " as expected." << std::endl;
  }
  result_effective = eph_status.effectiveBetween(abs_time_earlier, abs_time_during1);
  expected_effective = true;
  if (result_effective != expected_effective) {
    err() << "EphStatus::effectiveBetween(" << abs_time_earlier << ", " << abs_time_during1 << ") returned " <<
      result_effective << ", not " << expected_effective << " as expected." << std::endl;
  }
  result_effective = eph_status.effectiveBetween(abs_time_during1, abs_time_during2);
  expected_effective = true;
  if (result_effective != expected_effective) {
    err() << "EphStatus::effectiveBetween(" << abs_time_during1 << ", " << abs_time_during2 << ") returned " <<
      result_effective << ", not " << expected_effective << " as expected." << std::endl;
  }
  result_effective = eph_status.effectiveBetween(abs_time_during2, abs_time_later);
  expected_effective = true;
  if (result_effective != expected_effective) {
    err() << "EphStatus::effectiveBetween(" << abs_time_during2 << ", " << abs_time_later << ") returned " <<
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
  // Note: the following loads a set of remarks to the ephemeris computer as:
  //       1) a regular remark effective since abs_time_05 until abs_time_07,
  //       2) a regular remark effective since abs_time_09 until abs_time_10,
  //       3) a glitch at abs_time_03, resulting a remark effective since abs_time_03 until abs_time_05, and
  //       4) a glitch at abs_time_07, resulting a remark effective since abs_time_07 until abs_time_09.
  EphComputer eph_computer;
  eph_computer.loadEphRemark(EphStatus(abs_time_02, abs_time_03, Remarked, "Remark No. 1"));
  eph_computer.loadEphRemark(EphStatus(abs_time_03, abs_time_05, Remarked, "Remark No. 2"));
  eph_computer.loadEphRemark(EphStatus(abs_time_05, abs_time_07, Remarked, "Remark No. 3"));
  eph_computer.loadEphRemark(EphStatus(abs_time_07, abs_time_09, Remarked, "Remark No. 4"));
  eph_computer.loadEphRemark(EphStatus(abs_time_09, abs_time_10, Remarked, "Remark No. 5"));
  eph_computer.loadEphRemark(EphStatus(abs_time_03, abs_time_09, Remarked, "Remark No. 6"));
  const HighPrecisionEph::freq_type freq_empty(0);
  const HighPrecisionEph::wave_type wave_empty(0);
  HighPrecisionEph::glitch_type glitch_list(1);
  HighPrecisionEph::glitch_type::iterator glitch_itor = glitch_list.begin();
  glitch_itor->m_epoch = abs_time_02;
  eph_computer.loadPulsarEph(HighPrecisionEph("TDB", abs_time_01, abs_time_03, abs_time_02, 0., 0., 0., 0., 0., 0., abs_time_02,
    freq_empty, 1., wave_empty, wave_empty, glitch_list));
  glitch_itor->m_epoch = abs_time_03;
  eph_computer.loadPulsarEph(HighPrecisionEph("TDB", abs_time_01, abs_time_05, abs_time_02, 0., 0., 0., 0., 0., 0., abs_time_02,
    freq_empty, 1., wave_empty, wave_empty, glitch_list));
  glitch_itor->m_epoch = abs_time_05;
  eph_computer.loadPulsarEph(HighPrecisionEph("TDB", abs_time_01, abs_time_07, abs_time_02, 0., 0., 0., 0., 0., 0., abs_time_02,
    freq_empty, 1., wave_empty, wave_empty, glitch_list));
  glitch_list.resize(3);
  glitch_itor = glitch_list.begin();
  glitch_itor->m_epoch = abs_time_07;
  ++glitch_itor;
  glitch_itor->m_epoch = abs_time_09;
  ++glitch_itor;
  glitch_itor->m_epoch = abs_time_03;
  eph_computer.loadPulsarEph(HighPrecisionEph("TDB", abs_time_01, abs_time_10, abs_time_02, 0., 0., 0., 0., 0., 0., abs_time_02,
    freq_empty, 1., wave_empty, wave_empty, glitch_list));

  // Test selection of ephemeris remarks by an ephemeris computer.
  EphStatusCont expected_status_cont;
  expected_status_cont.push_back(EphStatus(abs_time_03, abs_time_05, Remarked, "Remark No. 2"));
  expected_status_cont.push_back(EphStatus(abs_time_05, abs_time_07, Remarked, "Remark No. 3"));
  expected_status_cont.push_back(EphStatus(abs_time_07, abs_time_09, Remarked, "Remark No. 4"));
  expected_status_cont.push_back(EphStatus(abs_time_03, abs_time_09, Remarked, "Remark No. 6"));
  expected_status_cont.push_back(EphStatus(abs_time_03, abs_time_05, Remarked, "Glitch No. 2"));
  expected_status_cont.push_back(EphStatus(abs_time_05, abs_time_07, Remarked, "Glitch No. 3"));
  expected_status_cont.push_back(EphStatus(abs_time_07, abs_time_10, Remarked, "Glitch No. 4"));
  expected_status_cont.push_back(EphStatus(abs_time_03, abs_time_10, Remarked, "Glitch No. 6"));
  eph_computer.getEphRemark(abs_time_04, abs_time_08, eph_status_cont);
  if (eph_status_cont.size() != expected_status_cont.size()) {
    err() << "EphComputer::getEphRemark method returned " << eph_status_cont.size() << " ephemeris remark(s), not " <<
      expected_status_cont.size() << " as expected" << std::endl;
  } else {
    EphStatusCont::const_iterator result_itor = eph_status_cont.begin();
    EphStatusCont::const_iterator expected_itor = expected_status_cont.begin();
    int remark_number = 1;
    for (; result_itor != eph_status_cont.end() && expected_itor != expected_status_cont.end(); ++result_itor, ++expected_itor,
      ++remark_number) {
      // Compare status code.
      const EphStatusCodeType & result_code = result_itor->getStatusCode();
      const EphStatusCodeType & expected_code = expected_itor->getStatusCode();
      if (result_code != expected_code) {
        err() << "EphComputer::getEphRemark method returned an EphStatus object with status codes of " <<
          result_code << ", not " << expected_code << ", for " << expected_itor->getDescription() << "." << std::endl;
      }

      // Compare effective-since time.
      const AbsoluteTime & result_since = result_itor->getEffectiveSince();
      const AbsoluteTime & expected_since = expected_itor->getEffectiveSince();
      if (!result_since.equivalentTo(expected_since, tolerance)) {
        err() << "EphComputer::getEphRemark method returned an EphStatus object effective since " <<
          result_since << ", not " << expected_since << ", for " << expected_itor->getDescription() << "." << std::endl;
      }

      // Compare effective-until time.
      const AbsoluteTime & result_until = result_itor->getEffectiveUntil();
      const AbsoluteTime & expected_until = expected_itor->getEffectiveUntil();
      if (!result_until.equivalentTo(expected_until, tolerance)) {
        err() << "EphComputer::getEphRemark method returned an EphStatus object effective until " <<
          result_until << ", not " << expected_until << ", for " << expected_itor->getDescription() << "." << std::endl;
      }
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
      ofs << prependDataPath("psrdb_spin_per.txt") << std::endl;
      ofs << prependDataPath("psrdb_spin_hp.txt") << std::endl;
      ofs << prependDataPath("psrdb_binary.txt") << std::endl;
      ofs << prependDataPath("psrdb_binary_bt.txt") << std::endl;
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
  test_name_cont.push_back("par16");
  test_name_cont.push_back("par17");

  // Prepare files to be used in the tests.
  std::string test_pulsardb = prependDataPath("testpsrdb_ephcomp.fits");
  std::string leap_file = prependDataPath("gtephem_leapsec.fits");
  std::string remark_file = prependDataPath("gtephem_remark.txt");
  std::string glitch_file = prependDataPath("gtephem_glitch.txt");
  std::string summary_file = getMethod() + "_summary.txt";
  std::ofstream ofs_summary(summary_file.c_str());
  ofs_summary << test_pulsardb << std::endl;
  ofs_summary << remark_file << std::endl;
  ofs_summary << glitch_file << std::endl;
  ofs_summary << prependDataPath("testpsrdb_crab.fits") << std::endl; // Not strictly necessary for this test.
  ofs_summary << prependDataPath("testpsrdb_text.fits") << std::endl; // Not strictly necessary for this test.
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
      pars["psrdbfile"] = test_pulsardb;
      pars["psrname"] = "PSR B0540-69";
      pars["reftime"] = 212380785.922;
      pars["timeformat"] = "FERMI";
      pars["timesys"] = "TDB";
      log_file_ref = prependOutrefPath(log_file);

    } else if ("par2" == test_name) {
      // Test strict=yes.
      pars["psrdbfile"] = test_pulsardb;
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
      pars["psrdbfile"] = test_pulsardb;
      pars["chatter"] = 3;
      pars["timeformat"] = "FERMI";
      pars["timesys"] = "TDB";
      log_file_ref = prependOutrefPath(log_file);

    } else if ("par4" == test_name) {
      // Test an orbital ephemeris.
      pars["psrname"] = "PSR J1834-0010";
      pars["reftime"] = 212502400.0;
      pars["psrdbfile"] = test_pulsardb;
      pars["timeformat"] = "FERMI";
      pars["timesys"] = "TDB";
      log_file_ref = prependOutrefPath(log_file);

    } else if ("par5" == test_name) {
      // Test chatter=3 with an orbital ephemeris.
      pars["psrname"] = "PSR J1834-0010";
      pars["reftime"] = 212502400.0;
      pars["psrdbfile"] = test_pulsardb;
      pars["chatter"] = 3;
      pars["timeformat"] = "FERMI";
      pars["timesys"] = "TDB";
      log_file_ref = prependOutrefPath(log_file);

    } else if ("par6" == test_name) {
      // Test a reference time during a leap second.
      pars["psrname"] = "PSR B0540-69";
      pars["reftime"] = 157766336.0;
      pars["psrdbfile"] = test_pulsardb;
      pars["timeformat"] = "FERMI";
      pars["timesys"] = "UTC";
      pars["leapsecfile"] = leap_file;
      log_file_ref = prependOutrefPath(log_file);

    } else if ("par7" == test_name) {
      // Create a reference output for testing ISO 8601 format.
      pars["psrname"] = "PSR B0540-69";
      pars["reftime"] = "54368.11169104166666666666"; // Need to be a string to preserve precision.
      pars["psrdbfile"] = test_pulsardb;
      pars["strict"] = "yes";
      pars["timeformat"] = "MJD";
      pars["timesys"] = "TDB";
      log_file_ref = prependOutrefPath(log_file);

    } else if ("par8" == test_name) {
      // Test ISO 8601 calendar format.
      pars["psrname"] = "PSR B0540-69";
      pars["reftime"] = "2007-09-25T02:40:50.106";
      pars["psrdbfile"] = test_pulsardb;
      pars["strict"] = "yes";
      pars["timeformat"] = "ISO";
      pars["timesys"] = "TDB";
      log_file_ref = prependOutrefPath(log_file);

    } else if ("par9" == test_name) {
      // Test ISO 8601 week-day format.
      pars["psrname"] = "PSR B0540-69";
      pars["reftime"] = "2007-W39-2T02:40:50.106";
      pars["psrdbfile"] = test_pulsardb;
      pars["strict"] = "yes";
      pars["timeformat"] = "ISO";
      pars["timesys"] = "TDB";
      log_file_ref = prependOutrefPath(log_file);

    } else if ("par10" == test_name) {
      // Test ISO 8601 ordinal-day format.
      pars["psrname"] = "PSR B0540-69";
      pars["reftime"] = "2007-268T02:40:50.106";
      pars["psrdbfile"] = test_pulsardb;
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
      pars["psrdbfile"] = test_pulsardb;
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
      pars["psrdbfile"] = test_pulsardb;
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
    table->setNumRecords(1);
    tip::Table::Iterator record_itor = table->begin();
    tip::Table::Record & record(*record_itor);
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
