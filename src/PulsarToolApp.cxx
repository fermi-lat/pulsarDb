/** \file PulsarToolApp.cxx
    \brief Implementation of base class for pulsar tool applications.
    \author Masaharu Hirayama, GSSC
            James Peachey, HEASARC/GSSC
*/
#include <algorithm>
#include <cctype>
#include <ctime>
#include <list>
#include <memory>
#include <set>
#include <stdexcept>
#include <utility>

#include "facilities/commonUtilities.h"

#include "hoops/hoops.h"
#include "hoops/hoops_exception.h"

#include "pulsarDb/EphChooser.h"
#include "pulsarDb/EphComputer.h"
#include "pulsarDb/FrequencyEph.h"
#include "pulsarDb/PeriodEph.h"
#include "pulsarDb/PulsarDb.h"
#include "pulsarDb/PulsarEph.h"
#include "pulsarDb/PulsarToolApp.h"
#include "pulsarDb/SimpleDdEph.h"

#include "st_app/AppParGroup.h"

#include "st_facilities/FileSys.h"

#include "st_stream/Stream.h"

#include "timeSystem/AbsoluteTime.h"
#include "timeSystem/BaryTimeComputer.h"
#include "timeSystem/CalendarFormat.h"
#include "timeSystem/Duration.h"
#include "timeSystem/ElapsedTime.h"
#include "timeSystem/EventTimeHandler.h"
#include "timeSystem/GlastTimeHandler.h"
#include "timeSystem/MjdFormat.h"
#include "timeSystem/TimeInterval.h"

#include "tip/Header.h"
#include "tip/IFileSvc.h"
#include "tip/Table.h"
#include "tip/TipException.h"

using namespace st_facilities;
using namespace timeSystem;

namespace pulsarDb {

  PulsarToolApp::PulsarToolApp(): m_event_handler_cont(), m_gti_handler_cont(), m_time_field(),
    m_gti_start_field(), m_gti_stop_field(), m_output_field_cont(),
    m_reference_handler(0), m_computer(0), m_tcmode_dict_bary(), m_tcmode_dict_bin(), m_tcmode_dict_pdot(),
    m_tcmode_bary(ALLOWED), m_tcmode_bin(ALLOWED), m_tcmode_pdot(ALLOWED),
    m_request_bary(false), m_demod_bin(false), m_cancel_pdot(false), m_vary_ra_dec(true),
    m_target_time_system(0), m_target_time_origin("TDB", 0, 0.),
    m_event_handler_itor(m_event_handler_cont.begin()), m_report_eph_status(false) {}

  PulsarToolApp::~PulsarToolApp() throw() {
    resetApp();
  }

  void PulsarToolApp::run() {
    // Clean up from any previous runs.
    resetApp();

    // Run the application.
    runApp();

    // Clean up (close files, reset variables, etc.).
    resetApp();
  }

  AbsoluteTime PulsarToolApp::parseTime(const std::string & time_format, const std::string & time_system,
    const std::string & time_value, std::string & parsed_time_format, std::string & parsed_time_system,
    const tip::Header * header) const {
    // Make upper case copies of input for case insensitive comparisons.
    std::string time_format_rat(time_format);
    for (std::string::iterator itor = time_format_rat.begin(); itor != time_format_rat.end(); ++itor) *itor = std::toupper(*itor);

    // Check whether the header is given or not when time format is specified as FILE.
    if (("FILE" == time_format_rat) && (0 == header)) {
      throw std::runtime_error("File was not given when time format is specified as FILE");
    }

    // Make upper case copies of input for case insensitive comparisons.
    std::string time_system_rat(time_system);
    for (std::string::iterator itor = time_system_rat.begin(); itor != time_system_rat.end(); ++itor) *itor = std::toupper(*itor);

    // First check whether time system should be read from the tip::Header.
    if ("FILE" == time_system_rat) {
      // Then check whether the header is given or not.
      if (header) {
        (*header)["TIMESYS"].get(time_system_rat);
        for (std::string::iterator itor = time_system_rat.begin(); itor != time_system_rat.end(); ++itor) *itor = std::toupper(*itor);
      } else {
        throw std::runtime_error("File was not given when time system is specified as FILE");
      }
    }

    // Replace arguments for time format and time system with rationalized ones.
    parsed_time_format = time_format_rat;
    parsed_time_system = time_system_rat;

    // Parse a time string and convert it to an absolute time.
    AbsoluteTime abs_time("TDB", 0, 0.);
    if ("FILE" == time_format_rat) {
      abs_time = m_reference_handler->parseTimeString(time_value, time_system_rat);

    } else {
      if ("GLAST" == time_format_rat) {
        // TODO: Should we use fitsGen/*/data/ft1.tpl instead?
        std::string tpl_file = facilities::commonUtilities::joinPath(facilities::commonUtilities::getDataPath("pulsarDb"),
          "timeformat_glast.tpl"); // Copied from fitsGen/*/data/ft1.tpl on May 25th, 2008.
        std::ostringstream oss;
        oss << "eventfile" << this << ".fits";
        tip::TipFile tip_file = tip::IFileSvc::instance().createMemFile(oss.str(), tpl_file);
        std::auto_ptr<EventTimeHandler> handler(GlastTimeHandler::createInstance(tip_file.getName(), "EVENTS"));
        if (handler.get()) {
          abs_time = handler->parseTimeString(time_value, time_system_rat);
        } else {
          throw std::runtime_error("Error occurred in parsing a time string \"" + time_value + "\"");
        }

      } else if ("MJD" == time_format_rat) {
        abs_time.set(time_system_rat, MjdFmt, time_value);

      } else if ("ISO" == time_format_rat) {
        bool abs_time_set = false;

        // Try Calendar Date and Time format.
        if (!abs_time_set) {
          try {
            abs_time.set(time_system_rat, CalendarFmt, time_value);
            abs_time_set = true;
          } catch (const std::exception &) {
            abs_time_set = false;
          }
        }

        // Try ISO Week Date and Time format.
        if (!abs_time_set) {
          try {
            abs_time.set(time_system_rat, IsoWeekFmt, time_value);
            abs_time_set = true;
          } catch (const std::exception &) {
            abs_time_set = false;
          }
        }

        // Try Ordinal Date and Time format.
        if (!abs_time_set) {
          try {
            abs_time.set(time_system_rat, OrdinalFmt, time_value);
            abs_time_set = true;
          } catch (const std::exception &) {
            abs_time_set = false;
          }
        }

        if (!abs_time_set) {
          throw std::runtime_error("Error in interpreting the given time in ISO format: \"" + time_value + "\"");
        }

      } else {
        throw std::runtime_error("Time format \"" + time_format + "\" is not supported for ephemeris time");
      }

    }

    return abs_time;
  }

  void PulsarToolApp::openEventFile(const st_app::AppParGroup & pars, bool read_only) {
    // Read parameters from the parameter file.
    std::string event_file = pars["evfile"];
    std::string event_extension = pars["evtable"];

    // Open the event table(s), either for reading or reading and writing.
    FileSys::FileNameCont file_name_cont = FileSys::expandFileList(event_file);
    for (FileSys::FileNameCont::const_iterator itor = file_name_cont.begin(); itor != file_name_cont.end(); ++itor) {
      std::string file_name = *itor;

      // Create and store an event time handler for EVENTS extension.
      EventTimeHandler * event_handler(IEventTimeHandlerFactory::createHandler(file_name, event_extension, read_only));
      m_event_handler_cont.push_back(event_handler);

      // Create and store an event time handler for GTI extension.
      EventTimeHandler * gti_handler(IEventTimeHandlerFactory::createHandler(file_name, "GTI", read_only));
      m_gti_handler_cont.push_back(gti_handler);
    }

    // Set names of TIME column in EVENTS exetension and START/STOP coulmns in GTI extensions.
    m_time_field = pars["timefield"].Value();
    m_gti_start_field = "START";
    m_gti_stop_field = "STOP";

    // Select and set reference handler.
    m_reference_handler = m_event_handler_cont.at(0);
  }

  void PulsarToolApp::reserveOutputField(const std::string & field_name, const std::string & field_format) {
    m_output_field_cont.push_back(std::make_pair(field_name, field_format));
  }

  void PulsarToolApp::defineTimeCorrectionMode(const std::string & mode_name, TimeCorrectionMode_e tcmode_bary,
    TimeCorrectionMode_e tcmode_bin, TimeCorrectionMode_e tcmode_pdot) {
    // Make mode_name argument case-insensitive.
    std::string mode_name_uc = mode_name;
    for (std::string::iterator itor = mode_name_uc.begin(); itor != mode_name_uc.end(); ++itor) *itor = std::toupper(*itor);

    // Store given modes to internal variables.
    m_tcmode_dict_bary[mode_name_uc] = tcmode_bary;
    m_tcmode_dict_bin[mode_name_uc] = tcmode_bin;
    m_tcmode_dict_pdot[mode_name_uc] = tcmode_pdot;
  }

  void PulsarToolApp::selectTimeCorrectionMode(const std::string & mode_name) {
    // Make mode_name argument case-insensitive.
    std::string mode_name_uc = mode_name;
    for (std::string::iterator itor = mode_name_uc.begin(); itor != mode_name_uc.end(); ++itor) *itor = std::toupper(*itor);

    // Select time correction mode.
    std::map<const std::string, TimeCorrectionMode_e>::iterator itor;
    itor = m_tcmode_dict_bary.find(mode_name_uc);
    if (itor != m_tcmode_dict_bary.end()) {
      m_tcmode_bary = m_tcmode_dict_bary[mode_name_uc];
      m_tcmode_bin = m_tcmode_dict_bin[mode_name_uc];
      m_tcmode_pdot = m_tcmode_dict_pdot[mode_name_uc];
    } else {
      throw std::runtime_error("Unknown time correction mode requested");
    }
  }

  void PulsarToolApp::selectTimeCorrectionMode(const st_app::AppParGroup & pars) {
    // Read tcorrect parameter.
    std::string t_correct = pars["tcorrect"];

    // Determine time correction mode.
    selectTimeCorrectionMode(t_correct);
  }

  void PulsarToolApp::initEphComputer(const st_app::AppParGroup & pars, const EphChooser & chooser, st_stream::OStream & os) {
    // Read ephstyle parameter.
    std::string eph_style = pars["ephstyle"];

    // Initialize EphComputer with given ephemeris style.
    initEphComputer(pars, chooser, eph_style, os);
  }

  void PulsarToolApp::initEphComputer(const st_app::AppParGroup & pars, const EphChooser & chooser, const std::string & eph_style,
    st_stream::OStream & os) {
    // Make eph_style argument case-insensitive.
    std::string eph_style_uc = eph_style;
    for (std::string::iterator itor = eph_style_uc.begin(); itor != eph_style_uc.end(); ++itor) *itor = std::toupper(*itor);

    if (eph_style_uc == "DB") {
      // Create ephemeris computer.
      m_computer = new EphComputer(chooser);

    } else {
      // Create ephemeris computer with the sloppy chooser, so that the spin ephemeris given by the user will be always chosen.
      SloppyEphChooser sloppy_chooser;
      m_computer = new EphComputer(sloppy_chooser);

      if (eph_style_uc != "NONE") {
        std::string epoch_time_format = pars["timeformat"];
        std::string epoch_time_sys = pars["timesys"];
        std::string epoch = pars["ephepoch"];
        tip::Header & header = m_reference_handler->getHeader();
        AbsoluteTime abs_epoch = parseTime(epoch_time_format, epoch_time_sys, epoch, epoch_time_format, epoch_time_sys, &header);

        // Read phi0 parameter.  If it doesn't exist in the parameter file, let phi0 = 0.
        double phi0 = 0.;
        try {
          // Try to read phi0 parameter from the parameter file.
          phi0 = pars["phi0"];
        } catch (const hoops::Hexception & x) {
          // Do nothing (i.e., leave phi0=0.) if phi0 is not found, otherwise re-throw the exception.
          if (hoops::PAR_NOT_FOUND != x.Code()) throw;
        }

        // Read RA and Dec from a parameter file.
        double ra = pars["ra"];
        double dec = pars["dec"];

        // Handle either period or frequency-style input.
        if (eph_style_uc == "FREQ") {
          double f0 = pars["f0"];
          double f1 = pars["f1"];
          double f2 = pars["f2"];

          if (0. >= f0) throw std::runtime_error("Non-positive frequency is given to f0 parameter (must be positive)");

          // Add the ephemeris the user provided.
          m_computer->loadPulsarEph(FrequencyEph(epoch_time_sys, abs_epoch, abs_epoch, abs_epoch, ra, dec, phi0, f0, f1, f2));
        } else if (eph_style_uc == "PER") {
          double p0 = pars["p0"];
          double p1 = pars["p1"];
          double p2 = pars["p2"];

          if (0. >= p0) throw std::runtime_error("Non-positive period is given to p0 parameter (must be positive)");

          // Add the ephemeris the user provided.
          m_computer->loadPulsarEph(PeriodEph(epoch_time_sys, abs_epoch, abs_epoch, abs_epoch, ra, dec, phi0, p0, p1, p2));
        } else {
          throw std::runtime_error("Unknown ephemeris style \"" + eph_style + "\" was specified");
        }
      }
    }

    // Determine whether ephemeris status should be reported or not.
    m_report_eph_status = pars["reportephstatus"];

    // Open pulsar database if necessary.
    if (eph_style_uc == "DB" || m_tcmode_bin != SUPPRESSED || m_report_eph_status) {
      // Create an empty pulsar ephemerides database, using template file.
      std::string tpl_file = facilities::commonUtilities::joinPath(facilities::commonUtilities::getDataPath("pulsarDb"), "PulsarDb.tpl");
      static const PulsarDb::TableCont::size_type default_spin_extension = 1;
      static const PulsarDb::TableCont::size_type default_orbital_extension = 2;
      PulsarDb database(tpl_file, default_spin_extension, default_orbital_extension);

      // Register PulsarEph and OrbitalEph subclasses for various ephemeris models.
      database.registerPulsarEph<FrequencyEph>("FREQ");
      database.registerOrbitalEph<SimpleDdEph>("DD");

      // Load the given ephemerides database(s).
      std::string psrdb_file = pars["psrdbfile"];
      FileSys::FileNameCont file_name_cont = FileSys::expandFileList(psrdb_file);
      for (FileSys::FileNameCont::const_iterator itor = file_name_cont.begin(); itor != file_name_cont.end(); ++itor) {
        database.load(*itor);
      }
      int num_eph_total = database.getNumEph();

      // Select only ephemerides for this pulsar.
      std::string psr_name = pars["psrname"];
      database.filterName(psr_name);
      int num_eph_psrname = database.getNumEph();

      // Select only ephemerides with the solar system ephemeris used for barycentering.
      std::string solar_eph = pars["solareph"];
      std::string match_solar_eph = pars["matchsolareph"];
      for (std::string::iterator itor = match_solar_eph.begin(); itor != match_solar_eph.end(); ++itor) *itor = std::toupper(*itor);
      bool solar_eph_must_match = (match_solar_eph == "PSRDB" || match_solar_eph == "ALL");
      if (solar_eph_must_match) database.filterSolarEph(solar_eph);
      int num_eph_solareph = database.getNumEph();

      // Report ephemeris loading summary.
      std::list<std::string> command_history;
      std::list<std::string> ancestry_record;
      std::string dashes(26, '=');
      database.getHistory(command_history, ancestry_record);
      if (command_history.size()) {
        os.prefix() << dashes << std::endl;
        os.prefix() << "Pulsar ephemerides are loaded and/or filtered as follows:" << std::endl;
        for (std::list<std::string>::const_iterator itor = command_history.begin(); itor != command_history.end(); ++itor) {
          os.prefix() << "   " << *itor << std::endl;
        }
        if (ancestry_record.size()) {
          os.prefix() << dashes << std::endl;
          os.prefix() << "Pulsar ephemerides FITS files loaded are:" << std::endl;
          int psrdb_index = 0;
          for (std::list<std::string>::const_iterator itor = ancestry_record.begin(); itor != ancestry_record.end(); ++itor) {
            std::string list_prefix;
            if (0 == itor->find("PULSARDB")) {
              ++psrdb_index;
              std::ostringstream oss;
              oss << "[" << psrdb_index << "] ";
              list_prefix = oss.str();
            } else {
              list_prefix = "   ";
            }
            os.prefix() << list_prefix << *itor << std::endl;
          }
        }
        os.prefix() << dashes << std::endl;
      }

      // Load the selected pulsar ephemerides.
      if (eph_style_uc == "DB") {
        // Thrown an exception if no ephemeris is left in the database.
        if (0 == database.getNumEph()) {
          std::ostringstream os;
          os << "No spin ephemeris is available for a requested condition. Brief summary of ephemeris selection is following." <<
            std::endl;
          os << num_eph_total << " spin ephemeri(de)s in the database." << std::endl;
          os << num_eph_psrname << " spin ephemeri(de)s for pulsar \"" << psr_name << "\" in the database." << std::endl;
          if (solar_eph_must_match) {
            os << num_eph_solareph << " spin ephemeri(de)s for pulsar \"" << psr_name << "\" with solar system ephemeris \"" <<
              solar_eph << "\" in the database.";
          } else {
            os << "(Solar system ephemeris in spin parameters was not requested to match \"" << solar_eph << "\" given by user.)";
          }
          throw std::runtime_error(os.str());
        }

        // Load spin parameters.
        m_computer->loadPulsarEph(database);
      }

      // Load orbital parameters.
      if (m_tcmode_bin != SUPPRESSED) m_computer->loadOrbitalEph(database);

      // Load ephemeris remarks.
      if (m_report_eph_status) m_computer->loadEphRemark(database);
    }
  }

  AbsoluteTime PulsarToolApp::computeTimeBoundary(bool request_start_time, bool request_time_correction) const {
    bool candidate_found = false;
    AbsoluteTime abs_candidate_time("TDB", 0, 0.);

    // First, look for requested time (start or stop) in the GTI.
    for (handler_cont_type::const_iterator itor = m_gti_handler_cont.begin(); itor != m_gti_handler_cont.end(); ++itor) {
      EventTimeHandler & gti_handler = **itor;

      // If possible, get tstart (or tstop) from first (or last) interval in GTI extension.
      gti_handler.setFirstRecord();
      if (!gti_handler.isEndOfTable()) {
        // Set field name and move to the last record if necessary.
        std::string field_name;
        if (request_start_time) {
          field_name = m_gti_start_field;
        } else {
          gti_handler.setLastRecord();
          field_name = m_gti_stop_field;
        }

        // Read GTI START (or STOP) column value as AbsoluteTime.
        AbsoluteTime abs_gti_time = readTimeColumn(gti_handler, field_name, request_time_correction);

        if (candidate_found) {
          // See if the time is "better" than the current candidate (i.e., earlier for start or later for stop).
          if (request_start_time == (abs_gti_time < abs_candidate_time)) abs_candidate_time = abs_gti_time;
        } else {
          // First candidate is the currently picked time.
          abs_candidate_time = abs_gti_time;
          candidate_found = true;
        }
      }
    }

    // In the unlikely event that there were no GTI files, no intervals in the GTI, and no event files, this is a
    // serious problem.
    if (!candidate_found) throw std::runtime_error("Could not determine start/stop of the observation interval");

    // Return the candidate.
    return abs_candidate_time;
  }

  void PulsarToolApp::initTimeCorrection(const st_app::AppParGroup & pars, bool vary_ra_dec, bool guess_pdot) {
    // Read timeorigin parameter.
    std::string str_origin = pars["timeorigin"];

    // Initialize time correction with the given time origin.
    initTimeCorrection(pars, vary_ra_dec, guess_pdot, str_origin);
  }

  void PulsarToolApp::initTimeCorrection(const st_app::AppParGroup & pars, bool vary_ra_dec, bool guess_pdot,
    const std::string & str_origin) {
    AbsoluteTime abs_origin("TDB", 0, 0.);

    // Make str_origin argument case-insensitive.
    std::string str_origin_uc = str_origin;
    for (std::string::iterator itor = str_origin_uc.begin(); itor != str_origin_uc.end(); ++itor) *itor = std::toupper(*itor);

    // Compute the time origin specified by str_origin argument.
    if (str_origin_uc == "START") {
      // Get the uncorrected start time of event list.
      abs_origin = computeTimeBoundary(true, false);

    } else if (str_origin_uc == "STOP") {
      // Get the uncorrected stop time of event list.
      abs_origin = computeTimeBoundary(false, false);

    } else if (str_origin_uc == "MIDDLE") {
      // Use the center of the observation as the time origin.
      AbsoluteTime abs_tstart = computeTimeBoundary(true, false);
      AbsoluteTime abs_tstop = computeTimeBoundary(false, false);

      // Get time system name from the header.
      std::string time_sys;
      tip::Header & header = m_reference_handler->getHeader();
      header["TIMESYS"].get(time_sys);

      // Compute mid-time difference between TSTART and TSTOP.
      double elapsed = (abs_tstop - abs_tstart).computeDuration(time_sys, "Sec");
      abs_origin = abs_tstart + ElapsedTime(time_sys, Duration(elapsed * 0.5, "Sec"));

    } else if (str_origin_uc == "USER") {
      // Get time of origin and its format and system from parameters.
      std::string origin_time = pars["usertime"];
      std::string origin_time_format = pars["userformat"];
      std::string origin_time_sys = pars["usersys"].Value();

      // Convert user-specified time into AbsoluteTime.
      std::string parsed_time_format;
      std::string parsed_time_sys;
      tip::Header & header = m_reference_handler->getHeader();
      abs_origin = parseTime(origin_time_format, origin_time_sys, origin_time, parsed_time_format, parsed_time_sys, &header);

    } else {
      throw std::runtime_error("Unsupported time origin " + str_origin);
    }

    // Initialize time correction with the time origin just computed.
    initTimeCorrection(pars, vary_ra_dec, guess_pdot, abs_origin);
  }

  void PulsarToolApp::initTimeCorrection(const st_app::AppParGroup & pars, bool vary_ra_dec, bool guess_pdot,
    const AbsoluteTime & abs_origin) {
    // Determine whether to perform binary demodulation.
    m_demod_bin = false;
    if ((m_tcmode_bin == REQUIRED) || (m_tcmode_bin == ALLOWED)) {
      // Check whether orbital parameters are available for binary demodulation.
      if (0 != m_computer->getNumOrbitalEph()) {
        m_demod_bin = true;
      } else if (m_tcmode_bin == REQUIRED) {
        throw std::runtime_error("Binary demodulation was required by user, but no orbital ephemeris was found");
      }
    }

    // Determine whether to cancel pdot.
    m_cancel_pdot = false;
    if ((m_tcmode_pdot == REQUIRED) || (m_tcmode_pdot == ALLOWED)) {
      if (!guess_pdot || 0 != m_computer->getNumPulsarEph()) {
        m_cancel_pdot = true;
      } else if (m_tcmode_pdot == REQUIRED) {
        throw std::runtime_error("Pdot cancellation was required by user, but no spin ephemeris was found");
      }
    }

    // Determine whether to request barycentric correction.
    if (m_tcmode_bary == SUPPRESSED && (m_demod_bin || m_cancel_pdot)) {
      throw std::runtime_error("Barycentric correction is suppressed when binary demodulation or pdot cancellation is requested");
    }
    m_request_bary = (m_tcmode_bary == REQUIRED || m_tcmode_bary == ALLOWED);

    // Initialize the time series to analyze.
    m_target_time_origin = abs_origin;
    std::string time_system_name("");
    bool time_system_set = false;
    if (!m_request_bary && !m_demod_bin && !m_cancel_pdot) {
      // When NO corrections are requested, the analysis will be performed in the time system written in event files,
      // requiring all event files have same time system.
      std::string this_time_system;
      for (handler_cont_type::const_iterator itor = m_event_handler_cont.begin(); itor != m_event_handler_cont.end(); ++itor) {
        const EventTimeHandler & event_handler = **itor;
        const tip::Header & header(event_handler.getHeader());
        header["TIMESYS"].get(this_time_system);
        if (!time_system_set) {
          time_system_name = this_time_system;
          time_system_set = true;
        } else if (this_time_system != time_system_name) {
          throw std::runtime_error("event files with different TIMESYS values cannot be combined");
        }
      }
    } else {
      // When ANY correction(s) are requested, the analysis will be performed in TDB system.
      time_system_name = "TDB";
      time_system_set = true;
    }
    m_target_time_system = &timeSystem::TimeSystem::getSystem(time_system_name);

    // Check whether time system is successfully set.
    if (!time_system_set) throw std::runtime_error("cannot determine time system for the time series to analyze");

    // Compute spin ephemeris to be used in pdot cancellation, and replace PulsarEph in EphComputer with it.
    if (m_cancel_pdot) {
      const int max_derivative = 2;
      if (guess_pdot) {
        // Compute an ephemeris at m_target_time_origin to use for pdot cancellation, up to the second time derivative.
        m_computer->setPdotCancelParameter(m_target_time_origin, max_derivative);

      } else {
        // Read parameters for pdot cancellation from pfile.
        std::string eph_style = pars["ephstyle"];
        for (std::string::iterator itor = eph_style.begin(); itor != eph_style.end(); ++itor) *itor = std::toupper(*itor);

        // Set dummy ra, dec, and phi0 (those are not used in pdot cancellation).
        double ra = 0.;
        double dec = 0.;
        double phi0 = 0.;
        if (eph_style == "FREQ") {
          std::vector<double> fdot_ratio(max_derivative, 0.);
          if (max_derivative > 0) fdot_ratio[0] = pars["f1f0ratio"];
          if (max_derivative > 1) fdot_ratio[1] = pars["f2f0ratio"];
          m_computer->setPdotCancelParameter(m_target_time_system->getName(), m_target_time_origin, fdot_ratio);
        } else if (eph_style == "PER") {
          double p0 = 1.;
          double p1 = pars["p1p0ratio"];
          double p2 = pars["p2p0ratio"];
          m_computer->setPdotCancelParameter(m_target_time_origin,
            PeriodEph(m_target_time_system->getName(), m_target_time_origin, m_target_time_origin, m_target_time_origin,
              ra, dec, phi0, p0, p1, p2), max_derivative);
        } else {
          throw std::runtime_error("Unknown ephemeris style given to ephstyle parameter (must be either FREQ or PER)");
        }
      }
    }

    // Initialize barycentric corrections.
    m_vary_ra_dec = vary_ra_dec;
    if (m_request_bary) {
      // Read parameters from the parameter file.
      std::string sc_file = pars["scfile"];
      std::string sc_extension = pars["sctable"];
      std::string solar_eph = pars["solareph"];
      double ang_tolerance = pars["angtol"];

      // Get the uncorrected start time of event list as a reference time for initial RA and Dec for barycentric corrections.
      AbsoluteTime abs_start = computeTimeBoundary(true, false);

      // Determine whether to check solar system ephemeris used for barycentered event files.
      std::string match_solar_eph = pars["matchsolareph"];
      for (std::string::iterator itor = match_solar_eph.begin(); itor != match_solar_eph.end(); ++itor) *itor = std::toupper(*itor);
      bool request_match_solar_eph = (match_solar_eph == "EVENT" || match_solar_eph == "ALL");

      // Get initial RA and Dec for barycentric corrections.
      double ra = 0.;
      double dec = 0.;
      if (!m_vary_ra_dec) {
        ra = pars["ra"];
        dec = pars["dec"];
      }

      // Initialize event extensions for barycentric corrections.
      for (handler_cont_type::const_iterator itor = m_event_handler_cont.begin(); itor != m_event_handler_cont.end(); ++itor) {
        EventTimeHandler & event_handler = **itor;
        event_handler.initTimeCorrection(sc_file, sc_extension, solar_eph, request_match_solar_eph, ang_tolerance);
        if (!m_vary_ra_dec) event_handler.setSourcePosition(ra, dec);
      }

      // Initialize GTI extensions for barycentric corrections.
      for (handler_cont_type::const_iterator itor = m_gti_handler_cont.begin(); itor != m_gti_handler_cont.end(); ++itor) {
        EventTimeHandler & gti_handler = **itor;
        gti_handler.initTimeCorrection(sc_file, sc_extension, solar_eph, request_match_solar_eph, ang_tolerance);
        if (!m_vary_ra_dec) gti_handler.setSourcePosition(ra, dec);
      }
    }
  }

  double PulsarToolApp::computeElapsedSecond(const AbsoluteTime & abs_time) const {
    return (abs_time - m_target_time_origin).computeDuration(m_target_time_system->getName(), "Sec");
  }

  AbsoluteTime PulsarToolApp::computeAbsoluteTime(double elapsed_time) const {
    return m_target_time_origin + ElapsedTime(m_target_time_system->getName(), Duration(elapsed_time, "Sec"));
  }

  void PulsarToolApp::setupCurrentEventTable() {
    if (m_event_handler_itor != m_event_handler_cont.end()) {
      // Set to the first event in the table.
      EventTimeHandler & handler = **m_event_handler_itor;
      handler.setFirstRecord();

      // Append requested output field(s), if missing.
      tip::Table & table = handler.getTable();
      for (std::vector<std::pair<std::string, std::string> >::const_iterator itor = m_output_field_cont.begin();
           itor != m_output_field_cont.end(); ++itor) {
        std::string field_name = itor->first;
        std::string field_format = itor->second;
        try {
          table.getFieldIndex(field_name);
        } catch (const tip::TipException &) {
          table.appendField(field_name, field_format);
        }
      }
    }
  }

  void PulsarToolApp::setFirstEvent() {
    // Set event table iterator.
    for (m_event_handler_itor = m_event_handler_cont.begin(); m_event_handler_itor != m_event_handler_cont.end();
         ++m_event_handler_itor) {
      // Setup current event table.
      setupCurrentEventTable();

      // Skip an empty event table.
      EventTimeHandler & handler = **m_event_handler_itor;
      if (!handler.isEndOfTable()) break;
    }
  }

  void PulsarToolApp::setNextEvent() {
    // Increment event record iterator.
    (*m_event_handler_itor)->setNextRecord();

    // Move on to the next file if it reaches to the end of table, until an event is found.
    while ((*m_event_handler_itor)->isEndOfTable()) {
      // Increment event table iterator.
      ++m_event_handler_itor;

      // Stop if the end of file list is reached.
      if (m_event_handler_itor == m_event_handler_cont.end()) break;

      // Setup new event table.
      setupCurrentEventTable();
    }
  }

  bool PulsarToolApp::isEndOfEventList() const {
    return (m_event_handler_itor == m_event_handler_cont.end());
  }

  AbsoluteTime PulsarToolApp::getEventTime() const {
    return readTimeColumn(**m_event_handler_itor, m_time_field, true);
  }

  AbsoluteTime PulsarToolApp::getStartTime() const {
    return computeTimeBoundary(true, true);
  }

  AbsoluteTime PulsarToolApp::getStopTime() const {
    return computeTimeBoundary(false, true);
  }

  EphComputer & PulsarToolApp::getEphComputer() const {
    return *m_computer;
  }

  void PulsarToolApp::reportEphStatus(st_stream::OStream & os, const std::set<EphStatusCodeType> & code_to_report) const {
    // Return if reporting is not requested.
    if (!m_report_eph_status) return;
    if (code_to_report.empty()) return;

    // Get time boundaries from GTI tables --- before time corrections applied.
    AbsoluteTime start_time_uc = computeTimeBoundary(true, false);
    AbsoluteTime stop_time_uc = computeTimeBoundary(false, false);

    // Try to compute time-corrected time boundaries. If unsuccessful, apply no corrections.
    AbsoluteTime start_time_tc = start_time_uc;
    try {
      start_time_tc = computeTimeBoundary(true, true);
    } catch (const std::exception &) {
      start_time_tc = start_time_uc;
    }
    AbsoluteTime stop_time_tc = start_time_uc;
    try {
      stop_time_tc = computeTimeBoundary(false, true);
    } catch (const std::exception &) {
      stop_time_tc = stop_time_uc;
    }

    // Find a time interval to investigate.
    const AbsoluteTime & start_time = std::min(start_time_uc, start_time_tc);
    const AbsoluteTime & stop_time = std::max(stop_time_uc, stop_time_tc);

    // Report ephemeris status with the time boundaries obtained.
    reportEphStatus(os, start_time, stop_time, code_to_report);
  }

  void PulsarToolApp::reportEphStatus(st_stream::OStream & os, const timeSystem::AbsoluteTime & start_time,
    const timeSystem::AbsoluteTime & stop_time, const std::set<EphStatusCodeType> & code_to_report) const {
    // Return if reporting is not requested.
    if (!m_report_eph_status) return;
    if (code_to_report.empty()) return;

    // Create a container for ephemeris status to report.
    typedef std::list<const EphStatus *> EphStatusList;
    EphStatusList status_to_report;

    // Collect ephemeris remark(s) that become effective during a Good Time Interval.
    EphStatusCont eph_remark_cont;
    if (code_to_report.find(Remarked) != code_to_report.end()) {
      m_computer->getEphRemark(start_time, stop_time, eph_remark_cont);
      for (EphStatusCont::const_iterator itor = eph_remark_cont.begin(); itor != eph_remark_cont.end(); ++itor) {
        const EphStatus & eph_status = *itor;
        const EphStatusCodeType & status_code = eph_status.getStatusCode();
        if (code_to_report.find(status_code) != code_to_report.end()) status_to_report.push_back(&eph_status);
      }
    }

    // Find ephemeris gaps.
    EphStatusCont eph_status_cont;
    if (code_to_report.find(Unavailable) != code_to_report.end() || code_to_report.find(Extrapolated) != code_to_report.end()) {
      m_computer->examinePulsarEph(start_time, stop_time, eph_status_cont);

      // Subselect ephemeris status by status code.
      EphStatusList eph_status_list;
      for (EphStatusCont::const_iterator itor = eph_status_cont.begin(); itor != eph_status_cont.end(); ++itor) {
        const EphStatus & eph_status = *itor;
        const EphStatusCodeType & status_code = eph_status.getStatusCode();
        if (code_to_report.find(status_code) != code_to_report.end()) eph_status_list.push_back(&eph_status);
      }

      // Flag ephemeris gaps that overlap at least one of GTI's.
      std::set<const EphStatus *> selected_eph_status;
      for (handler_cont_type::const_iterator gti_itor = m_gti_handler_cont.begin(); gti_itor != m_gti_handler_cont.end(); ++gti_itor) {
        EventTimeHandler & gti_handler = **gti_itor;
        for (gti_handler.setFirstRecord(); !gti_handler.isEndOfTable(); gti_handler.setNextRecord()) {
          // Read a GTI uncorrected.
          AbsoluteTime gti_start_uc = readTimeColumn(gti_handler, m_gti_start_field, false);
          AbsoluteTime gti_stop_uc = readTimeColumn(gti_handler, m_gti_stop_field, false);

          // Try to read a GTI time-corrected.
          AbsoluteTime gti_start_tc = gti_start_uc;
          try {
            gti_start_tc = readTimeColumn(gti_handler, m_gti_start_field, true);
          } catch (const std::exception &) {
            gti_start_tc = gti_start_uc;
          }
          AbsoluteTime gti_stop_tc = gti_stop_uc;
          try {
            gti_stop_tc = readTimeColumn(gti_handler, m_gti_stop_field, true);
          } catch (const std::exception &) {
            gti_stop_tc = gti_stop_uc;
          }

          // Take the widest time interval to examine ephemeris gaps.
          const AbsoluteTime & gti_start = std::min(gti_start_uc, gti_start_tc);
          const AbsoluteTime & gti_stop = std::max(gti_stop_uc, gti_stop_tc);

          // Select ephemeris status that ovelaps with this GTI.
          for (EphStatusList::const_iterator eph_itor = eph_status_list.begin(); eph_itor != eph_status_list.end(); ++eph_itor) {
            const EphStatus & eph_status = **eph_itor;
            if (eph_status.effectiveBetween(gti_start, gti_stop)) selected_eph_status.insert(&eph_status);
          }
        }
      }

      // Collect selected ephemeris gap(s) in the order of their original appearance.
      for (EphStatusList::const_iterator itor = eph_status_list.begin(); itor != eph_status_list.end(); ++itor) {
        const EphStatus & eph_status = **itor;
        if (selected_eph_status.find(&eph_status) != selected_eph_status.end()) status_to_report.push_back(&eph_status);
      }
    }

    // Return if nothing to report at this point.
    if (status_to_report.empty()) return;

    // Report collected ephemeris status following a header line.
    os.prefix() << "The following pulsar ephemeris status are reported." << std::endl;
    int status_number = 1;
    for (EphStatusList::const_iterator itor = status_to_report.begin(); itor != status_to_report.end(); ++itor, ++status_number) {
      const EphStatus & eph_status = **itor;
      std::string eph_status_string;
      try {
        eph_status_string = eph_status.report("TDB", MjdFmt);
      } catch (const std::exception &) {
        eph_status_string = eph_status.report("TDB", CalendarFmt);
      }
      os << "[" << status_number << "] " << eph_status_string << std::endl;
    }
  }

  void PulsarToolApp::writeParameter(const st_app::AppParGroup & pars) {
    // Prepare values for CREATOR and DATE keywords.
    const std::string creator_value = getName() + " " + getVersion();
    const time_t modification_time = time(0);

    for (handler_cont_type::iterator itor = m_event_handler_cont.begin(); itor != m_event_handler_cont.end(); ++itor) {
      EventTimeHandler & event_handler = **itor;
      tip::Header & header(event_handler.getHeader());

      // Format DATE keyword value.
      const std::string date_value(header.formatTime(modification_time));

      // Write out all the parameters into HISTORY keywords.
      header.addHistory("File modified by " + creator_value + " on " + date_value);
      for (hoops::ConstGenParItor par_itor = pars.begin(); par_itor != pars.end(); ++par_itor) {
        std::ostringstream oss_par;
        oss_par << getName() << ".par: " << **par_itor;
        header.addHistory(oss_par.str());
      }
    }
  }

  void PulsarToolApp::resetApp() {
    // Reset ephemeris status reporting flag.
    m_report_eph_status = false;

    // Reset target time settings.
    m_target_time_origin = AbsoluteTime("TDB", 0, 0.);
    m_target_time_system = 0;

    // Reset time correction flags.
    m_vary_ra_dec = true;
    m_cancel_pdot = false;
    m_demod_bin = false;
    m_request_bary = false;

    // Reset time correction modes.
    m_tcmode_pdot = ALLOWED;
    m_tcmode_bin = ALLOWED;
    m_tcmode_bary = ALLOWED;
    m_tcmode_dict_pdot.clear();
    m_tcmode_dict_bin.clear();
    m_tcmode_dict_bary.clear();

    // Get rid of current EphComputer.
    delete m_computer; m_computer = 0;

    // Reset reference handler.
    m_reference_handler = 0;

    // Reset columns to be modified.
    m_output_field_cont.clear();

    // Clear out any gti handlers.
    for (handler_cont_type::reverse_iterator itor = m_gti_handler_cont.rbegin(); itor != m_gti_handler_cont.rend(); ++itor) {
      delete *itor;
    }
    m_gti_handler_cont.clear();

    // Clear out any event handlers.
    for (handler_cont_type::reverse_iterator itor = m_event_handler_cont.rbegin(); itor != m_event_handler_cont.rend(); ++itor) {
      delete *itor;
    }
    m_event_handler_cont.clear();

    // Reset iterator so that it points to the beginning of the (empty) container.
    m_event_handler_itor = m_event_handler_cont.begin();
  }

  AbsoluteTime PulsarToolApp::readTimeColumn(EventTimeHandler & handler, const std::string & column_name,
    bool request_time_correction) const {

    // Read the original photon arrival time.
    AbsoluteTime abs_time = handler.readTime(column_name);

    // Apply selected corrections, if requested.
    if (request_time_correction) {
      // Apply barycentric correction, if requested.
      if (m_request_bary) {
        // Reset RA and Dec for the given arrival time, if requested.
        if (m_vary_ra_dec) {
          std::pair<double, double> ra_dec = m_computer->calcSkyPosition(abs_time);
          handler.setSourcePosition(ra_dec.first, ra_dec.second);
        }

        // Try barycentric correction with the RA and Dec.
        abs_time = handler.getBaryTime(column_name);
      }

      // Apply binary demodulation, if requested.
      if (m_demod_bin) m_computer->demodulateBinary(abs_time);

      // Apply pdot cancellation corrections, if requested.
      if (m_cancel_pdot) m_computer->cancelPdot(abs_time);
    }

    // Return the requested time.
    return abs_time;
  }
}
