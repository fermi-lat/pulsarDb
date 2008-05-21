/** \file PulsarToolApp.cxx
    \brief Implementation of base class for pulsar tool applications.
    \author Masaharu Hirayama, GSSC
            James Peachey, HEASARC/GSSC
*/
#include <cctype>
#include <memory>
#include <stdexcept>
#include <utility>

#include "facilities/commonUtilities.h"

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

#include "timeSystem/AbsoluteTime.h"
#include "timeSystem/BaryTimeComputer.h"
#include "timeSystem/EventTimeHandler.h"
#include "timeSystem/GlastMetRep.h"
#include "timeSystem/GlastTimeHandler.h"
#include "timeSystem/TimeRep.h"

#include "tip/Header.h"
#include "tip/Table.h"

using namespace st_facilities;
using namespace timeSystem;

namespace pulsarDb {

  PulsarToolApp::PulsarToolApp(): m_event_handler_cont(), m_gti_handler_cont(), m_time_field(),
    m_gti_start_field(), m_gti_stop_field(), m_output_field_cont(),
    m_reference_header(0), m_computer(0), m_tcmode_dict_bary(), m_tcmode_dict_bin(), m_tcmode_dict_pdot(),
    m_tcmode_bary(ALLOWED), m_tcmode_bin(ALLOWED), m_tcmode_pdot(ALLOWED),
    m_request_bary(false), m_demod_bin(false), m_cancel_pdot(false), m_target_time_rep(0),
    m_event_handler_itor(m_event_handler_cont.begin()) {}

  PulsarToolApp::~PulsarToolApp() throw() {
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

    // Create representation for this time format and time system.
    AbsoluteTime abs_time("TDB", 0, 0.);
    if ("FILE" == time_format_rat) {
      // Check TELESCOP keyword for supported missions.
      std::string telescope;
      (*header)["TELESCOP"].get(telescope);
      for (std::string::iterator itor = telescope.begin(); itor != telescope.end(); ++itor) *itor = std::toupper(*itor);
      if (telescope != "GLAST") throw std::runtime_error("Only GLAST supported for now");

      // Get the mjdref from the header, which is not as simple as just reading a single keyword.
      MjdRefDatabase mjd_ref_db;
      IntFracPair mjd_ref(mjd_ref_db(*header));
      MetRep time_rep(time_system_rat, mjd_ref, 0.);

      // Assign the time supplied by the user to the representation.
      time_rep.assign(time_value);
      abs_time = time_rep;

    } else {
      // Create representation for supported time formats.
      if ("GLAST" == time_format_rat) {
        GlastMetRep time_rep(time_system_rat, 0.);

        // Assign the time supplied by the user to the representation.
        time_rep.assign(time_value);
        abs_time = time_rep;

      } else if ("MJD" == time_format_rat) {
        abs_time.set(time_system_rat, "MJD", time_value);
      } else {
        throw std::runtime_error("Time format \"" + time_format + "\" is not supported for ephemeris time");
      }

    }

    return abs_time;
  }

  TimeRep * PulsarToolApp::createMetRep(const std::string & time_system, const AbsoluteTime & abs_reference) const {
    // Compute MJD of abs_reference (the origin of the time series to analyze), to be given as MJDREF of MetRep.
    // NOTE: MetRep should take AbsoluteTime for its MJDREF (Need refactor of AbsoluteTime for that).
    // TODO: Once MetRep is refactored, remove this method.
    Mjd mjd(0, 0.);
    abs_reference.get(time_system, mjd);

    // Create MetRep to represent the time series to analyze and return it.
    return new MetRep(time_system, mjd.m_int, mjd.m_frac, 0.);
  }

  void PulsarToolApp::openEventFile(const st_app::AppParGroup & pars, bool read_only) {
    // Read parameters from the parameter file.
    std::string event_file = pars["evfile"];
    std::string event_extension = pars["evtable"];
    std::string sc_file = pars["scfile"];
    std::string sc_extension = pars["sctable"];
    std::string solar_eph = pars["solareph"];
    double ang_tolerance = pars["angtol"];

    // Determine whether to check solar system ephemeris used for barycentered event files.
    std::string match_solar_eph = pars["matchsolareph"];
    for (std::string::iterator itor = match_solar_eph.begin(); itor != match_solar_eph.end(); ++itor) *itor = std::toupper(*itor);
    bool check_solar_eph = (match_solar_eph == "EVENT" || match_solar_eph == "ALL");

    // Open the event table(s), either for reading or reading and writing.
    FileSys::FileNameCont file_name_cont = FileSys::expandFileList(event_file);
    for (FileSys::FileNameCont::const_iterator itor = file_name_cont.begin(); itor != file_name_cont.end(); ++itor) {
      std::string file_name = *itor;

      // Create and store an event time handler for EVENTS extension.
      EventTimeHandler * event_handler(IEventTimeHandlerFactory::createHandler(file_name, event_extension, sc_file, sc_extension,
        ang_tolerance, read_only));
      m_event_handler_cont.push_back(event_handler);

      // Check solar system ephemeris for EVENTS extension.
      if (check_solar_eph) event_handler->checkSolarEph(solar_eph);

      // Create and store an event time handler for GTI extension.
      EventTimeHandler * gti_handler(IEventTimeHandlerFactory::createHandler(file_name, "GTI", sc_file, sc_extension,
        ang_tolerance, read_only));
      m_gti_handler_cont.push_back(gti_handler);

      // Check solar system ephemeris for GTI extension.
      if (check_solar_eph) gti_handler->checkSolarEph(solar_eph);
    }

    // Set names of TIME column in EVENTS exetension and START/STOP coulmns in GTI extensions.
    m_time_field = pars["timefield"].Value();
    m_gti_start_field = "START";
    m_gti_stop_field = "STOP";

    // Select and set reference header.
    EventTimeHandler * reference_handler = m_event_handler_cont.at(0);
    m_reference_header = &(reference_handler->getHeader());
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
      throw std::runtime_error("Unknown time correction mode requested.");
    }
  }

  void PulsarToolApp::selectTimeCorrectionMode(const st_app::AppParGroup & pars) {
    // Read tcorrect parameter.
    std::string t_correct = pars["tcorrect"];

    // Determine time correction mode.
    selectTimeCorrectionMode(t_correct);
  }

  void PulsarToolApp::initEphComputer(const st_app::AppParGroup & pars, const EphChooser & chooser) {
    // Read ephstyle parameter.
    std::string eph_style = pars["ephstyle"];

    // Initialize EphComputer with given ephemeris style.
    initEphComputer(pars, chooser, eph_style);
  }

  void PulsarToolApp::initEphComputer(const st_app::AppParGroup & pars, const EphChooser & chooser, const std::string & eph_style) {
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
        AbsoluteTime abs_epoch = parseTime(epoch_time_format, epoch_time_sys, epoch, epoch_time_format, epoch_time_sys,
          m_reference_header);

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

          if (0. >= f0) throw std::runtime_error("Frequency must be positive.");

          // Add the ephemeris the user provided.
          PulsarEphCont & ephemerides(m_computer->getPulsarEphCont());
          ephemerides.push_back(FrequencyEph(epoch_time_sys, abs_epoch, abs_epoch, abs_epoch, ra, dec, phi0, f0, f1, f2).clone());
        } else if (eph_style_uc == "PER") {
          double p0 = pars["p0"];
          double p1 = pars["p1"];
          double p2 = pars["p2"];

          if (0. >= p0) throw std::runtime_error("Period must be positive.");

          // Add the ephemeris the user provided.
          PulsarEphCont & ephemerides(m_computer->getPulsarEphCont());
          ephemerides.push_back(PeriodEph(epoch_time_sys, abs_epoch, abs_epoch, abs_epoch, ra, dec, phi0, p0, p1, p2).clone());
        } else {
          throw std::runtime_error("Unknown ephemeris style \"" + eph_style + "\" was specified.");
        }
      }
    }

    if (eph_style_uc == "DB" || m_tcmode_bin != SUPPRESSED) {
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

      // Select only ephemerides for this pulsar.
      std::string psr_name = pars["psrname"];
      database.filterName(psr_name);

      // Select only ephemerides with the solar system ephemeris used for barycentering.
      std::string match_solar_eph = pars["matchsolareph"];
      for (std::string::iterator itor = match_solar_eph.begin(); itor != match_solar_eph.end(); ++itor) *itor = std::toupper(*itor);
      if (match_solar_eph == "PSRDB" || match_solar_eph == "ALL") {
        std::string solar_eph = pars["solareph"];
        database.filterSolarEph(solar_eph);
      }

      // Load the selected ephemerides.
      if (eph_style_uc == "DB") m_computer->loadPulsarEph(database);
      if (m_tcmode_bin != SUPPRESSED) m_computer->loadOrbitalEph(database);
    }
  }

  AbsoluteTime PulsarToolApp::computeTimeBoundary(bool request_start_time, bool request_time_correction) {
    bool candidate_found = false;
    AbsoluteTime abs_candidate_time("TDB", 0, 0.);

    // First, look for requested time (start or stop) in the GTI.
    for (handler_cont_type::const_iterator itor = m_gti_handler_cont.begin(); itor != m_gti_handler_cont.end(); ++itor) {
      EventTimeHandler & gti_handler = *(*itor);

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
    if (!candidate_found) throw std::runtime_error("computeTimeBoundary: cannot determine start/stop of the observation interval");

    // Return the candidate.
    return abs_candidate_time;
  }

  void PulsarToolApp::initTimeCorrection(const st_app::AppParGroup & pars, bool guess_pdot) {
    // Read timeorigin parameter.
    std::string str_origin = pars["timeorigin"];

    // Initialize time correction with the given time origin.
    initTimeCorrection(pars, guess_pdot, str_origin);
  }

  void PulsarToolApp::initTimeCorrection(const st_app::AppParGroup & pars, bool guess_pdot, const std::string & str_origin) {
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

      std::string time_sys;
      (*m_reference_header)["TIMESYS"].get(time_sys);
      std::auto_ptr<TimeRep> time_rep(createMetRep(time_sys, abs_tstart));

      double tstart = 0.;
      double tstop = 0.;
      *time_rep = abs_tstart;
      time_rep->get("TIME", tstart);
      *time_rep = abs_tstop;
      time_rep->get("TIME", tstop);
      time_rep->set("TIME", .5 * (tstart + tstop));

      abs_origin = *time_rep;

    } else if (str_origin_uc == "USER") {
      // Get time of origin and its format and system from parameters.
      std::string origin_time = pars["usertime"];
      std::string origin_time_format = pars["userformat"];
      std::string origin_time_sys = pars["usersys"].Value();

      // Convert user-specified time into AbsoluteTime.
      std::string parsed_time_format;
      std::string parsed_time_sys;
      abs_origin = parseTime(origin_time_format, origin_time_sys, origin_time, parsed_time_format, parsed_time_sys,
        m_reference_header);

    } else {
      throw std::runtime_error("Unsupported time origin " + str_origin_uc);
    }

    // Initialize time correction with the time origin just computed.
    initTimeCorrection(pars, guess_pdot, abs_origin);
  }

  void PulsarToolApp::initTimeCorrection(const st_app::AppParGroup & pars, bool guess_pdot, const AbsoluteTime & abs_origin) {
    // Determine whether to perform binary demodulation.
    m_demod_bin = false;
    if ((m_tcmode_bin == REQUIRED) || (m_tcmode_bin == ALLOWED)) {
      // Check whether orbital parameters are available for binary demodulation.
      if (!m_computer->getOrbitalEphCont().empty()) {
        m_demod_bin = true;
      } else if (m_tcmode_bin == REQUIRED) {
        throw std::runtime_error("Binary demodulation was required by user, but no orbital ephemeris was found");
      }
    }

    // Determine whether to cancel pdot.
    m_cancel_pdot = false;
    if ((m_tcmode_pdot == REQUIRED) || (m_tcmode_pdot == ALLOWED)) {
      if (!guess_pdot || !m_computer->getPulsarEphCont().empty()) {
        m_cancel_pdot = true;
      } else if (m_tcmode_pdot == REQUIRED) {
        throw std::runtime_error("Pdot cancellation was required by user, but no spin ephemeris was found");
      }
    }

    // Determine whether to request barycentric correction.
    if (m_tcmode_bary == SUPPRESSED && (m_demod_bin || m_cancel_pdot)) {
      throw std::runtime_error("Barycentric correction is suppressed when binary demodulation or pdot cancellation is requested.");
    }
    m_request_bary = (m_tcmode_bary == REQUIRED || m_tcmode_bary == ALLOWED);

    // Initialize the time series to analyze.
    std::string target_time_sys;
    bool time_system_set = false;

    if (!m_request_bary && !m_demod_bin && !m_cancel_pdot) {
      // When NO corrections are requested, the analysis will be performed in the time system written in event files,
      // requiring all event files have same time system.
      std::string this_time_system;
      for (handler_cont_type::const_iterator itor = m_event_handler_cont.begin(); itor != m_event_handler_cont.end(); ++itor) {
        const EventTimeHandler & event_handler = *(*itor);
        const tip::Header & header(event_handler.getHeader());
        header["TIMESYS"].get(this_time_system);
        if (!time_system_set) {
          target_time_sys = this_time_system;
          time_system_set = true;
        } else if (this_time_system != target_time_sys) {
          throw std::runtime_error("event files with different TIMESYS values cannot be combined");
        }
      }
    } else {
      // When ANY correction(s) are requested, the analysis will be performed in TDB system.
      target_time_sys = "TDB";
      time_system_set = true;
    }

    // Check whether time system is successfully set.
    if (!time_system_set) throw std::runtime_error("cannot determine time system for the time series to analyze");

    // Set up target time representation, used to compute the time series to analyze.
    m_target_time_rep = createMetRep(target_time_sys, abs_origin);

    // Compute spin ephemeris to be used in pdot cancellation, and replace PulsarEph in EphComputer with it.
    if (m_cancel_pdot) {
      const int max_derivative = 2;
      if (guess_pdot) {
        // Compute an ephemeris at abs_origin to use for pdot cancellation, up to the second time derivative.
        m_computer->setPdotCancelParameter(abs_origin, max_derivative);

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
          m_computer->setPdotCancelParameter(target_time_sys, abs_origin, fdot_ratio);
        } else if (eph_style == "PER") {
          double p0 = 1.;
          double p1 = pars["p1p0ratio"];
          double p2 = pars["p2p0ratio"];
          m_computer->setPdotCancelParameter(abs_origin,
            PeriodEph(target_time_sys, abs_origin, abs_origin, abs_origin, ra, dec, phi0, p0, p1, p2), max_derivative);
        } else {
          throw std::runtime_error("Ephemeris style must be either FREQ or PER.");
        }
      }
    }

    // Initialize barycentric corrections.
    if (m_request_bary) {
      std::string solar_eph = pars["solareph"];
      BaryTimeComputer::getComputer().initialize(solar_eph);

      // Load RA and Dec to PulsarEph container for barycentric correction, if none is loaded.
      PulsarEphCont & ephemerides(m_computer->getPulsarEphCont());
      if (ephemerides.empty()) {
        double ra = pars["ra"];
        double dec = pars["dec"];
        
        // Set dummy phi0, f0, f1, and f2 (those are not used in barycentric correction.
        double phi0 = 0.;
        double f0 = 0.;
        double f1 = 0.;
        double f2 = 0.;
        ephemerides.push_back(FrequencyEph(target_time_sys, abs_origin, abs_origin, abs_origin, ra, dec, phi0, f0, f1, f2).clone());
      }
    }

  }

  double PulsarToolApp::computeElapsedSecond(const AbsoluteTime & abs_time) {
    double time_value = 0.;

    // Assign the absolute time to the time representation.
    *m_target_time_rep = abs_time;

    // Get value from the time representation.
    m_target_time_rep->get("TIME", time_value);

    return time_value;
  }

  AbsoluteTime PulsarToolApp::computeAbsoluteTime(double elapsed_time) {
    // Assign the elapsed time to the time representation.
    m_target_time_rep->set("TIME", elapsed_time);

    // Convert it to AbsoluteTime object and return it.
    return AbsoluteTime(*m_target_time_rep);
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
    m_event_handler_itor = m_event_handler_cont.begin();

    // Setup current event table.
    setupCurrentEventTable();
  }

  void PulsarToolApp::setNextEvent() {
    // Increment event record iterator.
    (*m_event_handler_itor)->setNextRecord();

    // Check whether it reaches to the end of table.
    if ((*m_event_handler_itor)->isEndOfTable()) {
      // Increment event table iterator.
      ++m_event_handler_itor;

      // Setup new event table.
      setupCurrentEventTable();
    }
  }

  bool PulsarToolApp::isEndOfEventList() const {
    return (m_event_handler_itor == m_event_handler_cont.end());
  }

  AbsoluteTime PulsarToolApp::getEventTime() {
    return readTimeColumn(**m_event_handler_itor, m_time_field, true);
  }

  AbsoluteTime PulsarToolApp::getStartTime() {
    return computeTimeBoundary(true, true);
  }

  AbsoluteTime PulsarToolApp::getStopTime() {
    return computeTimeBoundary(false, true);
  }

  EphComputer & PulsarToolApp::getEphComputer() const {
    return *m_computer;
  }

  void PulsarToolApp::resetApp() {
    // Destroy target TimeRep object.
    delete m_target_time_rep; m_target_time_rep = 0;

    // Reset time correction flags.
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

    // Reset reference header.
    m_reference_header = 0;

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
    bool request_time_correction) {

    // Read the original photon arrival time.
    AbsoluteTime abs_time = handler.readColumn(column_name);

    // Apply selected corrections, if requested.
    if (request_time_correction) {
      // Apply barycentric correction, if requested.
      if (m_request_bary) {
        // Get RA and Dec for the given arrival time.
        std::pair<double, double> ra_dec = m_computer->calcSkyPosition(abs_time);

        // Try barycentric correction with the RA and Dec.
        abs_time = handler.readColumn(column_name, ra_dec.first, ra_dec.second);
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
