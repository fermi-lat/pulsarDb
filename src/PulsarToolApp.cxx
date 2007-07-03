/** \file PulsarToolApp.cxx
    \brief Implementation of base class for pulsar tool applications.
    \author Masaharu Hirayama, GSSC
            James Peachey, HEASARC/GSSC
*/
#include <cctype>
#include <memory>
#include <stdexcept>

#include "pulsarDb/PulsarDb.h"
#include "pulsarDb/PulsarToolApp.h"

#include "st_app/AppParGroup.h"

#include "timeSystem/AbsoluteTime.h"
#include "timeSystem/GlastMetRep.h"
#include "timeSystem/TimeRep.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "st_facilities/Env.h"

using namespace timeSystem;

namespace pulsarDb {

  PulsarToolApp::PulsarToolApp(): m_time_field(""), m_gti_start_field(""), m_gti_stop_field(""), m_reference_header(0), m_computer(0),
    m_request_bary(false), m_demod_bin(false), m_cancel_pdot(false), m_target_time_rep(0) {}

  PulsarToolApp::~PulsarToolApp() throw() {
    if (m_computer) delete m_computer;
    if (m_target_time_rep) delete m_target_time_rep;
    for (std::map<const tip::Table *, TimeRep *>::iterator itor = m_time_rep_dict.begin(); itor != m_time_rep_dict.end(); ++itor) {
      delete itor->second;
    }
  }

  TimeRep * PulsarToolApp::createTimeRep(const std::string & time_format, const std::string & time_system,
    const std::string & time_value) {
    TimeRep * time_rep(0);

    // Make upper case copies of input for case insensitive comparisons.
    std::string time_format_uc(time_format);
    for (std::string::iterator itor = time_format_uc.begin(); itor != time_format_uc.end(); ++itor) *itor = std::toupper(*itor);

    std::string time_system_uc(time_system);
    for (std::string::iterator itor = time_system_uc.begin(); itor != time_system_uc.end(); ++itor) *itor = std::toupper(*itor);

    // Create representation for this time format and time system.
    if ("GLAST" == time_format_uc) {
      time_rep = new GlastMetRep(time_system, 0.);
    } else if ("MJD" == time_format_uc) {
      time_rep = new MjdRep(time_system, 0, 0.);
    } else {
      throw std::runtime_error("Time format \"" + time_format + "\" is not supported for ephemeris time");
    }

    // Assign the ephtime supplied by the user to the representation.
    time_rep->assign(time_value);

    return time_rep;
  }

  TimeRep * PulsarToolApp::createTimeRep(const std::string & time_format, const std::string & time_system,
    const std::string & time_value, const tip::Header & header) {
    TimeRep * time_rep(0);

    // Make upper case copies of input for case insensitive comparisons.
    std::string time_format_uc(time_format);
    for (std::string::iterator itor = time_format_uc.begin(); itor != time_format_uc.end(); ++itor) *itor = std::toupper(*itor);

    std::string time_system_uc(time_system);
    for (std::string::iterator itor = time_system_uc.begin(); itor != time_system_uc.end(); ++itor) *itor = std::toupper(*itor);

    // Make a local modifiable copy to hold the rationalized time system.
    std::string time_system_rat(time_system);

    // First check whether time system should be read from the tip::Header.
    if ("FILE" == time_system_uc) header["TIMESYS"].get(time_system_rat);

    // Create representation for this time format and time system.
    if ("FILE" == time_format_uc) {
      // Check TELESCOP keyword for supported missions.
      std::string telescope;
      header["TELESCOP"].get(telescope);
      for (std::string::iterator itor = telescope.begin(); itor != telescope.end(); ++itor) *itor = std::toupper(*itor);
      if (telescope != "GLAST") throw std::runtime_error("Only GLAST supported for now");

      // Get the mjdref from the header, which is not as simple as just reading a single keyword.
      MjdRefDatabase mjd_ref_db;
      IntFracPair mjd_ref(mjd_ref_db(header));
      time_rep = new MetRep(time_system_rat, mjd_ref, 0.);
    } else {
      // Delegate to overload that does not use tip.
      return createTimeRep(time_format, time_system_rat, time_value);
    }

    // Assign the time supplied by the user to the representation.
    time_rep->assign(time_value);

    return time_rep;
  }

  TimeRep * PulsarToolApp::createMetRep(const std::string & time_system, const AbsoluteTime & abs_reference) {
    // Compute MJD of abs_reference (the origin of the time series to analyze), to be given as MJDREF of MetRep.
    // NOTE: MetRep should take AbsoluteTime for its MJDREF (Need refactor of AbsoluteTime for that).
    // TODO: Once MetRep is refactored, remove this method.
    std::auto_ptr<TimeRep> mjd_rep(new MjdRep(time_system, 0, 0.));
    *mjd_rep = abs_reference;
    long mjd_int = 0;
    double mjd_frac = 0.;
    mjd_rep->get("MJDI", mjd_int);
    mjd_rep->get("MJDF", mjd_frac);

    // Create MetRep to represent the time series to analyze and return it.
    return new MetRep(time_system, mjd_int, mjd_frac, 0.);
  }

  void PulsarToolApp::openEventFile(const st_app::AppParGroup & pars) {
    std::string event_file = pars["evfile"];
    std::string event_extension = pars["evtable"];

    // Clear out any gti tables already in gti_table_cont.
    for (table_cont_type::reverse_iterator itor = m_gti_table_cont.rbegin(); itor != m_gti_table_cont.rend(); ++itor) {
      delete *itor;
    }
    m_gti_table_cont.clear();

    // Clear out any event tables already in event_table_cont.
    for (table_cont_type::reverse_iterator itor = m_event_table_cont.rbegin(); itor != m_event_table_cont.rend(); ++itor) {
      delete *itor;
    }
    m_event_table_cont.clear();

    // List all event tables and GTI tables.
    table_cont_type all_table_cont;

    // Open the event table.
    const tip::Table * event_table(tip::IFileSvc::instance().readTable(event_file, event_extension));

    // Add the table to the container.
    m_event_table_cont.push_back(event_table);
    all_table_cont.push_back(event_table);

    // Open the GTI table.
    const tip::Table * gti_table(tip::IFileSvc::instance().readTable(event_file, "GTI"));

    // Add the table to the container.
    m_gti_table_cont.push_back(gti_table);
    all_table_cont.push_back(gti_table);

    // Set names of TIME column in EVENTS exetension and START/STOP coulmns in GTI extensions.
    if (m_time_field.empty()) m_time_field = pars["timefield"].Value();
    m_gti_start_field = "START";
    m_gti_stop_field = "STOP";

    // Analyze all event tables and GTI tables, and set the results to internal variables.
    for (table_cont_type::const_iterator itor = all_table_cont.begin(); itor != all_table_cont.end(); ++itor) {
      const tip::Table * table = *itor;
      const tip::Header & header(table->getHeader());
      // Create and store TimeRep for this table.
      m_time_rep_dict[table] = createTimeRep("FILE", "FILE", "0.", header);

      // Check whether this table needs barycentering or not, and store it.
      std::string time_ref;
      header["TIMEREF"].get(time_ref);
      m_need_bary_dict[table] = ("SOLARSYSTEM" != time_ref);
    }

    // Select and set reference header.
    const tip::Table * reference_table = m_event_table_cont.at(0);
    m_reference_header = &(reference_table->getHeader());
  }

  void PulsarToolApp::initEphComputer(const st_app::AppParGroup & pars, const TimingModel & model,
    const EphChooser & chooser) {

    // Create ephemeris computer.
    m_computer = new EphComputer(model, chooser);

    std::string eph_style = pars["ephstyle"];
    for (std::string::iterator itor = eph_style.begin(); itor != eph_style.end(); ++itor) *itor = std::toupper(*itor);

    // Determine the time system used for the ephemeris epoch.
    std::string epoch_time_sys;
    if (eph_style == "DB") epoch_time_sys = "TDB";
    else epoch_time_sys = pars["timesys"].Value();

    // Ignored but needed for timing model.
    double phi0 = 0.;

    std::string psr_name = pars["psrname"];

    if (eph_style != "DB") {
      std::string epoch_time_format = pars["timeformat"];
      std::string epoch = pars["ephepoch"];
      std::auto_ptr<TimeRep> time_rep(createTimeRep(epoch_time_format, epoch_time_sys, epoch, *m_reference_header));
      AbsoluteTime abs_epoch(*time_rep);

      // Handle either period or frequency-style input.
      if (eph_style == "FREQ") {
        double f0 = pars["f0"];
        double f1 = pars["f1"];
        double f2 = pars["f2"];

        if (0. >= f0) throw std::runtime_error("Frequency must be positive.");

        // Override any ephemerides which may have been found in the database with the ephemeris the user provided.
        PulsarEphCont & ephemerides(m_computer->getPulsarEphCont());
        ephemerides.push_back(FrequencyEph(epoch_time_sys, abs_epoch, abs_epoch, abs_epoch, phi0, f0, f1, f2).clone());
      } else if (eph_style == "PER") {
        double p0 = pars["p0"];
        double p1 = pars["p1"];
        double p2 = pars["p2"];

        if (0. >= p0) throw std::runtime_error("Period must be positive.");

        // Override any ephemerides which may have been found in the database with the ephemeris the user provided.
        PulsarEphCont & ephemerides(m_computer->getPulsarEphCont());
        ephemerides.push_back(PeriodEph(epoch_time_sys, abs_epoch, abs_epoch, abs_epoch, phi0, p0, p1, p2).clone());
      }
    }

    std::string demod_bin_string = pars["demodbin"];
    for (std::string::iterator itor = demod_bin_string.begin(); itor != demod_bin_string.end(); ++itor) *itor = std::toupper(*itor);
    if (eph_style == "DB" || demod_bin_string != "NO") {
      // Find the pulsar database.
      std::string psrdb_file = pars["psrdbfile"];
      std::string psrdb_file_uc = psrdb_file;
      for (std::string::iterator itor = psrdb_file_uc.begin(); itor != psrdb_file_uc.end(); ++itor) *itor = std::toupper(*itor);
      if ("DEFAULT" == psrdb_file_uc) {
        using namespace st_facilities;
        psrdb_file = Env::appendFileName(Env::getDataDir("periodSearch"), "master_pulsardb.fits");
      }

      // Open the database.
      PulsarDb database(psrdb_file);

      // Select only ephemerides for this pulsar.
      database.filterName(psr_name);

      // Load the selected ephemerides.
      if (eph_style == "DB") m_computer->loadPulsarEph(database);
      m_computer->loadOrbitalEph(database);
    }
  }

  void PulsarToolApp::computeTimeBoundary(AbsoluteTime & abs_tstart, AbsoluteTime & abs_tstop) {
    bool candidate_found = false;

    // First, look for first and last times in the GTI.
    for (table_cont_type::const_iterator itor = m_gti_table_cont.begin(); itor != m_gti_table_cont.end(); ++itor) {
      const tip::Table & gti_table = *(*itor);

      // If possible, get tstart and tstop from first and last interval in GTI extension.
      tip::Table::ConstIterator gti_itor = gti_table.begin();
      if (gti_itor != gti_table.end()) {
        // Get start of the first interval and stop of last interval in GTI.
        AbsoluteTime abs_gti_start = readTimeColumn(gti_table, *gti_itor, m_gti_start_field);
        gti_itor = gti_table.end();
        --gti_itor;
        AbsoluteTime abs_gti_stop = readTimeColumn(gti_table, *gti_itor, m_gti_stop_field);

        if (candidate_found) {
          // See if current interval extends the observation.
          if (abs_gti_start < abs_tstart) abs_tstart = abs_gti_start;
          if (abs_gti_stop > abs_tstop) abs_tstop = abs_gti_stop;
        } else {
          // First candidate is the current interval.
          abs_tstart = abs_gti_start;
          abs_tstop = abs_gti_stop;
          candidate_found = true;
        }
      }
    }

    // In the unlikely event that there were no GTI files, no intervals in the GTI, and no event files, this is a
    // serious problem.
    if (!candidate_found) throw std::runtime_error("computeTimeBoundary: cannot determine start/stop of the observation interval");
  }

  AbsoluteTime PulsarToolApp::initTargetTime(const st_app::AppParGroup & pars, const AbsoluteTime & abs_tstart,
    const AbsoluteTime & abs_tstop) {
    std::string target_time_sys;
    bool time_system_set = false;

    // TODO: the following if-statement will become 'if (tcorrect == "NONE")' when tcorrect is introduced.
    if (!m_request_bary && !m_demod_bin && !m_cancel_pdot) {
      // When NO corrections are requested, the analysis will be performed in the time system written in event files,
      // requiring all event files have same time system.
      std::string this_time_system;
      for (table_cont_type::const_iterator itor = m_event_table_cont.begin(); itor != m_event_table_cont.end(); ++itor) {
        const tip::Table * event_table = *itor;
        const tip::Header & header(event_table->getHeader());
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

    // Handle styles of origin input.
    std::string origin_style = pars["timeorigin"];
    for (std::string::iterator itor = origin_style.begin(); itor != origin_style.end(); ++itor) *itor = std::toupper(*itor);
    AbsoluteTime abs_origin("TDB", Duration(0, 0.), Duration(0, 0.));
    if (origin_style == "START") {
      // Get time of origin and its time system from event file.
      abs_origin = abs_tstart;

    } else if (origin_style == "STOP") {
      // Get time of origin and its time system from event file.
      abs_origin = abs_tstop;

    } else if (origin_style == "MIDDLE") {
      // Use the center of the observation as the time origin.
      double tstart = computeElapsedSecond(abs_tstart);
      double tstop = computeElapsedSecond(abs_tstop);
      std::auto_ptr<TimeRep> time_rep(createMetRep(target_time_sys, abs_tstart));
      time_rep->set("TIME", .5 * (tstart + tstop));
      abs_origin = *time_rep;

    } else if (origin_style == "USER") {
      // Get time of origin and its format and system from parameters.
      std::string origin_time = pars["usertime"];
      std::string origin_time_format = pars["userformat"];
      std::string origin_time_sys = pars["usersys"].Value();

      // Convert user-specified time into AbsoluteTime.
      std::auto_ptr<TimeRep> time_rep(createTimeRep(origin_time_format, origin_time_sys, origin_time, *m_reference_header));
      abs_origin = *time_rep;

    } else {
      throw std::runtime_error("Unsupported origin style " + origin_style);
    }

    // Set up target time representation, used to compute the time series to analyze.
    m_target_time_rep = createMetRep(target_time_sys, abs_origin);

    return abs_origin;
  }

  PulsarEph & PulsarToolApp::updateEphComputer(const AbsoluteTime & abs_origin) {
    // Compute an ephemeris at abs_origin to use for the test.
    PulsarEph * eph(m_computer->calcPulsarEph(abs_origin).clone());

    // Reset computer to contain only the corrected ephemeris which was just computed.
    PulsarEphCont & ephemerides(m_computer->getPulsarEphCont());
    ephemerides.clear();
    ephemerides.push_back(eph);

    // Return reference to the computed ephemeris.
    return *eph;
  }

  void PulsarToolApp::initTimeCorrection(const st_app::AppParGroup & pars) {
    // Determine whether to request barycentric correction.
    // TODO: Read tcorrect parameter and set m_request_bary unless tcorrect == NONE.
    m_request_bary = false;

    std::string demod_bin_string = pars["demodbin"];
    for (std::string::iterator itor = demod_bin_string.begin(); itor != demod_bin_string.end(); ++itor) *itor = std::toupper(*itor);

    // Determine whether to perform binary demodulation.
    m_demod_bin = false;
    if (demod_bin_string != "NO") {
      // User selected not "no", so attempt to perform demodulation
      if (!m_computer->getOrbitalEphCont().empty()) {
        m_demod_bin = true;
      } else if (demod_bin_string == "YES") {
        throw std::runtime_error("Binary demodulation was required by user, but no orbital ephemeris was found");
      }
    }

    // Determine whether to cancel pdot.
    m_cancel_pdot = bool(pars["cancelpdot"]);

  }

  double PulsarToolApp::computeElapsedSecond(const AbsoluteTime & abs_time) {
    double time_value = 0.;

    // Assign the absolute time to the time representation.
    *m_target_time_rep = abs_time;

    // Get value from the time representation.
    m_target_time_rep->get("TIME", time_value);

    return time_value;
  }

  void PulsarToolApp::setFirstEvent() {
    // Set event table iterator.
    m_table_itor = m_event_table_cont.begin();

    // Setup current event table.
    if (m_table_itor != m_event_table_cont.end()) m_event_itor = (*m_table_itor)->begin();
  }

  void PulsarToolApp::setNextEvent() {
    // Increment event record iterator.
    ++m_event_itor;

    // Check whether it reaches to the end of table.
    if (m_event_itor == (*m_table_itor)->end()) {
      // Increment event table iterator.
      ++m_table_itor;

      // Setup new event table.
      if (m_table_itor != m_event_table_cont.end()) m_event_itor = (*m_table_itor)->begin();
    }
  }

  bool PulsarToolApp::isEndOfEventList() {
    return (m_table_itor == m_event_table_cont.end());
  }

  AbsoluteTime PulsarToolApp::getEventTime() {
    return AbsoluteTime(readTimeColumn(**m_table_itor, *m_event_itor, m_time_field));
  }

  AbsoluteTime PulsarToolApp::readTimeColumn(const tip::Table & table, tip::ConstTableRecord & record, const std::string & column_name) {
    // Get time value from given record.
    double time_value = record[column_name].get();

    // Get TimeRep for this table.
    TimeRep * time_rep(m_time_rep_dict[&table]);

    // Assign the value to the time representation.
    time_rep->set("TIME", time_value);

    // Convert TimeRep into AbsoluteTime so that computer can perform the necessary corrections.
    AbsoluteTime abs_time(*time_rep);

    // Apply selected corrections.
    bool correct_bary = m_request_bary && m_need_bary_dict[&table];
    if (correct_bary) throw std::runtime_error("Automatic barycentric correction not implemented.");
    if (m_demod_bin) m_computer->demodulateBinary(abs_time);
    if (m_cancel_pdot) m_computer->cancelPdot(abs_time);

    return abs_time;
  }
}
