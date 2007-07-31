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

  EphComputer2::EphComputer2(): EphComputer(), m_pdot_pars(0), m_model2(new TimingModel) {}

  EphComputer2::EphComputer2(const TimingModel & model, const EphChooser & chooser): EphComputer(model, chooser), m_pdot_pars(0),
    m_model2(model.clone()) {}

  EphComputer2::~EphComputer2() {
    if (m_pdot_pars) delete m_pdot_pars;
    delete m_model2;
  }

  void EphComputer2::setPdotCancelParameter(const PulsarEph & pdot_pars) {
    m_pdot_pars = pdot_pars.clone();
  }

  void EphComputer2::cancelPdot(timeSystem::AbsoluteTime & ev_time) const {
    if (m_pdot_pars) {
      m_model2->cancelPdot(*m_pdot_pars, ev_time);
    } else {
      throw std::runtime_error("Parameters for pdot cancellation are not set");
    }
  }

  PulsarToolApp::PulsarToolApp(): m_time_field(""), m_gti_start_field(""), m_gti_stop_field(""), m_reference_header(0), m_computer(0),
    m_tcmode_bary(ALLOWED), m_tcmode_bin(ALLOWED), m_tcmode_pdot(ALLOWED),
    m_request_bary(false), m_demod_bin(false), m_cancel_pdot(false), m_target_time_rep(0) {}

  PulsarToolApp::~PulsarToolApp() throw() {
    if (m_computer) delete m_computer;
    if (m_target_time_rep) delete m_target_time_rep;
    for (std::map<const tip::Table *, TimeRep *>::iterator itor = m_time_rep_dict.begin(); itor != m_time_rep_dict.end(); ++itor) {
      delete itor->second;
    }
    for (table_cont_type::reverse_iterator itor = m_gti_table_cont.rbegin(); itor != m_gti_table_cont.rend(); ++itor) {
      delete *itor;
    }
    for (table_cont_type::reverse_iterator itor = m_event_table_cont.rbegin(); itor != m_event_table_cont.rend(); ++itor) {
      delete *itor;
    }
  }

  TimeRep * PulsarToolApp::createTimeRep(const std::string & time_format, const std::string & time_system,
    const std::string & time_value) const {
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
    const std::string & time_value, const tip::Header & header) const {
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

  TimeRep * PulsarToolApp::createMetRep(const std::string & time_system, const AbsoluteTime & abs_reference) const {
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

  void PulsarToolApp::openEventFile(const st_app::AppParGroup & pars, bool read_only) {
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

    // Open the event table, either for reading or reading and writing.
    // Note: for convenience, read-only and read-write tables are stored as const Table pointers
    // in the container. In the few cases where writing to the tables is necessary, a const_cast will
    // be necessary.
    const tip::Table * event_table = 0;
    if (read_only) event_table = tip::IFileSvc::instance().readTable(event_file, event_extension);
    else event_table = tip::IFileSvc::instance().editTable(event_file, event_extension);

    // Add the table to the container.
    m_event_table_cont.push_back(event_table);
    all_table_cont.push_back(event_table);

    // Open the GTI table.
    // Note: At present, GTI is never modified, so no need to open it read-write.
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

// TODO: Uncomment the line below when tcorrect parameter replaces cancelpdot and demodbin parameters.
//#define tcorrect_parameter_supported
#ifdef tcorrect_parameter_supported
  void PulsarToolApp::selectTimeCorrectionMode(const st_app::AppParGroup & pars) {
    // Read tcorrect parameter.
    std::string t_correct = pars["tcorrect"];

    // Determine time correction mode.
    selectTimeCorrectionMode(t_correct);
  }
#else
  void PulsarToolApp::selectTimeCorrectionMode(const st_app::AppParGroup & pars) {
    // Determine time correction mode for barycentric correction.
    // TODO: Change this when barycentering-on-the-fly is implemented.
    // NOTE: Barycentering-on-the-fly is NOT implemented.
    m_tcmode_bary = SUPPRESSED;

    // Determine time correction mode for binary demodulation.
    std::string demod_bin_string = pars["demodbin"];
    for (std::string::iterator itor = demod_bin_string.begin(); itor != demod_bin_string.end(); ++itor) *itor = std::toupper(*itor);
    if (demod_bin_string == "YES") {
      m_tcmode_bin = REQUIRED;
    } else if (demod_bin_string == "NO") {
      m_tcmode_bin = SUPPRESSED;
    } else if (demod_bin_string == "AUTO") {
      m_tcmode_bin = ALLOWED;
    }

    // Determine time correction mode for pdot cancellation.
    try {
      bool cancel_pdot_bool = pars["cancelpdot"];
      if (cancel_pdot_bool) {
        m_tcmode_pdot = REQUIRED;
      } else {
        m_tcmode_pdot = SUPPRESSED;
      }
    } catch (...) {
      m_tcmode_pdot = SUPPRESSED;
    }
  }
#endif

  void PulsarToolApp::initEphComputer(const st_app::AppParGroup & pars, const TimingModel & model,
    const EphChooser & chooser) {
    // Read ephstyle parameter.
    std::string eph_style = pars["ephstyle"];

    // Initialize EphComputer with given ephemeris style.
    initEphComputer(pars, model, chooser, eph_style);
  }

  void PulsarToolApp::initEphComputer(const st_app::AppParGroup & pars, const TimingModel & model,
    const EphChooser & chooser, const std::string & eph_style) {
    // Create ephemeris computer.
    m_computer = new EphComputer2(model, chooser);

    // Make eph_style argument case-insensitive.
    std::string eph_style_uc = eph_style;
    for (std::string::iterator itor = eph_style_uc.begin(); itor != eph_style_uc.end(); ++itor) *itor = std::toupper(*itor);

    // Determine the time system used for the ephemeris epoch.
    std::string epoch_time_sys;
    if (eph_style_uc == "DB") epoch_time_sys = "TDB";
    else epoch_time_sys = pars["timesys"].Value();

    // Ignored but needed for timing model.
    double phi0 = 0.;

    std::string psr_name = pars["psrname"];

    if (eph_style_uc != "DB") {
      std::string epoch_time_format = pars["timeformat"];
      std::string epoch = pars["ephepoch"];
      std::auto_ptr<TimeRep> time_rep(createTimeRep(epoch_time_format, epoch_time_sys, epoch, *m_reference_header));
      AbsoluteTime abs_epoch(*time_rep);

      // Handle either period or frequency-style input.
      if (eph_style_uc == "FREQ") {
        double f0 = pars["f0"];
        double f1 = pars["f1"];
        double f2 = pars["f2"];

        if (0. >= f0) throw std::runtime_error("Frequency must be positive.");

        // Override any ephemerides which may have been found in the database with the ephemeris the user provided.
        PulsarEphCont & ephemerides(m_computer->getPulsarEphCont());
        ephemerides.push_back(FrequencyEph(epoch_time_sys, abs_epoch, abs_epoch, abs_epoch, phi0, f0, f1, f2).clone());
      } else if (eph_style_uc == "PER") {
        double p0 = pars["p0"];
        double p1 = pars["p1"];
        double p2 = pars["p2"];

        if (0. >= p0) throw std::runtime_error("Period must be positive.");

        // Override any ephemerides which may have been found in the database with the ephemeris the user provided.
        PulsarEphCont & ephemerides(m_computer->getPulsarEphCont());
        ephemerides.push_back(PeriodEph(epoch_time_sys, abs_epoch, abs_epoch, abs_epoch, phi0, p0, p1, p2).clone());
      } else {
        throw std::runtime_error("Ephemeris style must be either FREQ or PER.");
      }
    }

    if (eph_style_uc == "DB" || m_tcmode_bin != SUPPRESSED) {
      // Find the pulsar database.
      std::string psrdb_file = pars["psrdbfile"];
      std::string psrdb_file_uc = psrdb_file;
      for (std::string::iterator itor = psrdb_file_uc.begin(); itor != psrdb_file_uc.end(); ++itor) *itor = std::toupper(*itor);
      if ("DEFAULT" == psrdb_file_uc) {
        using namespace st_facilities;
        // TODO: Change the folowing directory/file names of master pulsar database to the "official" one.
        psrdb_file = Env::appendFileName(Env::getDataDir("pulsarDb"), "master_pulsardb.fits");
      }

      // Open the database.
      PulsarDb database(psrdb_file);

      // Select only ephemerides for this pulsar.
      database.filterName(psr_name);

      // Load the selected ephemerides.
      if (eph_style_uc == "DB") m_computer->loadPulsarEph(database);
      m_computer->loadOrbitalEph(database);
    }
  }

  timeSystem::AbsoluteTime PulsarToolApp::computeTimeBoundary(bool request_start_time, bool request_time_correction) {
    bool candidate_found = false;
    AbsoluteTime abs_candidate_time("TDB", Duration(0, 0.), Duration(0, 0.));

    // First, look for requested time (start or stop) in the GTI.
    for (table_cont_type::const_iterator itor = m_gti_table_cont.begin(); itor != m_gti_table_cont.end(); ++itor) {
      const tip::Table & gti_table = *(*itor);

      // If possible, get tstart (or tstop) from first (or last) interval in GTI extension.
      tip::Table::ConstIterator gti_itor;
      tip::Table::ConstIterator gti_begin = gti_table.begin();
      tip::Table::ConstIterator gti_end = gti_table.end();
      std::string field_name;
      if (gti_begin != gti_end) {
        if (request_start_time) {
          // Get start of the first interval in GTI.
          gti_itor = gti_begin;
          field_name = m_gti_start_field;
        } else {
          // Get stop of the last interval in GTI.
          gti_itor = gti_end;
          --gti_itor;
          field_name = m_gti_stop_field;
        }

        // Read GTI START (or STOP) column value as AbsoluteTime.
        AbsoluteTime abs_gti_time = readTimeColumn(gti_table, *gti_itor, field_name, request_time_correction);

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
    // Determine whether to request barycentric correction.
    // TODO: Change this when barycentering-on-the-fly is implemented.
    m_request_bary = false;

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

    // Initialize the time series to analyze.
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

    // Set up target time representation, used to compute the time series to analyze.
    AbsoluteTime abs_origin = getTimeOrigin(pars);
    m_target_time_rep = createMetRep(target_time_sys, abs_origin);

    // Compute spin ephemeris to be used in pdot cancellation, and replace PulsarEph in EphComputer with it.
    if (m_cancel_pdot) {
      if (guess_pdot) {
        // Compute an ephemeris at abs_origin to use for pdot cancellation.
        m_computer->setPdotCancelParameter(m_computer->calcPulsarEph(abs_origin));
      } else {
        // Read parameters for pdot cancellation from pfile.
        std::string eph_style = pars["ephstyle"];
        for (std::string::iterator itor = eph_style.begin(); itor != eph_style.end(); ++itor) *itor = std::toupper(*itor);
        double phi0 = 0.;
        if (eph_style == "FREQ") {
          double f0 = 1.;
          double f1 = pars["f1f0ratio"];;
          double f2 = pars["f2f0ratio"];
          m_computer->setPdotCancelParameter(FrequencyEph(target_time_sys, abs_origin, abs_origin, abs_origin, phi0, f0, f1, f2));
        } else if (eph_style == "PER") {
          double p0 = 1.;
          double p1 = pars["p1p0ratio"];
          double p2 = pars["p2p0ratio"];
          m_computer->setPdotCancelParameter(PeriodEph(target_time_sys, abs_origin, abs_origin, abs_origin, phi0, p0, p1, p2));
        } else {
          throw std::runtime_error("Ephemeris style must be either FREQ or PER.");
        }
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

  void PulsarToolApp::setupEventTable(const tip::Table & table) {
    m_event_itor = table.begin();

    // Add output field if missing.
    for (std::vector<std::pair<std::string, std::string> >::const_iterator itor = m_output_field_cont.begin();
         itor != m_output_field_cont.end(); ++itor) {
      std::string field_name = itor->first;
      std::string field_format = itor->second;
      try {
        table.getFieldIndex(field_name);
      } catch (const tip::TipException &) {
        // Container stores tables as "const Table *" whether the table was really opened read-only or not,
        // so need to convert back to non-const. Note that if the file actually is opened read-only,
        // appendField will fail.
        tip::Table & table_nc = const_cast<tip::Table &>(table);
        table_nc.appendField(field_name, field_format);
      }
    }
  }

  void PulsarToolApp::setFirstEvent() {
    // Set event table iterator.
    m_table_itor = m_event_table_cont.begin();

    // Setup current event table.
    if (m_table_itor != m_event_table_cont.end()) setupEventTable(**m_table_itor);
  }

  void PulsarToolApp::setNextEvent() {
    // Increment event record iterator.
    ++m_event_itor;

    // Check whether it reaches to the end of table.
    if (m_event_itor == (*m_table_itor)->end()) {
      // Increment event table iterator.
      ++m_table_itor;

      // Setup new event table.
      if (m_table_itor != m_event_table_cont.end()) setupEventTable(**m_table_itor);
    }
  }

  bool PulsarToolApp::isEndOfEventList() const {
    return (m_table_itor == m_event_table_cont.end());
  }

  AbsoluteTime PulsarToolApp::getEventTime() {
    return readTimeColumn(**m_table_itor, *m_event_itor, m_time_field, true);
  }

  void PulsarToolApp::setFieldValue(const std::string & field_name, double field_value) {
    // Container stores tables as "const Table *" whether the table was really opened read-only or not,
    // so need to convert back to non-const. Note that if the file actually is opened read-only,
    // TableRecord::set method will fail.
    // TODO: Is there some other way to do this? The static_cast gives us the creeps!
    tip::TableRecord & record = static_cast<tip::TableRecord &>(*m_event_itor);
    record[field_name].set(field_value);
  }

  AbsoluteTime PulsarToolApp::getStartTime() {
    return computeTimeBoundary(true, true);
  }

  AbsoluteTime PulsarToolApp::getStopTime() {
    return computeTimeBoundary(false, true);
  }

  AbsoluteTime PulsarToolApp::getTimeOrigin(const st_app::AppParGroup & pars) {
    AbsoluteTime abs_origin("TDB", Duration(0, 0.), Duration(0, 0.));

    // Handle styles of origin input.
    std::string origin_style = pars["timeorigin"];
    for (std::string::iterator itor = origin_style.begin(); itor != origin_style.end(); ++itor) *itor = std::toupper(*itor);
    if (origin_style == "START") {
      // Get the uncorrected start time of event list.
      abs_origin = computeTimeBoundary(true, false);

    } else if (origin_style == "STOP") {
      // Get the uncorrected stop time of event list.
      abs_origin = computeTimeBoundary(false, false);

    } else if (origin_style == "MIDDLE") {
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

    return abs_origin;
  }

  EphComputer2 & PulsarToolApp::getEphComputer() const {
    return *m_computer;
  }

  AbsoluteTime PulsarToolApp::readTimeColumn(const tip::Table & table, tip::ConstTableRecord & record, const std::string & column_name,
    bool request_time_correction) {
    // Get time value from given record.
    double time_value = record[column_name].get();

    // Get TimeRep for this table.
    TimeRep * time_rep(m_time_rep_dict[&table]);

    // Assign the value to the time representation.
    time_rep->set("TIME", time_value);

    // Convert TimeRep into AbsoluteTime so that computer can perform the necessary corrections.
    AbsoluteTime abs_time(*time_rep);

    // Apply selected corrections, if requested.
    if (request_time_correction) {
      bool correct_bary = m_request_bary && m_need_bary_dict[&table];
      if (correct_bary) throw std::runtime_error("Automatic barycentric correction not implemented.");
      if (m_demod_bin) m_computer->demodulateBinary(abs_time);
      if (m_cancel_pdot) m_computer->cancelPdot(abs_time);
    }

    return abs_time;
  }

}
