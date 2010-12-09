/** \file EphComputerApp.cxx
    \brief Implementation of the EphComputerApp class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#include "pulsarDb/EphComputerApp.h"

#include "pulsarDb/EphChooser.h"
#include "pulsarDb/EphComputer.h"
#include "pulsarDb/EphStatus.h"

#include "st_app/AppParGroup.h"

#include "timeSystem/AbsoluteTime.h"
#include "timeSystem/CalendarFormat.h"
#include "timeSystem/Duration.h"
#include "timeSystem/ElapsedTime.h"
#include "timeSystem/IntFracUtility.h"
#include "timeSystem/MjdFormat.h"
#include "timeSystem/SourcePosition.h"
#include "timeSystem/TimeSystem.h"

#include <cctype>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>

static const std::string s_cvs_id("$Name:  $");

namespace pulsarDb {

  EphComputerApp::EphComputerApp(): m_os("EphComputerApp", "", 2) {
    setName("gtephem");
    setVersion(s_cvs_id);
  }

  EphComputerApp::~EphComputerApp() throw() {}

  void EphComputerApp::runApp() {
    m_os.setMethod("runApp()");

    using namespace st_app;
    using namespace st_stream;

    // Suppress 'INFO' in the prefix (cosmetic).
    m_os.info().setPrefix(m_os.out().getPrefix());

    // Get parameters.
    AppParGroup & pars(getParGroup());

    // Prompt and save.
    pars.Prompt();
    pars.Save();

    // Get parameters.
    std::string ref_time = pars["reftime"];
    std::string time_format = pars["timeformat"];
    std::string time_sys = pars["timesys"];
    bool strict = pars["strict"];

    // Handle leap seconds.
    std::string leap_sec_file = pars["leapsecfile"];
    timeSystem::TimeSystem::setDefaultLeapSecFileName(leap_sec_file);

    // Create the correct time representation for this time system and format,
    // and set the time of the representation to be the given reference time.
    std::string time_format_parsed;
    std::string time_sys_parsed;
    timeSystem::AbsoluteTime abs_ref_time = parseTime(time_format, time_sys, ref_time, time_format_parsed, time_sys_parsed);

    // Set up EphComputer for ephemeris computations.
    std::auto_ptr<EphChooser> chooser(0);
    if (strict) {
      chooser.reset(new StrictEphChooser);
    } else {
      chooser.reset(new SloppyEphChooser);
    }
    initEphComputer(pars, *chooser, "DB", m_os.info(4));
    EphComputer & computer(getEphComputer());

    // Set off the optional output.
    std::string dashes(26, '-');
    m_os.info(3) << prefix << dashes << std::endl;

    // Choose the best ephemeris for the given time.
    const PulsarEph * chosen_pulsar_eph(0);
    try {
      chosen_pulsar_eph = &(computer.choosePulsarEph(abs_ref_time));
    } catch (const std::exception &) {
      chosen_pulsar_eph = 0;
    }

    // Choose the best ephemeris for the given time.
    const OrbitalEph * chosen_orbital_eph(0);
    try {
      chosen_orbital_eph = &(computer.chooseOrbitalEph(abs_ref_time));
    } catch (const std::exception &) {
      chosen_orbital_eph = 0;
    }

    // Report the chosen spin and orbital ephemerides.
    std::string no_eph_available("   None");
    m_os.info(3) << prefix << "Spin ephemeris chosen from database is:" << std::endl;
    if (chosen_pulsar_eph) {
      m_os.info(3) << *chosen_pulsar_eph << std::endl;
    } else {
      m_os.info(3) << prefix << no_eph_available << std::endl;
    }
    m_os.info(3) << prefix << "Orbital ephemeris chosen from database is:" << std::endl;
    if (chosen_orbital_eph) {
      m_os.info(3) << *chosen_orbital_eph << std::endl;
    } else {
      m_os.info(3) << prefix << no_eph_available << std::endl;
    }

    // Set off the optional output.
    m_os.info(3) << prefix << dashes << std::endl;

    // Compose a string expression of the given time.
    std::string ref_time_string;
    try {
      ref_time_string = abs_ref_time.represent(time_sys_parsed, timeSystem::MjdFmt);
    } catch (const std::exception &) {
      ref_time_string = abs_ref_time.represent(time_sys_parsed, timeSystem::CalendarFmt);
    }

    // Calculate spin ephemeris for the given reference time, and report the result.
    m_os.out() << prefix << "Spin ephemeris estimated is:" << std::endl;
    m_os.out().prefix().width(30); m_os.out() << "Reference Time : " << ref_time_string << std::endl;
    if (0 == chosen_pulsar_eph) {
      // Report no spin ephemeris is found.
      m_os.out().prefix().width(30); m_os.out() << "Ephemeris Data : " << "Not Available" << std::endl;

    } else {
      // Compute extrapolated ephemeris.
      timeSystem::SourcePosition src_pos(0., 0.);
      double phi0 = 0.;
      double f0 = 0.;
      double f1 = 0.;
      double f2 = 0.;
      bool computed_ok = false;
      try {
        src_pos = chosen_pulsar_eph->calcPosition(abs_ref_time);
        phi0 = chosen_pulsar_eph->calcPulsePhase(abs_ref_time);
        f0 = chosen_pulsar_eph->calcFrequency(abs_ref_time, 0);
        f1 = chosen_pulsar_eph->calcFrequency(abs_ref_time, 1);
        f2 = chosen_pulsar_eph->calcFrequency(abs_ref_time, 2);
        computed_ok = true;
      } catch (const std::exception & x) {
        m_os.err() << prefix << "Unexpected error occurred in ephemeris computations." << std::endl << x.what() << std::endl;
      }

      // Compute RA and Dec.
      static const double degree_per_radian = 180. / M_PI;
      const std::vector<double> & src_dir = src_pos.getDirection();
      double dec = std::asin(src_dir[2]) * degree_per_radian;
      double ra = 0.;
      if (src_dir[0] != 0. || src_dir[1] != 0.) ra = std::atan2(src_dir[1], src_dir[0]) * degree_per_radian;
      if (ra < 0.) ra += 360.;

      // Print computed ephemeris.
      if (computed_ok) {
        m_os.out().precision(std::numeric_limits<double>::digits10);
        m_os.out().prefix().width(30); m_os.out() << "Right Ascension (degree) : " << ra << std::endl;
        m_os.out().prefix().width(30); m_os.out() << "Declination (degree) : " << dec << std::endl;
        m_os.out().prefix().width(30); m_os.out() << "Distance (light-second) : ";
        if (src_pos.hasDistance()) m_os.out() << src_pos.getDistance();
        else m_os.out() << "Unknown";
        m_os.out() << std::endl;
        m_os.out().prefix().width(30); m_os.out() << "Pulse Phase : " << phi0 << std::endl;
        m_os.out().prefix().width(30); m_os.out() << "Pulse Frequency (Hz) : " << f0 << std::endl;
        m_os.out().prefix().width(30); m_os.out() << "1st Derivative (Hz/s) : " << f1 << std::endl;
        m_os.out().prefix().width(30); m_os.out() << "2nd Derivative (Hz/s/s) : " << f2 << std::endl;
      }
    }


    // Calculate the orbital position of the pulsar for the given reference time, and report the result.
    if (chosen_orbital_eph) {
      m_os.out() << prefix << "Orbital position estimated is:" << std::endl;
      m_os.out().prefix().width(30); m_os.out() << "Reference Time : " << ref_time_string << std::endl;

      // Compute the orbital position and related quantities.
      double orbital_phase = 0.;
      timeSystem::ElapsedTime orbital_delay("TDB", timeSystem::Duration(0., "Sec"));
      timeSystem::AbsoluteTime abs_mod_time = abs_ref_time;
      bool computed_ok = false;
      try {
        orbital_phase = chosen_orbital_eph->calcOrbitalPhase(abs_ref_time);
        orbital_delay = chosen_orbital_eph->calcOrbitalDelay(abs_ref_time);
        chosen_orbital_eph->modulateBinary(abs_mod_time);
        computed_ok = true;
      } catch (const std::exception & x) {
        m_os.err() << prefix << "Unexpected error occurred in computing the orbital position." << std::endl << x.what() << std::endl;
      }

      // Compose a string expression of the modulated time.
      std::string mod_time_string;
      try {
        mod_time_string = abs_mod_time.represent(time_sys_parsed, timeSystem::MjdFmt);
      } catch (const std::exception &) {
        mod_time_string = abs_mod_time.represent(time_sys_parsed, timeSystem::CalendarFmt);
      }

      // Parse the orbital delay and compose a string expression of it.
      const std::string orbital_delay_unit = orbital_delay.getSystem().getName() + " second";
      long orbital_delay_int = 0;
      double orbital_delay_frac = 0.;
      orbital_delay.getDuration("Sec", orbital_delay_int, orbital_delay_frac);
      const timeSystem::IntFracUtility & utility(timeSystem::IntFracUtility::getUtility());
      std::string orbital_delay_second = utility.format(orbital_delay_int, orbital_delay_frac);

      // Print computed position.
      if (computed_ok) {
        m_os.out().precision(std::numeric_limits<double>::digits10);
        m_os.out().prefix().width(30); m_os.out() << "Orbital Phase : " << orbital_phase << std::endl;
        m_os.out().prefix().width(30); m_os.out() << "Orbital Delay (" + orbital_delay_unit + ") : " << orbital_delay_second <<
          std::endl;
        m_os.out().prefix().width(30); m_os.out() << "Barycenter Arrival Time : " << mod_time_string << std::endl;
      }
    }

    // Report ephemeris status between the given reference time and the reference epoch of the chosen pulsar ephemeris.
    if (chosen_pulsar_eph) {
      std::set<EphStatusCodeType> code_to_report;
      code_to_report.insert(Remarked);
      reportEphStatus(m_os.warn(), abs_ref_time, chosen_pulsar_eph->getEpoch(), code_to_report);
    }
  }

}
