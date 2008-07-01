/** \file EphComputerApp.cxx
    \brief Implementation of the EphComputerApp class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#include "pulsarDb/EphChooser.h"
#include "pulsarDb/EphComputer.h"
#include "pulsarDb/EphComputerApp.h"

#include "st_app/AppParGroup.h"

#include "timeSystem/AbsoluteTime.h"
#include "timeSystem/CalendarFormat.h"
#include "timeSystem/MjdFormat.h"
#include "timeSystem/TimeSystem.h"

#include <cctype>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

static const std::string s_cvs_id("$Name:  $");

using namespace timeSystem;

namespace pulsarDb {

  EphComputerApp::EphComputerApp(): m_os("EphComputerApp", "", 2) {
    setName("gtephem");
    setVersion(s_cvs_id);
  }

  EphComputerApp::~EphComputerApp() throw() {}

  void EphComputerApp::run() {
    // Clean up from any previous runs.
    resetApp();

    using namespace st_app;
    using namespace st_stream;

    m_os.setMethod("run()");

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
    AbsoluteTime abs_ref_time = parseTime(time_format, time_sys, ref_time, time_format_parsed, time_sys_parsed);

    // Set up EphComputer for ephemeris computations.
    std::auto_ptr<EphChooser> chooser(0);
    if (strict) {
      chooser.reset(new StrictEphChooser);
    } else {
      chooser.reset(new SloppyEphChooser);
    }
    initEphComputer(pars, *chooser, "DB");
    EphComputer & computer(getEphComputer());

    // Cosmetic: suppress info.
    m_os.info().setPrefix(m_os.out().getPrefix());

    // Set off the optional output.
    std::string dashes(26, '-');
    m_os.info(3) << prefix << dashes << std::endl;

    // Report the best spin ephemeris.
    bool found_pulsar_eph = false;
    try {
      const PulsarEph & eph(computer.choosePulsarEph(abs_ref_time));
      m_os.info(3) << prefix << "Spin ephemeris chosen from database is:" << std::endl << eph << std::endl;
      found_pulsar_eph = true;
    } catch (const std::exception & x) {
      m_os.out() << prefix << "No spin ephemeris was found in database:" << std::endl << x.what() << std::endl;
    }

    // Report the best binary ephemeris.
    try {
      const OrbitalEph & eph(computer.chooseOrbitalEph(abs_ref_time));
      m_os.info(3) << prefix << "Orbital ephemeris chosen from database is:" << std::endl << eph << std::endl;
    } catch (const std::exception & x) {
      m_os.info(3) << prefix << "No orbital ephemeris was found in database:" << std::endl << x.what() << std::endl;
    }

    // Set off the optional output.
    m_os.info(3) << prefix << dashes << std::endl;

    // Report the calculated spin ephemeris, provided at least a spin ephemeris was found above.
    if (found_pulsar_eph) {
      // Choose the best ephemeris for the given time.
      const PulsarEph & chosen_eph = computer.choosePulsarEph(abs_ref_time);

      // Compute extrapolated ephemeris.
      std::pair<double, double> ra_dec;
      double phi0 = 0.;
      double f0 = 0.;
      double f1 = 0.;
      double f2 = 0.;
      bool computed_ok = false;
      try {
        ra_dec = chosen_eph.calcSkyPosition(abs_ref_time);
        phi0 = chosen_eph.calcPulsePhase(abs_ref_time);
        f0 = chosen_eph.calcFrequency(abs_ref_time, 0);
        f1 = chosen_eph.calcFrequency(abs_ref_time, 1);
        f2 = chosen_eph.calcFrequency(abs_ref_time, 2);
        computed_ok = true;
      } catch (const std::exception & x) {
        m_os.err() << prefix << "Unexpected problem computing ephemeris." << std::endl << x.what() << std::endl;
      }

      // Compose a string expression of the given time.
      std::string time_string;
      try {
        time_string = abs_ref_time.represent(time_sys_parsed, MjdFmt);
      } catch (const std::exception &) {
        time_string = abs_ref_time.represent(time_sys_parsed, CalendarFmt);
      }

      // Print computed ephemeris.
      if (computed_ok) {
        m_os.out() << prefix << "Spin ephemeris estimated is:" << std::endl;
        m_os.out().precision(std::numeric_limits<double>::digits10);
        m_os.out().prefix().width(30); m_os.out() << "Reference Time : " << time_string << std::endl;
        m_os.out().prefix().width(30); m_os.out() << "Right Ascension (degree) : " << ra_dec.first << std::endl;
        m_os.out().prefix().width(30); m_os.out() << "Declination (degree) : " << ra_dec.second << std::endl;
        m_os.out().prefix().width(30); m_os.out() << "Pulse Phase : " << phi0 << std::endl;
        m_os.out().prefix().width(30); m_os.out() << "Pulse Frequency (Hz) : " << f0 << std::endl;
        m_os.out().prefix().width(30); m_os.out() << "1st Derivative (Hz/s) : " << f1 << std::endl;
        m_os.out().prefix().width(30); m_os.out() << "2nd Derivative (Hz/s/s) : " << f2 << std::endl;
      }
    }
  }
}
