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

using namespace timeSystem;

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
    AbsoluteTime abs_ref_time = parseTime(time_format, time_sys, ref_time, time_format_parsed, time_sys_parsed);

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

    // Report the best spin ephemeris.
    bool found_pulsar_eph = false;
    m_os.info(3) << prefix << "Spin ephemeris chosen from database is:" << std::endl;
    try {
      const PulsarEph & eph(computer.choosePulsarEph(abs_ref_time));
      m_os.info(3) << eph << std::endl;
      found_pulsar_eph = true;
    } catch (const std::exception &) {
      m_os.info(3) << prefix << "   None" << std::endl;
    }

    // Report the best binary ephemeris.
    m_os.info(3) << prefix << "Orbital ephemeris chosen from database is:" << std::endl;
    try {
      const OrbitalEph & eph(computer.chooseOrbitalEph(abs_ref_time));
      m_os.info(3) << eph << std::endl;
    } catch (const std::exception &) {
      m_os.info(3) << prefix << "   None" << std::endl;
    }

    // Set off the optional output.
    m_os.info(3) << prefix << dashes << std::endl;

    // Compose a string expression of the given time.
    std::string time_string;
    try {
      time_string = abs_ref_time.represent(time_sys_parsed, MjdFmt);
    } catch (const std::exception &) {
      time_string = abs_ref_time.represent(time_sys_parsed, CalendarFmt);
    }

    // Calculate spin ephemeris for the given reference time, if at least a spin ephemeris was found above, and report the result.
    m_os.out() << prefix << "Spin ephemeris estimated is:" << std::endl;
    m_os.out().prefix().width(30); m_os.out() << "Reference Time : " << time_string << std::endl;
    if (!found_pulsar_eph) {
      // Report no spin ephemeris is found.
      m_os.out().prefix().width(30); m_os.out() << "Ephemeris Data : " << "Not Available" << std::endl;

    } else {
      // Choose the best ephemeris for the given time.
      const PulsarEph & chosen_eph = computer.choosePulsarEph(abs_ref_time);

      // Compute extrapolated ephemeris.
      SourcePosition src_pos(0., 0.);
      double phi0 = 0.;
      double f0 = 0.;
      double f1 = 0.;
      double f2 = 0.;
      bool computed_ok = false;
      try {
        src_pos = chosen_eph.calcPosition(abs_ref_time);
        phi0 = chosen_eph.calcPulsePhase(abs_ref_time);
        f0 = chosen_eph.calcFrequency(abs_ref_time, 0);
        f1 = chosen_eph.calcFrequency(abs_ref_time, 1);
        f2 = chosen_eph.calcFrequency(abs_ref_time, 2);
        computed_ok = true;
      } catch (const std::exception & x) {
        m_os.err() << prefix << "Unexpected problem computing ephemeris." << std::endl << x.what() << std::endl;
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
#if 0
        // TODO: Uncomment this section once an ephemeris with a finite distance is implemented.
        m_os.out().prefix().width(30); m_os.out() << "Distance (light-second) : ";
        if (src_pos.hasDistance()) m_os.out() << src_pos.getDistance();
        else m_os.out() << "Unknown";
        m_os.out() << std::endl;
#endif
        m_os.out().prefix().width(30); m_os.out() << "Pulse Phase : " << phi0 << std::endl;
        m_os.out().prefix().width(30); m_os.out() << "Pulse Frequency (Hz) : " << f0 << std::endl;
        m_os.out().prefix().width(30); m_os.out() << "1st Derivative (Hz/s) : " << f1 << std::endl;
        m_os.out().prefix().width(30); m_os.out() << "2nd Derivative (Hz/s/s) : " << f2 << std::endl;
      }

      // Report ephemeris status between reference time and reference epoch of chosen ephemeris.
      std::set<EphStatusCodeType> code_to_report;
      code_to_report.insert(Remarked);
      reportEphStatus(m_os.warn(), abs_ref_time, chosen_eph.getEpoch(), code_to_report);
    }
  }
}
