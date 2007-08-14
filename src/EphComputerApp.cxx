/** \file EphComputerApp.cxx
    \brief Implementation of the EphComputerApp class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#include "pulsarDb/EphChooser.h"
#include "pulsarDb/EphComputer.h"
#include "pulsarDb/EphComputerApp.h"
#include "pulsarDb/PulsarDb.h"
#include "pulsarDb/TimingModel.h"

#include "st_app/AppParGroup.h"

#include "st_facilities/Env.h"

#include "timeSystem/AbsoluteTime.h"
#include "timeSystem/GlastMetRep.h"
#include "timeSystem/TimeRep.h"
#include "timeSystem/TimeSystem.h"

#include <cctype>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

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
    std::auto_ptr<TimeRep> time_rep(createTimeRep(time_format, time_sys, ref_time));
    AbsoluteTime abs_ref_time(*time_rep);

    // Set up EphComputer for ephemeris computations.
    TimingModel model;
    std::auto_ptr<EphChooser> chooser(0);
    if (strict) {
      chooser.reset(new StrictEphChooser);
    } else {
      chooser.reset(new SloppyEphChooser);
    }
    initEphComputer(pars, model, *chooser, "DB");
    EphComputer & computer(getEphComputer());

    m_os.out() << prefix << "User supplied time " << *time_rep << std::endl;

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
      try {
        FrequencyEph eph(computer.calcPulsarEph(abs_ref_time));
        m_os.out() << prefix << "Spin ephemeris estimated at the user supplied time is:" << std::endl << eph << std::endl;
      } catch (const std::exception & x) {
        m_os.err() << prefix << "Unexpected problem computing ephemeris." << std::endl << x.what() << std::endl;
      }
    }
  }
}
