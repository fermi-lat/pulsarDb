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

#include <cctype>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

static const std::string s_cvs_id("$Name:  $");

using namespace timeSystem;

namespace pulsarDb {

  EphComputerApp::EphComputerApp(): m_os("EphComputerApp", "", 2) {
    setName("gtephcomp");
    setVersion(s_cvs_id);
  }

  EphComputerApp::~EphComputerApp() throw() {}

  void EphComputerApp::run() {
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
    for (std::string::iterator itor = time_format.begin(); itor != time_format.end(); ++itor) *itor = toupper(*itor);
    std::string time_sys = pars["timesys"];
    for (std::string::iterator itor = time_sys.begin(); itor != time_sys.end(); ++itor) *itor = toupper(*itor);
    bool strict = pars["strict"];

    std::auto_ptr<TimeRep> time_rep(0);

    // Create the correct time representation for this time system and format.
    if (time_format == "GLAST") {
      time_rep.reset(new GlastMetRep("TDB", 0.));
    } else if (time_format == "MJD") {
      time_rep.reset(new MjdRep("TDB", 0, 0.));
    } else {
      throw std::runtime_error("Time format \"" + time_format + "\" is not supported");
    }

    // Set the time of the representation to be the given reference time.
    time_rep->assign(ref_time);

    // Find the pulsar database.
    std::string psrdb_file = pars["psrdbfile"];
    std::string psrdb_file_uc = psrdb_file;
    for (std::string::iterator itor = psrdb_file_uc.begin(); itor != psrdb_file_uc.end(); ++itor) *itor = toupper(*itor);
    if ("DEFAULT" == psrdb_file_uc) {
      using namespace st_facilities;
      psrdb_file = Env::appendFileName(Env::getDataDir("pulsarDb"), "master_pulsardb.fits");
    }
    std::string psr_name = pars["psrname"];

    // TODO Remove the following line when correct conversion from TT is implemented.
    if (time_sys != "TDB") throw std::runtime_error("Only TDB time system is supported for now");

    AbsoluteTime abs_ref_time(*time_rep);

    std::auto_ptr<EphComputer> computer(0);
    if (strict) {
      computer.reset(new EphComputer);
    } else {
      TimingModel model;
      SloppyEphChooser chooser;
      computer.reset(new EphComputer(model, chooser));
    }

    // Open the database.
    PulsarDb database(psrdb_file);

    // Select only ephemerides for this pulsar.
    database.filterName(psr_name);

    // Load the selected ephemerides.
    computer->load(database);

    m_os.out() << prefix << "User supplied time " << *time_rep << std::endl;

    // Cosmetic: suppress info.
    m_os.info().setPrefix(m_os.out().getPrefix());

    // Set off the optional output.
    std::string dashes(26, '-');
    m_os.info(3) << prefix << dashes << std::endl;

    // Report the best spin ephemeris.
    bool found_pulsar_eph = false;
    try {
      const PulsarEph & eph(computer->choosePulsarEph(abs_ref_time));
      m_os.info(3) << prefix << "Spin ephemeris chosen from database is:" << std::endl << eph << std::endl;
      found_pulsar_eph = true;
    } catch (const std::exception & x) {
      m_os.out() << prefix << "No spin ephemeris was found in database." << std::endl;
    }

    // Report the best binary ephemeris.
    try {
      const OrbitalEph & eph(computer->chooseOrbitalEph(abs_ref_time));
      m_os.info(3) << prefix << "Orbital ephemeris chosen from database is:" << std::endl << eph << std::endl;
    } catch (const std::exception & x) {
      m_os.info(3) << prefix << "No orbital ephemeris was found in database." << std::endl;
    }

    // Set off the optional output.
    m_os.info(3) << prefix << dashes << std::endl;

    // Report the calculated spin ephemeris, provided at least a spin ephemeris was found above.
    if (found_pulsar_eph) {
      try {
        FrequencyEph eph(computer->calcPulsarEph(abs_ref_time));
        m_os.out() << prefix << "Spin ephemeris estimated at the user supplied time is:" << std::endl << eph << std::endl;
      } catch (const std::exception & x) {
        m_os.err() << prefix << "Unexpected problem computing ephemeris." << std::endl << x.what() << std::endl;
      }
    }
  }
}
