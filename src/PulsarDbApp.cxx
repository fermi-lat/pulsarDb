/** \file PulsarDbApp.cxx
    \brief Implementation of the PulsarDbApp class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#include <cctype>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>

#include "facilities/commonUtilities.h"

#include "pulsarDb/PulsarDb.h"
#include "pulsarDb/PulsarDbApp.h"

#include "st_app/AppParGroup.h"

#include "st_facilities/FileSys.h"

#include "tip/IFileSvc.h"

static const std::string s_cvs_id("$Name:  $");

namespace pulsarDb {

  PulsarDbApp::PulsarDbApp(): m_os("PulsarDbApp", "PulsarDbApp()", 2) {
    setName("gtpulsardb");
    setVersion(s_cvs_id);
  }

  PulsarDbApp::~PulsarDbApp() throw() {}

  void PulsarDbApp::run() {
    m_os.setMethod("run()");
    using namespace st_app;

    // Get parameters.
    AppParGroup & pars(getParGroup());

    // Prompt and save.
    pars.Prompt("psrdbfile");
    pars.Prompt("outfile");
    pars.Prompt("filter");
    std::string filter = pars["filter"];
    for (std::string::iterator itor = filter.begin(); itor != filter.end(); ++itor) *itor = tolower(*itor);
    if (filter == "name") {
      pars.Prompt("psrname");
    } else if (filter == "time") {
      pars.Prompt("tstart");
      pars.Prompt("tstop");
    } else if (filter == "solareph") {
      pars.Prompt("solareph");
    }
    pars.Prompt("author");
    pars.Prompt("chatter");
    pars.Prompt("clobber");
    pars.Prompt("debug");
    pars.Prompt("gui");
    pars.Prompt("mode");
    pars.Save();

    // Handle leap seconds.
    std::string leap_sec_file = pars["leapsecfile"];
    timeSystem::TimeSystem::setDefaultLeapSecFileName(leap_sec_file);

    // Check whether the output file already exists, when clobber is set to no.
    std::string out_file = pars["outfile"];
    bool clobber = pars["clobber"];
    if (!clobber && tip::IFileSvc::instance().fileExists(out_file)) {
      throw std::runtime_error("Output file \"" + out_file + "\" already exists");
    }

    // Get contents of input file, which may be a list of files.
    std::string in_file = pars["psrdbfile"];
    st_facilities::FileSys::FileNameCont file_names = st_facilities::FileSys::expandFileList(in_file);

    // Display an error message if no input files were found.
    if (file_names.empty()) throw std::runtime_error("No files were found matching input file \"" + in_file + "\"");

    // Create an empty pulsar ephemerides database, using the template file.
    std::string tpl_file = facilities::commonUtilities::joinPath(facilities::commonUtilities::getDataPath("pulsarDb"), "PulsarDb.tpl");
    PulsarDb data_base(tpl_file);

    // Load input ephemerides.
    for (st_facilities::FileSys::FileNameCont::const_iterator itor = file_names.begin(); itor != file_names.end(); ++itor) {
      data_base.load(*itor);
    }

    // Filter ephemerides.
    if (filter == "name") {
      // Filter on pulsar name.
      std::string psr_name = pars["psrname"];
      data_base.filterName(psr_name);
      
    } else if (filter == "time") {
      // Filter on time.
      double t_start = pars["tstart"];
      double t_stop = pars["tstop"];
      data_base.filterInterval(t_start, t_stop);

    } else if (filter == "solareph") {
      // Filter on solar system ephemeris.
      std::string solar_eph = pars["solareph"];
      data_base.filterSolarEph(solar_eph);
    }

    // Display a warning message if no ephemerides are left after filtering.
    if (0 >= data_base.getNumEph()) m_os.warn(1).prefix() << "No matching ephemerides were found." << std::endl;

    // Write output.
    std::string author = pars["author"];
    data_base.save(out_file, getName() + " " + getVersion(), author, clobber);
  }

}
