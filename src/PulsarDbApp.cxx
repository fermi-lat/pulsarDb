/** \file PulsarDbApp.cxx
    \brief Implementation of the PulsarDbApp class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>

#include "pulsarDb/PulsarDb.h"
#include "pulsarDb/PulsarDbApp.h"
#include "pulsarDb/TextPulsarDb.h"

#include "st_app/AppParGroup.h"
#include "st_facilities/Env.h"
#include "st_facilities/FileSys.h"

#include "tip/IFileSvc.h"
#include "tip/TipException.h"

static const std::string s_cvs_id("$Name$");

namespace pulsarDb {

  PulsarDbApp::PulsarDbApp(): m_tpl_file() {
    setName("gtpulsardb");
    setVersion(s_cvs_id);
  }

  PulsarDbApp::~PulsarDbApp() throw() {}

  void PulsarDbApp::run() {
    using namespace st_app;

    // Get parameters.
    AppParGroup & pars(getParGroup("gtpulsardb"));

    // Prompt and save.
    pars.Prompt();
    pars.Save();

    // Interpret pars.
    std::string in_file = pars["psrdbfile"].Value();
    std::string expression = pars["filter"];
    std::string out_file = pars["outfile"];
    std::string psr_name = pars["psrname"];
    double t_start = pars["tstart"];
    double t_stop = pars["tstop"];

    // Find template file.
    m_tpl_file = st_facilities::Env::appendFileName(st_facilities::Env::getDataDir("pulsarDb"), "PulsarEph.tpl");

    // Create output file, respecting clobber.
    tip::IFileSvc::instance().createFile(out_file, m_tpl_file, pars["clobber"]);

    // Get contents of input file, which may be a list of files.
    st_facilities::FileSys::FileNameCont file_names = st_facilities::FileSys::expandFileList(in_file);

    if (file_names.empty()) throw std::runtime_error("No files were found matching input file \"" + in_file + "\"");

    for (st_facilities::FileSys::FileNameCont::const_iterator itor = file_names.begin(); itor != file_names.end(); ++itor) {
      // Create data base representation.
      std::auto_ptr<PulsarDb> data_base(openDbFile(*itor));

      // Write output.
      data_base->save(out_file, m_tpl_file);
    }

    // Re-open the output file as the data_base object.
    std::auto_ptr<PulsarDb> data_base(openDbFile(out_file, true));

    // Filter using user's expression.
    data_base->filter(expression);

    // Filter on pulsar name.
    data_base->filterName(psr_name);

    // Filter on time.
    data_base->filterInterval(t_start, t_stop);
  }

  PulsarDb * PulsarDbApp::openDbFile(const std::string & in_file, bool edit_in_place) {
    PulsarDb * data_base = 0;

    // First try opening a tip file containing the ephemerides.
    try {
      data_base = new PulsarDb(in_file, edit_in_place);
    } catch (const tip::TipException &) {
      // Ignore an error; hope that maybe it is a text file.
    }

    // Next try opening it as a text file if there was a problem.
    if (0 == data_base) {
      try {
        data_base = new TextPulsarDb(in_file, m_tpl_file);
      } catch (const std::exception &) {
        throw std::runtime_error("Could not open file " + in_file + " as either a FITS file or a text data file");
      }
    }
    return data_base;
  }
}
