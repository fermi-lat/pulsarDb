/** \file PulsarDbApp.cxx
    \brief Implementation of the PulsarDbApp class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#include <memory>
#include <string>

#include "pulsarDb/PulsarDb.h"
#include "pulsarDb/PulsarDbApp.h"

#include "st_app/AppParGroup.h"

namespace pulsarDb {

  PulsarDbApp::PulsarDbApp() { setName("gtpulsardb"); }

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

    // Create data base representation.
    std::auto_ptr<PulsarDb> data_base(new PulsarDb(in_file));

    // Filter using user's expression.
    data_base->filter(expression);

    // Filter on pulsar name.
    data_base->filterName(psr_name);

    // Filter on time.
    data_base->filterInterval(t_start, t_stop);

    // Write output.
    data_base->save(out_file);
  }

}
