/** \file PulsarEph.cxx
    \brief Implementation of the PulsarEph class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#include <iomanip>
#include <iostream>
#include <limits>

#include "pulsarDb/PulsarEph.h"

#include "timeSystem/TimeRep.h"
#include "timeSystem/TimeSystem.h"

using namespace timeSystem;

namespace pulsarDb {

  st_stream::OStream & PulsarEph::write(st_stream::OStream & os) const {
    // Save the original settings and set the prefered formats.
    std::ios::fmtflags orig_flags = os.flags();
    int orig_prec = os.precision(std::numeric_limits<double>::digits10);
    os << std::right;

    // Prepare for MJD expression of time.
    std::string time_system_name = getSystem().getName();
    MjdRep mjd_rep(time_system_name, 0, 0.);

    // Write validity window.
    mjd_rep = getValidSince();
    os << format("Valid Since", mjd_rep, " : ") << std::endl;
    mjd_rep = getValidUntil();
    os << format("Valid Until", mjd_rep, " : ") << std::endl;

    // Write subclass-specific parameters (delegated to subclass).
    writeModelParameter(os);

    // Restore the saved settings.
    os.flags(orig_flags);
    os.precision(orig_prec);
    return os;
  }
}
