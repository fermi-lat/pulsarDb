/** \file PulsarEph.cxx
    \brief Implementation of the PulsarEph class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#include <iomanip>
#include <iostream>
#include <limits>

#include "pulsarDb/PulsarEph.h"

#include "timeSystem/AbsoluteTime.h"
#include "timeSystem/CalendarFormat.h"
#include "timeSystem/MjdFormat.h"
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
    std::string time_string;

    // Write the start time of the validity window.
    try {
      time_string = getValidSince().represent(time_system_name, "MJD");
    } catch (const std::exception &) {
      time_string = getValidSince().represent(time_system_name, "Calendar");
    }
    os << format("Valid Since", time_string, " : ") << std::endl;

    // Write the stop time of the validity window.
    try {
      time_string = getValidUntil().represent(time_system_name, "MJD");
    } catch (const std::exception &) {
      time_string = getValidUntil().represent(time_system_name, "Calendar");
    }
    os << format("Valid Until", time_string, " : ") << std::endl;

    // Write subclass-specific parameters (delegated to subclass).
    writeModelParameter(os);

    // Restore the saved settings.
    os.flags(orig_flags);
    os.precision(orig_prec);
    return os;
  }
}
