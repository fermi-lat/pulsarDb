/** \file EphChooser.cxx
    \brief Implementation of the EphChooser class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#include <limits>
#include <sstream>
#include <stdexcept>

#include "pulsarDb/AbsoluteTime.h"
#include "pulsarDb/EphChooser.h"

namespace pulsarDb {

  const PulsarEph & EphChooser::choose(const PulsarEphCont & ephemerides, const AbsoluteTime & t) const {
    PulsarEphCont::const_iterator candidate = ephemerides.end();

    for (PulsarEphCont::const_iterator itor = ephemerides.begin(); itor != ephemerides.end(); ++itor) {
      // See if this ephemeris contains the time.
      if ((*itor)->valid_since() <= t && t < (*itor)->valid_until()) {

        // See if this is the first candidate, which is automatically accepted.
        if (ephemerides.end() == candidate) {
          candidate = itor;
        // Otherwise, prefer the eph which starts later.
        } else if ((*itor)->valid_since() > (*candidate)->valid_since()) {
          candidate = itor;
        } else if ((*itor)->valid_since() == (*candidate)->valid_since()) {
          // The two start at the same time, so break the tie based on which one is valid longer.
          // Note that in a tie here, the one selected is the one appearing last in the sequence.
          if ((*itor)->valid_until() >= (*candidate)->valid_until())
            candidate = itor;
        }
      }
    }

    // If no candidate was found, throw an exception.
    if (ephemerides.end() == candidate) {
      std::ostringstream os;
      os << "EphChooser::choose could not find a spin ephemeris for time " << t;
      throw std::runtime_error(os.str());
    }

    return *(*candidate);
  }

  const OrbitalEph & EphChooser::choose(const OrbitalEphCont & ephemerides, const AbsoluteTime & t) const {
    if (ephemerides.empty()) throw std::runtime_error("EphChooser::choose was passed empty container of ephemerides");
    // Start with minimum = maximum value.
    double min_time_diff = std::numeric_limits<double>::max();

    OrbitalEphCont::const_iterator candidate = ephemerides.end();
    
    // Find the closest ephemeris time to the given time.
    for (OrbitalEphCont::const_iterator itor = ephemerides.begin(); itor != ephemerides.end(); ++itor) {
      double time_diff = fabs(((*itor)->t0() - t).sec());
      if (time_diff < min_time_diff) {
        candidate = itor;
        min_time_diff = time_diff;
      }
    }

    // If no candidate was found, throw an exception.
    if (ephemerides.end() == candidate) {
      std::ostringstream os;
      os << "EphChooser::choose could not find an orbital ephemeris for time " << t;
      throw std::runtime_error(os.str());
    }

    return *(*candidate);
  }

  EphChooser * EphChooser::clone() const {
    return new EphChooser(*this);
  }

  const PulsarEph & SloppyEphChooser::choose(const PulsarEphCont & ephemerides, const AbsoluteTime & t) const {
    // First try to get a strictly correct choice.
    try {
      return EphChooser::choose(ephemerides, t);
    } catch (const std::exception &) {
      // Ignore this exception because the sloppy code below will then try.
    }

    PulsarEphCont::const_iterator candidate = ephemerides.begin();

    Duration diff = std::min((t - (*candidate)->valid_since()).op(fabs), (t - (*candidate)->valid_until()).op(fabs));
    
    for (PulsarEphCont::const_iterator itor = ephemerides.begin(); itor != ephemerides.end(); ++itor) {
      Duration new_diff = std::min((t - (*candidate)->valid_since()).op(fabs), (t - (*candidate)->valid_until()).op(fabs));
      if (new_diff < diff) {
        candidate = itor;
        diff = new_diff;
      }
    }

    // If no candidate was found, throw an exception.
    if (ephemerides.end() == candidate) {
      std::ostringstream os;
      os << "SloppyEphChooser::choose could not find an ephemeris for time " << t;
      throw std::runtime_error(os.str());
    }

    return *(*candidate);
  }

  EphChooser * SloppyEphChooser::clone() const {
    return new SloppyEphChooser(*this);
  }

}
