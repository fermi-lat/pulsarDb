/** \file EphChooser.cxx
    \brief Implementation of the EphChooser class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#include <limits>
#include <sstream>
#include <stdexcept>

#include "pulsarDb/EphChooser.h"

#include "timeSystem/AbsoluteTime.h"
#include "timeSystem/Duration.h"
#include "timeSystem/ElapsedTime.h"
#include "timeSystem/IntFracPair.h"
#include "timeSystem/TimeInterval.h"

using namespace timeSystem;

namespace pulsarDb {

  const PulsarEph & EphChooser::choose(const PulsarEphCont & ephemerides, const timeSystem::AbsoluteTime & t) const {
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
//        } else if ((*itor)->valid_since() == (*candidate)->valid_since()) {
        } else if ((*itor)->valid_since() >= (*candidate)->valid_since()) {
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

  const OrbitalEph & EphChooser::choose(const OrbitalEphCont & ephemerides, const timeSystem::AbsoluteTime & t) const {
    if (ephemerides.empty()) throw std::runtime_error("EphChooser::choose was passed empty container of ephemerides");
    // Start with minimum = maximum value.
    double min_time_diff = std::numeric_limits<double>::max();

    OrbitalEphCont::const_iterator candidate = ephemerides.end();
    
    // Find the closest ephemeris time to the given time.
    for (OrbitalEphCont::const_iterator itor = ephemerides.begin(); itor != ephemerides.end(); ++itor) {
      IntFracPair diff_pair = ((*itor)->t0() - t).computeElapsedTime("TDB").getTime().getValue(Sec);
      double time_diff = std::fabs(diff_pair.getIntegerPart() + diff_pair.getFractionalPart());
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

  const PulsarEph & SloppyEphChooser::choose(const PulsarEphCont & ephemerides, const timeSystem::AbsoluteTime & t) const {
    if (ephemerides.empty()) throw std::runtime_error("SloppyEphChooser::choose was passed empty container of ephemerides");
    using std::fabs;

    // First try to get a strictly correct choice.
    try {
      return EphChooser::choose(ephemerides, t);
    } catch (const std::exception &) {
      // Ignore this exception because the sloppy code below will then try.
    }

    PulsarEphCont::const_iterator candidate = ephemerides.begin();

    IntFracPair diff_since = (t - (*candidate)->valid_since()).computeElapsedTime("TDB").getTime().getValue(Sec);
    IntFracPair diff_until = (t - (*candidate)->valid_until()).computeElapsedTime("TDB").getTime().getValue(Sec);
    double diff = std::min(std::fabs(diff_since.getIntegerPart() + diff_since.getFractionalPart()),
      std::fabs(diff_until.getIntegerPart() + diff_until.getFractionalPart()));
    
    for (PulsarEphCont::const_iterator itor = ephemerides.begin(); itor != ephemerides.end(); ++itor) {
      diff_since = (t - (*candidate)->valid_since()).computeElapsedTime("TDB").getTime().getValue(Sec);
      diff_until = (t - (*candidate)->valid_until()).computeElapsedTime("TDB").getTime().getValue(Sec);
      double new_diff = std::min(std::fabs(diff_since.getIntegerPart() + diff_since.getFractionalPart()),
        std::fabs(diff_until.getIntegerPart() + diff_until.getFractionalPart()));
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
