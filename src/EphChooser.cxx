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
#include "timeSystem/IntFracPair.h"
#include "timeSystem/TimeInterval.h"

using namespace timeSystem;

namespace pulsarDb {

  const PulsarEph & EphChooser::findClosest(const PulsarEphCont & ephemerides, const timeSystem::AbsoluteTime & t) const {
    if (ephemerides.empty()) throw std::runtime_error("EphChooser::findClosest was passed empty container of ephemerides");
    using std::fabs;

    PulsarEphCont::const_iterator candidate = ephemerides.begin();

    // Seed the current difference with a value which will be larger than any other.
    double diff = std::numeric_limits<double>::max();
    
    for (PulsarEphCont::const_iterator itor = ephemerides.begin(); itor != ephemerides.end(); ++itor) {
      double diff_since = (*itor)->dt(t, (*itor)->valid_since());
      double diff_until = (*itor)->dt(t, (*itor)->valid_until());
      double new_diff = std::min(std::fabs(diff_since), std::fabs(diff_until));
      // Found a better candidate if the new difference is smaller than the previous difference.
      if (new_diff <= diff) {
        candidate = itor;
        diff = new_diff;
      }
    }

    // If no candidate was found, throw an exception.
    if (ephemerides.end() == candidate) {
      std::ostringstream os;
      os << "EphChooser::findClosest could not find an ephemeris for time " << t;
      throw std::runtime_error(os.str());
    }

    return *(*candidate);
  }

  const OrbitalEph & EphChooser::findClosest(const OrbitalEphCont & ephemerides, const timeSystem::AbsoluteTime & t) const {
    if (ephemerides.empty()) throw std::runtime_error("EphChooser::findClosest was passed empty container of ephemerides");
    // Start with minimum = maximum value.
    double min_time_diff = std::numeric_limits<double>::max();

    OrbitalEphCont::const_iterator candidate = ephemerides.end();
    
    // Find the closest ephemeris time to the given time.
    for (OrbitalEphCont::const_iterator itor = ephemerides.begin(); itor != ephemerides.end(); ++itor) {
      double time_diff = std::fabs((*itor)->dt(t));
      if (time_diff <= min_time_diff) {
        candidate = itor;
        min_time_diff = time_diff;
      }
    }

    // If no candidate was found, throw an exception.
    if (ephemerides.end() == candidate) {
      std::ostringstream os;
      os << "EphChooser::findClosest could not find an orbital ephemeris for time " << t;
      throw std::runtime_error(os.str());
    }

    return *(*candidate);
  }

  StrictEphChooser::StrictEphChooser(const timeSystem::ElapsedTime & tolerance): m_tolerance(tolerance) {}

  const PulsarEph & StrictEphChooser::choose(const PulsarEphCont & ephemerides, const timeSystem::AbsoluteTime & t) const {
    PulsarEphCont::const_iterator candidate = ephemerides.end();

    for (PulsarEphCont::const_iterator itor = ephemerides.begin(); itor != ephemerides.end(); ++itor) {
      // See if this ephemeris contains the time.
      if ((*itor)->valid_since() <= t && t < (*itor)->valid_until()) {

        // See if this is the first candidate, which is automatically accepted.
        if (ephemerides.end() == candidate) {
          candidate = itor;
        } else if ((*itor)->valid_since().equivalentTo((*candidate)->valid_since(), m_tolerance)) {
          // The two start at the same time, so break the tie based on which one is valid longer.
          // Note that in a tie here, the one selected is the one appearing last in the sequence.
          if ((*itor)->valid_until() > (*candidate)->valid_until() ||
            (*itor)->valid_until().equivalentTo((*candidate)->valid_until(), m_tolerance))
            candidate = itor;
        // Otherwise, prefer the eph which starts later.
        } else if ((*itor)->valid_since() > (*candidate)->valid_since()) {
          candidate = itor;
        }
      }
    }

    // If no candidate was found, throw an exception.
    if (ephemerides.end() == candidate) {
      std::ostringstream os;
      os << "StrictEphChooser::choose could not find a spin ephemeris for time " << t;
      throw std::runtime_error(os.str());
    }

    return *(*candidate);
  }

  const OrbitalEph & StrictEphChooser::choose(const OrbitalEphCont & ephemerides, const timeSystem::AbsoluteTime & t) const {
    return findClosest(ephemerides, t);
  }

  EphChooser * StrictEphChooser::clone() const {
    return new StrictEphChooser(*this);
  }

  SloppyEphChooser::SloppyEphChooser(): m_strict_chooser() {}

  SloppyEphChooser::SloppyEphChooser(const timeSystem::ElapsedTime & tolerance): m_strict_chooser(tolerance) {}

  const PulsarEph & SloppyEphChooser::choose(const PulsarEphCont & ephemerides, const timeSystem::AbsoluteTime & t) const {
    if (ephemerides.empty()) throw std::runtime_error("SloppyEphChooser::choose was passed empty container of ephemerides");
    using std::fabs;

    // First try to get a strictly correct choice.
    try {
      return m_strict_chooser.choose(ephemerides, t);
    } catch (const std::exception &) {
      // Ignore this exception because the sloppy code below will then try.
    }

    return findClosest(ephemerides, t);
  }

  const OrbitalEph & SloppyEphChooser::choose(const OrbitalEphCont & ephemerides, const timeSystem::AbsoluteTime & t) const {
    return findClosest(ephemerides, t);
  }

  EphChooser * SloppyEphChooser::clone() const {
    return new SloppyEphChooser(*this);
  }
}
