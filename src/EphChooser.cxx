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
#include "timeSystem/TimeInterval.h"

using namespace timeSystem;

namespace pulsarDb {

  const PulsarEph & EphChooser::findClosest(const PulsarEphCont & ephemerides, const timeSystem::AbsoluteTime & t) const {
    if (ephemerides.empty()) throw std::runtime_error("EphChooser::findClosest was passed empty container of ephemerides");
    PulsarEphCont::const_iterator candidate = ephemerides.begin();

    // Seed the current difference with a value which will be larger than any other.
    double diff = std::numeric_limits<double>::max();
    
    for (PulsarEphCont::const_iterator itor = ephemerides.begin(); itor != ephemerides.end(); ++itor) {
      double diff_since = measureTimeSeparation(t, (*itor)->getValidSince());
      double diff_until = measureTimeSeparation(t, (*itor)->getValidUntil());
      double new_diff = std::min(diff_since, diff_until);
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
      double time_diff = measureTimeSeparation(t, (*itor)->t0());
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

  double EphChooser::measureTimeSeparation(const AbsoluteTime & at1, const AbsoluteTime & at2) const {
    // Use TDB because: 1) we must choose *some* system, and 2) TDB is "steadier" than TT, TAI or UTC.
    return std::fabs((at1 - at2).computeDuration("TDB", "Day"));
  }

  StrictEphChooser::StrictEphChooser(const timeSystem::ElapsedTime & tolerance): m_tolerance(tolerance) {}

  const PulsarEph & StrictEphChooser::choose(const PulsarEphCont & ephemerides, const timeSystem::AbsoluteTime & t) const {
    PulsarEphCont::const_iterator candidate = ephemerides.end();

    for (PulsarEphCont::const_iterator itor = ephemerides.begin(); itor != ephemerides.end(); ++itor) {
      // See if this ephemeris contains the time.
      if ((*itor)->getValidSince() <= t && t < (*itor)->getValidUntil()) {

        // See if this is the first candidate, which is automatically accepted.
        if (ephemerides.end() == candidate) {
          candidate = itor;
        } else if ((*itor)->getValidSince().equivalentTo((*candidate)->getValidSince(), m_tolerance)) {
          // The two start at the same time, so break the tie based on which one is valid longer.
          // Note that in a tie here, the one selected is the one appearing last in the sequence.
          if ((*itor)->getValidUntil() > (*candidate)->getValidUntil() ||
            (*itor)->getValidUntil().equivalentTo((*candidate)->getValidUntil(), m_tolerance))
            candidate = itor;
        // Otherwise, prefer the eph which starts later.
        } else if ((*itor)->getValidSince() > (*candidate)->getValidSince()) {
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
