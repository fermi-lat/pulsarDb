/** \file EphChooser.cxx
    \brief Implementation of the EphChooser class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#include <algorithm>
#include <limits>
#include <sstream>
#include <stdexcept>

#include "pulsarDb/EphChooser.h"

#include "timeSystem/AbsoluteTime.h"
#include "timeSystem/TimeInterval.h"

using namespace timeSystem;

namespace pulsarDb {

  const PulsarEph & EphChooser::findClosest(const PulsarEphCont & ephemerides, const timeSystem::AbsoluteTime & abs_time) const {
    if (ephemerides.empty()) throw std::runtime_error("EphChooser::findClosest was passed empty container of ephemerides");
    PulsarEphCont::const_iterator candidate = ephemerides.begin();

    // Seed the current difference with a value which will be larger than any other.
    double diff = std::numeric_limits<double>::max();

    for (PulsarEphCont::const_iterator itor = ephemerides.begin(); itor != ephemerides.end(); ++itor) {
      double diff_since = measureTimeSeparation(abs_time, (*itor)->getValidSince());
      double diff_until = measureTimeSeparation(abs_time, (*itor)->getValidUntil());
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
      os << "EphChooser::findClosest could not find an ephemeris for time " << abs_time;
      throw std::runtime_error(os.str());
    }

    return *(*candidate);
  }

  const OrbitalEph & EphChooser::findClosest(const OrbitalEphCont & ephemerides, const timeSystem::AbsoluteTime & abs_time) const {
    if (ephemerides.empty()) throw std::runtime_error("EphChooser::findClosest was passed empty container of ephemerides");
    // Start with minimum = maximum value.
    double min_time_diff = std::numeric_limits<double>::max();

    OrbitalEphCont::const_iterator candidate = ephemerides.end();
    
    // Find the closest ephemeris time to the given time.
    for (OrbitalEphCont::const_iterator itor = ephemerides.begin(); itor != ephemerides.end(); ++itor) {
      double time_diff = measureTimeSeparation(abs_time, (*itor)->t0());
      if (time_diff <= min_time_diff) {
        candidate = itor;
        min_time_diff = time_diff;
      }
    }

    // If no candidate was found, throw an exception.
    if (ephemerides.end() == candidate) {
      std::ostringstream os;
      os << "EphChooser::findClosest could not find an orbital ephemeris for time " << abs_time;
      throw std::runtime_error(os.str());
    }

    return *(*candidate);
  }

  double EphChooser::measureTimeSeparation(const AbsoluteTime & at1, const AbsoluteTime & at2) const {
    // Use TDB because: 1) we must choose *some* system, and 2) TDB is "steadier" than TT, TAI or UTC.
    return std::fabs((at1 - at2).computeDuration("TDB", "Day"));
  }

  StrictEphChooser::StrictEphChooser(const timeSystem::ElapsedTime & tolerance): m_tolerance(tolerance) {}

  const PulsarEph & StrictEphChooser::choose(const PulsarEphCont & ephemerides, const timeSystem::AbsoluteTime & abs_time) const {
    PulsarEphCont::const_iterator candidate = ephemerides.end();

    for (PulsarEphCont::const_iterator itor = ephemerides.begin(); itor != ephemerides.end(); ++itor) {
      // See if this ephemeris contains the time.
      if ((*itor)->getValidSince() <= abs_time && abs_time < (*itor)->getValidUntil()) {

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
      os << "StrictEphChooser::choose could not find a spin ephemeris for time " << abs_time;
      throw std::runtime_error(os.str());
    }

    return *(*candidate);
  }

  const OrbitalEph & StrictEphChooser::choose(const OrbitalEphCont & ephemerides, const timeSystem::AbsoluteTime & abs_time) const {
    return findClosest(ephemerides, abs_time);
  }

  void StrictEphChooser::examine(const PulsarEphCont & ephemerides, const timeSystem::AbsoluteTime & start_time,
    const timeSystem::AbsoluteTime & stop_time, EphStatusCont & eph_status) const {
    // Check the time order of the time interval arguments.
    if (start_time > stop_time) {
      std::ostringstream os;
      os << "StrictEphChooser::examine is given a negative time interval: " << start_time << " - " << stop_time;
      throw std::runtime_error(os.str());
    }

    // Clear the contents of the returning container.
    eph_status.clear();

    // Check ephemeris availability.
    if (0 == ephemerides.size()) {
      // Put an ephemeris status if no ephemeris is available.
      eph_status.push_back(EphStatus(start_time, stop_time, Unavailable, ""));

    } else {
      // Prepare variables to count up the number of ephemerides that are valid at the edges of the given time interval.
      int num_eph_start = 0;
      int num_eph_stop = 0;
      typedef std::multimap<AbsoluteTime, int> ValidityEdgeCont;
      ValidityEdgeCont validity_edge;

      // Collect edges of ephemeris validity windows.
      for (PulsarEphCont::const_iterator itor = ephemerides.begin(); itor != ephemerides.end(); ++itor) {
        const PulsarEph & this_eph = **itor;
        const AbsoluteTime & valid_since = this_eph.getValidSince();
        const AbsoluteTime & valid_until = this_eph.getValidUntil();

        // Select ephemerides whose validity window has a positive time duration and overlaps with the given time interval,
        // assuming the end of a validity window is considered outside of the window.
        if (valid_since < valid_until && valid_since <= stop_time && start_time < valid_until) {

          // Check where the beginning of this validity window is.
          if (valid_since <= start_time) {
            // This ephemeris is valid at the start time of the given time interval.
            ++num_eph_start;
          } else {
            // This ephemeris becomes valid during the given time interval.
            validity_edge.insert(ValidityEdgeCont::value_type(valid_since, +1));
          }

          // Check where the end of this validity window is.
          if (stop_time < valid_until) {
            // This ephemeris is valid at the stop time of thegiven time interval.
            ++num_eph_stop;
          } else {
            // This ephemeris becomes invalid during the given time interval.
            validity_edge.insert(ValidityEdgeCont::value_type(valid_until, -1));
          }
        }
      }

      // Find ephemeris gaps in the given time interval, where the number of valid ephemerides becomes zero (0).
      AbsoluteTime gap_start = start_time;
      int num_eph_current = num_eph_start;
      for (ValidityEdgeCont::const_iterator itor = validity_edge.begin(); itor != validity_edge.end(); ++itor) {
        const AbsoluteTime & edge_time = itor->first;
        const int & num_eph_change = itor->second;
        if (1 == num_eph_current && -1 == num_eph_change) {
          // The beginning of an ephemeris gap.
          gap_start = edge_time;

        } else if (0 == num_eph_current && +1 == num_eph_change) {
          // The end of the current ephemeris gap.
          eph_status.push_back(EphStatus(gap_start, edge_time, Unavailable, ""));
        }

        // Update the number of ephemerides that are currently valid.
        num_eph_current += num_eph_change;
      }

      // Close an ephemeris gap that has not ended within the given time interval.
      if (0 == num_eph_current) eph_status.push_back(EphStatus(gap_start, stop_time, Unavailable, ""));

      // Sanity check on the number of valid ephemerides at the end of the given time interval.
      if (num_eph_current != num_eph_stop) {
        std::ostringstream os;
        os << "StrictEphChooser::examine found " << num_eph_current << " ephemeri(de)s, not " << num_eph_stop <<
          ", at the end of the given time interval " << stop_time;
        throw std::runtime_error(os.str());
      }
    }
  }

  EphChooser * StrictEphChooser::clone() const {
    return new StrictEphChooser(*this);
  }

  SloppyEphChooser::SloppyEphChooser(): m_strict_chooser() {}

  SloppyEphChooser::SloppyEphChooser(const timeSystem::ElapsedTime & tolerance): m_strict_chooser(tolerance) {}

  const PulsarEph & SloppyEphChooser::choose(const PulsarEphCont & ephemerides, const timeSystem::AbsoluteTime & abs_time) const {
    if (ephemerides.empty()) throw std::runtime_error("SloppyEphChooser::choose was passed empty container of ephemerides");
    // First try to get a strictly correct choice.
    try {
      return m_strict_chooser.choose(ephemerides, abs_time);
    } catch (const std::exception &) {
      // Ignore this exception because the sloppy code below will then try.
    }

    return findClosest(ephemerides, abs_time);
  }

  const OrbitalEph & SloppyEphChooser::choose(const OrbitalEphCont & ephemerides, const timeSystem::AbsoluteTime & abs_time) const {
    return findClosest(ephemerides, abs_time);
  }

  void SloppyEphChooser::examine(const PulsarEphCont & ephemerides, const timeSystem::AbsoluteTime & start_time,
    const timeSystem::AbsoluteTime & stop_time, EphStatusCont & eph_status) const {
    // Check the time order of the time interval arguments.
    if (start_time > stop_time) {
      std::ostringstream os;
      os << "StrictEphChooser::examine is given a negative time interval: " << start_time << " - " << stop_time;
      throw std::runtime_error(os.str());
    }

    // Clear the contents of the returning container.
    eph_status.clear();

    // Check ephemeris availability.
    if (0 == ephemerides.size()) {
      // Put an ephemeris status if no ephemeris is available.
      eph_status.push_back(EphStatus(start_time, stop_time, Unavailable, ""));

    } else {
      // Use StrictEphChooser to examine ephemerides.
      m_strict_chooser.examine(ephemerides, start_time, stop_time, eph_status);

      // Replace "Unavailable" status with "Extrapolated" status.
      for (EphStatusCont::iterator itor = eph_status.begin(); itor != eph_status.end(); ++itor) {
        const EphStatusCodeType & status_code = itor->getStatusCode();
        if (Unavailable == status_code) {
          const AbsoluteTime & since = itor->getEffectiveSince();
          const AbsoluteTime & until = itor->getEffectiveUntil();
          const std::string & description = itor->getDescription();
          *itor = EphStatus(since, until, Extrapolated, description);
        }
      }
    }
  }

  EphChooser * SloppyEphChooser::clone() const {
    return new SloppyEphChooser(*this);
  }
}
