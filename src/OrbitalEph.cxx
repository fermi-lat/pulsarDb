/** \file OrbitalEph.cxx
    \brief Implementation of the OrbitalEph class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#include <limits>
#include <iomanip>
#include <iostream>

#include "st_stream/Stream.h"

#include "pulsarDb/OrbitalEph.h"

#include "timeSystem/AbsoluteTime.h"
#include "timeSystem/ElapsedTime.h"

namespace pulsarDb {

  OrbitalEph::OrbitalEph(const timeSystem::ElapsedTime & tolerance, int max_iteration): m_tolerance(tolerance),
    m_max_iteration(max_iteration) {}

  void OrbitalEph::modulateBinary(timeSystem::AbsoluteTime & ev_time) const {
    ev_time += calcOrbitalDelay(ev_time);
  }

  void OrbitalEph::demodulateBinary(timeSystem::AbsoluteTime & ev_time) const {
    // Save target arrival time (ev_time) in orig_time.
    timeSystem::AbsoluteTime orig_time = ev_time;

    // Initial guess of orbital delay.
    timeSystem::ElapsedTime delay = calcOrbitalDelay(ev_time);

    // Iterative approximation of demodulated time.
    int ii;
    for (ii=0; ii<m_max_iteration; ii++) {

      // Compute next candidate of demodulated time.
      ev_time = orig_time - delay;

      // Compute orbital delay at ev_time.
      delay = calcOrbitalDelay(ev_time);

      // Compare time difference between candidate demodulated time
      // (ev_time) and target arrival time (orig_time) with the
      // estimated orbital delay based on the binary model (delay).
      if (orig_time.equivalentTo(ev_time + delay, m_tolerance)) break;
    }

    // Check for non-convergence.
    if (ii == m_max_iteration) throw std::runtime_error("Binary demodulation did not converge");

  }

  st_stream::OStream & OrbitalEph::write(st_stream::OStream & os) const {
    // Save the original settings and set the prefered formats.
    std::ios::fmtflags orig_flags = os.flags();
    int orig_prec = os.precision(std::numeric_limits<double>::digits10);
    os << std::right;

    // Write subclass-specific parameters (delegated to subclass).
    writeModelParameter(os);

    // Restore the saved settings.
    os.flags(orig_flags);
    os.precision(orig_prec);
    return os;
  }
}
