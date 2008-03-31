/** \file PeriodEph.cxx
    \brief Implementation of the PeriodEph class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#include <iostream>

#include "pulsarDb/PulsarEph.h"
#include "pulsarDb/PeriodEph.h"

#include "timeSystem/ElapsedTime.h"
#include "timeSystem/TimeInterval.h"
#include "timeSystem/TimeRep.h"
#include "timeSystem/TimeSystem.h"

using namespace timeSystem;

namespace pulsarDb {

  double PeriodEph::calcElapsedSecond(const AbsoluteTime & at1, const AbsoluteTime & at2) const {
    return (at1 - at2).computeElapsedTime(m_system->getName()).getTime().getValue(Sec).getDouble();
  }

  void PeriodEph::writeModelParameter(st_stream::OStream & os) const {
    MjdRep mjd_rep(m_system->getName(), 0, 0.);
    mjd_rep = m_epoch;
    os << format("Epoch", mjd_rep) << std::endl;
    os << format("RA",    m_ra)    << std::endl;
    os << format("Dec",   m_dec)   << std::endl;
    os << format("Phi0",  m_phi0)  << std::endl;
    os << format("P0",    m_p0)    << std::endl;
    os << format("P1",    m_p1)    << std::endl;
    os << format("P2",    m_p2);
  }

  double PeriodEph::calcPulsePhase(const AbsoluteTime & ev_time, double phase_offset) const {
    // TODO: Compute pulse phase directly from p0, p1, and p2, using indefinite integral of ax^2+bx+c.
    double f0 = 1. / m_p0;
    double f1 = - m_p1 / (m_p0 * m_p0);
    double p0sq = m_p0 * m_p0;
    double f2 = 2. * m_p1 * m_p1 / (m_p0 * p0sq) - m_p2 / p0sq;

    // Compute a pulse phase value.
    double dt = calcElapsedSecond(ev_time);
    double dt_squared = dt * dt;
    double phase = m_phi0 + f0 * dt + f1/2.0 * dt_squared + f2/6.0 * dt * dt_squared;
    return phase;

    // Express phase as a value between 0. and 1., after adding a global phase offset.
    return trimPhaseValue(phase, phase_offset);
  }

  double PeriodEph::calcFrequency(const AbsoluteTime & ev_time, int derivative_order) const {
    double return_value = 0.;

    // TODO: Compute frequency (or its derivative) directly from p0, p1, and p2, using generic formula (Bell's polynomial?).
    double f0 = 1. / m_p0;
    double f1 = - m_p1 / (m_p0 * m_p0);
    double p0sq = m_p0 * m_p0;
    double f2 = 2. * m_p1 * m_p1 / (m_p0 * p0sq) - m_p2 / p0sq;
    if (0 == derivative_order) {
      double dt = calcElapsedSecond(ev_time);
      return_value = f0 + f1 * dt + 0.5 * f2 * dt * dt;
    } else if (1 == derivative_order) {
      double dt = calcElapsedSecond(ev_time);
      return_value = f1 + f2 * dt;
    } else if (2 == derivative_order) {
      return_value = f2;
    } else {
      return_value = 0.;
    }
    return return_value;
  }

  std::pair<double, double> PeriodEph::calcSkyPosition(const AbsoluteTime & /* ev_time */) const {
    return std::make_pair(m_ra, m_dec);
  }
}
