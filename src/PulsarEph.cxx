/** \file PulsarEph.cxx
    \brief Implementation of the PulsarEph class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#include <iomanip>

#include "pulsarDb/PulsarEph.h"

#include "timeSystem/ElapsedTime.h"
#include "timeSystem/IntFracPair.h"
#include "timeSystem/TimeInterval.h"
#include "timeSystem/TimeRep.h"

using namespace timeSystem;

namespace pulsarDb {

  double PulsarEph::calcPulsePhase(const timeSystem::AbsoluteTime & ev_time, double phase_offset) const {
    // Compute the cycle count (which includes the iteger part of pulse phase).
    double cycle_count = calcCycleCount(ev_time);

    // Express phase as a value between 0. and 1., after adding a global phase offset.
    double int_part; // ignored, needed for modf.
    double phase = std::modf(phase_offset + cycle_count, &int_part);
    if (phase < 0.) ++phase;
    return phase;
  }

  double FrequencyEph::calcElapsedSecond(const AbsoluteTime & at1, const AbsoluteTime & at2) const {
    return (at1 - at2).computeElapsedTime(m_system->getName()).getTime().getValue(Sec).getDouble();
  }

  double PeriodEph::calcElapsedSecond(const AbsoluteTime & at1, const AbsoluteTime & at2) const {
    return (at1 - at2).computeElapsedTime(m_system->getName()).getTime().getValue(Sec).getDouble();
  }

  // TODO: Unify/generalize this method with PeriodEph::write method.
  st_stream::OStream & FrequencyEph::write(st_stream::OStream & os) const {
    std::ios::fmtflags orig_flags = os.flags();
    int orig_prec = os.precision(15);
    os << std::right;
    os.prefix().width(14);
    std::string time_system_name = getSystem().getName();
    MjdRep mjd_rep(time_system_name, 0, 0.);

    // Note: below, break into two consecutive strings so that width applies to first part only.
    if (!getValidSince().equivalentTo(getValidUntil(), ElapsedTime(time_system_name, Duration(0, 1.e-9)))) {
      mjd_rep = getValidSince();
      os << "Validity : " << "in range " << "[" << mjd_rep << ", ";
      mjd_rep = getValidUntil();
      os << mjd_rep << ")" << std::endl;
    } else {
      mjd_rep = getValidSince();
      os << "Validity : " << "only at time " << mjd_rep << std::endl;
    }
    mjd_rep = getEpoch();
    os.prefix().width(14); os << "Epoch = " << mjd_rep << std::endl;
    os.prefix().width(14); os << "RA = " << m_ra << std::endl;
    os.prefix().width(14); os << "Dec = " << m_dec << std::endl;
    os.prefix().width(14); os << "Phi0 = " << m_phi0 << std::endl;
    os.prefix().width(14); os << "F0 = " << m_f0 << std::endl;
    os.prefix().width(14); os << "F1 = " << m_f1 << std::endl;
    os.prefix().width(14); os << "F2 = " << m_f2;
    os.flags(orig_flags);
    os.precision(orig_prec);
    return os;
  }

  // TODO: Unify/generalize this method with FrequencyEph::write method.
  st_stream::OStream & PeriodEph::write(st_stream::OStream & os) const {
    std::ios::fmtflags orig_flags = os.flags();
    int orig_prec = os.precision(15);
    os << std::right;
    os.prefix().width(14);
    std::string time_system_name = getSystem().getName();
    MjdRep mjd_rep(time_system_name, 0, 0.);

    // Note: below, break into two consecutive strings so that width applies to first part only.
    if (!getValidSince().equivalentTo(getValidUntil(), ElapsedTime(time_system_name, Duration(0, 1.e-9)))) {
      mjd_rep = getValidSince();
      os << "Validity : " << "in range " << "[" << mjd_rep << ", ";
      mjd_rep = getValidUntil();
      os << mjd_rep << ")" << std::endl;
    } else {
      mjd_rep = getValidSince();
      os << "Validity : " << "only at time " << mjd_rep << std::endl;
    }
    mjd_rep = getEpoch();
    os.prefix().width(14); os << "Epoch = " << mjd_rep << std::endl;
    os.prefix().width(14); os << "RA = " << m_ra << std::endl;
    os.prefix().width(14); os << "Dec = " << m_dec << std::endl;
    os.prefix().width(14); os << "Phi0 = " << m_phi0 << std::endl;
    os.prefix().width(14); os << "P0 = " << m_p0 << std::endl;
    os.prefix().width(14); os << "P1 = " << m_p1 << std::endl;
    os.prefix().width(14); os << "P2 = " << m_p2;
    os.flags(orig_flags);
    os.precision(orig_prec);
    return os;
  }

  double FrequencyEph::calcCycleCount(const timeSystem::AbsoluteTime & ev_time) const {
    double dt = calcElapsedSecond(ev_time);
    double dt_squared = dt * dt;
    double cycle_count = m_phi0 + m_f0 * dt + m_f1/2.0 * dt_squared + m_f2/6.0 * dt * dt_squared;
    return cycle_count;
  }

  double PeriodEph::calcCycleCount(const timeSystem::AbsoluteTime & ev_time) const {
    // TODO: Compute pulse phase directly from p0, p1, and p2, using indefinite integral of ax^2+bx+c.
    double f0 = 1. / m_p0;
    double f1 = - m_p1 / (m_p0 * m_p0);
    double p0sq = m_p0 * m_p0;
    double f2 = 2. * m_p1 * m_p1 / (m_p0 * p0sq) - m_p2 / p0sq;

    double dt = calcElapsedSecond(ev_time);
    double dt_squared = dt * dt;
    double cycle_count = m_phi0 + f0 * dt + f1/2.0 * dt_squared + f2/6.0 * dt * dt_squared;
    return cycle_count;
  }

  double FrequencyEph::calcFrequency(const timeSystem::AbsoluteTime & ev_time, int derivative_order) const {
    double return_value = 0.;
    if (0 == derivative_order) {
      double dt = calcElapsedSecond(ev_time);
      return_value = m_f0 + m_f1 * dt + 0.5 * m_f2 * dt * dt;
    } else if (1 == derivative_order) {
      double dt = calcElapsedSecond(ev_time);
      return_value = m_f1 + m_f2 * dt;
    } else if (2 == derivative_order) {
      return_value = m_f2;
    } else {
      return_value = 0.;
    }
    return return_value;
  }

  double PeriodEph::calcFrequency(const timeSystem::AbsoluteTime & ev_time, int derivative_order) const {
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

  std::pair<double, double> FrequencyEph::calcSkyPosition(const timeSystem::AbsoluteTime & /* ev_time */) const {
    return std::make_pair(m_ra, m_dec);
  }

  std::pair<double, double> PeriodEph::calcSkyPosition(const timeSystem::AbsoluteTime & /* ev_time */) const {
    return std::make_pair(m_ra, m_dec);
  }
}
