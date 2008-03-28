/** \file PulsarEph.cxx
    \brief Implementation of the PulsarEph class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#include <iomanip>
#include <limits>

#include "pulsarDb/PulsarEph.h"

#include "timeSystem/ElapsedTime.h"
#include "timeSystem/IntFracPair.h"
#include "timeSystem/TimeInterval.h"
#include "timeSystem/TimeRep.h"

#include "tip/Header.h"

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

  FrequencyEph::FrequencyEph(const tip::Table::ConstRecord & record, const tip::Header & /* header */):
    m_system(&TimeSystem::getSystem("TDB")), m_since("TDB", Duration(0, 0.), Duration(0, 0.)),
    m_until("TDB", Duration(0, 0.), Duration(0, 0.)), m_epoch("TDB", Duration(0, 0.), Duration(0, 0.)), m_ra(0.), m_dec(0.),
    m_phi0(0.), m_f0(0.), m_f1(0.), m_f2(0.) {
    // Epoch and toa are split into int and frac parts.
    long epoch_int = 0;
    double epoch_frac = 0.;
    long toa_int = 0;
    double toa_frac = 0.;

    // Read reference epoch and pulse TOA (integer parts required, fractional parts optional).
    read(record, "EPOCH_INT", epoch_int);
    read(record, "EPOCH_FRAC", epoch_frac, 0.);
    read(record, "TOABARY_INT", toa_int);
    read(record, "TOABARY_FRAC", toa_frac, 0.);

    // Combine separate parts of epoch and toa to get single values.
    m_epoch = AbsoluteTime(MjdRep("TDB", epoch_int, epoch_frac));
    AbsoluteTime toa(MjdRep("TDB", toa_int, toa_frac));

    // Read the start time of validity window (required).
    long valid_since_date = 0;
    read(record, "VALID_SINCE", valid_since_date);
    m_since = AbsoluteTime(MjdRep("TDB", valid_since_date, 0.));

    // Read the end time of validity window (required).
    // Note: One is added to the endpoint because the "VALID_UNTIL" field in the file expires at the end of that day,
    // whereas the valid_until argument to the ephemeris object is the absolute cutoff.
    long valid_until_date = 0;
    read(record, "VALID_UNTIL", valid_until_date);
    m_until = AbsoluteTime(MjdRep("TDB", valid_until_date + 1, 0.));

    // Read the sky position and frequency coefficients (RA, Dec, F0: required, F1, F2: optional).
    read(record, "RA",  m_ra);
    read(record, "Dec", m_dec);
    read(record, "F0",  m_f0);
    read(record, "F1",  m_f1 , 0.);
    read(record, "F2",  m_f2 , 0.);

    // Create temporary copy of this ephemeris with phi0 == 0.
    FrequencyEph tmp("TDB", m_since, m_until, m_epoch, m_ra, m_dec, 0., m_f0, m_f1, m_f2);

    // Use the temporary ephemeris to compute the phase from the negative of the toa field.
    m_phi0 = - tmp.calcPulsePhase(toa);

    // Make sure it is in the range [0, 1). calcPulsePhase is bounded in this way.
    if (0. > m_phi0) m_phi0 += 1.;
  }

  double FrequencyEph::calcElapsedSecond(const AbsoluteTime & at1, const AbsoluteTime & at2) const {
    return (at1 - at2).computeElapsedTime(m_system->getName()).getTime().getValue(Sec).getDouble();
  }

  double PeriodEph::calcElapsedSecond(const AbsoluteTime & at1, const AbsoluteTime & at2) const {
    return (at1 - at2).computeElapsedTime(m_system->getName()).getTime().getValue(Sec).getDouble();
  }

  void FrequencyEph::writeModelParameter(st_stream::OStream & os) const {
    MjdRep mjd_rep(m_system->getName(), 0, 0.);
    mjd_rep = m_epoch;
    os << format("Epoch", mjd_rep) << std::endl;
    os << format("RA",    m_ra)    << std::endl;
    os << format("Dec",   m_dec)   << std::endl;
    os << format("Phi0",  m_phi0)  << std::endl;
    os << format("F0",    m_f0)    << std::endl;
    os << format("F1",    m_f1)    << std::endl;
    os << format("F2",    m_f2);
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

  double FrequencyEph::calcPulsePhase(const AbsoluteTime & ev_time, double phase_offset) const {
    // Compute a pulse phase value.
    double dt = calcElapsedSecond(ev_time);
    double dt_squared = dt * dt;
    double phase = m_phi0 + m_f0 * dt + m_f1/2.0 * dt_squared + m_f2/6.0 * dt * dt_squared;

    // Express phase as a value between 0. and 1., after adding a global phase offset.
    return trimPhaseValue(phase, phase_offset);
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

  double FrequencyEph::calcFrequency(const AbsoluteTime & ev_time, int derivative_order) const {
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

  std::pair<double, double> FrequencyEph::calcSkyPosition(const AbsoluteTime & /* ev_time */) const {
    return std::make_pair(m_ra, m_dec);
  }

  std::pair<double, double> PeriodEph::calcSkyPosition(const AbsoluteTime & /* ev_time */) const {
    return std::make_pair(m_ra, m_dec);
  }
}
