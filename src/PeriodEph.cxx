/** \file PeriodEph.cxx
    \brief Implementation of the PeriodEph class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "pulsarDb/PulsarEph.h"
#include "pulsarDb/PeriodEph.h"

#include "timeSystem/CalendarFormat.h"
#include "timeSystem/MjdFormat.h"
#include "timeSystem/SourcePosition.h"

using namespace timeSystem;

namespace {
  /** \brief Helper function to compute coefficients to compute a frequency derivative of an arbitrary degree
             of a frequency trend expressed by a polinomial of period derivatives.
      \param em Parameter for a coefficient to compute, which corresponds to the degree of frequency derivative to compute.
      \param el Parameter for a coefficient to compute, which corresponds to the degree of period derivative to compute from.
  */
  int computeCoeff(unsigned int em, unsigned int el) {
    if (el == 0) return 1;
    else if (em == 2*el - 1) return 0;
    else return computeCoeff(em-1, el) + (em - 2*el + 1)*computeCoeff(em-1, el-1);
  }
}

namespace pulsarDb {

  PeriodEph::PeriodEph(const tip::Table::ConstRecord & record, const tip::Header & /* header */):
    m_system(&TimeSystem::getSystem("TDB")), m_since("TDB", 0, 0.), m_until("TDB", 0, 0.), m_epoch("TDB", 0, 0.),
    m_ra(0.), m_dec(0.), m_phi0(0.), m_p0(0.), m_p1(0.), m_p2(0.) {
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
    m_epoch = AbsoluteTime("TDB", Mjd(epoch_int, epoch_frac));
    AbsoluteTime toa("TDB", Mjd(toa_int, toa_frac));

    // Read the start time of validity window (required).
    long valid_since_date = 0;
    read(record, "VALID_SINCE", valid_since_date);
    m_since = AbsoluteTime("TDB", Mjd(valid_since_date, 0.));

    // Read the end time of validity window (required).
    // Note: One is added to the endpoint because the "VALID_UNTIL" field in the file expires at the end of that day,
    // whereas the valid_until argument to the ephemeris object is the absolute cutoff.
    long valid_until_date = 0;
    read(record, "VALID_UNTIL", valid_until_date);
    m_until = AbsoluteTime("TDB", Mjd(valid_until_date + 1, 0.));

    // Read the sky position and frequency coefficients (RA, DEC, F0: required, F1, F2: optional).
    read(record, "RA",  m_ra);
    read(record, "DEC", m_dec);
    read(record, "P0",  m_p0);
    read(record, "P1",  m_p1 , 0.);
    read(record, "P2",  m_p2 , 0.);

    // Create temporary copy of this ephemeris with phi0 == 0.
    PeriodEph tmp("TDB", m_since, m_until, m_epoch, m_ra, m_dec, 0., m_p0, m_p1, m_p2);

    // Use the temporary ephemeris to compute the phase from the negative of the toa field.
    m_phi0 = - tmp.calcPulsePhase(toa);

    // Make sure it is in the range [0, 1). calcPulsePhase is bounded in this way.
    if (0. > m_phi0) m_phi0 += 1.;
  }

  void PeriodEph::writeModelParameter(st_stream::OStream & os) const {
    std::string epoch_string;
    try {
      epoch_string = m_epoch.represent(m_system->getName(), MjdFmt);
    } catch (const std::exception &) {
      epoch_string = m_epoch.represent(m_system->getName(), CalendarFmt);
    }
    os << format("Epoch", epoch_string, "") << std::endl;
    os << format("RA",    m_ra,         "degrees") << std::endl;
    os << format("Dec",   m_dec,        "degrees") << std::endl;
    os << format("Phi0",  m_phi0,       "")        << std::endl;
    os << format("P0",    m_p0,         "s")       << std::endl;
    os << format("P1",    m_p1,         "")        << std::endl;
    os << format("P2",    m_p2,         "s**(-1)");
  }

  double PeriodEph::calcPulsePhase(const AbsoluteTime & ev_time, double phase_offset) const {
    // Set a pulse phase at the reference epoch.
    double phase = m_phi0;

    // Compute an elapsed time in seconds.
    double dt = calcElapsedSecond(ev_time);

    // Set the error message for problems in integration.
    static const std::string s_integral_error("PeriodEph: pulse period predicted to be zero at a time between the reference epoch and the time for which a pulse phase is to be computed");

    // Compute contribution from p0, p1, and p2.
    bool p0_is_zero = (0. == m_p0);
    bool p1_is_zero = (0. == m_p1);
    bool p2_is_zero = (0. == m_p2);
    if (p0_is_zero && p1_is_zero && p2_is_zero) {
      // Throw an exception for p0 == p1 == p2 == 0.
      throw std::runtime_error("Unphysical ephemeris is given: all coefficients are zeros (p0 = p1 = p2 = 0.)");

    } else if (p1_is_zero && p2_is_zero) {
      // Compute pulse phase for p0 != 0 and p1 == p2 == 0.
      phase += dt / m_p0;

    } else if (p2_is_zero) {
      // Compute pulse phase for any p0, p1 != 0, and p2 == 0.
      if (m_p0 * (m_p0 + m_p1 * dt) > 0.) {
        // Note: In this branch, it is guaranteed that m_p0 is not zero.
        phase += std::log(m_p1 * dt / m_p0 + 1.) / m_p1;

      } else {
        throw std::runtime_error(s_integral_error);
      }

    } else {
      // Compute pulse phase for any p0, any p1, and p2 != 0.
      double two_p0p2_minus_p1sq = 2. * m_p0 * m_p2 - m_p1 * m_p1;
      if (two_p0p2_minus_p1sq > 0.) {
        double sqrt_term = std::sqrt(two_p0p2_minus_p1sq);
        double u0 = m_p1 / sqrt_term;
        double u1 = (m_p1 + m_p2 * dt) / sqrt_term;
        if ((u0 > 1. && u1 > 1.) || (u0 < -1.) && (u1 < -1.)) {
          phase += 2. / sqrt_term * std::atan(sqrt_term * dt / (2. * m_p0 + m_p1 * dt));

        } else {
          phase += 2. / sqrt_term * (std::atan(u1) - std::atan(u0));
        }

      } else if (two_p0p2_minus_p1sq == 0.) {
        double denominator = m_p1 * (m_p1 + m_p2 * dt);
        if (denominator > 0.) {
          phase += 2 * m_p2 * dt / denominator;

        } else {
          throw std::runtime_error(s_integral_error);
        }

      } else {
        double sqrt_term = std::sqrt(-two_p0p2_minus_p1sq);
        double x_plus = -(m_p1 + sqrt_term) / m_p2;
        double x_minus = -(m_p1 - sqrt_term) / m_p2;
        if (x_plus * (x_plus - dt) > 0. && x_minus * (x_minus - dt) > 0.) {
          double numerator = 2. * m_p0 + (m_p1 + sqrt_term)* dt;
          double denominator = 2. * m_p0 + (m_p1 - sqrt_term)* dt;
          phase += std::log(std::fabs(numerator/denominator)) / sqrt_term;

        } else {
          throw std::runtime_error(s_integral_error);
        }
      }
    }

    // Express phase as a value between 0. and 1., after adding a global phase offset.
    return trimPhaseValue(phase, phase_offset);
  }

  double PeriodEph::calcFrequency(const AbsoluteTime & ev_time, int derivative_order) const {
    double return_value = 0.;

    // Compute period and its derivatives at the requested time.
    double dt = calcElapsedSecond(ev_time);
    double q0 = m_p0 + m_p1 * dt + m_p2 / 2. * dt * dt;
    double q1 = m_p1 + m_p2 * dt;

    // Compute the requested value.
    if (derivative_order == 0) {
      return_value = 1. / q0;

    } else if (derivative_order > 0) {
      // Precompute frequently-used values.
      int kk_max = derivative_order / 2;
      double q1sq = q1 * q1;

      // Compute the product of factorial, q0, and q1 for the maximum index of p2 (kk_max).
      int factorial = 1;
      double q0q1_component = (derivative_order % 2 ? q1 : 1.) / q0;
      for (int ii = 1; ii <= derivative_order - kk_max; ++ii) {
        factorial *= ii;
        q0q1_component /= -q0;
      }

      // Compute products of factorial, q0, and q1 for all terms, and store them in an array.
      std::vector<double> precomputed(kk_max+1, 1.);
      int jj = derivative_order - kk_max + 1;
      for (std::vector<double>::reverse_iterator itor = precomputed.rbegin(); itor != precomputed.rend(); ++itor, ++jj) {
        // Compute and store the product for this index.
        *itor = factorial * q0q1_component;

        // Update the components for the next iteration.
        q0q1_component *= q1sq / (-q0);
        factorial *= jj;
      }

      // Compute the powers of p2 and the integer coefficient, multiply them to the stored products, and sum them up.
      return_value = 0.;
      double p2_component = 1.;
      int kk = 0;
      for (std::vector<double>::iterator itor = precomputed.begin(); itor != precomputed.end(); ++itor, ++kk) {
        return_value += *itor * p2_component * computeCoeff(derivative_order, kk);

        // Update the components for the next iteration.
        p2_component *= m_p2;
      }

    } else {
      std::ostringstream os;
      os << "Negative order of period derivative is given: " << derivative_order;
      throw std::runtime_error(os.str());
    }

    // Return the computed value.
    return return_value;
  }

  SourcePosition PeriodEph::calcPosition(const AbsoluteTime & /* ev_time */) const {
    return SourcePosition(m_ra, m_dec);
  }

}
