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
#include "timeSystem/TimeSystem.h"

using namespace timeSystem;

namespace {
  int computeCoeff(unsigned int em, unsigned int el) {
    if (el == 0) return 1;
    else if (em == 2*el - 1) return 0;
    else return computeCoeff(em-1, el) + (em - 2*el + 1)*computeCoeff(em-1, el-1);
  }
}

namespace pulsarDb {

  double PeriodEph::calcElapsedSecond(const AbsoluteTime & at1, const AbsoluteTime & at2) const {
    return (at1 - at2).computeDuration(m_system->getName(), "Sec");
  }

  void PeriodEph::writeModelParameter(st_stream::OStream & os) const {
    std::string epoch_string = m_epoch.represent(m_system->getName(), "MJD");
    os << format("Epoch", epoch_string) << std::endl;
    os << format("RA",    m_ra)         << std::endl;
    os << format("Dec",   m_dec)        << std::endl;
    os << format("Phi0",  m_phi0)       << std::endl;
    os << format("P0",    m_p0)         << std::endl;
    os << format("P1",    m_p1)         << std::endl;
    os << format("P2",    m_p2);
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
      throw std::runtime_error("PeriodEph: all coefficients are zeros (p0 = p1 = p2 = 0.)");

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
      throw std::runtime_error("PeriodEph is given a negative order of derivative");
    }

    // Return the computed value.
    return return_value;
  }

  std::pair<double, double> PeriodEph::calcSkyPosition(const AbsoluteTime & /* ev_time */) const {
    return std::make_pair(m_ra, m_dec);
  }
}
