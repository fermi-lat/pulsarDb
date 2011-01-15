/** \file Ell1ModelEph.cxx
    \brief Implementation of the Ell1ModelEph class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#include <cmath>
#include <limits>
#include <iomanip>
#include <stdexcept>

#include "pulsarDb/OrbitalEph.h"
#include "pulsarDb/Ell1ModelEph.h"

#include "timeSystem/AbsoluteTime.h"
#include "timeSystem/CalendarFormat.h"
#include "timeSystem/Duration.h"
#include "timeSystem/ElapsedTime.h"
#include "timeSystem/MjdFormat.h"

using namespace timeSystem;

namespace tip {
  class Header;
}

namespace pulsarDb {

  const double Ell1ModelEph::s_two_pi = 2. * M_PI;
  const double Ell1ModelEph::s_sec_per_microsec = 1.e-6;

  Ell1ModelEph::Ell1ModelEph(const std::string & time_system_name, double pb, double pb_dot, double a1, double x_dot,
    double eps1, double eps1_dot, double eps2, double eps2_dot, const AbsoluteTime & tasc, double shapiro_r, double shapiro_s):
    OrbitalEph(ElapsedTime(time_system_name, Duration(10.e-9, "Sec")), 100), m_system(&TimeSystem::getSystem(time_system_name)),
    m_pb(pb), m_pb_dot(pb_dot), m_a1(a1), m_x_dot(x_dot), m_eps1(eps1), m_eps1_dot(eps1_dot), m_eps2(eps2), m_eps2_dot(eps2_dot),
    m_tasc(tasc), m_shapiro_r(shapiro_r * s_sec_per_microsec), m_shapiro_s(shapiro_s) {}

  Ell1ModelEph::Ell1ModelEph(const tip::Table::ConstRecord & record, const tip::Header & /* header */):
    OrbitalEph(ElapsedTime("TDB", Duration(10.e-9, "Sec")), 100), m_system(&TimeSystem::getSystem("TDB")), m_tasc("TDB", 0, 0.) {
    // Get parameters from record.
    // Required fields: PB, A1, EPS1, EPS2, TASC.
    // Optional fields: PBDOT, XDOT, EPS1DOT, EPS2DOT, SHAPIRO_R, SHAPIRO_S.
    read(record, "PB",        m_pb);
    read(record, "PBDOT",     m_pb_dot, 0.);
    read(record, "A1",        m_a1);
    read(record, "XDOT",      m_x_dot, 0.);
    read(record, "EPS1",      m_eps1);
    read(record, "EPS1DOT",   m_eps1_dot, 0.);
    read(record, "EPS2",      m_eps2);
    read(record, "EPS2DOT",   m_eps2_dot, 0.);
    double dbl_tasc = 0.;
    read(record, "TASC",      dbl_tasc);
    read(record, "SHAPIRO_R", m_shapiro_r, 0.);
    read(record, "SHAPIRO_S", m_shapiro_s, 0.);

    // Create an AbsoluteTime object from the value of "TASC" column.
    m_tasc = AbsoluteTime("TDB", Mjd1(dbl_tasc));

    // Adjust units.
    m_shapiro_r *= s_sec_per_microsec;
  }

  void Ell1ModelEph::writeModelParameter(st_stream::OStream & os) const {
    os << format("PB",        m_pb,      "s")         << std::endl;
    os << format("PBDOT",     m_pb_dot,  "")          << std::endl;
    os << format("A1",        m_a1,      "lt-s")      << std::endl;
    os << format("XDOT",      m_x_dot,   "lt-s/s")    << std::endl;
    os << format("EPS1",      m_eps1,     "")         << std::endl;
    os << format("EPS1DOT",   m_eps1_dot, "s**(-1)")  << std::endl;
    os << format("EPS2",      m_eps2,     "")         << std::endl;
    os << format("EPS2DOT",   m_eps2_dot, "s**(-1)")  << std::endl;
    std::string tasc_string;
    try {
      tasc_string = m_tasc.represent(m_system->getName(), MjdFmt);
    } catch (const std::exception &) {
      tasc_string = m_tasc.represent(m_system->getName(), CalendarFmt);
    } 
    os << format("TASC",      tasc_string, "")  << std::endl;
    os << format("SHAPIRO_R", m_shapiro_r, "s") << std::endl;
    os << format("SHAPIRO_S", m_shapiro_s, "");
  }

  OrbitalEph * Ell1ModelEph::clone() const { return new Ell1ModelEph(*this); }

  double Ell1ModelEph::calcOrbitalPhase(const AbsoluteTime & ev_time, double phase_offset) const {
    // Compute elapsed time from epoch of periastron in seconds.
    double elapsed_second = calcElapsedSecond(ev_time);

    // Compute orbital phase multiplied by two pi.
    double large_phi = calcLargePhi(elapsed_second);

    // Compute the complete phase.
    double phase = large_phi / Ell1ModelEph::s_two_pi;

    // Express phase as a value between 0. and 1., after adding a global phase offset.
    return trimPhaseValue(phase, phase_offset);
  }

  ElapsedTime Ell1ModelEph::calcOrbitalDelay(const AbsoluteTime & ev_time) const {
    // Compute elapsed time from epoch of periastron in seconds.
    double elapsed_second = calcElapsedSecond(ev_time);

    // Compute projected semimajor axis.
    // --- Eq. A10 in Lange, et al., MNRAS 326, 274 (2001)
    double semiax = m_a1 + m_x_dot * elapsed_second;

    // Compute the first Laplace-Lagrange parameter, eta, at the given event time.
    // --- Eq. A10 in Lange, et al., MNRAS 326, 274 (2001)
    double eta = m_eps1 + m_eps1_dot * elapsed_second;

    // Compute the second Laplace-Lagrange parameter, kappa, at the given event time.
    // --- Eq. A10 in Lange, et al., MNRAS 326, 274 (2001)
    double kappa = m_eps2 + m_eps2_dot * elapsed_second;

    // Compute orbital phase multiplied by two pi.
    double large_phi = calcLargePhi(elapsed_second);

    // Compute time delays due to orbital motion.
    // --- Eqs. A6 and A16 in Lange, et al., MNRAS 326, 274 (2001)
    double roemer = semiax * (std::sin(large_phi) + kappa / 2.0 * std::sin(2.0*large_phi) - eta / 2.0 * std::cos(2.0*large_phi));
    double shapiro = - 2.0 * m_shapiro_r * std::log(1.0 - m_shapiro_s*std::sin(large_phi));

    // Return total delay.
    // --- Eq. A1 in Lange, et al., MNRAS 326, 274 (2001), with the Shapiro delay added
    return ElapsedTime(m_system->getName(), Duration(roemer + shapiro, "Sec"));
  }

  double Ell1ModelEph::calcLargePhi(double elapsed_second) const {
    // Rename the Laplace-Lagrange parameters, to follow notations in Lange, et al., MNRAS 326, 274 (2001).
    const double & eta_zero = m_eps1;
    const double & eta_dot  = m_eps1_dot;
    const double & kappa_zero = m_eps2;
    const double & kappa_dot  = m_eps2_dot;

    // Compute the periastron longitude at the epoch, omega_0, and its time derivative, omega_dot.
    // --- Drived from Eq. A8 in Lange, et al., MNRAS 326, 274 (2001)
    double omega_zero = 0.;
    double omega_dot = 0.;
    if (eta_zero != 0. || kappa_zero != 0.) {
      // Compute the periastron longitude at the epoch, omega_0, and put it in range [0, 2*pi).
      omega_zero = std::atan2(eta_zero, kappa_zero);
      if (omega_zero < 0.0) omega_zero += Ell1ModelEph::s_two_pi;

      // Compute the time derivative of the periastron longitude from the Laplace-Lagrange parameters.
      omega_dot = (eta_dot * kappa_zero - kappa_dot * eta_zero) / (eta_zero * eta_zero + kappa_zero * kappa_zero);

    } else {
      // Note: In this case, eta_zero = kappa_zero = 0.0, hence eta and kappa at an arbitrary time lie
      //       on a straight line passing through the origin of (eta, kappa)-plane. Mathematically, this
      //       makes omega remains constant over time, with the exception of a sudden increase/decrease
      //       by 180 degrees when point (eta, kappa) passes through the origin. Physically, the sudden
      //       change in omega is unphysical and simply an artifact of this particular parametrization.
      //       In practical situations, omega stays constant throughout times of our interest, which is
      //       assumed in the computations below.
      omega_dot = 0.;

      if (eta_dot != 0. || kappa_dot != 0.) {
        // Compute the periastron longitude at the given time, put it in range [0, 2*pi), and take it as
        // the value at the epoch (see the note above).
        omega_zero = std::atan2(eta_dot, kappa_dot);
        if (omega_zero < 0.0) omega_zero += Ell1ModelEph::s_two_pi;

      } else {
        // In this case, eccentricity is zero (0) at all times, and the periastron longitude cannot
        // be uniquely defined. Take a zero (0) as an arbitrary choice here.
        omega_zero = 0.;
      }
    }

    // Compute the original orbital angular frequency and its time derivative.
    // Note: An orbital angular frequency is the reciprocal of an orbital period, multiplied by two pi.
    // --- Eq. A4 in Lange, et al., MNRAS 326, 274 (2001).
    double nb = Ell1ModelEph::s_two_pi / m_pb;
    double nb_dot = - m_pb_dot / m_pb * nb;

    // Compute the apparent orbital angular frequency, and its time derivative.
    // --- Eqs. A11, A12, and A13 in Lange, et al., MNRAS 326, 274 (2001)
    double nb_bar = nb + omega_dot - nb_dot * omega_zero / (nb + omega_dot);
    double nb_bar_dot = nb_dot;

    // Compute the orbital phase measured from the time of the ascending node, multiplied by two pi.
    // --- Eq. A9 in Lange, et al., MNRAS 326, 274 (2001)
    double large_phi = (nb_bar + nb_bar_dot / 2.0 * elapsed_second) * elapsed_second;

    // Return the computed result.
    return large_phi;
  }

}
