/** \file Ell1ModelEph.h
    \brief Interface for Ell1ModelEph class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_Ell1ModelEph_h
#define pulsarDb_Ell1ModelEph_h

#include <vector>

#include "tip/Table.h"

#include "pulsarDb/OrbitalEph.h"

#include "timeSystem/AbsoluteTime.h"
#include "timeSystem/TimeInterval.h"

namespace st_stream {
  class OStream;
}

namespace timeSystem {
  class ElapsedTime;
}

namespace tip {
  class Header;
}

namespace pulsarDb {

  /** \class Ell1ModelEph
      \brief Class that represents an orbital ephemeris based on the ELL1 model by Ch. Lange, F. Camilo, N. Wex,
             M. Kramer, D. C. Backer, A. G. Lyne, and O.Doroshenko, MNRAS 326, 274 (2001).
  */
  class Ell1ModelEph : public OrbitalEph {
    public:
      typedef std::vector<double>::size_type size_type;

      static const double s_two_pi;
      static const double s_sec_per_microsec;

      /** \brief Construct a Ell1ModelEph object from parameter values.
          \param time_system_name Name of time system in which this orbital ephemeris is defined.
          \param pb Orbital period in seconds.
          \param pb_dot First time derivative of the orbital period (dimensionless).
          \param a1 Projected semi-major axis in light-seconds.
          \param x_dot First time derivative of the projected semi-major axis in light-seconds per second.
          \param eps1 Eccentricity multiplied by the sine of the periastron longitude (dimensionless).
          \param eps1_dot First time derivative of eps1 in inverse of seconds.
          \param eps2 Eccentricity multiplied by the cosine of the periastron longitude (dimensionless).
          \param eps2_dot First time derivative of eps2 in inverse of seconds.
          \param tasc Time of ascending node.
          \param shapiro_r Range parameter of Shapiro delay in binary system in microseconds.
          \param shapiro_s Shape parameter of Shapiro delay in binary system (dimensionless).
      */
      Ell1ModelEph(const std::string & time_system_name, double pb, double pb_dot, double a1, double x_dot, double eps1,
        double eps1_dot, double eps2, double eps2_dot, const timeSystem::AbsoluteTime & tasc, double shapiro_r, double shapiro_s);

      /** \brief Construct a Ell1ModelEph object from a FITS row.
          \param record FITS row from which orbital parameters are to be read.
          \param header FITS header to read other information if necessary (not used).
      */
      Ell1ModelEph(const tip::Table::ConstRecord & record, const tip::Header & header);

      /// \brief Destruct this Ell1ModelEph object.
      virtual ~Ell1ModelEph() {}

      /// \brief Return a time system in which binary demodulation is performed.
      virtual const timeSystem::TimeSystem & getSystem() const { return *m_system; }

      /// \brief Return the Tasc parameter value, which is the time of ascending node.
      virtual const timeSystem::AbsoluteTime & t0() const { return m_tasc; }

      /** \brief Compute an orbital phase of a given time, and return it.
          \param ev_time Absolute time for which an orbital phase is to be computed.
          \param phase_offset Value to be added to an orbital phase. This value is added to the computed orbital phase
                 before truncated to a value in range [0, 1).
      */
      virtual double calcOrbitalPhase(const timeSystem::AbsoluteTime & ev_time, double phase_offset = 0.) const;

      /** \brief Compute a propagation delay in a binary system, in reference to the center of gravity of the system.
                 Note that the propagation delay includes not only a light-travel time, but also gravitational delays.
          \param ev_time Absolute time for which a propagation delay is to be computed.
      */
      virtual timeSystem::ElapsedTime calcOrbitalDelay(const timeSystem::AbsoluteTime & ev_time) const;

      /// \brief Create a copy of this object, and return a pointer to it.
      virtual OrbitalEph * clone() const;

    protected:
      /** \brief Write a text representation of this object to an output stream.
          \param os Output stream to write a text representation of this object to.
      */
      virtual void writeModelParameter(st_stream::OStream & os) const;

    private:
      /** \brief Compute the number of elapsed seconds since the time of ascending node (Tasc parameter), and return it.
          \param at Absolute time for which the number of elapsed seconds is to be computed.
      */
      inline double calcElapsedSecond(const timeSystem::AbsoluteTime & at) const {
        return (at - m_tasc).computeDuration(m_system->getName(), "Sec");
      }

      /** \brief Compute the orbital phase measured from the time of the ascending node, multiplied by two pi, and
                 return it in the unit of radians. The return value is the quantity denoted by a capital phi in
                 Eq. A9 in Lange, et al., MNRAS 326, 274 (2001).
          \param elapsed_second The number of elapsed seconds for which the orbital parameters are to be computed.
      */
      double calcLargePhi(double elapsed_second) const;

      const timeSystem::TimeSystem * m_system;
      double m_pb, m_pb_dot;
      double m_a1, m_x_dot;
      double m_eps1, m_eps1_dot;
      double m_eps2, m_eps2_dot;
      timeSystem::AbsoluteTime m_tasc;
      double m_shapiro_r, m_shapiro_s;
  };

}

#endif
