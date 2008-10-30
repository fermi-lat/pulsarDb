/** \file SimpleDdEph.h
    \brief Interface for SimpleDdEph class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_SimpleDdEph_h
#define pulsarDb_SimpleDdEph_h

#include <vector>

#include "tip/Table.h"

#include "pulsarDb/OrbitalEph.h"

#include "timeSystem/AbsoluteTime.h"

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

  /** \class SimpleDdEph
      \brief Class that represents an orbital ephemeris based on the simplified DD model.
  */
  class SimpleDdEph : public OrbitalEph {
    public:
      typedef std::vector<double>::size_type size_type;

      static const double s_one_pi;
      static const double s_two_pi;
      static const double s_rad_per_deg;
      static const double s_sec_per_day;
      static const double s_sec_per_year;
      static const double s_rad_year_per_deg_sec;
      static const double s_sec_per_microsec;

      /** \brief Construct a SimpleDdEph object from parameter values.
          \param time_system_name Name of time system in which this orbital ephemeris is defined.
          \param pb Orbital period in seconds.
          \param pb_dot First time derivative of the orbital period (dimension-less)
          \param a1 Projected semi-major axis in light seconds.
          \param x_dot First time derivative of the projected semi-major axis in light-seconds per second.
          \param ecc Orbital eccentricity (dimension-less)
          \param ecc_dot First time derivative of eccentricity in inverse of seconds.
          \param om Longitude of periastron in degrees.
          \param om_dot First Time derivative of periastron longitude in degrees per year.
          \param t0 Barycentric time of periastron.
          \param gamma Time-dilation and gravitational redshift parameter (dimension-less)
          \param shapiro_r Range parameter of Shapiro delay in binary system in micro-seconds.
          \param shapiro_s Shape parameter of Shapiro delay in binary system (dimension-less)
      */
      SimpleDdEph(const std::string & time_system_name, double pb, double pb_dot, double a1, double x_dot, double ecc, double ecc_dot,
        double om, double om_dot, const timeSystem::AbsoluteTime & t0, double gamma, double shapiro_r, double shapiro_s);

      /** \brief Construct a SimpleDdEph object from a FITS row.
          \param record FITS row from which orbital parameters are to be read.
          \param header FITS header to read other information if necessary (not used).
      */
      SimpleDdEph(const tip::Table::ConstRecord & record, const tip::Header & header);

      /// \brief Destruct this SimpleDdEph object.
      virtual ~SimpleDdEph();

      /// \brief Return a time system in which binary demodulation is performed.
      virtual const timeSystem::TimeSystem & getSystem() const { return *m_system; }

      /// \brief Return the T0 parameter value, which is the barycentric time of periastron.
      virtual const timeSystem::AbsoluteTime & t0() const { return m_t0; }

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
      /** \brief Compute the number of elapsed seconds since the barycentric time of periastron (T0 parameter), and return it.
          \param at Absolute time for which the number of elapsed seconds is to be computed.
      */
      virtual double calcElapsedSecond(const timeSystem::AbsoluteTime & at) const;

      const timeSystem::TimeSystem * m_system;
      double m_pb, m_pb_dot;
      double m_a1, m_x_dot;
      double m_ecc, m_ecc_dot;
      double m_om, m_om_dot;
      timeSystem::AbsoluteTime m_t0;
      double m_gamma;
      double m_shapiro_r, m_shapiro_s;
  };
}

#endif
