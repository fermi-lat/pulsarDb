/** \file HighPrecisionEph.h
    \brief Interface for HighPrecisionEph class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_HighPrecisionEph_h
#define pulsarDb_HighPrecisionEph_h

#include <list>
#include <string>
#include <utility>
#include <vector>

#include "pulsarDb/EphStatus.h"
#include "pulsarDb/FormattedEph.h"
#include "pulsarDb/PulsarEph.h"

#include "tip/Table.h"

#include "timeSystem/AbsoluteTime.h"
#include "timeSystem/TimeSystem.h"

namespace tip {
  class Header;
}

namespace st_stream {
  class OStream;
}

namespace pulsarDb {

  /** \class HighPrecisionEph
      \brief Class representing a single pulsar ephemeris for high-precision pulsar timing.
  */
  class HighPrecisionEph : public PulsarEph {
    public:
      typedef std::vector<double> freq_type;
      typedef std::vector<double> wave_type;
      typedef std::vector<double> jump_type;
      typedef std::list<std::pair<double, double> > decay_type;

      static const double s_rad_per_deg; // Radians per degree.
      static const double s_rad_per_mas; // Radians per milliarcsecond.
      static const double s_rad_day_per_mas_year; // "Radians per day" per "milliarcseconds per Julian year (365.25 days)"
      static const double s_km_per_lts; // Kilometers per light-second (i.e., one light-second in kilometers).
      static const double s_lts_per_au; // Light-seconds per astronomical unit (i.e., one astronomical unit in light-seconds).

      /** \class GlitchParameters
          \brief Helper class to store all parameters for one pulsar glitch, to be used by HighPrecisionEph class.
      */
      struct GlitchParameter {
        timeSystem::AbsoluteTime m_epoch;
        jump_type m_perm_jump;
        decay_type m_decay_comp;
        GlitchParameter(): m_epoch("TDB", 0, 0.), m_perm_jump(0), m_decay_comp(0) {}
      };
      typedef std::list<GlitchParameter> glitch_type;

      /** \brief Create a pulsar ephemeris object with the given parameters.
          \param time_system_name Name of time system in which this pulsar ephemeris is defined.
          \param valid_since Beginning of time period during which this ephemeris is considered valid.
          \param valid_until End of time period during which this ephemeris is considered valid.
          \param pos_epoch Reference epoch of the proper motion.
          \param ra Right Ascension of the pulsar in degrees at a time of reference (pos_epoch).
          \param dec Declination of the pulsar in degrees at a time of reference (pos_epoch).
          \param ra_velocity Angular change in Right Ascension in milliarcseconds per Julian year (365.25 days).
          \param dec_velocity Angular change in Declination in milliarcseconds per Julian year (365.25 days).
          \param radial_velocity Velocity in the direction of the line of sight in kilometers per second.
          \param parallax Annual parallax in milliarcseconds. Give a non-positive value if unknown.
          \param freq_epoch Reference epoch of frequency parameters (freq).
          \param freq_pars Pulse phase (dimensionless), pulse frequency in the units of s^(-1), and its time derivatives
                 in the units of s^(-2), s^(-3), etc., at the given epoch (freq_epoch).
          \param wave_omega Fundamental frequency of sinusoidal timing residual in the units of radians per day.
          \param wave_sine List of the amplitudes of the sine components in seconds, whose i-th element is
                 the coefficient of the sine term for the i-th sinusoid.
          \param wave_cosine List of the amplitudes of the cosine components in seconds, whose i-th element is
                 the coefficient of the sine term for the i-th sinusoid.
          \param glitch_list List of pulsar glitches. Each element of this list is a GlitchParameter object,
                 whose public data members must contain the following glitch parameters.
                   1) m_epoch: Reference epoch of the glitch, after which glitch parameters will be in effect.
                   2) m_perm_jump: Permanent increment in pulse phase, pulse frequency, and its time derivatives, where
                      m_perm_jump[0] is a permanent increment in pulse phase at the glitch epoch (dimensionless),
                      m_perm_jump[1] is a permanent increment in pulse frequency at the glitch epoch in the units of s^(-1),
                      m_perm_jump[2] that in its first time derivative in the units of s^(-2), and so on.
                   3) m_decay_comp: List of parameters of decaying components of pulse frequency increments. Each element of
                      this argument is a pair of the amplitude and the decay time of a decaying increment in pulse frequency,
                      where the amplitude is in the units of s^(-1) and the decay time in days.
      */
      HighPrecisionEph(const std::string & time_system_name, const timeSystem::AbsoluteTime & valid_since,
        const timeSystem::AbsoluteTime & valid_until, const timeSystem::AbsoluteTime & pos_epoch, double ra, double dec,
        double ra_velocity, double dec_velocity, double radial_velocity, double parallax,
        const timeSystem::AbsoluteTime & freq_epoch, const freq_type & freq_pars, double wave_omega,
        const wave_type & wave_sine, const wave_type & wave_cosine, const glitch_type & glitch_list):
        m_system(&timeSystem::TimeSystem::getSystem(time_system_name)), m_since(valid_since), m_until(valid_until),
        m_pos_epoch(pos_epoch), m_ra(ra), m_dec(dec), m_ra_vel(ra_velocity), m_dec_vel(dec_velocity),
        m_radial_vel(radial_velocity), m_parallax(parallax), m_freq_epoch(freq_epoch), m_freq_pars(freq_pars),
        m_wave_omega(wave_omega), m_wave_sine(wave_sine), m_wave_cosine(wave_cosine), m_glitch_list(glitch_list),
        m_remark_cont() { setRemark(); }

      /** \brief Create a pulsar ephemeris object with the parameters stored in tip record.
          \param record FITS row from which all parameters for an ephemeris being created are to be read.
          \param header FITS header to read other information if necessary (not used).
      */
      HighPrecisionEph(const tip::Table::ConstRecord & record, const tip::Header & header);

      /// \brief Destruct this HighPrecisionEph object.
      virtual ~HighPrecisionEph() {}

      /// \brief Return a time system to be used to interpret return values of calcFrequency method.
      virtual const timeSystem::TimeSystem & getSystem() const { return *m_system; }

      /** \brief Return a start time of a time interval, during which this ephemeris is considered valid.
                 Note: The ephemeris is also considered valid on the start time itself.
      */
      virtual const timeSystem::AbsoluteTime & getValidSince() const { return m_since; }

      /** \brief Return an end time of a time interval, during which this ephemeris is considered valid.
                 Note: The ephemeris is also considered valid on the end time itself.
      */
      virtual const timeSystem::AbsoluteTime & getValidUntil() const { return m_until; }

      /// \brief Return a reference epoch of this ephemeris.
      virtual const timeSystem::AbsoluteTime & getEpoch() const { return m_freq_epoch; }

      /// \brief Return the container of ephemeris remarks.
      virtual const EphStatusCont & getRemark() const { return m_remark_cont; }

      /// \brief Create a copy of this object, and return a pointer to it.
      virtual PulsarEph * clone() const { return new HighPrecisionEph(*this); }

      /** \brief Compute a pulse phase of a given time, and return it.
                 Note: validity of the ephemeris (valid since and valid until) are not checked.
          \param ev_time Absolute time for which a pulse phase is to be computed.
          \param phase_offset Value to be added to a pulse phase. This value is added to the computed pulse phase
                 before truncated to a value in range [0, 1).
      */
      virtual double calcPulsePhase(const timeSystem::AbsoluteTime & ev_time, double phase_offset = 0.) const;

      /** \brief Compute a pulse frequency or its time derivative at a given time in the time system given to
                 this object upon its construction. Call getSystem method to obtain the time system to interpret
                 the return value of this method. The unit of frequency and its derivatives must be derived from
                 the time unit of seconds, i.e., Hz (sE-1), Hz/s (sE-2), and so on.
                 Note: validity of the ephemeris (valid since and valid until) are not checked. 
          \param ev_time Absolute time at which a pulse frequency is to be computed.
          \param derivative_order Order of frequency derivative to be computed. Set zero (0) to this argument
                 to compute a pulse frequency, one (1) the first time derivative of pulse frequency, and so on.
      */
      virtual double calcFrequency(const timeSystem::AbsoluteTime & ev_time, int derivative_order = 0) const;

      /** \brief Compute the pulsar position at a given time, and return it.
                 Note: validity of the ephemeris (valid since and valid until) are not checked.
          \param ev_time Absolute time at which the pulsar position is to be computed.
      */
      virtual timeSystem::SourcePosition calcPosition(const timeSystem::AbsoluteTime & ev_time) const;

    protected:
      /// \brief Output text expression of subclass-specific parameters of this object to a given output stream.
      virtual void writeModelParameter(st_stream::OStream & os) const;

    private:
      const timeSystem::TimeSystem * m_system;
      timeSystem::AbsoluteTime m_since;
      timeSystem::AbsoluteTime m_until;
      timeSystem::AbsoluteTime m_pos_epoch;
      double m_ra;
      double m_dec;
      double m_ra_vel;
      double m_dec_vel;
      double m_radial_vel;
      double m_parallax;
      timeSystem::AbsoluteTime m_freq_epoch;
      freq_type m_freq_pars;
      double m_wave_omega;
      wave_type m_wave_sine;
      wave_type m_wave_cosine;
      glitch_type m_glitch_list;
      EphStatusCont m_remark_cont;

      /// \brief Set a list of ephemeris remarks to a data member.
      void setRemark();
  };

}

#endif
