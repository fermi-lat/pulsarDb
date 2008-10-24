/** \file PulsarEph.h
    \brief Interface for PulsarEph class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_PulsarEph_h
#define pulsarDb_PulsarEph_h

#include <utility>
#include <vector>

#include "pulsarDb/FormattedEph.h"

namespace st_stream {
  class OStream;
}

namespace timeSystem {
  class AbsoluteTime;
  class TimeSystem;
}

namespace tip {
  class Header;
}

namespace pulsarDb {

  /** \class PulsarEph
      \brief Base class representing a single pulsar ephemeris.
  */
  class PulsarEph: public FormattedEph {
    public:
      /// \brief Destruct this PulsarEph object.
      virtual ~PulsarEph() {}

      /// \brief Return a time system to be used to interpret return values of calcFrequency method.
      virtual const timeSystem::TimeSystem & getSystem() const = 0;

      /** \brief Return a start time of a time interval, during which this ephemeris is considered valid.
                 Note: The ephemeris is also considered valid on the start time itself.
      */
      virtual const timeSystem::AbsoluteTime & getValidSince() const = 0;

      /** \brief Return an end time of a time interval, during which this ephemeris is considered valid.
                 Note: The ephemeris is also considered valid on the end time itself.
      */
      virtual const timeSystem::AbsoluteTime & getValidUntil() const = 0;

      /// \brief Return a reference epoch of this ephemeris.
      virtual const timeSystem::AbsoluteTime & getEpoch() const = 0;

      /** \brief Compute the spin phase of the given time.
                 Note: validity of the ephemeris (valid since and valid until) are not checked.
          \param ev_time Absolute time for which a pulse phase is to be compted.
          \param phase_offset Value to be added to a pulse phase. This value is added to the computed pulse phase
                 before truncated to a value in range [0, 1).
      */
      virtual double calcPulsePhase(const timeSystem::AbsoluteTime & ev_time, double phase_offset = 0.) const = 0;

      /** \brief Compute the pulse frequency at a given time in the time system given to this object
                 upon its construction. Call getSytem method to obtain the time system to interpret
                 the return value of this method. The unit of frequency and its derivatives must be
                 based on seconds, i.e., Hz (sE-1), Hz/s (sE-2), and so on.
                 Note: validity of the ephemeris (valid since and valid until) are not checked. 
          \param ev_time Time of the event.
          \param derivative_order Order of time derivative of frequency to compute. Set 0 (zero) to obtain
                 frequency, 1 (one) the first time derivative of frequency, and so on.
      */
      virtual double calcFrequency(const timeSystem::AbsoluteTime & ev_time, int derivative_order = 0) const = 0;

      /** \brief Compute a sky position at a given time, and return it. The returned value is a pair of Right Ascension
                 and Declination of the position.
                 Note: validity of the ephemeris (valid since and valid until) are not checked.
          \param ev_time Absolute time at which a sky position is to be computed.
      */
      virtual std::pair<double, double> calcSkyPosition(const timeSystem::AbsoluteTime & ev_time) const = 0;

      /// \brief Create a copy of this object, and return a pointer to it.
      virtual PulsarEph * clone() const = 0;

      /// \brief Output text expression of this PulsarEph to a given output stream.
      virtual st_stream::OStream & write(st_stream::OStream & os) const;

    protected:
      /// \brief Construct a PulsarEph object.
      PulsarEph() {};

      /// \brief Output text expression of subclass-specific parameters of this object to a given output stream.
      virtual void writeModelParameter(st_stream::OStream & os) const = 0;
  };

  inline st_stream::OStream & operator <<(st_stream::OStream & os, const PulsarEph & eph) { return eph.write(os); }

  typedef std::vector<PulsarEph *> PulsarEphCont;
}

#endif
