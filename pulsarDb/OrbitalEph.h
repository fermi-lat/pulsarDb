/** \file OrbitalEph.h
    \brief Interface for OrbitalEph class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_OrbitalEph_h
#define pulsarDb_OrbitalEph_h

#include <vector>

#include "tip/Table.h"

#include "pulsarDb/FormattedEph.h"

#include "timeSystem/ElapsedTime.h"

namespace st_stream {
  class OStream;
}

namespace timeSystem {
  class AbsoluteTime;
}

namespace tip {
  class Header;
}

namespace pulsarDb {

  class OrbitalEph: public FormattedEph {
    public:
      /// \brief Destruct this OrbitalEph object.
      virtual ~OrbitalEph() {}

      /// \brief Return a time system in which binary demodulation is performed.
      virtual const timeSystem::TimeSystem & getSystem() const = 0;

      /** \brief Compute an orbital delay for a photon emission time, add it to the time, and set the result to the time.
          \param ev_time Time at which a photo leaves a star in a binary system. This argument is to be updated with
                 the result of the computation.
      */
      void modulateBinary(timeSystem::AbsoluteTime & ev_time) const;

      /** \brief Estimate a photon emission time from a photon arrival time, and set the result to the argument as a return value.
          \param ev_time Time at which a photon arrives at an observer outside of a binary system. This argument is to be updated
                 with the result of the computation.
      */
      void demodulateBinary(timeSystem::AbsoluteTime & ev_time) const;

      /// \brief Return the T0 parameter value, which is the barycentric time of periastron.
      virtual const timeSystem::AbsoluteTime & t0() const = 0;

      /** \brief Compute an orbital phase of a given time, and return it.
          \param ev_time Absolute time for which an orbital phase is to be computed.
          \param phase_offset Value to be added to an orbital phase. This value is added to the computed orbital phase
                 before truncated to a value in range [0, 1).
      */
      virtual double calcOrbitalPhase(const timeSystem::AbsoluteTime & ev_time, double phase_offset = 0.) const = 0;

      /** \brief Compute a propagation delay in a binary system, in reference to the center of gravity of the system.
                 Note that the propagation delay includes not only a light-travel time, but also gravitational delays.
          \param ev_time Absolute time for which a propagation delay is to be computed.
      */
      virtual timeSystem::ElapsedTime calcOrbitalDelay(const timeSystem::AbsoluteTime & ev_time) const = 0;

      /// \brief Create a copy of this object.
      virtual OrbitalEph * clone() const = 0;

      /// \brief Output text expression of this PulsarEph to a given output stream.
      virtual st_stream::OStream & write(st_stream::OStream & os) const;

    protected:
      /** \brief Construct an OrbitalEph object.
          \param tolerance Maximum time difference that must be reached for a convergence in binary demodulation.
          \param max_iteration Maximum number of iterations to be tried to reach a convergence in binary demodulation.
      */
      OrbitalEph(const timeSystem::ElapsedTime & tolerance, int max_iteration);

      /// \brief Output text expression of subclass-specific parameters of this PulsarEph to a given output stream.
      virtual void writeModelParameter(st_stream::OStream & os) const = 0;

    private:
      timeSystem::ElapsedTime m_tolerance;
      int m_max_iteration;
  };

  typedef std::vector<OrbitalEph *> OrbitalEphCont;

  inline st_stream::OStream & operator <<(st_stream::OStream & os, const OrbitalEph & eph) { return eph.write(os); }
}

#endif
