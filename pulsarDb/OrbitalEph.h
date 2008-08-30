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
      virtual ~OrbitalEph() {}

      void modulateBinary(timeSystem::AbsoluteTime & ev_time) const;

      void demodulateBinary(timeSystem::AbsoluteTime & ev_time) const;

      virtual const timeSystem::AbsoluteTime & t0() const = 0;

      virtual double calcOrbitalPhase(const timeSystem::AbsoluteTime & ev_time, double phase_offset = 0.) const = 0;

      virtual timeSystem::ElapsedTime calcOrbitalDelay(const timeSystem::AbsoluteTime & ev_time) const = 0;

      /// \brief Create a copy of this object.
      virtual OrbitalEph * clone() const = 0;

      /// \brief Output text expression of this PulsarEph to a given output stream.
      virtual st_stream::OStream & write(st_stream::OStream & os) const;

    protected:
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
