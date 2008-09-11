/** \file EphComputer.h
    \brief Interface for EphComputer class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_EphComputer_h
#define pulsarDb_EphComputer_h

#include <utility>

#include "pulsarDb/OrbitalEph.h"
#include "pulsarDb/PulsarEph.h"

namespace timeSystem {
  class AbsoluteTime;
}

namespace pulsarDb {
  class EphChooser;
  class PdotCanceler;
  class PulsarDb;

  /** \class EphComputer
      \brief Abstraction providing high-level access to all main functions of pulsarDb package.
  */
  class EphComputer {
    public:
      EphComputer();

      EphComputer(const EphChooser & chooser);

      ~EphComputer();

      void load(const PulsarDb & database);

      void loadPulsarEph(const PulsarDb & database);

      void loadPulsarEph(const PulsarEph & pulsar_eph);

      void loadOrbitalEph(const PulsarDb & database);

      void loadOrbitalEph(const OrbitalEph & orbital_eph);

      void loadEphRemark(const PulsarDb & database);

      void loadEphRemark(const EphStatus & eph_remark);

      void setPdotCancelParameter(const std::string & time_system_name, const timeSystem::AbsoluteTime & time_origin,
        const std::vector<double> & fdot_ratio);

      void setPdotCancelParameter(const timeSystem::AbsoluteTime & time_origin, const PulsarEph & pulsar_eph, int max_derivative);

      void setPdotCancelParameter(const timeSystem::AbsoluteTime & time_origin, int max_derivative);

      void cancelPdot(timeSystem::AbsoluteTime & ev_time) const;

      double calcPulsePhase(const timeSystem::AbsoluteTime & ev_time, double phase_offset = 0.) const;

      std::pair<double, double> calcSkyPosition(const timeSystem::AbsoluteTime & ev_time) const;

      double calcOrbitalPhase(const timeSystem::AbsoluteTime & ev_time, double phase_offset = 0.) const;

      void modulateBinary(timeSystem::AbsoluteTime & emission_time) const;

      void demodulateBinary(timeSystem::AbsoluteTime & arrival_time) const;

      PulsarEphCont::size_type getNumPulsarEph() const;

      OrbitalEphCont::size_type getNumOrbitalEph() const;

      EphStatusCont::size_type getNumEphRemark() const;

      const PulsarEph & choosePulsarEph(const timeSystem::AbsoluteTime & ev_time) const;

      const OrbitalEph & chooseOrbitalEph(const timeSystem::AbsoluteTime & ev_time) const;

      void examinePulsarEph(const timeSystem::AbsoluteTime & start_time, const timeSystem::AbsoluteTime & stop_time,
        EphStatusCont & eph_status_cont) const;

      void getEphRemark(const timeSystem::AbsoluteTime & start_time, const timeSystem::AbsoluteTime & stop_time,
        EphStatusCont & eph_status_cont) const;

    private:
      EphComputer(const EphComputer &);
      EphComputer & operator =(const EphComputer &);

      PulsarEphCont m_pulsar_eph_cont;
      OrbitalEphCont m_orbital_eph_cont;
      EphStatusCont m_eph_remark_cont;
      PdotCanceler * m_pdot_canceler;
      EphChooser * m_chooser;
  };

}

#endif
