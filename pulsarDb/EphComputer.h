/** \file EphComputer.h
    \brief Interface for EphComputer class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_EphComputer_h
#define pulsarDb_EphComputer_h

#include "pulsarDb/OrbitalEph.h"
#include "pulsarDb/PulsarEph.h"

namespace timeSystem {
  class AbsoluteTime;
}

namespace pulsarDb {
  class EphChooser;
  class PulsarDb;
  class TimingModel;

  /** \class EphComputer
      \brief Abstraction providing high-level access to all main functions of pulsarDb package.
  */
  class EphComputer {
    public:
      EphComputer();

      EphComputer(const TimingModel & model, const EphChooser & chooser);

      ~EphComputer();

      void load(const PulsarDb & database);

      void loadPulsarEph(const PulsarDb & database);

      void loadOrbitalEph(const PulsarDb & database);

      void setPdotCancelParameter(const PulsarEph & pdot_pars);

      FrequencyEph calcPulsarEph(const timeSystem::AbsoluteTime & ev_time) const;

      void cancelPdot(timeSystem::AbsoluteTime & ev_time) const;

      double calcPulsePhase(const timeSystem::AbsoluteTime & ev_time) const;

      double calcOrbitalPhase(const timeSystem::AbsoluteTime & ev_time) const;

      void modulateBinary(timeSystem::AbsoluteTime & emission_time) const;

      void demodulateBinary(timeSystem::AbsoluteTime & arrival_time) const;

      PulsarEphCont & getPulsarEphCont();

      const PulsarEphCont & getPulsarEphCont() const;

      OrbitalEphCont & getOrbitalEphCont();

      const OrbitalEphCont & getOrbitalEphCont() const;

      const PulsarEph & choosePulsarEph(const timeSystem::AbsoluteTime & ev_time) const;

      const OrbitalEph & chooseOrbitalEph(const timeSystem::AbsoluteTime & ev_time) const;

    private:
      EphComputer(const EphComputer &);
      EphComputer & operator =(const EphComputer &);

      PulsarEphCont m_pulsar_eph_cont;
      OrbitalEphCont m_orbital_eph_cont;
      PulsarEph * m_pdot_pars;
      TimingModel * m_model;
      EphChooser * m_chooser;
  };

}

#endif
