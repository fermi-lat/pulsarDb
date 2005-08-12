/** \file EphComputer.h
    \brief Interface for EphComputer class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_EphComputer_h
#define pulsarDb_EphComputer_h

#include "pulsarDb/OrbitalEph.h"
#include "pulsarDb/PulsarEph.h"

namespace pulsarDb {

  class AbsoluteTime;
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

      FrequencyEph calcPulsarEph(const AbsoluteTime & ev_time) const;

      void cancelPdot(AbsoluteTime & ev_time) const;

      double calcPulsePhase(const AbsoluteTime & ev_time) const;

      double calcOrbitalPhase(const AbsoluteTime & ev_time) const;

      void modulateBinary(AbsoluteTime & emission_time) const;

      void demodulateBinary(AbsoluteTime & arrival_time) const;

      PulsarEphCont & getPulsarEphCont();

      const PulsarEphCont & getPulsarEphCont() const;

      OrbitalEphCont & getOrbitalEphCont();

      const OrbitalEphCont & getOrbitalEphCont() const;

      const PulsarEph & choosePulsarEph(const AbsoluteTime & ev_time) const;

      const OrbitalEph & chooseOrbitalEph(const AbsoluteTime & ev_time) const;

    private:
      PulsarEphCont m_pulsar_eph_cont;
      OrbitalEphCont m_orbital_eph_cont;
      TimingModel * m_model;
      EphChooser * m_chooser;
  };

}

#endif
