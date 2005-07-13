/** \file EphComputer.cxx
    \brief Implementation for EphComputer class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#include "pulsarDb/EphChooser.h"
#include "pulsarDb/EphComputer.h"
#include "pulsarDb/PulsarDb.h"
#include "pulsarDb/TimingModel.h"

namespace pulsarDb {

  EphComputer::EphComputer(): m_pulsar_eph_cont(), m_orbital_eph_cont(), m_model(new TimingModel), m_chooser(new EphChooser) {
  }

  EphComputer::EphComputer(const TimingModel & model, const EphChooser & chooser): m_pulsar_eph_cont(), m_orbital_eph_cont(),
    m_model(model.clone()), m_chooser(chooser.clone()) {
  }

  EphComputer::~EphComputer() {
    delete m_chooser;
    delete m_model;
  }

  void EphComputer::load(const PulsarDb & pulsar_db) {
    pulsar_db.getEph(m_pulsar_eph_cont);
    pulsar_db.getEph(m_orbital_eph_cont);
  }

  FrequencyEph EphComputer::calcPulsarEph(const AbsoluteTime & ev_time) const {
    const PulsarEph & eph(m_chooser->choose(m_pulsar_eph_cont, ev_time));
    return m_model->calcEphemeris(eph, ev_time);
  }

  void EphComputer::cancelPdot(AbsoluteTime & ev_time) const {
    const PulsarEph & eph(m_chooser->choose(m_pulsar_eph_cont, ev_time));
    m_model->cancelPdot(eph, ev_time);
  }

  double EphComputer::calcPulsePhase(const AbsoluteTime & ev_time) const {
    const PulsarEph & eph(m_chooser->choose(m_pulsar_eph_cont, ev_time));
    return m_model->calcPulsePhase(eph, ev_time);
  }

  void EphComputer::modulateBinary(AbsoluteTime & emission_time) const {
    const OrbitalEph & eph(m_chooser->choose(m_orbital_eph_cont, emission_time));
    m_model->modulateBinary(eph, emission_time);
  }

  void EphComputer::demodulateBinary(AbsoluteTime & arrival_time) const {
    const OrbitalEph & eph(m_chooser->choose(m_orbital_eph_cont, arrival_time));
    m_model->demodulateBinary(eph, arrival_time);
  }

  PulsarEphCont & EphComputer::getPulsarEphCont() {
    return m_pulsar_eph_cont;
  }

  const PulsarEphCont & EphComputer::getPulsarEphCont() const {
    return m_pulsar_eph_cont;
  }

  OrbitalEphCont & EphComputer::getOrbitalEphCont() {
    return m_orbital_eph_cont;
  }

  const OrbitalEphCont & EphComputer::getOrbitalEphCont() const {
    return m_orbital_eph_cont;
  }

}
