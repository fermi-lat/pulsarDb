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

  EphComputer::EphComputer(): m_pulsar_eph_cont(), m_orbital_eph_cont(), m_pdot_pars(0), m_model(new TimingModel),
    m_chooser(new StrictEphChooser) {
  }

  EphComputer::EphComputer(const TimingModel & model, const EphChooser & chooser): m_pulsar_eph_cont(), m_orbital_eph_cont(),
    m_pdot_pars(0), m_model(model.clone()), m_chooser(chooser.clone()) {
  }

  EphComputer::~EphComputer() {
    delete m_pdot_pars;
    delete m_chooser;
    delete m_model;
  }

  void EphComputer::load(const PulsarDb & database) {
    loadPulsarEph(database);
    loadOrbitalEph(database);
  }

  void EphComputer::loadPulsarEph(const PulsarDb & database) {
    database.getEph(m_pulsar_eph_cont);
  }

  void EphComputer::loadOrbitalEph(const PulsarDb & database) {
    database.getEph(m_orbital_eph_cont);
  }

  void EphComputer::setPdotCancelParameter(const PulsarEph & pdot_pars) {
    m_pdot_pars = pdot_pars.clone();
  }

  FrequencyEph EphComputer::calcPulsarEph(const timeSystem::AbsoluteTime & ev_time) const {
    const PulsarEph & eph(m_chooser->choose(m_pulsar_eph_cont, ev_time));
    return m_model->calcEphemeris(eph, ev_time);
  }

  void EphComputer::cancelPdot(timeSystem::AbsoluteTime & ev_time) const {
    if (m_pdot_pars) {
      m_model->cancelPdot(*m_pdot_pars, ev_time);
    } else {
      throw std::runtime_error("Parameters for pdot cancellation are not set");
    }
  }

  double EphComputer::calcPulsePhase(const timeSystem::AbsoluteTime & ev_time, double phase_offset) const {
    const PulsarEph & eph(m_chooser->choose(m_pulsar_eph_cont, ev_time));
    return m_model->calcPulsePhase(eph, ev_time, phase_offset);
  }

  double EphComputer::calcOrbitalPhase(const timeSystem::AbsoluteTime & ev_time, double phase_offset) const {
    const OrbitalEph & eph(m_chooser->choose(m_orbital_eph_cont, ev_time));
    return m_model->calcOrbitalPhase(eph, ev_time, phase_offset);
  }

  void EphComputer::modulateBinary(timeSystem::AbsoluteTime & emission_time) const {
    const OrbitalEph & eph(m_chooser->choose(m_orbital_eph_cont, emission_time));
    m_model->modulateBinary(eph, emission_time);
  }

  void EphComputer::demodulateBinary(timeSystem::AbsoluteTime & arrival_time) const {
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

  const PulsarEph & EphComputer::choosePulsarEph(const timeSystem::AbsoluteTime & ev_time) const {
    return m_chooser->choose(m_pulsar_eph_cont, ev_time);
  }

  const OrbitalEph & EphComputer::chooseOrbitalEph(const timeSystem::AbsoluteTime & ev_time) const {
    return m_chooser->choose(m_orbital_eph_cont, ev_time);
  }

}
