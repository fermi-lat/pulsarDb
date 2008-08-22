/** \file EphComputer.cxx
    \brief Implementation for EphComputer class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#include "pulsarDb/EphChooser.h"
#include "pulsarDb/EphComputer.h"
#include "pulsarDb/PdotCanceler.h"
#include "pulsarDb/PulsarDb.h"

namespace pulsarDb {

  EphComputer::EphComputer(): m_pulsar_eph_cont(), m_orbital_eph_cont(), m_pdot_canceler(0), m_chooser(new StrictEphChooser) {}

  EphComputer::EphComputer(const EphChooser & chooser): m_pulsar_eph_cont(), m_orbital_eph_cont(),
    m_pdot_canceler(0), m_chooser(chooser.clone()) {
  }

  EphComputer::~EphComputer() {
    delete m_chooser;
    delete m_pdot_canceler;
    for (OrbitalEphCont::reverse_iterator itor = m_orbital_eph_cont.rbegin(); itor != m_orbital_eph_cont.rend(); ++itor) delete *itor;
    for (PulsarEphCont::reverse_iterator itor = m_pulsar_eph_cont.rbegin(); itor != m_pulsar_eph_cont.rend(); ++itor) delete *itor;
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

  void EphComputer::setPdotCancelParameter(const std::string & time_system_name, const timeSystem::AbsoluteTime & time_origin,
    const std::vector<double> & fdot_ratio) {
    delete m_pdot_canceler;
    m_pdot_canceler = new PdotCanceler(time_system_name, time_origin, fdot_ratio);
  }

  void EphComputer::setPdotCancelParameter(const timeSystem::AbsoluteTime & time_origin, const PulsarEph & pulsar_eph,
    int max_derivative) {
    delete m_pdot_canceler;
    m_pdot_canceler = new PdotCanceler(time_origin, pulsar_eph, max_derivative);
  }

  void EphComputer::setPdotCancelParameter(const timeSystem::AbsoluteTime & time_origin, int max_derivative) {
    const PulsarEph & eph(m_chooser->choose(m_pulsar_eph_cont, time_origin));
    delete m_pdot_canceler;
    m_pdot_canceler = new PdotCanceler(time_origin, eph, max_derivative);
  }

  void EphComputer::cancelPdot(timeSystem::AbsoluteTime & ev_time) const {
    if (m_pdot_canceler) {
      m_pdot_canceler->cancelPdot(ev_time);
    } else {
      throw std::runtime_error("Parameters for pdot cancellation are not set");
    }
  }

  double EphComputer::calcPulsePhase(const timeSystem::AbsoluteTime & ev_time, double phase_offset) const {
    const PulsarEph & eph(m_chooser->choose(m_pulsar_eph_cont, ev_time));
    return eph.calcPulsePhase(ev_time, phase_offset);
  }

  std::pair<double, double> EphComputer::calcSkyPosition(const timeSystem::AbsoluteTime & ev_time) const {
    const PulsarEph & eph(m_chooser->choose(m_pulsar_eph_cont, ev_time));
    return eph.calcSkyPosition(ev_time);
  }

  double EphComputer::calcOrbitalPhase(const timeSystem::AbsoluteTime & ev_time, double phase_offset) const {
    const OrbitalEph & eph(m_chooser->choose(m_orbital_eph_cont, ev_time));
    return eph.calcOrbitalPhase(ev_time, phase_offset);
  }

  void EphComputer::modulateBinary(timeSystem::AbsoluteTime & emission_time) const {
    const OrbitalEph & eph(m_chooser->choose(m_orbital_eph_cont, emission_time));
    eph.modulateBinary(emission_time);
  }

  void EphComputer::demodulateBinary(timeSystem::AbsoluteTime & arrival_time) const {
    const OrbitalEph & eph(m_chooser->choose(m_orbital_eph_cont, arrival_time));
    eph.demodulateBinary(arrival_time);
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
