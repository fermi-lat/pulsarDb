/** \file EphComputer.cxx
    \brief Implementation for EphComputer class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#include <iostream>
#include <map>
#include <sstream>

#include "pulsarDb/EphChooser.h"
#include "pulsarDb/EphComputer.h"
#include "pulsarDb/PdotCanceler.h"
#include "pulsarDb/PulsarDb.h"

namespace pulsarDb {

  EphComputer::EphComputer(): m_pulsar_eph_cont(), m_orbital_eph_cont(), m_eph_remark_cont(), m_pdot_canceler(0),
    m_chooser(new StrictEphChooser) {}

  EphComputer::EphComputer(const EphChooser & chooser): m_pulsar_eph_cont(), m_orbital_eph_cont(), m_eph_remark_cont(),
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
    loadEphRemark(database);
  }

  void EphComputer::loadPulsarEph(const PulsarDb & database) {
    database.getEph(m_pulsar_eph_cont);
  }

  void EphComputer::loadOrbitalEph(const PulsarDb & database) {
    database.getEph(m_orbital_eph_cont);
  }

  void EphComputer::loadEphRemark(const PulsarDb & database) {
    database.getRemark(m_eph_remark_cont);
  }

  void EphComputer::loadPulsarEph(const PulsarEph & pulsar_eph) {
    m_pulsar_eph_cont.push_back(pulsar_eph.clone());
  }

  void EphComputer::loadOrbitalEph(const OrbitalEph & orbital_eph) {
    m_orbital_eph_cont.push_back(orbital_eph.clone());
  }

  void EphComputer::loadEphRemark(const EphStatus & eph_remark) {
    if (Remarked == eph_remark.getStatusCode()) m_eph_remark_cont.push_back(eph_remark);
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

  PulsarEphCont::size_type EphComputer::getNumPulsarEph() const {
    return m_pulsar_eph_cont.size();
  }

  OrbitalEphCont::size_type EphComputer::getNumOrbitalEph() const {
    return m_orbital_eph_cont.size();
  }

  EphStatusCont::size_type EphComputer::getNumEphRemark() const {
    return m_eph_remark_cont.size();
  }

  const PulsarEph & EphComputer::choosePulsarEph(const timeSystem::AbsoluteTime & ev_time) const {
    return m_chooser->choose(m_pulsar_eph_cont, ev_time);
  }

  const OrbitalEph & EphComputer::chooseOrbitalEph(const timeSystem::AbsoluteTime & ev_time) const {
    return m_chooser->choose(m_orbital_eph_cont, ev_time);
  }

  void EphComputer::examinePulsarEph(const timeSystem::AbsoluteTime & start_time, const timeSystem::AbsoluteTime & stop_time,
    EphStatusCont & eph_status_cont) const {
    m_chooser->examine(m_pulsar_eph_cont, start_time, stop_time, eph_status_cont);
  }

  void EphComputer::getEphRemark(const timeSystem::AbsoluteTime & start_time, const timeSystem::AbsoluteTime & stop_time,
    EphStatusCont & eph_status_cont) const {
    // Clear the output container.
    eph_status_cont.clear();

    // Subselect ephemeris remarks by the given time interval, and copy selected remarks to the output.
    for (EphStatusCont::const_iterator itor = m_eph_remark_cont.begin(); itor != m_eph_remark_cont.end(); ++itor) {
      if (itor->effectiveBetween(start_time, stop_time)) eph_status_cont.push_back(*itor);
    }
  }

  void EphComputer::summarizeTimeSystem(std::string & spin_summary, std::string & orbital_summary, std::string & pdot_summary) const {
    typedef std::map<std::string, long> map_type;

    // Collect time system names of spin ephemerides.
    map_type num_system_spin;
    for (PulsarEphCont::const_iterator spin_itor = m_pulsar_eph_cont.begin(); spin_itor != m_pulsar_eph_cont.end(); ++spin_itor) {
      std::string time_system_name = (*spin_itor)->getSystem().getName();
      if (num_system_spin.find(time_system_name) == num_system_spin.end()) num_system_spin[time_system_name] = 1;
      else num_system_spin[time_system_name] += 1;
    }

    // Summarize the time systems of spin ephemerides.
    spin_summary = "";
    for (map_type::const_iterator map_itor = num_system_spin.begin(); map_itor != num_system_spin.end(); ++map_itor) {
      std::ostringstream os_str;
      os_str << map_itor->first << "(" << map_itor->second << ") ";
      spin_summary += os_str.str();
    }
    if (spin_summary.empty()) spin_summary = "None";
    else spin_summary.erase(spin_summary.size() - 1);

    // Collect time system names of orbital ephemerides.
    map_type num_system_orbital;
    for (OrbitalEphCont::const_iterator orb_itor = m_orbital_eph_cont.begin(); orb_itor != m_orbital_eph_cont.end(); ++orb_itor) {
      std::string time_system_name = (*orb_itor)->getSystem().getName();
      if (num_system_orbital.find(time_system_name) == num_system_orbital.end()) num_system_orbital[time_system_name] = 1;
      else num_system_orbital[time_system_name] += 1;
    }

    // Summarize the time systems of orbital ephemerides.
    orbital_summary = "";
    for (map_type::const_iterator map_itor = num_system_orbital.begin(); map_itor != num_system_orbital.end(); ++map_itor) {
      std::ostringstream os_str;
      os_str << map_itor->first << "(" << map_itor->second << ") ";
      orbital_summary += os_str.str();
    }
    if (orbital_summary.empty()) orbital_summary = "None";
    else orbital_summary.erase(orbital_summary.size() - 1);

    // Summarize the time system for pdot cancellation.
    if (m_pdot_canceler) pdot_summary = m_pdot_canceler->getSystem().getName();
    else pdot_summary = "None";
  }

}
