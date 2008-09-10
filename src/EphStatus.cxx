/** \file EphStatus.cxx
    \brief Implementation of the EphStatus class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#include "pulsarDb/EphStatus.h"

namespace pulsarDb {

  EphStatus::EphStatus(const timeSystem::AbsoluteTime & effective_since, const timeSystem::AbsoluteTime & effective_until,
    const EphStatusCodeType & status_code, const std::string & description): m_since(effective_since), m_until(effective_until),
    m_code(status_code), m_description(description) {}

  const timeSystem::AbsoluteTime & EphStatus::getEffectiveSince() const {
    return m_since;
  }

  const timeSystem::AbsoluteTime & EphStatus::getEffectiveUntil() const {
    return m_until;
  }

  const EphStatusCodeType & EphStatus::getStatusCode() const {
    return m_code;
  }

  const std::string & EphStatus::getDescription() const {
    return m_description;
  }

  bool EphStatus::effectiveBetween(const timeSystem::AbsoluteTime & at1, const timeSystem::AbsoluteTime & at2) const {
    // Return false if this status never becomes effective.
    if (m_since > m_until) return false;

    // Time-order start_time and stop_time arguments, w/o using std::min/max.
    bool in_order = (at1 <= at2);
    const timeSystem::AbsoluteTime & start_time = (in_order ? at1 : at2);
    const timeSystem::AbsoluteTime & stop_time  = (in_order ? at2 : at1);

    // Return 
    return (start_time <= m_until && stop_time >= m_since);
  }

}
