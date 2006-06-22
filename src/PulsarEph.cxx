/** \file PulsarEph.cxx
    \brief Implementation of the PulsarEph class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#include <iomanip>

#include "pulsarDb/PulsarEph.h"
#include "pulsarDb/TimingModel.h"

#include "timeSystem/Duration.h"
#include "timeSystem/ElapsedTime.h"
#include "timeSystem/IntFracPair.h"
#include "timeSystem/TimeInterval.h"

using namespace timeSystem;

namespace pulsarDb {

//  double PulsarEph::dt(const AbsoluteTime & at) const {
//    IntFracPair numerator = (at - m_epoch).computeElapsedTime(m_system->getName()).getTime().getValue(Sec);
//    return (numerator.getIntegerPart() + numerator.getFractionalPart()) / m_unit_time;
//  }

  double PulsarEph::dt(const AbsoluteTime & at1, const AbsoluteTime & at2) const {
    IntFracPair numerator = (at1 - at2).computeElapsedTime(m_system->getName()).getTime().getValue(Sec);
    return (numerator.getIntegerPart() + numerator.getFractionalPart()) / m_unit_time;
  }

  st_stream::OStream & operator <<(st_stream::OStream & os, const PulsarEph & eph) {
    std::ios::fmtflags orig_flags = os.flags();
    int orig_prec = os.precision(15);
    os << std::right;
    os.prefix().width(14);
    // Note: below, break into two consecutive strings so that width applies to first part only.
    if (eph.valid_since() < eph.valid_until())
      os << "Validity : " << "in range " << "[" << eph.valid_since() << ", " << eph.valid_until() << ")" << std::endl;
    else
      os << "Validity : " << "only at epoch" << std::endl;
    os.prefix().width(14); os << "Epoch = " << eph.epoch() << std::endl;
    os.prefix().width(14); os << "Phi0 = " << eph.phi0() << std::endl;
    os.prefix().width(14); os << "F0 = " << eph.f0() << std::endl;
    os.prefix().width(14); os << "F1 = " << eph.f1() << std::endl;
    os.prefix().width(14); os << "F2 = " << eph.f2();
    os.flags(orig_flags);
    os.precision(orig_prec);
    return os;
  }

}
