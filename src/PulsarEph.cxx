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
#include "timeSystem/TimeRep.h"

using namespace timeSystem;

namespace pulsarDb {

  double PulsarEph::dt(const AbsoluteTime & at1, const AbsoluteTime & at2) const {
    IntFracPair numerator = (at1 - at2).computeElapsedTime(m_system->getName()).getTime().getValue(Sec);
    return numerator.getDouble() / m_unit_time;
  }

  st_stream::OStream & operator <<(st_stream::OStream & os, const PulsarEph & eph) {
    std::ios::fmtflags orig_flags = os.flags();
    int orig_prec = os.precision(15);
    os << std::right;
    os.prefix().width(14);
    std::string time_system_name = eph.getSystem().getName();
    MjdRep mjd_rep(time_system_name, 0, 0.);

    // Note: below, break into two consecutive strings so that width applies to first part only.
    if (!eph.valid_since().equivalentTo(eph.valid_until(), ElapsedTime(time_system_name, Duration(0, 1.e-9)))) {
      eph.valid_since().getTime(mjd_rep);
      os << "Validity : " << "in range " << "[" << mjd_rep << ", ";
      eph.valid_until().getTime(mjd_rep);
      os << mjd_rep << ")" << std::endl;
    } else {
      eph.valid_since().getTime(mjd_rep);
      os << "Validity : " << "only at time " << mjd_rep << std::endl;
    }
    eph.epoch().getTime(mjd_rep);
    os.prefix().width(14); os << "Epoch = " << mjd_rep << std::endl;
    os.prefix().width(14); os << "Phi0 = " << eph.phi0() << std::endl;
    os.prefix().width(14); os << "F0 = " << eph.f0() << std::endl;
    os.prefix().width(14); os << "F1 = " << eph.f1() << std::endl;
    os.prefix().width(14); os << "F2 = " << eph.f2();
    os.flags(orig_flags);
    os.precision(orig_prec);
    return os;
  }

}
