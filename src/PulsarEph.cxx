/** \file PulsarEph.cxx
    \brief Implementation of the PulsarEph class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#include <iomanip>

#include "pulsarDb/PulsarEph.h"
#include "pulsarDb/TimingModel.h"

namespace pulsarDb {

  st_stream::OStream & operator <<(st_stream::OStream & os, const PulsarEph & eph) {
    std::ios::fmtflags orig_flags = os.flags();
    int orig_prec = os.precision(15);
    os << std::right;
    os.prefix().width(14); os << "Valid range = " << "[" << eph.valid_since() << ", " << eph.valid_until() << "]" << std::endl;
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
