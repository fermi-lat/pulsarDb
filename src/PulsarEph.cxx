/** \file PulsarEph.cxx
    \brief Implementation of the PulsarEph class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#include "pulsarDb/PulsarEph.h"

namespace pulsarDb {

  PulsarEph::PulsarEph(double valid_since, double valid_until, long epoch_int, double epoch_frac, long toa_geo_int,
    double toa_geo_frac, double phi0, double f0, double f1, double f2):
    m_epoch((long double)(epoch_int) + epoch_frac),
    m_toa((long double)(toa_geo_int) + toa_geo_frac),
    m_since(valid_since), m_until(valid_until), m_phi0(phi0), m_f0(f0), m_f1(f1), m_f2(f2) {}

  PulsarEph::PulsarEph(double valid_since, double valid_until, long double epoch, long double toa_geo, double phi0,
    double f0, double f1, double f2): m_epoch(epoch), m_toa(toa_geo), m_since(valid_since), m_until(valid_until),
    m_phi0(phi0), m_f0(f0), m_f1(f1), m_f2(f2) {}

  PulsarEph::~PulsarEph() {}

}
