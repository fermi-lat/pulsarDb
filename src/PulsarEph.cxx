/** \file PulsarEph.cxx
    \brief Implementation of the PulsarEph class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#include "pulsarDb/PulsarEph.h"

namespace pulsarDb {

  PulsarEph::PulsarEph(double valid_since, double valid_until, long epoch_int, double epoch_frac, long t0geo_int, double t0geo_frac,
    double f0, double f1, double f2):
    m_epoch((long double)(epoch_int) + epoch_frac),
    m_t0((long double)(t0geo_int) + t0geo_frac),
    m_since(valid_since), m_until(valid_until), m_f0(f0), m_f1(f1), m_f2(f2) {}

  PulsarEph::PulsarEph(double valid_since, double valid_until, long double epoch, long double t0geo, double f0, double f1,
    double f2): m_epoch(epoch), m_t0(t0geo), m_since(valid_since), m_until(valid_until), m_f0(f0), m_f1(f1), m_f2(f2) {}

  PulsarEph::~PulsarEph() {}

}
