/** \file PulsarEph.cxx
    \brief Implementation of the PulsarEph class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#include "pulsarDb/PulsarEph.h"
#include "pulsarDb/TimingModel.h"

namespace pulsarDb {

  DatabaseEph::DatabaseEph(const TimingModel & model, long double valid_since, long double valid_until, long double epoch,
    long double toa, double f0, double f1, double f2): PulsarEph(valid_since, valid_until, epoch), m_model(&model),
    m_toa(toa), m_f0(f0), m_f1(f1), m_f2(f2) {}

  DatabaseEph::~DatabaseEph() {}

  double DatabaseEph::phi0() const {
    // The phase is computed from the phase of the toa.
    return -m_model->calcPhase(*this, m_toa);
  }

}
