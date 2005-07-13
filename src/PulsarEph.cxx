/** \file PulsarEph.cxx
    \brief Implementation of the PulsarEph class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#include "pulsarDb/PulsarEph.h"
#include "pulsarDb/TimingModel.h"

namespace pulsarDb {

  DatabaseEph::DatabaseEph(const AbsoluteTime & valid_since, const AbsoluteTime & valid_until,
    const AbsoluteTime & epoch, const AbsoluteTime & toa, double f0, double f1, double f2):
    PulsarEph(valid_since, valid_until, epoch), m_toa(toa.clone()), m_f0(f0), m_f1(f1), m_f2(f2) {}

  DatabaseEph::~DatabaseEph() { delete m_toa; }

  double DatabaseEph::phi0() const {
    TimingModel model;

    // Create temporary copy of this ephemeris with phi0 == 0.
    FrequencyEph tmp(*m_since, *m_until, *m_epoch, 0., m_f0, m_f1, m_f2);

    // Use the timing model and temporary ephemeris to compute the phase from the negative of the toa field.
    double r = - model.calcPulsePhase(tmp, *m_toa);

    // Make sure it is in the range [0, 1). calcPulsePhase is bounded in this way.
    if (0. > r) r += 1.;

    return r;
  }

}
