/** \file PulsarEph.h
    \brief Interface for PulsarEph class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_PulsarEph_h
#define pulsarDb_PulsarEph_h

#include <vector>

namespace pulsarDb {

  class TimingModel;

  /** \class PulsarEph
      \brief Class representing a single pulsar ephemeris.
  */
  class PulsarEph {
    public:
      PulsarEph(long double valid_since, long double valid_until, long double epoch): m_since(valid_since), m_until(valid_until),
        m_epoch(epoch) {}

      virtual ~PulsarEph() {}

      virtual long double valid_since() const { return m_since; }
      virtual long double valid_until() const { return m_until; }
      virtual long double epoch() const { return m_epoch; }
      virtual double phi0() const = 0;
      virtual double f0() const = 0;
      virtual double f1() const = 0;
      virtual double f2() const = 0;
      virtual PulsarEph * clone() const = 0;

    protected:
      long double m_since;
      long double m_until;
      long double m_epoch;
  };

  /** \class DatabaseEph
      \brief Class representing a single pulsar ephemeris.
  */
  class DatabaseEph : public PulsarEph {
    public:
      /** \brief Create pulsar ephemeris with the given properties.
          \param valid_since The start time of the valid interval.
          \param valid_until The stop time of the valid interval.
          \param epoch_int The integer portion of the epoch (time origin).
          \param epoch_frac The fractional portion of the epoch (time origin).
          \param toa_int The integer portion of the pulse arrival time.
          \param toa_frac The fractional portion of the pulse arrival time.
          \param f0 The frequency at the epoch (time origin).
          \param f1 The first time derivative of the frequency at the epoch (time origin).
          \param f2 The second time derivative of the frequency at the epoch (time origin).
      */
      DatabaseEph(const TimingModel & model, long double valid_since, long double valid_until, long epoch_int, double epoch_frac,
        long toa_int, double toa_frac, double f0, double f1, double f2);

      /** \brief Create pulsar ephemeris with the given properties.
          \param valid_since The start time of the valid interval.
          \param valid_until The stop time of the valid interval.
          \param epoch The epoch (time origin).
          \param toa_int The pulse arrival time.
          \param f0 The frequency at the epoch (time origin).
          \param f1 The first time derivative of the frequency at the epoch (time origin).
          \param f2 The second time derivative of the frequency at the epoch (time origin).
      */
      DatabaseEph(const TimingModel & model, long double valid_since, long double valid_until, long double epoch, long double toa,
        double f0, double f1, double f2);

      virtual ~DatabaseEph();

      virtual double phi0() const;
      virtual double f0() const { return m_f0; }
      virtual double f1() const { return m_f1; }
      virtual double f2() const { return m_f2; }

      long double toa() const { return m_toa; }
      virtual PulsarEph * clone() const { return new DatabaseEph(*this); }

    private:
      const TimingModel * m_model;
      long double m_toa;
      double m_f0;
      double m_f1;
      double m_f2;
  };

  /** \class FrequencyEph
      \brief Class representing a single pulsar ephemeris.
  */
  class FrequencyEph : public PulsarEph {
    public:
      /** \brief Create pulsar ephemeris with the given properties.
          \param epoch The epoch (time origin).
          \param f0 The frequency at the epoch (time origin).
          \param f1 The first time derivative of the frequency at the epoch (time origin).
          \param f2 The second time derivative of the frequency at the epoch (time origin).
      */
      FrequencyEph(long double valid_since, long double valid_until, long double epoch, double phi0, double f0,
        double f1, double f2): PulsarEph(valid_since, valid_until, epoch), m_phi0(phi0), m_f0(f0), m_f1(f1), m_f2(f2) {}

      virtual ~FrequencyEph() {}

      virtual double phi0() const { return m_phi0; }
      virtual double f0() const { return m_f0; }
      virtual double f1() const { return m_f1; }
      virtual double f2() const { return m_f2; }
      virtual PulsarEph * clone() const { return new FrequencyEph(*this); }

    private:
      double m_phi0;
      double m_f0;
      double m_f1;
      double m_f2;
  };

  /** \class PeriodEph
      \brief Class representing a single pulsar ephemeris.
  */
  class PeriodEph : public PulsarEph {
    public:
      /** \brief Create pulsar ephemeris with the given properties.
          \param epoch The epoch (time origin).
          \param p0 The period at the epoch (time origin).
          \param p1 The first time derivative of the period at the epoch (time origin).
          \param p2 The second time derivative of the period at the epoch (time origin).
      */
      PeriodEph(long double valid_since, long double valid_until, long double epoch, double phi0, double p0, double p1, double p2):
        PulsarEph(valid_since, valid_until, epoch), m_phi0(phi0), m_p0(p0), m_p1(p1), m_p2(p2) {}

      virtual ~PeriodEph() {}

      virtual double phi0() const { return m_phi0; }
      virtual double f0() const { return 1. / m_p0; }
      virtual double f1() const { return - m_p1 / (m_p0 * m_p0); }
      virtual double f2() const { double p0sq = m_p0 * m_p0; return 2. * m_p1 * m_p1 / (m_p0 * p0sq) - m_p2 / p0sq; }
      virtual PulsarEph * clone() const { return new PeriodEph(*this); }

    private:
      double m_phi0;
      double m_p0;
      double m_p1;
      double m_p2;
  };

}

#endif
