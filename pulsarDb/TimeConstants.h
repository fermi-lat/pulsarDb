/** \file TimeConstants
    \brief Declaration of globally accessible static inline functions to assist in various conversions.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_TimeConstants_h
#define pulsarDb_TimeConstants_h

namespace pulsarDb {

  // Declarations.

  // Unit conversion constants.
  long double DayPerSec();

  long double SecPerDay();

  // Time system conversion constants.
  long double TaiMinusTtDay();

  long double TtMinusTaiDay();

  // Definitions.

  // Unit conversion constants.
  inline long double DayPerSec() { static long double r = 1.L / SecPerDay(); return r; }

  inline long double SecPerDay() { static long double r = 86400.L; return r; }

  // Time system conversion constants.
  inline long double TaiMinusTtDay() { static long double r = -TtMinusTaiDay(); return r; }

  inline long double TtMinusTaiDay() { static long double r = 32.184L * DayPerSec(); return r;}

}

#endif
