/** \file HighPrecisionEph.cxx
    \brief Implementation of the HighPrecisionEph class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <utility>

#include "pulsarDb/PulsarEph.h"
#include "pulsarDb/HighPrecisionEph.h"

#include "timeSystem/CalendarFormat.h"
#include "timeSystem/MjdFormat.h"
#include "timeSystem/SourcePosition.h"
#include "timeSystem/TimeConstant.h"
#include "timeSystem/TimeInterval.h"

using namespace timeSystem;

namespace tip {
  class Header;
}

namespace pulsarDb {

  const double HighPrecisionEph::s_rad_per_deg  = M_PI / 180.;
  const double HighPrecisionEph::s_rad_per_mas  = HighPrecisionEph::s_rad_per_deg / 3600. / 1000.;
  const double HighPrecisionEph::s_rad_day_per_mas_year = HighPrecisionEph::s_rad_per_mas / 365.25;
  const double HighPrecisionEph::s_km_per_lts  = 2.99792458e+5;
  const double HighPrecisionEph::s_lts_per_au = 1.49597870691e+8 / HighPrecisionEph::s_km_per_lts; // From JPL DE405.

  HighPrecisionEph::HighPrecisionEph(const tip::Table::ConstRecord & record, const tip::Header & /* header */):
    m_system(&TimeSystem::getSystem("TDB")), m_since("TDB", 0, 0.), m_until("TDB", 0, 0.), m_pos_epoch("TDB", 0, 0.),
    m_ra(0.), m_dec(0.), m_ra_vel(0.), m_dec_vel(0.), m_radial_vel(0.), m_parallax(-1.), m_freq_epoch("TDB", 0, 0.),
    m_freq_pars(0), m_wave_omega(-1.), m_wave_sine(0), m_wave_cosine(0), m_glitch_list(0), m_remark_cont() {
    // Read the start time of validity window (required).
    long valid_since_date = 0;
    read(record, "VALID_SINCE", valid_since_date);
    m_since = AbsoluteTime("TDB", Mjd(valid_since_date, 0.));

    // Read the end time of validity window (required).
    // Note: One is added to the endpoint because the "VALID_UNTIL" field in the file expires at the end of that day,
    // whereas the valid_until argument to the ephemeris object is the absolute cutoff.
    long valid_until_date = 0;
    read(record, "VALID_UNTIL", valid_until_date);
    m_until = AbsoluteTime("TDB", Mjd(valid_until_date + 1, 0.));

    // Read the position epoch (integer parts required, fractional parts optional).
    long pos_epoch_int = 0;
    double pos_epoch_frac = 0.;
    read(record, "POS_EPOCH_INT", pos_epoch_int);
    read(record, "POS_EPOCH_FRAC", pos_epoch_frac, 0.);
    m_pos_epoch = AbsoluteTime("TDB", Mjd(pos_epoch_int, pos_epoch_frac));

    // Read the position parameters (RA, DEC: required, RA_VELOCITY, DEC_VELOCITY, RADIAL_VELOCITY, PARALLAX: optional).
    read(record, "RA",  m_ra);
    read(record, "DEC", m_dec);
    read(record, "RA_VELOCITY",  m_ra_vel, 0.);
    read(record, "DEC_VELOCITY", m_dec_vel, 0.);
    read(record, "RADIAL_VELOCITY", m_radial_vel, 0.);
    read(record, "PARALLAX", m_parallax, -1.); // To indicate parallax not known.

    // Read the frequency epoch (integer parts required, fractional parts optional).
    long freq_epoch_int = 0;
    double freq_epoch_frac = 0.;
    read(record, "FREQ_EPOCH_INT", freq_epoch_int);
    read(record, "FREQ_EPOCH_FRAC", freq_epoch_frac, 0.);
    m_freq_epoch = AbsoluteTime("TDB", Mjd(freq_epoch_int, freq_epoch_frac));

    // Read the pulse TOA (integer parts required, fractional parts optional).
    long toa_int = 0;
    double toa_frac = 0.;
    read(record, "TOABARY_INT", toa_int);
    read(record, "TOABARY_FRAC", toa_frac, 0.);
    AbsoluteTime toa("TDB", Mjd(toa_int, toa_frac));

    // Read the frequency parameters (one element required, others optional).
    std::vector<double> double_array;
    std::vector<bool> null_array;
    read(record, "FREQ_PARAMETERS", double_array, 0., null_array);
    // Note: Allow an empty array for frequency parameters (i.e., F0 = F1 = F2 = ... = 0.0). Technically it is possible
    //       that a positive frequency is resulted from other terms such as wave parameters and glitches.
    m_freq_pars.resize(double_array.size() + 1, 0.);
    for (std::size_t ii = 0; ii != double_array.size(); ++ii) m_freq_pars[ii+1] = double_array[ii];

    // Read the wave parameters (optional).
    bool undefined = read(record, "WAVE_OMEGA", m_wave_omega, 1.);
    if (!undefined) {
      read(record, "WAVE_SINE", m_wave_sine, 0., null_array);
      read(record, "WAVE_COSINE", m_wave_cosine, 0., null_array);
    }

    // Read the glitch parameters (optional).
    read(record, "GLITCH_PARAMETERS", double_array, 0., null_array);
    if (double_array.size() > 0) {
      // Read the dimensions (required in this scope).
      std::vector<int> int_array;
      read(record, "GLITCH_DIMENSIONS", int_array);

      // Check the dimensions and the total number of parameters.
      std::vector<double>::size_type num_elem = 0;
      for (std::vector<int>::const_iterator itor = int_array.begin(); itor != int_array.end(); ++itor) {
        if (*itor < 0) throw std::runtime_error("Negative value found in GLITCH_DIMENSIONS (all must be non-negative)");
        num_elem += *itor;
      }
      if (double_array.size() != num_elem) {
        std::ostringstream oss;
        oss << "The number of elements in GLITCH_PARAMETERS (" << double_array.size() <<
          ") does not match the total of GLITCH_DIMENSIONS elements (" << num_elem << ")";
        throw std::runtime_error(oss.str());
      }

      // Decode GLITCH_PARAMETERS values.
      std::vector<double>::size_type first_index = 0;
      std::vector<double>::size_type last_index = double_array.size();
      for (std::vector<int>::const_iterator num_itor = int_array.begin(); num_itor != int_array.end(); ++num_itor) {
        // Append a new GlitchParameter object.
        m_glitch_list.push_back(GlitchParameter());
        GlitchParameter & this_glitch(m_glitch_list.back());

        // Read a glitch epoch (required in this scope).
        if (1 != *num_itor) throw std::runtime_error("More than one (1) field is assigned to a glitch epoch");
        if (null_array[first_index]) throw std::runtime_error("Glitch epoch is not given for a glitch");
        this_glitch.m_epoch = AbsoluteTime("TDB", Mjd1(double_array[first_index]));
        ++first_index;

        // Read permanent jump parameters (optional).
        ++num_itor;
        if (num_itor != int_array.end()) {
          last_index = first_index + *num_itor;
          this_glitch.m_perm_jump.insert(this_glitch.m_perm_jump.end(), double_array.begin() + first_index,
            double_array.begin() + last_index);
          first_index = last_index;

          // Read parameters for decaying components (optional).
          ++num_itor;
          if (num_itor != int_array.end()) {
            last_index = first_index + *num_itor;
            for (; first_index != last_index; ++first_index) {
              // Read the amplitude (optional).
              double amplitude = double_array[first_index];

              // Read the decay time (required in this scope).
              ++first_index;
              if (first_index == double_array.size()) throw std::runtime_error("Decay time is not given for a decaying component");
              if (null_array[first_index]) throw std::runtime_error("Decay time for a decaying component is undefined");
              double decay_time = double_array[first_index];

              // Append a decaying component.
              this_glitch.m_decay_comp.push_back(std::make_pair(amplitude, decay_time));
            }
          }
        }
      }
    }

    // Create temporary copy of this ephemeris with m_freq_pars[0] == 0 (where m_freq_pars[0] stores PHI0).
    HighPrecisionEph tmp("TDB", m_since, m_until, m_pos_epoch, m_ra, m_dec, m_ra_vel, m_dec_vel, m_radial_vel, m_parallax,
      m_freq_epoch, m_freq_pars, m_wave_omega, m_wave_sine, m_wave_cosine, m_glitch_list);

    // Use the temporary ephemeris to compute the phase from the negative of the toa field.
    double phi0 = - tmp.calcPulsePhase(toa);

    // Make sure it is in the range [0, 1). calcPulsePhase is bounded in this way.
    if (0. > phi0) phi0 += 1.;

    // Put the computed phi0.
    m_freq_pars[0] = phi0;

    // Set ephemeris remarks to a data member.
    setRemark();
  }

  void HighPrecisionEph::writeModelParameter(st_stream::OStream & os) const {
    // Write the astrometric parameters.
    std::string pos_epoch_string;
    try {
      pos_epoch_string = m_pos_epoch.represent(m_system->getName(), MjdFmt);
    } catch (const std::exception &) {
      pos_epoch_string = m_pos_epoch.represent(m_system->getName(), CalendarFmt);
    } 
    os << format("Position Epoch",  pos_epoch_string, "") << std::endl;
    os << format("RA",              m_ra,             "degrees")            << std::endl;
    os << format("Dec",             m_dec,            "degrees")            << std::endl;
    os << format("RA Velocity",     m_ra_vel,         "milliarcseconds/yr") << std::endl;
    os << format("Dec Velocity",    m_dec_vel,        "milliarcseconds/yr") << std::endl;
    os << format("Radial Velocity", m_radial_vel,     "km/s")               << std::endl;
    if (m_parallax > 0.) os << format("Annual Parallax", m_parallax, "milliarcseconds") << std::endl;

    // Write the spin parameters.
    std::string freq_epoch_string;
    try {
      freq_epoch_string = m_freq_epoch.represent(m_system->getName(), MjdFmt);
    } catch (const std::exception &) {
      freq_epoch_string = m_freq_epoch.represent(m_system->getName(), CalendarFmt);
    } 
    os << format("Frequency Epoch", freq_epoch_string, "");
    int deriv_order = -1;
    for (freq_type::const_iterator freq_itor = m_freq_pars.begin(); freq_itor != m_freq_pars.end(); ++freq_itor, ++deriv_order) {
      // Create the parameter name.
      std::string par_name = (deriv_order < 0 ? "Phi" : "F");
      std::ostringstream oss_name;
      oss_name << std::max(deriv_order, 0);
      par_name += oss_name.str();

      // Create the physical unit.
      std::string par_unit = "";
      if (deriv_order >= 0) {
        std::ostringstream oss_unit;
        oss_unit << "s**(-" << deriv_order + 1 << ")";
        par_unit += oss_unit.str();
      }

      // Write the parameter.
      os << std::endl << format(par_name, *freq_itor, par_unit);
    }

    // Write the wave parameters.
    wave_type::size_type num_elem_sine = m_wave_sine.size();
    wave_type::size_type num_elem_cosine = m_wave_cosine.size();
    if (num_elem_sine > 0 || num_elem_cosine > 0) {
      os << std::endl << format("Wave Frequency", m_wave_omega, "radians/day");
      for (wave_type::size_type ii = 0; ii < std::max(num_elem_sine, num_elem_cosine); ++ii) {
        std::ostringstream oss;
        oss << ii + 1;
        if (ii < num_elem_sine) os << std::endl << format("Sin" + oss.str(), m_wave_sine[ii], "s");
        if (ii < num_elem_cosine) os << std::endl << format("Cos" + oss.str(), m_wave_cosine[ii], "s");
      }
    }

    // Write the glitch parameters.
    for (glitch_type::const_iterator glitch_itor = m_glitch_list.begin(); glitch_itor != m_glitch_list.end(); ++glitch_itor) {
      // Write the glitch epoch.
      std::string epoch_string;
      try {
        epoch_string = glitch_itor->m_epoch.represent(m_system->getName(), MjdFmt);
      } catch (const std::exception &) {
        epoch_string = glitch_itor->m_epoch.represent(m_system->getName(), CalendarFmt);
      }
      os << std::endl << format("Glitch Epoch", epoch_string, "");

      // Write the permanent jumps.
      int deriv_order = -1;
      const jump_type & perm_jump(glitch_itor->m_perm_jump);
      for (jump_type::const_iterator jump_itor = perm_jump.begin(); jump_itor != perm_jump.end(); ++jump_itor, ++deriv_order) {
        // Create the parameter name.
        std::string par_name = (deriv_order < 0 ? "dPhi" : "dF");
        std::ostringstream oss_name;
        oss_name << std::max(deriv_order, 0);
        par_name += oss_name.str();

        // Create the physical unit.
        std::string par_unit = "";
        if (deriv_order >= 0) {
          std::ostringstream oss_unit;
          oss_unit << "s**(-" << deriv_order + 1 << ")";
          par_unit += oss_unit.str();
        }

        // Write the parameter.
        os << std::endl << format(par_name, *jump_itor, par_unit);
      }

      // Write the decaying components.
      const decay_type & decay_comp(glitch_itor->m_decay_comp);
      for (decay_type::const_iterator decay_itor = decay_comp.begin(); decay_itor != decay_comp.end(); ++decay_itor) {
        os << std::endl << format("Amplitude", decay_itor->first, "s**(-1)");
        os << std::endl << format("Decay Time", decay_itor->second, "days");
      }
    }
  }

  double HighPrecisionEph::calcPulsePhase(const AbsoluteTime & ev_time, double phase_offset) const {
    // Prepare a return value.
    double phase = 0.;

    // Compute a phase value that comes from frequency parameters.
    if (m_freq_pars.size() > 0) {
      // Add the constant term.
      phase += m_freq_pars[0];

      // Add the frequency term, if necessary.
      if (m_freq_pars.size() > 1) {
        // Split the frequency into an integral part and a fractional.
        // Note: This is to give an extra cushon after the least significant digit of the significand of the frequency,
        //       and protect it from computational errors at the lower edge of floating-point expression.
        double freq_int = 0.;
        double freq_frac = std::modf(m_freq_pars[1], &freq_int);

        // Compute an elapsed time in seconds as a pair of an integer part and a fractional.
        long dt_int_long = 0;
        double dt_frac = 0.;
        (ev_time - m_freq_epoch).computeDuration(m_system->getName(), "Sec", dt_int_long, dt_frac);
        double dt_int = dt_int_long;

        // Compute the frequency term of pulse phase.
        phase = trimPhaseValue(freq_int *dt_frac, phase);
        phase = trimPhaseValue(freq_frac*dt_int,  phase);
        phase = trimPhaseValue(freq_frac*dt_frac, phase);

        // Add the frequency derivative terms, if necessary.
        if (m_freq_pars.size() > 2) {
          // Compute an elapsed time in seconds as a double precision floating-point number.
          // Note: Typically, the precision of a double variable exceedes measurement accuracy of frequency derivatives,
          //       and it is precise enough to compute these terms using regular floating-point computations with double.
          double dt = (ev_time - m_freq_epoch).computeDuration(m_system->getName(), "Sec");
          double factor = dt;
          for (std::vector<double>::size_type exponent = 2; exponent < m_freq_pars.size(); ++exponent) {
            factor *= dt / exponent;
            phase = trimPhaseValue(m_freq_pars[exponent]*factor, phase);
          }
        }
      }
    }

    // Compute a phase value that comes from wave parameters.
    if (m_wave_sine.size() > 0 || m_wave_cosine.size() > 0) {
      // Require the frequency parameter.
      // Note: Because m_wave_sine and m_wave_cosine represent a time residual, not a phase residual,
      //       an approximate pulse frequency is needed to convert a time residual to a phase residual.
      if (m_freq_pars.size() < 2) throw std::runtime_error("Pulse frequency must be given when wave parameters are given");
      const double & pulse_frequency = m_freq_pars[1];

      // Compute an elapsed time in seconds.
      double dt_wave = (ev_time - m_freq_epoch).computeDuration(m_system->getName(), "Day");

      // Compute sine and cosine terms and update the phase value to return.
      for (wave_type::size_type wave_index = 0; wave_index < std::max(m_wave_sine.size(), m_wave_cosine.size()); ++wave_index) {
        double argument = m_wave_omega * (wave_index + 1) * dt_wave;

        // Add sine components.
        if (wave_index < m_wave_sine.size()) {
          double sine_term = m_wave_sine[wave_index] * pulse_frequency * std::sin(argument);
          phase = trimPhaseValue(sine_term, phase);
        }

        // Add cosine components.
        if (wave_index < m_wave_cosine.size()) {
          double cosine_term = m_wave_cosine[wave_index] * pulse_frequency * std::cos(argument);
          phase = trimPhaseValue(cosine_term, phase);
        }
      }
    }

    // Compute a phase value that comes from glitches.
    for (glitch_type::const_iterator glitch_itor = m_glitch_list.begin(); glitch_itor != m_glitch_list.end(); ++glitch_itor) {
      if (ev_time >= glitch_itor->m_epoch) {
        const jump_type & perm_jump(glitch_itor->m_perm_jump);
        const decay_type & decay_comp(glitch_itor->m_decay_comp);

        // Add permanent jumps.
        if (perm_jump.size() > 0) {
          // Add the phase jump.
          phase = trimPhaseValue(perm_jump[0], phase);

          // Add the frequency jumps.
          if (perm_jump.size() > 1) {
            double dt_jump = (ev_time - glitch_itor->m_epoch).computeDuration(m_system->getName(), "Sec");
            double factor_jump = 1.;
            for (jump_type::size_type jump_order = 1; jump_order < perm_jump.size(); ++jump_order) {
              factor_jump *= dt_jump / jump_order;
              phase = trimPhaseValue(perm_jump[jump_order]*factor_jump, phase);
            }
          }
        }

        // Add decaying components.
        if (glitch_itor->m_decay_comp.size() > 0) {
          double dt_decay = (ev_time - glitch_itor->m_epoch).computeDuration(m_system->getName(), "Day");
          for (decay_type::const_iterator decay_itor = decay_comp.begin(); decay_itor != decay_comp.end(); ++decay_itor) {
            const double & amplitude(decay_itor->first);
            const double & decay_time(decay_itor->second);
            double decay_term = amplitude * decay_time * SecPerDay() * (1. - std::exp(-dt_decay/decay_time));
            phase = trimPhaseValue(decay_term, phase);
          }
        }
      }
    }

    // Express phase as a value between 0. and 1., after adding a global phase offset.
    return trimPhaseValue(phase, phase_offset);
  }

  double HighPrecisionEph::calcFrequency(const AbsoluteTime & ev_time, int derivative_order) const {
    // Throw an exception for a negative order of frequency derivative.
    if (0 > derivative_order) {
      std::ostringstream os;
      os << "Negative order of frequency derivative is given: " << derivative_order;
      throw std::runtime_error(os.str());
    }

    // Prepare a return value.
    double return_value = 0.;

    // Compute contributions from frequency parameters.
    if (m_freq_pars.size() > static_cast<freq_type::size_type>(derivative_order + 1)) {
      return_value += m_freq_pars[derivative_order + 1];

      // Add contributions from higher derivatives.
      if (m_freq_pars.size() > static_cast<freq_type::size_type>(derivative_order + 2)) {
        double dt = (ev_time - m_freq_epoch).computeDuration(m_system->getName(), "Sec");
        double factor = 1.;
        for (freq_type::size_type exponent = 1; exponent + 1 + derivative_order < m_freq_pars.size(); ++exponent) {
          factor *= dt / exponent;
          return_value += m_freq_pars[exponent + 1 + derivative_order] * factor;
        }
      }
    }

    // Compute contributions from wave parameters.
    if (m_wave_sine.size() > 0 || m_wave_cosine.size() > 0) {
      // Require the pulse_frequency parameter.
      if (m_freq_pars.size() < 2) throw std::runtime_error("Pulse frequency must be given when wave parameters are given");
      const double & pulse_frequency = m_freq_pars[1];

      // Compute an elapsed time in seconds.
      double dt = (ev_time - m_freq_epoch).computeDuration(m_system->getName(), "Sec");

      // Compute sine and cosine terms and update the phase value to return.
      for (wave_type::size_type wave_index = 0; wave_index < std::max(m_wave_sine.size(), m_wave_cosine.size()); ++wave_index) {
        double wave_frequency = m_wave_omega / SecPerDay() * (wave_index + 1);
        double factor_wave = pulse_frequency * wave_frequency;
        for (int ii = 0; ii < derivative_order; ++ii) factor_wave *= wave_frequency;
        double argument = wave_frequency * dt;

        // Add sine components.
        if (wave_index < m_wave_sine.size()) {
          double sine_term = m_wave_sine[wave_index] * factor_wave;
          if ((derivative_order + 1) / 2 % 2) sine_term *= -1.;
          sine_term *= (derivative_order % 2 ? std::sin(argument) : std::cos(argument));
          return_value += sine_term;
        }

        // Add cosine components.
        if (wave_index < m_wave_cosine.size()) {
          double cosine_term = m_wave_cosine[wave_index] * factor_wave;
          if ((derivative_order + 2) / 2 % 2) cosine_term *= -1.;
          cosine_term *= (derivative_order % 2 ? std::cos(argument) : std::sin(argument));
          return_value += cosine_term;
        }
      }
    }

    // Compute contributions from glitches.
    for (glitch_type::const_iterator glitch_itor = m_glitch_list.begin(); glitch_itor != m_glitch_list.end(); ++glitch_itor) {
      if (ev_time >= glitch_itor->m_epoch) {
        const jump_type & perm_jump(glitch_itor->m_perm_jump);
        const decay_type & decay_comp(glitch_itor->m_decay_comp);

        // Add permanent jumps.
        if (perm_jump.size() > static_cast<freq_type::size_type>(derivative_order + 1)) {
          return_value += perm_jump[derivative_order + 1];

          // Add contributions from higher derivatives.
          if (perm_jump.size() > static_cast<freq_type::size_type>(derivative_order + 2)) {
            double dt_jump = (ev_time - glitch_itor->m_epoch).computeDuration(m_system->getName(), "Sec");
            double factor_jump = 1.;
            for (jump_type::size_type jump_order = derivative_order + 2; jump_order < perm_jump.size(); ++jump_order) {
              factor_jump *= dt_jump / (jump_order - derivative_order - 1);
              return_value += perm_jump[jump_order] * factor_jump;
            }
          }
        }

        // Add decaying components.
        if (glitch_itor->m_decay_comp.size() > 0) {
          double dt_decay = (ev_time - glitch_itor->m_epoch).computeDuration(m_system->getName(), "Day");
          for (decay_type::const_iterator decay_itor = decay_comp.begin(); decay_itor != decay_comp.end(); ++decay_itor) {
            const double & amplitude(decay_itor->first);
            const double & decay_time(decay_itor->second);
            double factor_decay = amplitude;
            for (int ii = 0; ii < derivative_order; ++ii) factor_decay /= -decay_time * SecPerDay();
            return_value += factor_decay * std::exp(-dt_decay/decay_time);
          }
        }
      }
    }

    // Return the computed value.
    return return_value;
  }

  SourcePosition HighPrecisionEph::calcPosition(const AbsoluteTime & ev_time) const {
    // Compute the original source position.
    SourcePosition original_srcpos(m_ra, m_dec);

    // Compute the parallax distance, or set a negative value if unknown.
    double distance = -1.;
    if (m_parallax > 0.) distance = s_lts_per_au / (m_parallax * s_rad_per_mas);

    // Check whether RA and Dec change in time.
    if (m_ra_vel == 0. && m_dec_vel == 0.) {
      // Check whether the parallax distance is known.
      if (distance > 0.) {
        // Correct the distance for the radial velocity.
        if (m_radial_vel != 0.) {
          double dt_sec = (ev_time - m_pos_epoch).computeDuration(m_system->getName(), "Sec");
          distance += m_radial_vel / s_km_per_lts * dt_sec;
        }

        // Return the source position with the original direction and the computed distance.
        return SourcePosition(original_srcpos.getDirection(), distance);

      } else {
        // Return the source position at the epoch.
        // Note: Need to ignore the radial velocity because the parallax distance is unknown.
        return original_srcpos;
      }

    } else {
      // Convert physical units.
      double ra_rad = m_ra * s_rad_per_deg; // From degrees to radians.
      double dec_rad = m_dec * s_rad_per_deg; // From degrees to radians.
      double ra_vel_rad = m_ra_vel * s_rad_day_per_mas_year; // From milliarcseconds per Julian year (365.25 days) to radians per day.
      double dec_vel_rad = m_dec_vel * s_rad_day_per_mas_year; // From milliarcseconds per Julian year (365.25 days) to radians per day.

      // Compute the proper motion vector, in the units of radian per day.
      std::vector<double> vel_vector(3);
      vel_vector[0] = -std::sin(ra_rad) * std::cos(dec_rad) * ra_vel_rad - std::cos(ra_rad) * std::sin(dec_rad) * dec_vel_rad;
      vel_vector[1] = +std::cos(ra_rad) * std::cos(dec_rad) * ra_vel_rad - std::sin(ra_rad) * std::sin(dec_rad) * dec_vel_rad;
      vel_vector[2] = +std::cos(dec_rad) * dec_vel_rad;

      // Compute the elapsed time in days.
      double dt_day = (ev_time - m_pos_epoch).computeDuration(m_system->getName(), "Day");

      // Check whether the parallax distance is known.
      if (distance > 0.) {
        // Compute the source position vector at the epoch, in the units of light-second.
        std::vector<double> pos_vector = original_srcpos.getDirection();
        for (std::vector<double>::size_type ii = 0; ii < 3; ++ii) pos_vector[ii] *= distance;

        // Convert the unit of the proper motion vector from radian per day to light-second per day.
        for (std::vector<double>::size_type ii = 0; ii < 3; ++ii) vel_vector[ii] *= distance;

        // Correct the proper motion vector for the radial velocity.
        if (m_radial_vel != 0.) {
          const std::vector<double> & radial_vel_dir = original_srcpos.getDirection();
          for (std::vector<double>::size_type ii = 0; ii < 3; ++ii) {
            vel_vector[ii] += radial_vel_dir[ii] * m_radial_vel / s_km_per_lts * SecPerDay();
          }
        }

        // Correct the source position vector for the proper motion.
        for (std::vector<double>::size_type ii = 0; ii < 3; ++ii) pos_vector[ii] += vel_vector[ii] * dt_day;

        // Return the source position expressed by the position vector.
        double distance_squared = 0.;
        for (std::vector<double>::size_type ii = 0; ii < 3; ++ii) distance_squared += pos_vector[ii] * pos_vector[ii];
        return SourcePosition(pos_vector, std::sqrt(distance_squared));

      } else {
        // Correct the source direction for the proper motion.
        // Note: Need to ignore the radial velocity because the parallax distance is unknown.
        std::vector<double> src_direction = original_srcpos.getDirection();
        for (std::vector<double>::size_type ii = 0; ii < 3; ++ii) src_direction[ii] += vel_vector[ii] * dt_day;

        // Return the corrected source direction.
        return SourcePosition(src_direction);
      }
    }

    // Throw an exception, because this should not be reached.
    std::ostringstream oss;
    oss << "Logic error in computing the source position for " << ev_time;
    throw std::runtime_error(oss.str());
  }

  void HighPrecisionEph::setRemark() {
    // Add glitches to the internal container of ephemeris remarks.
    for (glitch_type::const_iterator glitch_itor = m_glitch_list.begin(); glitch_itor != m_glitch_list.end(); ++glitch_itor) {
      std::ostringstream oss;
      oss << "Glitch observed at " << glitch_itor->m_epoch.represent("TDB", MjdFmt);
      m_remark_cont.push_back(EphStatus(glitch_itor->m_epoch, m_until, Remarked, oss.str()));
    }
  }

}
