/** \file PulsarToolApp.h
    \brief Declaration of base class for pulsar tool applications.
    \author Masaharu Hirayama, GSSC
            James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_PulsarToolApp_h
#define pulsarDb_PulsarToolApp_h

#include <map>
#include <string>
#include <vector>

#include "pulsarDb/EphChooser.h"
#include "pulsarDb/EphComputer.h"
#include "pulsarDb/PulsarEph.h"
#include "pulsarDb/TimingModel.h"

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"

#include "timeSystem/AbsoluteTime.h"
#include "timeSystem/TimeRep.h"

#include "tip/Header.h"
#include "tip/Table.h"

namespace pulsarDb {

  // TODO: Merge EphComputer2 into EphComputer once all the pulsar tools use PulsarToolApp.
  class EphComputer2 : public EphComputer {
    public:
      EphComputer2();

      EphComputer2(const TimingModel & model, const EphChooser & chooser);

      ~EphComputer2();

      void setPdotCancelParameter(const PulsarEph & pdot_pars);

      void setPdotCancelParameter(const std::string & time_system_name, const timeSystem::AbsoluteTime & time_origin,
        double f1f0ratio, double f2f0ratio = 0.);

      void cancelPdot(timeSystem::AbsoluteTime & ev_time) const;

    private:
      PulsarEph * m_pdot_pars;
      TimingModel * m_model2; // Not needed when merged to EphComputer.
  };

  class PulsarToolApp : public st_app::StApp {
    public:
      typedef std::vector<const tip::Table *> table_cont_type;

      PulsarToolApp();
      virtual ~PulsarToolApp() throw();
      virtual void run() = 0;

      virtual timeSystem::TimeRep * createTimeRep(const std::string & time_format, const std::string & time_system,
        const std::string & time_value) const;

      virtual timeSystem::TimeRep * createTimeRep(const std::string & time_format, const std::string & time_system,
        const std::string & time_value, const tip::Header & header) const;

      void openEventFile(const st_app::AppParGroup & pars);

      void initEphComputer(const st_app::AppParGroup & pars, const TimingModel & model, const EphChooser & chooser);

      void initTimeCorrection(const st_app::AppParGroup & pars);

      double computeElapsedSecond(const timeSystem::AbsoluteTime & abs_time);

      void setFirstEvent();

      void setNextEvent();

      bool isEndOfEventList() const;

      timeSystem::AbsoluteTime getEventTime();

      timeSystem::AbsoluteTime getStartTime();

      timeSystem::AbsoluteTime getStopTime();

      timeSystem::AbsoluteTime getTimeOrigin(const st_app::AppParGroup & pars);

      EphComputer2 & getEphComputer() const;

    private:
      table_cont_type m_event_table_cont;
      table_cont_type m_gti_table_cont;
      std::string m_time_field;
      std::string m_gti_start_field;
      std::string m_gti_stop_field;
      std::map<const tip::Table *, timeSystem::TimeRep *> m_time_rep_dict;
      std::map<const tip::Table *, bool> m_need_bary_dict;
      const tip::Header * m_reference_header;
      EphComputer2 * m_computer;
      bool m_request_bary;
      bool m_demod_bin;
      bool m_cancel_pdot;
      timeSystem::TimeRep * m_target_time_rep;

      table_cont_type::iterator m_table_itor;
      tip::Table::ConstIterator m_event_itor;

      // TODO: refactor MetRep to include functionality of this method and remove this method.
      timeSystem::TimeRep * createMetRep(const std::string & time_system, const timeSystem::AbsoluteTime & abs_reference) const;

      timeSystem::AbsoluteTime readTimeColumn(const tip::Table & table, tip::ConstTableRecord & record,
        const std::string & column_name, bool request_time_correction = false);

      timeSystem::AbsoluteTime computeTimeBoundary(bool request_start_time, bool request_time_correction);
  };

}

#endif
