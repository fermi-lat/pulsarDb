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

  class PulsarToolApp : public st_app::StApp {
    public:
      typedef std::vector<const tip::Table *> table_cont_type;

      PulsarToolApp();
      virtual ~PulsarToolApp() throw();
      virtual void run() = 0;

      virtual timeSystem::TimeRep * createTimeRep(const std::string & time_format, const std::string & time_system,
        const std::string & time_value);

      virtual timeSystem::TimeRep * createTimeRep(const std::string & time_format, const std::string & time_system,
        const std::string & time_value, const tip::Header & header);

      void openEventFile(const st_app::AppParGroup & pars);

      void initEphComputer(const st_app::AppParGroup & pars, const TimingModel & model, const EphChooser & chooser);

      void computeTimeBoundary(timeSystem::AbsoluteTime & abs_tstart, timeSystem::AbsoluteTime & abs_tstop);

      timeSystem::AbsoluteTime initTargetTime(const st_app::AppParGroup & pars,
        const timeSystem::AbsoluteTime & abs_tstart, const timeSystem::AbsoluteTime & abs_tstop);

      PulsarEph & updateEphComputer(const timeSystem::AbsoluteTime & abs_time);

      void initTimeCorrection(const st_app::AppParGroup & pars);

      double computeElapsedSecond(const timeSystem::AbsoluteTime & abs_time);

      void setFirstEvent(const st_app::AppParGroup & pars);

      void setNextEvent();

      bool isEndOfEventList();

      timeSystem::AbsoluteTime getEventTime();

    private:
      table_cont_type m_event_table_cont;
      table_cont_type m_gti_table_cont;
      const tip::Header * m_reference_header;
      std::map<const tip::Table *, timeSystem::TimeRep *> m_time_rep_dict;
      std::map<const tip::Table *, bool> m_need_bary_dict;
      EphComputer * m_computer;
      bool m_request_bary;
      bool m_demod_bin;
      bool m_cancel_pdot;
      timeSystem::TimeRep * m_target_time_rep;
      std::string m_time_field;

      table_cont_type::iterator m_table_itor;
      tip::Table::ConstIterator m_event_itor;

      // TODO: refactor MetRep to include functionality of this method and remove this method.
      timeSystem::TimeRep * createMetRep(const std::string & time_system, const timeSystem::AbsoluteTime & abs_reference);

      timeSystem::AbsoluteTime readTimeColumn(const tip::Table & table, tip::ConstTableRecord & record,
        const std::string & column_name);
  };

}

#endif