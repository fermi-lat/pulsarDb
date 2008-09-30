/** \file PulsarDbApp.h
    \brief Interface for PulsarDbApp class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_PulsarDbApp_h
#define pulsarDb_PulsarDbApp_h

#include "st_app/StApp.h"

#include "st_stream/StreamFormatter.h"

namespace pulsarDb {

  class PulsarDb;

  /** \class PulsarDbApp
      \brief Main application class for pulsar database access.
  */
  class PulsarDbApp: public st_app::StApp {
    public:
      PulsarDbApp();
      virtual ~PulsarDbApp() throw();
      virtual void run();

    private:
      st_stream::StreamFormatter m_os;
  };

}
#endif
