/** \file gtpulsardb.cxx
    \brief Factory for gtpulsardb application.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#include "st_app/StAppFactory.h"
#include "pulsarDb/PulsarDbApp.h"

st_app::StAppFactory<pulsarDb::PulsarDbApp> g_factory("gtpulsardb");
