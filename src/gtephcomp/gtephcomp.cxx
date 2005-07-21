/** \file gtephcomp.cxx
    \brief Factory for gtephcomp application.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#include "st_app/StAppFactory.h"
#include "pulsarDb/EphComputerApp.h"

st_app::StAppFactory<pulsarDb::EphComputerApp> g_factory("gtephcomp");
