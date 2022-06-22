// Thorn      : NRPyEOS
// File       : NRPyEOS_Tabulated_headers.hh
// Author(s)  : Leo Werneck (wernecklr@gmail.com)
// Description: This file contains all the variables, structs, and
//              function prototypes we need in the rest of the thorn.
//              It also provides an interface with EOS_Omni variables.

#ifndef NRPyEOS_TABULATED_HEADERS_H
#define NRPyEOS_TABULATED_HEADERS_H

#include "NRPyEOS.h"

// Table reader
void NRPyEOS_readtable_set_EOS_params(const char *nuceos_table_name, NRPyEOS_params_tabulated *restrict eos_params);

// Free all memory allocated for the table
void NRPyEOS_free_memory();

// ------------------------------------------------------
// ------------- New general interpolators --------------
// ------------------------------------------------------
void NRPyEOS_from_rho_Ye_T_interpolate_n_quantities( const NRPyEOS_params_tabulated *restrict eos_params,
                                                     const int  n,
                                                     const double rho,
                                                     const double Ye,
                                                     const double T,
                                                     const int *restrict tablevars_keys,
                                                     double *restrict tablevars,
                                                     NRPyEOS_error_report *restrict report );

void NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities( const NRPyEOS_params_tabulated *restrict eos_params,
                                                                  const int  n,
                                                                  const double prec,
                                                                  const double rho,
                                                                  const double Ye,
                                                                  const double tablevar_in,
                                                                  const int  tablevar_in_key,
                                                                  const int *restrict tablevars_keys,
                                                                  double *restrict tablevars,
                                                                  double *restrict T,
                                                                  NRPyEOS_error_report *restrict report );

// ------------------------------------------------------
// ------ Functions where the temperature is known ------
// ------------------------------------------------------
void NRPyEOS_P_and_eps_from_rho_Ye_T( const NRPyEOS_params_tabulated *restrict eos_params,
                                      const double rho,
                                      const double Ye,
                                      const double T,
                                      double *restrict P,
                                      double *restrict eps );

void NRPyEOS_P_eps_and_S_from_rho_Ye_T( const NRPyEOS_params_tabulated *restrict eos_params,
                                        const double rho,
                                        const double Ye,
                                        const double T,
                                        double *restrict P,
                                        double *restrict eps,
                                        double *restrict S );

void NRPyEOS_P_eps_S_and_cs2_from_rho_Ye_T( const NRPyEOS_params_tabulated *restrict eos_params,
                                            const double rho,
                                            const double Ye,
                                            const double T,
                                            double *restrict P,
                                            double *restrict eps,
                                            double *restrict S,
                                            double *restrict cs2 );

void NRPyEOS_P_eps_and_depsdT_from_rho_Ye_T( const NRPyEOS_params_tabulated *restrict eos_params,
                                             const double rho,
                                             const double Ye,
                                             const double T,
                                             double *restrict P,
                                             double *restrict eps,
                                             double *restrict depsdT );

void NRPyEOS_P_eps_dPdrho_dPdT_depsdrho_and_depsdT_from_rho_Ye_T( const NRPyEOS_params_tabulated *restrict eos_params,
                                                                  const double rho,
                                                                  const double Ye,
                                                                  const double T,
                                                                  double *restrict P,
                                                                  double *restrict eps,
                                                                  double *restrict dPdrho,
                                                                  double *restrict dPdT,
                                                                  double *restrict depsdrho,
                                                                  double *restrict depsdT );

void NRPyEOS_mue_mup_mun_muhat_Xn_and_Xp_from_rho_Ye_T( const NRPyEOS_params_tabulated *restrict eos_params,
                                                        const double rho,
                                                        const double Ye,
                                                        const double T,
                                                        double *restrict mu_e,
                                                        double *restrict mu_p,
                                                        double *restrict mu_n,
                                                        double *restrict muhat,
                                                        double *restrict X_p,
                                                        double *restrict X_n );

// ------------------------------------------------------
// ---- Functions where the temperature is not known ----
// ------------------------------------------------------
void NRPyEOS_P_and_T_from_rho_Ye_eps( const NRPyEOS_params_tabulated *restrict eos_params,
                                      const double rho,
                                      const double Ye,
                                      const double eps,
                                      double *restrict P,
                                      double *restrict T );

void NRPyEOS_P_S_and_T_from_rho_Ye_eps( const NRPyEOS_params_tabulated *restrict eos_params,
                                        const double rho,
                                        const double Ye,
                                        const double eps,
                                        double *restrict P,
                                        double *restrict S,
                                        double *restrict T );

void NRPyEOS_eps_S_and_T_from_rho_Ye_P( const NRPyEOS_params_tabulated *restrict eos_params,
                                        const double rho,
                                        const double Ye,
                                        const double P,
                                        double *restrict eps,
                                        double *restrict S,
                                        double *restrict T );

void NRPyEOS_P_eps_and_T_from_rho_Ye_S( const NRPyEOS_params_tabulated *restrict eos_params,
                                        const double rho,
                                        const double Ye,
                                        const double S,
                                        double *restrict P,
                                        double *restrict eps,
                                        double *restrict T );

// ----------------------------------------

#endif // NRPyEOS_TABULATED_HEADERS_H
