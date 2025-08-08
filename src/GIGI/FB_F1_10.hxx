/*
(***********************************************************************************)
(*                                                                                 *)
(* The GIGI project                                                               *)
(*                                                                                 *)
(* Copyright (c) 2025, Mattia Piazza                                               *)
(*                                                                                 *)
(*    Mattia Piazza                                                                *)
(*    Department of Industrial Engineering                                         *)
(*    University of Trento                                                         *)
(*    e-mail: mattia.piazza@unitn.it                                               *)
(*                                                                                 *)
(***********************************************************************************)
*/

#ifndef FB_F1_10_HH
#define FB_F1_10_HH

#include "FWBW.hxx"
#include "types.hxx"
#include <functional>

namespace GG
{

  struct ggv_spline_data
  {
    std::vector<real> ay; // Lateral acceleration [m/s²]
    std::vector<real> vx; // Longitudinal velocity [m/s]
    std::vector<real> ax_max; // Maximum longitudinal acceleration [m/s²]
    std::vector<real> ax_min; // Minimum longitudinal acceleration [m/s²]
    std::vector<real> ay_max; // Maximum lateral acceleration [m/s²]
    std::vector<real> ay_min; // Minimum lateral acceleration [m/s²]
  };

  class FB_F1_10 : public FWBW
  {
  private:
    /* data */
    ggv_spline_data sp_data; // Spline data
    GG::BilinearInterpolator ax_max;
    GG::BilinearInterpolator ax_min;
    GG::LinearInterpolatorSet ay_lim;
    //
    std::function<real(real, real)> GG_UP = nullptr; // Upper bound function
    std::function<real(real, real)> GG_LO = nullptr; // Lower bound function
    std::function<real(real)> GG_R_MIN = nullptr; // Minimum lateral acceleration function
    std::function<real(real)> GG_R_MAX = nullptr; // Maximum lateral acceleration function
  public:
    // constructors

    FB_F1_10( const std::vector<real> & ay_spline,
              const std::vector<real> & vx_spline,
              const std::vector<real> & ax_bispline_max,
              const std::vector<real> & ax_bispline_min,
              const std::vector<real> & ay_spline_max,
              const std::vector<real> & ay_spline_min);

    explicit FB_F1_10( const ggv_spline_data & spline_data); 

  private:
    void setup_splines();

  };

} // namespace GG

#endif // FB_F1_10_HH
