#include "GIGI/FB_F1_10.hxx"
#include "GIGI/gg_utils.hxx"
#include "json.hpp"
#include <functional>

#define DEBUG_FB 0

#include <iostream>

using namespace GG;


FB_F1_10::FB_F1_10( 
  const std::vector<real> & ay_spline,
  const std::vector<real> & vx_spline,
  const std::vector<real> & ax_bispline_max,
  const std::vector<real> & ax_bispline_min,
  const std::vector<real> & ay_spline_max,
  const std::vector<real> & ay_spline_min)
: FWBW(), sp_data({ay_spline, vx_spline, 
    ax_bispline_max, ax_bispline_min, 
    ay_spline_max, ay_spline_min})
{
  setup_splines();
}

// --------------------------------------------------------------------------------------------

FB_F1_10::FB_F1_10( 
  const ggv_spline_data & spline_data)
: FWBW(), sp_data(spline_data) 
{
  setup_splines();
}

// --------------------------------------------------------------------------------------------

void 
FB_F1_10::setup_splines()
{
  // create the ax_max and ax_min interpolators
  this->ax_max.setup(sp_data.ay, sp_data.vx, sp_data.ax_max);
  this->ax_min.setup(sp_data.ay, sp_data.vx, sp_data.ax_min);
  // create the ay_lim interpolator
  this->ay_lim.setup(sp_data.vx, {sp_data.ay_max, sp_data.ay_min},{"ay_max", "ay_min"}); 
  this->GG_UP    = [this](real ay, real v) -> real { return this->ax_max.eval(ay, v);};
  this->GG_LO    = [this](real ay, real v) -> real { return this->ax_min.eval(ay, v);};
  this->GG_R_MIN = [this](real v)          -> real { return this->ay_lim.eval("ay_min", v);};
  this->GG_R_MAX = [this](real v)          -> real { return this->ay_lim.eval("ay_max", v);};
  this->setup_functions(GG_UP, GG_LO, {GG_R_MIN, GG_R_MAX});
}

