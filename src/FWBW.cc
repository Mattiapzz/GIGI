#include "GIGI/FWBW.hxx"
#include "GIGI/gg_utils.hxx"

#define DEBUG_FB 0

#include <iostream>

using namespace GG;

constexpr size_t DEFAULT_SIZE{100};

// --------------------------------------------------------------------------------------------

FWBW::FWBW(
  const std::function<real(real, real)> &gg_Upper,
  const std::function<real(real, real)> &gg_Lower,
  const gg_range_max_min &gg_range
)
    : gg_Upper(gg_Upper), gg_Lower(gg_Lower), gg_range(gg_range)
{
  this->Segments.reserve(DEFAULT_SIZE);
  this->Vmax_vec.reserve(DEFAULT_SIZE);
  this->dump_seg_id.reserve(DEFAULT_SIZE);
}

FWBW::FWBW()
  : gg_Upper(nullptr), gg_Lower(nullptr), gg_range({nullptr, nullptr})
{
  this->Segments.reserve(DEFAULT_SIZE);
  this->Vmax_vec.reserve(DEFAULT_SIZE);
  this->dump_seg_id.reserve(DEFAULT_SIZE);
}

void 
FWBW::setup_functions(
  const std::function<real(real, real)> &gg_Upper,
  const std::function<real(real, real)> &gg_Lower,
  const gg_range_max_min &gg_range
)
{
  this->gg_Upper = gg_Upper;
  this->gg_Lower = gg_Lower;
  this->gg_range = gg_range;
}

// --------------------------------------------------------------------------------------------

bool FWBW::is_in_range(const real ax, const real ay, const real v) const
{
  // Check if ax is within the upper or lower bounds
  return ax <= this->gg_Upper(ay, v) && ax >= this->gg_Lower(ay, v) && ay >= this->gg_range.min(v)
         && ay <= this->gg_range.max(v);
}

// --------------------------------------------------------------------------------------------

real FWBW::signed_distance(const real ax, const real ay, const real v) const
{
  real ayclip = clip(ay, this->gg_range.min(v), this->gg_range.max(v));
  return GG::signed_distance(
    ax,
    this->gg_Lower(ayclip, v),
    this->gg_Upper(ayclip, v),
    ay,
    this->gg_range.min(v),
    this->gg_range.max(v)
  );
}

// --------------------------------------------------------------------------------------------

void FWBW::compute_Vmax()
{
  // extimate a small value of curvature.
  constexpr real k_small     = 1e-4;
  constexpr real v_top_speed = 130;
  real vmax                  = QUIET_NAN;
  // chose bracketing interval
  constexpr real v_b = v_top_speed;
  // instantiate vector same size of SS
  this->Vmax_vec.clear();
  this->Vmax_vec.resize(this->S_vec.size());
  // fill the vector
  for (int i = 0; i < this->S_vec.size(); i++)
  {
    constexpr real v_a = 0;
    real k = this->K_vec[i];
    if (std::abs(k) < k_small)
    {
      vmax = v_top_speed;
    }
    else if (k >= k_small)
    {
      auto F2solve = [this, k](const real v) -> real { return k * v * v - this->gg_range.max(v); };
      bool ok = this->BD.solve(F2solve, v_a, v_b, vmax);
    }
    else
    {
      auto F2solve = [this, k](const real v) -> real { return k * v * v - this->gg_range.min(v); };
      bool ok = this->BD.solve(F2solve, v_a, v_b, vmax);
    }
    vmax              = std::isnan(vmax) ? v_top_speed : vmax;
    this->Vmax_vec[i] = vmax;
  }
}

// --------------------------------------------------------------------------------------------

real FWBW::compute(std::vector<real> const &SS, std::vector<real> const &KK, const real v0)
{
  this->S_vec = SS;
  this->K_vec = KK;
  this->v_I   = v0;
  this->dump_seg_id.clear();
  this->compute_Vmax();
  this->FW();
  this->BW();
  return this->compute_time();
}


// --------------------------------------------------------------------------------------------

real FWBW::compute_time() 
{
  real T = 0;
  for (auto &seg : this->Segments)
  {
    seg.set_times(T);
    T += seg.getT();
  }
  return T;
}

// --------------------------------------------------------------------------------------------

void FWBW::FW()
{
  real v0 = this->v_I;
  this->Vmax_vec[0] = v0;
  this->Segments.clear();
  this->Segments.resize(this->S_vec.size() - 1);
  for (int i = 0; i < this->S_vec.size() - 1; i++)
  {
    // extract data for convenience of notation
    const real k0 = this->K_vec[i];
    const real k1 = this->K_vec[i + 1];
    const real S0 = this->S_vec[i];
    const real S1 = this->S_vec[i + 1];
    const real L0 = S1 - S0;
    this->Segments[i] = segment(S0, L0, v0, k0, k1);
    const real ay0     = this->Segments[i].AY0();
    const real ay0clip = clip(ay0, this->gg_range.min(v0), this->gg_range.max(v0));
    const real axmax0  = this->gg_Upper(ay0clip, v0);
    const real axmin0  = this->gg_Lower(ay0clip, v0);
    const real distance_amax  = this->signed_distance( axmax0, this->Segments[i].AYA(axmax0), this->Segments[i].VA(axmax0));
    const real distance_azero = this->signed_distance(0.0, this->Segments[i].AYA(0.0), this->Segments[i].VA(0.0));
    bool ok = false;
    real ax0 = 0.0;
    if (distance_amax <= 0)
    {
      ax0 = axmax0;
      ok  = true;
    }
    else
    {
      const real a_solver = std::min( 
        this->gg_Lower(this->gg_range.min(v0), v0), 
        this->gg_Lower(this->gg_range.max(v0), v0)
      );
      ok = this->BD.solve(
        [this, &i, &v0](const real ax) -> real
        { return this->signed_distance(ax, this->Segments[i].AYA(ax), this->Segments[i].VA(ax)); },
        a_solver,
        axmax0,
        ax0
      );
    }
    if(ok)
    {
      this->Segments[i].set_a(ax0);
      this->Segments[i].set_type(FORWARD);
      v0 = std::max(0.0, std::min(this->Vmax_vec[i + 1], this->Segments[i].VF()));
      this->Vmax_vec[i + 1] = v0;
    }
    else 
    {
      this->Segments[i].set_a(QUIET_NAN);
      this->Segments[i].set_type(FORWARD_NAN);
      v0 = std::max(0.0, std::min(this->Vmax_vec[i + 1], v0));
      this->Vmax_vec[i + 1] = v0;
    }
  }
}

// --------------------------------------------------------------------------------------------

void FWBW::BW()
{
  real v1 = this->Vmax_vec.back();
  if (this->Segments.back().type() == FORWARD)
  {
    v1 = this->Segments.back().VF();
  }
  // explore in reverse order
  for (auto i = static_cast<integer>(this->Segments.size() - 1); i >= 0; i--)
  {
    const real v0  = this->Segments[i].v0();
    const real kf  = this->Segments[i].k1();
    const real ayf = kf * v1 * v1;
    const real ayfclip = clip(ayf, this->gg_range.min(v1), this->gg_range.max(v1));
    const real axmaxf  = this->gg_Upper(ayfclip, v1);
    const real axminf  = this->gg_Lower(ayfclip, v1);
    const real amean = (v1 * v1 - v0 * v0) / (2 * this->Segments[i].L());
    const real v0_reach_max = std::min(this->Vmax_vec[i], this->Segments[i].VB(axminf, v1));
    const real v0_reach_min = std::max(0.0, this->Segments[i].VB(axmaxf, v1));
    const bool is_v0_reachable = (v0 >= v0_reach_min && v0 <= v0_reach_max);
    const bool is_amean_candidate = (amean >= axminf && amean <= axmaxf);
    const bool is_amean_valid = (this->signed_distance(amean, this->Segments[i].AYB(amean, v1), this->Segments[i].VB(amean, v1)) <= this->solver_p.tolerance);
    const bool is_valid_forward = ((this->Segments[i].type() == FORWARD) &&  (std::abs(this->Segments[i].VF() - v1) <= this->solver_p.tolerance));
    if(is_valid_forward)
    {
      v1 = v0; 
      continue;
    }
    if (is_amean_candidate && is_v0_reachable && is_amean_valid )
    {
      this->Segments[i].set_a(amean);
      this->Segments[i].set_type(TRANSITION);
      v1 = this->Segments[i].VB(amean, v1);
      continue;
    }
    else 
    {
      real distance_amax = this->signed_distance(axmaxf, this->Segments[i].AYB(axmaxf, v1), this->Segments[i].VB(axmaxf, v1));
      real distance_azero = this->signed_distance(0.0, this->Segments[i].AYB(0.0, v1), this->Segments[i].VB(0.0, v1));
      real distance_amin = this->signed_distance(axminf, this->Segments[i].AYB(axminf, v1), this->Segments[i].VB(axminf, v1));
      //
      real ax0 = 0.0;
      bool ok = false;
      if (distance_amin <= 0)
      {
        ax0 = axminf;
        ok = true;
      }
      else 
      {
        const real a_solver = std::max(
          this->gg_Upper(this->gg_range.min(v1), v1),
          this->gg_Upper(this->gg_range.max(v1), v1)
        );
        ok = this->BD.solve(
          [this, &i, &v1](const real ax) -> real {
            return this->signed_distance(ax, this->Segments[i].AYB(ax, v1), this->Segments[i].VB(ax, v1));
          },
          axminf,
          a_solver,
          ax0
        );
      }
      if (ok)
      {
        this->Segments[i].set_a(ax0);
        this->Segments[i].set_v0(this->Segments[i].VB(ax0, v1));
        this->Segments[i].set_type(BACKWARD);
        v1 = this->Segments[i].VB(ax0, v1);
      }
      else
      {
        std::cout << "FWBW::BW() >> No solution found for segment " << i << "\n";
        this->dump_seg_id.push_back(i);
        this->Segments[i].set_a(QUIET_NAN);
        this->Segments[i].set_type(BACKWARD_NAN);
        v1 = std::max(0.0, std::min(this->Vmax_vec[i], v1));
      }
    }
  }
}

// --------------------------------------------------------------------------------------------

void FWBW::evaluate(
  std::vector<real> const &SS, std::vector<real> &AX, std::vector<real> &AY, std::vector<real> &V
)
{
  for (integer i = 0; i < SS.size(); i++)
  {
    real s = SS[i];
    // Find the segment that contains s
    auto segment_it = std::find_if(
      this->Segments.begin(),
      this->Segments.end(),
      [s](const segment &seg) { return s >= seg.s0() && s <= seg.s1(); }
    );
    if (segment_it != this->Segments.end())
    {
      // Populate AX, AY, and V using the found segment
      AX[i] = segment_it->AX(s);
      AY[i] = segment_it->AY(s);
      V[i]  = segment_it->V(s);
    }
    else
    {
      // Handle the case where no segment is found
      std::cerr << "Error: No segment found for s = " << s << "\n";
      AX[i] = AY[i] = V[i] = 0.0;
    }
  }
}

// --------------------------------------------------------------------------------------------

integer FWBW::get_seg_idx(const real s) const
{
  integer seg_idx = -1;
  for (integer i = 0; i < this->Segments.size(); i++)
  {
    if (s >= this->Segments[i].s0() && s <= this->Segments[i].s1())
    {
      seg_idx = i;
      break;
    }
  }
  if (seg_idx == -1)
  {
    std::cerr << "Error: No segment found for s = " << s << "\n";
    return -1;
  }
  return seg_idx;
}

// --------------------------------------------------------------------------------------------

integer FWBW::get_seg_idx_t(const real t) const
{
  integer seg_idx = -1;
  for (integer i = 0; i < this->Segments.size(); i++)
  {
    if (t >= this->Segments[i].getT0() && t <= this->Segments[i].getT1())
    {
      seg_idx = i;
      break;
    }
  }
  if (seg_idx == -1)
  {
    std::cerr << "Error: No segment found for t = " << t << "\n";
    return -1;
  }
  return seg_idx;
}

// --------------------------------------------------------------------------------------------

real FWBW::evalS(const real t) const
{
  const integer seg_idx = this->get_seg_idx_t(t);
  return this->Segments[seg_idx].S(t);
}

// --------------------------------------------------------------------------------------------

real FWBW::evalV(const real s) const
{
  const integer seg_idx = this->get_seg_idx(s);
  return this->Segments[seg_idx].V(s - this->Segments[seg_idx].s0());
}

// --------------------------------------------------------------------------------------------

real FWBW::evalV_t(const real t) const
{
  const integer seg_idx = this->get_seg_idx_t(t);
  auto const s = this->Segments[seg_idx].S(t);
  return this->Segments[seg_idx].V(s - this->Segments[seg_idx].s0());
}

// --------------------------------------------------------------------------------------------

real FWBW::evalAx(const real s) const
{
  const integer seg_idx = this->get_seg_idx(s);
  return this->Segments[seg_idx].AX(s - this->Segments[seg_idx].s0());
}

// --------------------------------------------------------------------------------------------

real FWBW::evalAx_t(const real t) const
{
  const integer seg_idx = this->get_seg_idx_t(t);
  auto const s = this->Segments[seg_idx].S(t);
  return this->Segments[seg_idx].AX(s - this->Segments[seg_idx].s0());
}

// --------------------------------------------------------------------------------------------

real FWBW::evalAy(const real s) const
{
  const integer seg_idx = this->get_seg_idx(s);
  return this->Segments[seg_idx].AY(s - this->Segments[seg_idx].s0());
}

// --------------------------------------------------------------------------------------------

real FWBW::evalAy_t(const real t) const
{
  const integer seg_idx = this->get_seg_idx_t(t);
  auto const s = this->Segments[seg_idx].S(t);
  return this->Segments[seg_idx].AY(s - this->Segments[seg_idx].s0());
}

// --------------------------------------------------------------------------------------------

real FWBW::evalVmax(const real s) const
{
  const integer seg_idx = this->get_seg_idx(s);
  const real s_norm = (s - this->Segments[seg_idx].s0()) / this->Segments[seg_idx].L();
  return this->Vmax_vec[seg_idx] * (1-s_norm) + this->Vmax_vec[seg_idx + 1] * s_norm;
}

// --------------------------------------------------------------------------------------------
