#include <cmath>
#include <iostream>
#include <fstream>

#include "GIGI/FWBW.hxx"
#include "GIGI/gg_utils.hxx"
#include "GIGI/types.hxx"

int main() {
  std::cout << "FBGA - Forward Backward Generic Acceleration constraints\n";
  std::cout << " > Running test basic\n";
  // 1. Define track data
  std::vector<GG::real> SS_vec = {0.0,50.0, 100.0, 150.0, 200.0};  // Arc length [m]
  std::vector<GG::real> KK_vec = {0.001, 0.01, 0.005, -0.01, 0.0};    // Curvature [1/m]
  GG::real v_initial = 5.0;  // Initial velocity [m/s]

  std::cout << " > Track data:\n";
  std::cout << " >   SS_vec: " << SS_vec.size() << " points\n";
  std::cout << " >   KK_vec: " << KK_vec.size() << " points\n";
  std::cout << " >   v_initial: " << v_initial << " m/s\n";
  for (size_t i = 0; i < SS_vec.size(); ++i) {
    std::cout << " >   SS[" << i << "]: " << SS_vec[i] 
              << ", KK[" << i << "]: " << KK_vec[i] << "\n";
  }


  GG::LinearInterpolatorSet spline_1D(
    SS_vec,
    { KK_vec },
    {"kappa"});

  GG::integer numpt_eval = std::ceil(SS_vec.back());

  std::vector<GG::real> SS_eval_vec(numpt_eval);
  std::vector<GG::real> KK_eval_vec(numpt_eval);
  for (GG::integer i = 0; i < std::ceil(SS_vec.back()); ++i) {
    SS_eval_vec[i] = static_cast<GG::real>(i);
    KK_eval_vec[i] = spline_1D.eval("kappa", static_cast<GG::real>(i));
  }

  
  // 2. Define splines data
  const std::vector<GG::real> AY_spline = {-5.0,-4.0,-3.0,-2.0,-1.0,0.0,1.0,2.0,3.0,4.0,5.0}; // Lateral acceleration [m/s²]
  const std::vector<GG::real> VX_spline = {0.0, 5.0, 10.0, 15.0, 20.0}; // Longitudinal velocity [m/s]
  const std::vector<GG::real> AX_spline_max = {
    0,0,0,0,0,
    0.360000000000000,0.337500000000000,0.270000000000000,0.157500000000000,0,
    0.640000000000000,0.600000000000000,0.480000000000000,0.280000000000000,0,
    0.840000000000000,0.787500000000000,0.630000000000000,0.367500000000000,0,
    0.960000000000000,0.900000000000000,0.720000000000000,0.420000000000000,0,
    1,0.937500000000000,0.750000000000000,0.437500000000000,0,
    0.960000000000000,0.900000000000000,0.720000000000000,0.420000000000000,0,
    0.840000000000000,0.787500000000000,0.630000000000000,0.367500000000000,0,
    0.640000000000000,0.600000000000000,0.480000000000000,0.280000000000000,0,
    0.360000000000000,0.337500000000000,0.270000000000000,0.157500000000000,0,
    0,0,0,0,0};

    const std::vector<GG::real> AX_spline_min = {
      0,0,0,0,0,
      -0.720000000000000,-0.675000000000000,-0.540000000000000,-0.315000000000000,0,
      -1.28000000000000,-1.20000000000000,-0.960000000000000,-0.560000000000000,0,
      -1.68000000000000,-1.57500000000000,-1.26000000000000,-0.735000000000000,0,
      -1.92000000000000,-1.80000000000000,-1.44000000000000,-0.840000000000000,0,
      -2,-1.87500000000000,-1.50000000000000,-0.875000000000000,0,
      -1.92000000000000,-1.80000000000000,-1.44000000000000,-0.840000000000000,0,
      -1.68000000000000,-1.57500000000000,-1.26000000000000,-0.735000000000000,0,
      -1.28000000000000,-1.20000000000000,-0.960000000000000,-0.560000000000000,0,
      -0.720000000000000,-0.675000000000000,-0.540000000000000,-0.315000000000000,0,
      0,0,0,0,0
    };

  const std::vector<GG::real> AY_spline_max = {5.0, 4.5, 4.5, 4.0, 3.5};
  const std::vector<GG::real> AY_spline_min = {-5.0, -5, -4, -4.0, -3.0};

  GG::BilinearInterpolator ax_max(AY_spline, VX_spline, AX_spline_max);
  GG::BilinearInterpolator ax_min(AY_spline, VX_spline, AX_spline_min);

  GG::LinearInterpolatorSet ay_lim(VX_spline,{AY_spline_max, AY_spline_min},
    {"ay_max", "ay_min"});

    
  // open a file to write the values of 
  std::ofstream output_file("./data/output_fb_spline.txt");
  // 
  GG::integer numpt_enve = 100;
  std::vector<GG::real> VX_linspace(numpt_enve);
  std::vector<GG::real> AY_linspace(numpt_enve);
  for (GG::integer i = 0; i < numpt_enve; ++i) {
    auto stepsize = static_cast<GG::real>(i)  / (numpt_enve - 1);
    VX_linspace[i] = stepsize * (VX_spline.back() - VX_spline.front()) + VX_spline.front();
    AY_linspace[i] = stepsize * (AY_spline.back() - AY_spline.front()) + AY_spline.front();
  }

  // write the header of the file
  output_file << "vx, ay, ay_max, ay_min, ay_clip, ax_max, ax_min\n";

  for(auto V: VX_linspace) {
    const GG::real ay_max = ay_lim.eval("ay_max", V);
    const GG::real ay_min = ay_lim.eval("ay_min", V);
    for(auto AY: AY_linspace) {
      auto ay_clip = std::clamp(AY, ay_min, ay_max);
      GG::real ax_max_val = ax_max.eval(ay_clip, V);
      GG::real ax_min_val = ax_min.eval(ay_clip, V);
      // write the values to the file
      output_file << V << ", " << AY << ", " << ay_max << ", " << ay_min << ", "
                  << ay_clip << ", " << ax_max_val << ", " << ax_min_val << "\n";
    }
  }

  output_file.close();
  
  // 3. Define constraint functions
  auto GG_shape1 = [&ax_max](GG::real ay, GG::real v) -> GG::real {
    return ax_max.eval(ay, v);
  };
  
  auto GG_shape2 = [&ax_min](GG::real ay, GG::real v) -> GG::real {
    return ax_min.eval(ay, v);
  };
  
  auto gg_range_min = [&ay_lim](GG::real v) -> GG::real {
    return ay_lim.eval("ay_min", v);
  };
  
  auto gg_range_max = [&ay_lim](GG::real v) -> GG::real {
    return ay_lim.eval("ay_max", v);
  };
  
  GG::gg_range_max_min gg_range = {gg_range_min, gg_range_max};
  
  // 4. Create FWBW object
  GG::FWBW fbga(GG_shape1, GG_shape2, gg_range);
  
  // 5. Compute optimal velocity profile
  GG::real total_time = fbga.compute(SS_eval_vec, KK_eval_vec, v_initial);
  
  // 6. Extract results
  std::vector<GG::real> VX_eval_vec(SS_eval_vec.size());
  std::vector<GG::real> AX_eval_vec(SS_eval_vec.size());
  std::vector<GG::real> AY_eval_vec(SS_eval_vec.size());
  
  for (size_t i = 0; i < SS_eval_vec.size(); ++i) {
    VX_eval_vec[i] = fbga.evalV(SS_eval_vec[i]);   // Velocity [m/s]
    AX_eval_vec[i] = fbga.evalAx(SS_eval_vec[i]);  // Longitudinal acceleration [m/s²]
    AY_eval_vec[i] = fbga.evalAy(SS_eval_vec[i]);  // Lateral acceleration [m/s²]
  }
  
  std::cout << "Total lap time: " << total_time << " seconds" << std::endl;

  std::cout << "Results:\n";
  std::cout << "SS: ";
  for (const auto& s : SS_eval_vec) {
    std::cout << s << " ";
  }
  std::cout << "\nKK: ";
  for (const auto& k : KK_eval_vec) {
    std::cout << k << " ";
  }
  std::cout << "\nVX: ";
  for (const auto& v : VX_eval_vec) {
    std::cout << v << " ";
  }
  std::cout << "\nAX: ";
  for (const auto& a : AX_eval_vec) {
    std::cout << a << " ";
  }
  std::cout << "\nAY: ";
  for (const auto& a : AY_eval_vec) {
    std::cout << a << " ";
  }
  std::cout << std::endl;
  return 0;
}