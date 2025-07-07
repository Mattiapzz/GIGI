#include <chrono>
#include <fstream>
#include <iostream>

#include "GIGI/FWBW.hxx"
#include "GIGI/gg_utils.hxx"
#include "GIGI/types.hxx"

// rapidcsv
#include "rapidcsv.h"

#include "json.hpp"

using Colour = struct Colour {
  uint8_t r{0};
  uint8_t g{0};
  uint8_t b{0};
};


#include "cxxopts.hpp"


struct AeroData{
  GG::real b = 0.73; // wheelbase
  GG::real L__W = 1.5; // length of the wheel
  GG::real h = 0.69; // height of the center of mass
  GG::real mu__X = 1.30; // longitudinal friction coefficient
  GG::real mu__Y = 1.40; // lateral friction coefficient
  GG::real c__a__0 = 0.00; // aerodynamic coefficient
  GG::real c__a__1 = 0.00; // aerodynamic coefficient
  GG::real c__a__2 = 0.5*1.2*0.25/(250.0*9.81); // aerodynamic coefficient
  GG::real h__a = 0.51; // height of the aerodynamic center
  GG::real g = 9.81; // gravitational acceleration
};

GG::real ax_max_engine(GG::real v, const AeroData& aero_data) {
  const auto b       = aero_data.b       ;
  const auto L__W    = aero_data.L__W    ;
  const auto h       = aero_data.h       ; 
  const auto mu__X   = aero_data.mu__X   ;
  const auto mu__Y   = aero_data.mu__Y   ;
  const auto c__a__0 = aero_data.c__a__0 ; 
  const auto c__a__1 = aero_data.c__a__1 ; 
  const auto c__a__2 = aero_data.c__a__2 ; 
  const auto h__a    = aero_data.h__a    ;
  const auto g       = aero_data.g       ;

  return  145.0 * 1000.0 / (250 * v) - (c__a__2*v*v+c__a__1*v+c__a__0)*g;
}

GG::real ax_min_adherence(GG::real ayg, GG::real v, const AeroData& aero_data) {
  const auto b       = aero_data.b       ;
  const auto L__W    = aero_data.L__W    ;
  const auto h       = aero_data.h       ; 
  const auto mu__X   = aero_data.mu__X   ;
  const auto mu__Y   = aero_data.mu__Y   ;
  const auto c__a__0 = aero_data.c__a__0 ; 
  const auto c__a__1 = aero_data.c__a__1 ; 
  const auto c__a__2 = aero_data.c__a__2 ; 
  const auto h__a    = aero_data.h__a    ;
  const auto g       = aero_data.g       ;
  //
  const auto ay = ayg / g;
  //
  return g*((-c__a__2 * mu__Y * v *v - c__a__1 * mu__Y * v - c__a__0 * mu__Y - sqrt(-std::pow(ay * mu__X ,2)+ std::pow(mu__Y * mu__X ,2))) / mu__Y);
}

GG::real ax_min_stoppie(GG::real ayg, GG::real v, const AeroData& aero_data) {
  const auto b       = aero_data.b       ;
  const auto L__W    = aero_data.L__W    ;
  const auto h       = aero_data.h       ; 
  const auto mu__X   = aero_data.mu__X   ;
  const auto mu__Y   = aero_data.mu__Y   ;
  const auto c__a__0 = aero_data.c__a__0 ; 
  const auto c__a__1 = aero_data.c__a__1 ; 
  const auto c__a__2 = aero_data.c__a__2 ; 
  const auto h__a    = aero_data.h__a    ;
  const auto g       = aero_data.g       ;
  //
  const auto ay = ayg / g;
  //
  return g*(-(-c__a__2 * h__a * v * v + c__a__2 * h * v * v - c__a__1 * h__a * v + c__a__1 * h * v + L__W * sqrt((ay * ay + 1)) - b * sqrt((ay * ay + 1)) - c__a__0 * h__a + c__a__0 * h) / h);
}

GG::real ax_max_wheeling(GG::real ayg, GG::real v, const AeroData& aero_data) {
  const auto b       = aero_data.b       ;
  const auto L__W    = aero_data.L__W    ;
  const auto h       = aero_data.h       ; 
  const auto mu__X   = aero_data.mu__X   ;
  const auto mu__Y   = aero_data.mu__Y   ;
  const auto c__a__0 = aero_data.c__a__0 ; 
  const auto c__a__1 = aero_data.c__a__1 ; 
  const auto c__a__2 = aero_data.c__a__2 ; 
  const auto h__a    = aero_data.h__a    ;
  const auto g       = aero_data.g       ;
  //
  const auto ay = ayg / g;
  //
  return g*(c__a__2 * h__a * v * v - c__a__2 * h * v * v + c__a__1 * h__a * v - c__a__1 * h * v + b * sqrt((ay * ay + 1)) + c__a__0 * h__a - c__a__0 * h) / h;
}

GG::real ax_max_adherence(GG::real ayg, GG::real v, const AeroData& aero_data) {
  const auto b       = aero_data.b       ;
  const auto L__W    = aero_data.L__W    ;
  const auto h       = aero_data.h       ; 
  const auto mu__X   = aero_data.mu__X   ;
  const auto mu__Y   = aero_data.mu__Y   ;
  const auto c__a__0 = aero_data.c__a__0 ; 
  const auto c__a__1 = aero_data.c__a__1 ; 
  const auto c__a__2 = aero_data.c__a__2 ; 
  const auto h__a    = aero_data.h__a    ;
  const auto g       = aero_data.g       ;
  //
  const auto ay = ayg / g;
  //
  return g*(h *std::pow(mu__X,2) * (ay - mu__Y) * (ay + mu__Y) * (b - L__W) * sqrt((std::pow(ay,2) + 1)) - (std::pow(L__W,2) * (std::pow(mu__Y,2)) + h *std::pow(mu__X,2) * (h - h__a)) * (c__a__2 * std::pow(v,2) + c__a__1 * v + c__a__0) * (std::pow(ay,2)) + (-std::pow(L__W,2) + h *std::pow(mu__X,2) * (h - h__a)) * (c__a__2 * std::pow(v,2) + c__a__1 * v + c__a__0) * (std::pow(mu__Y,2)) + sqrt(-(ay + mu__Y) * std::pow(L__W,2) * (std::pow(ay,2) + 1) * (0.2e1 * h__a * (c__a__2 * std::pow(v,2) + c__a__1 * v + c__a__0) * (b - L__W) * sqrt((std::pow(ay,2) + 1)) + std::pow(b - L__W,2) * (std::pow(ay,2)) + b * b - 0.2e1 * L__W * b + std::pow(L__W,2) + std::pow(c__a__2 * std::pow(v,2) + c__a__1 * v + c__a__0,2) * h__a * h__a) *std::pow(mu__X,2) * (ay - mu__Y) * (std::pow(mu__Y,2)))) / ((std::pow(L__W,2) * (std::pow(mu__Y,2)) + std::pow(h,2) *std::pow(mu__X,2)) * (std::pow(ay,2)) + (-std::pow(h,2) *std::pow(mu__X,2) + std::pow(L__W,2)) * (std::pow(mu__Y,2)));
}




int main(int argc, char *argv[])
{
  // Parse command line arguments
  cxxopts::Options options("GIGI");
  options.add_options()
    ("h,help", "Print help")
    ("c,circuit", "Circuit name", cxxopts::value<std::string>()->default_value("Catalunya"));
  
  auto result = options.parse(argc, argv); 
  // Print help
  if (result.count("help"))
  {
    std::cout << options.help() << std::endl;
    exit(0);
  }


  // ██████╗  ██████╗  █████╗ ██████╗
  // ██╔══██╗██╔═══██╗██╔══██╗██╔══██╗
  // ██████╔╝██║   ██║███████║██║  ██║
  // ██╔══██╗██║   ██║██╔══██║██║  ██║
  // ██║  ██║╚██████╔╝██║  ██║██████╔╝
  // ╚═╝  ╚═╝ ╚═════╝ ╚═╝  ╚═╝╚═════╝

  const std::string circuit_name = result["circuit"].as<std::string>();

  const std::string sol_file_path = "./data/paper/ard_2d_moto_" + circuit_name + ".txt";
  const std::string comp_file_path = "./data/paper/ard_2d_fwbw_moto_" + circuit_name + ".txt";


  rapidcsv::Document doc(
    sol_file_path, 
    rapidcsv::LabelParams(0, -1), // Ensure the first non-comment line is treated as headers
    rapidcsv::SeparatorParams('\t',true),
    rapidcsv::ConverterParams(),
    rapidcsv::LineReaderParams(true, '#')); // Skip only comment lines 

    // extract zeta, curv_trans, v__x, a__x, a__y

  std::vector<GG::real> SS_vec_MLT = doc.GetColumn<GG::real>("zeta");
  std::vector<GG::real> KK_vec_MLT = doc.GetColumn<GG::real>("curv_trans");
  std::vector<GG::real> VX_vec_MLT = doc.GetColumn<GG::real>("v__x");
  std::vector<GG::real> VY_vec_MLT = doc.GetColumn<GG::real>("v__y");
  std::vector<GG::real> AX_vec_MLT = doc.GetColumn<GG::real>("a__x");
  std::vector<GG::real> AY_vec_MLT = doc.GetColumn<GG::real>("ay_2d");

  GG::LinearInterpolatorSet traj_spline_MLT(
    SS_vec_MLT, 
    { KK_vec_MLT, VX_vec_MLT, VY_vec_MLT, AX_vec_MLT, AY_vec_MLT},
    {"kappa", "v_x", "v_y", "a_x", "a_y"});
  

  GG::real mesh_size_MLT = 0.0;
  for (size_t i = 0; i < SS_vec_MLT.size() - 1; ++i) {
    mesh_size_MLT += (SS_vec_MLT[i + 1] - SS_vec_MLT[i]) / static_cast<GG::real>(SS_vec_MLT.size() - 1);
  }


  rapidcsv::Document doc_1D(
    comp_file_path, 
    rapidcsv::LabelParams(0, -1), // Ensure the first non-comment line is treated as headers
    rapidcsv::SeparatorParams('\t',true),
    rapidcsv::ConverterParams(),
    rapidcsv::LineReaderParams(true, '#')); // Skip only comment lines 

  // extract
  std::vector<GG::real> SS_vec_1D = doc_1D.GetColumn<GG::real>("zeta");
  std::vector<GG::real> KK_vec_1D = doc_1D.GetColumn<GG::real>("curvature");
  std::vector<GG::real> VX_vec_1D = doc_1D.GetColumn<GG::real>("v__x");
  std::vector<GG::real> AX_vec_1D = doc_1D.GetColumn<GG::real>("a__x");
  // ay is curvature times v^2
  std::vector<GG::real> AY_vec_1D;
  AY_vec_1D.resize(VX_vec_1D.size());
  for (size_t i = 0; i < VX_vec_1D.size(); ++i) {
    AY_vec_1D[i] = KK_vec_1D[i] * VX_vec_1D[i] * VX_vec_1D[i];
  }

  GG::LinearInterpolatorSet traj_spline_1D(
    SS_vec_1D, 
    { KK_vec_1D, VX_vec_1D, AX_vec_1D, AY_vec_1D},
    {"kappa", "v_x", "a_x", "a_y"});
  GG::real mesh_size_1D = 0.0;
  for (size_t i = 0; i < SS_vec_1D.size() - 1; ++i) {
    mesh_size_1D += (SS_vec_1D[i + 1] - SS_vec_1D[i]) / static_cast<GG::real>(SS_vec_1D.size() - 1);
  }



  // ███████╗███╗   ██╗██╗   ██╗███████╗██╗      ██████╗ ██████╗ ███████╗
  // ██╔════╝████╗  ██║██║   ██║██╔════╝██║     ██╔═══██╗██╔══██╗██╔════╝
  // █████╗  ██╔██╗ ██║██║   ██║█████╗  ██║     ██║   ██║██████╔╝█████╗
  // ██╔══╝  ██║╚██╗██║╚██╗ ██╔╝██╔══╝  ██║     ██║   ██║██╔═══╝ ██╔══╝
  // ███████╗██║ ╚████║ ╚████╔╝ ███████╗███████╗╚██████╔╝██║     ███████╗
  // ╚══════╝╚═╝  ╚═══╝  ╚═══╝  ╚══════╝╚══════╝ ╚═════╝ ╚═╝     ╚══════╝


  AeroData aero_data;

  auto GG_shape1 = [&aero_data](GG::real ay, GG::real v) -> GG::real {
    return std::min(
      std::min(
        ax_max_wheeling(ay, v, aero_data),
        ax_max_adherence(ay, v, aero_data)),
      ax_max_engine(v, aero_data));
  };
  auto GG_shape2 = [&aero_data](GG::real ay, GG::real v) -> GG::real {
    return std::max(
      ax_min_stoppie(ay, v, aero_data),
      ax_min_adherence(ay, v, aero_data));
  };
  auto gg_range_min = [&aero_data](GG::real v) -> GG::real {
    return -aero_data.mu__Y * aero_data.g * 0.999;
  };
  auto gg_range_max = [&aero_data](GG::real v) -> GG::real {
    return +aero_data.mu__Y * aero_data.g * 0.999;
  };

  GG::gg_range_max_min gg_range = {gg_range_min, gg_range_max};

  GG::real Vminenv = 5.0;
  GG::real Vmaxenv = 100.0;

                                    
// ███████╗██╗    ██╗██████╗ ██╗    ██╗
// ██╔════╝██║    ██║██╔══██╗██║    ██║
// █████╗  ██║ █╗ ██║██████╔╝██║ █╗ ██║
// ██╔══╝  ██║███╗██║██╔══██╗██║███╗██║
// ██║     ╚███╔███╔╝██████╔╝╚███╔███╔╝
// ╚═╝      ╚══╝╚══╝ ╚═════╝  ╚══╝╚══╝ 
                                            

  GG::FWBW fwbw(GG_shape1, GG_shape2, gg_range);

  std::vector<GG::real> SS_gg_vec = SS_vec_MLT;
  std::vector<GG::real> KK_gg_vec = KK_vec_MLT;
  GG::integer numpts_gg = SS_gg_vec.size();

  GG::real v_initial = VX_vec_MLT.front();

  auto start = std::chrono::steady_clock::now();
  GG::real T = fwbw.compute(SS_gg_vec, KK_gg_vec, v_initial);
  auto end = std::chrono::steady_clock::now();

  auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  auto elapsed_microsec = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  std::cout << "CPU time: " << elapsed_microsec.count() / 1000.0 << " milliseconds\n";
  std::cout << "Total lap time: " << T << "\n";
  
  std::cout << "Execution time: " << elapsed.count() << " miliseconds" << "\n";
  std::cout << "Num segments: " << numpts_gg << "\n";
  std::cout << "Average time per segment: " << (double) elapsed_microsec.count() / (double) numpts_gg << " microseconds\n";

  // get the solution of the forward backward in terms of VX, AX, AY
  std::vector<GG::real> VX_gg_vec(numpts_gg);
  std::vector<GG::real> AX_gg_vec(numpts_gg);
  std::vector<GG::real> AY_gg_vec(numpts_gg);
  for (GG::integer i = 0; i < numpts_gg; ++i) {
    VX_gg_vec[i] = fwbw.evalV(SS_gg_vec[i]);
    AX_gg_vec[i] = fwbw.evalAx(SS_gg_vec[i]);
    AY_gg_vec[i] = fwbw.evalAy(SS_gg_vec[i]);
  }

  // write an output csv file with the solution
  std::ofstream output_file("./data/output_" + circuit_name + "_fwbw_moto.txt");
  if (output_file.is_open()) {
    output_file << "# circuit_name = " << circuit_name << "\n";
    output_file << "# v_initial = " << v_initial << "\n";
    output_file << "# FWBW solution\n";
    output_file << "# \n";
    output_file << "# cpu_time = " << elapsed_microsec.count() / 1000.0 << "ms \n";
    output_file << "# iterations = " << numpts_gg << "\n";
    output_file << "# time\n";
    output_file << "# headers\n";
    output_file << "zeta\tcurvature\tv_x\ta_x\ta_y\n";
    for (GG::integer i = 0; i < numpts_gg; ++i) {
      output_file << SS_gg_vec[i] << "\t" 
                  << KK_gg_vec[i] << "\t" 
                  << VX_gg_vec[i] << "\t" 
                  << AX_gg_vec[i] << "\t" 
                  << AY_gg_vec[i] << "\n";
    }
    output_file.close();
  } else {
    std::cerr << "Unable to open file for writing: ./data/FWBW/" + circuit_name + "_fwbw.txt\n";
  }
  

  return 0;
}
