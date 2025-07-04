#include <chrono>
#include <iostream>

#include "GIGI/gg_utils.hxx"
#include "GIGI/types.hxx"

#include "GIGI/FWBW.hxx"

// rapidcsv
#include "rapidcsv.h"

using Colour = struct Colour
{
  uint8_t r{0};
  uint8_t g{0};
  uint8_t b{0};
};


int main()
{

  const std::string road_file_path = "./data/circuit/road_Mugello_3D.txt";

  rapidcsv::Document doc(
    road_file_path, 
    rapidcsv::LabelParams(0, -1),
                          rapidcsv::SeparatorParams('\t')
  );

  std::vector<GG::real> abscissa  = doc.GetColumn<GG::real>("abscissa");

  GG::LinearInterpolatorSet RD;
  RD.setup(
    doc.GetColumn<GG::real>("abscissa"),
    {
      doc.GetColumn<GG::real>("X"),
      doc.GetColumn<GG::real>("Y"),
      doc.GetColumn<GG::real>("elevation"),
      doc.GetColumn<GG::real>("theta"),
      doc.GetColumn<GG::real>("curvature"),
      doc.GetColumn<GG::real>("banking"),
      doc.GetColumn<GG::real>("slope"),
      doc.GetColumn<GG::real>("upsilon"),
      doc.GetColumn<GG::real>("torsion"),
      doc.GetColumn<GG::real>("width_L"),
      doc.GetColumn<GG::real>("width_R")
    },
    {
      "X", "Y", "elevation", "theta", "curvature", "banking", "slope", "upsilon", "torsion",
      "width_L", "width_R"
    }
  );

  
  auto GG_shapeUP = [](GG::real ay, GG::real v) -> GG::real { return +20.0*sqrt(20.0*20.0 - ay*ay)*std::max(-0.01*(v-100.0),0.0)/20.0; };
  auto GG_shapeDW = [](GG::real ay, GG::real v) -> GG::real { return -20.0*sqrt(20.0*20.0 - ay*ay)*std::max(+0.0001*v*v+1.0,0.0)/20.0; };
  
  auto gg_range_min = [](GG::real v) -> GG::real { return -18.0+0*v; };
  auto gg_range_max = [](GG::real v) -> GG::real { return +18.0+0*v; };
  
  GG::gg_range_max_min gg_range = { gg_range_min, gg_range_max };
  
  GG::FWBW fwbw(GG_shapeUP, GG_shapeDW, gg_range);
  
  
  
  // for loop to extract 1 point for each metre.
  GG::integer numpoints = ceil(abscissa.back()-abscissa.front());
 std::vector<GG::real> SS_gg_vec(numpoints);
  std::vector<GG::real> KK_gg_vec(numpoints);
  for (GG::integer i = 0; i < numpoints; i++)
  {
    GG::real s = abscissa.front() + static_cast<GG::real>(i);
    SS_gg_vec[i] = s;
    KK_gg_vec[i] = RD.eval("curvature", s);
  }

  GG::real v_initial = 20.0;

  auto start = std::chrono::high_resolution_clock::now();
  GG::real T = fwbw.compute(SS_gg_vec, KK_gg_vec, v_initial);
  auto end = std::chrono::high_resolution_clock::now();

  auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  std::cout << "Elapsed time:     " << elapsed.count() / 1000.0 << " milliseconds\n";
  std::cout << "Total lap time:   " << T << "\n";
  std::cout << "Num segments:     " << numpoints << "\n";
  std::cout << "Avg time/segment: " << (GG::real) elapsed.count() / (GG::real) numpoints << " microseconds\n";


  return 0;
}


