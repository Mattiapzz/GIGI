#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>

#include "GIGI/FB_F1_10.hxx"
#include "GIGI/FWBW.hxx"
#include "GIGI/types.hxx"

namespace py = pybind11;
using namespace GG;

PYBIND11_MODULE(pygigi, m) {
    m.doc() = "Python bindings for GIGI library";

    // Bind the gg_range_max_min struct
    py::class_<gg_range_max_min>(m, "GGRangeMaxMin")
        .def(py::init<>())
        .def_readwrite("min", &gg_range_max_min::min, "Minimum function")
        .def_readwrite("max", &gg_range_max_min::max, "Maximum function");

    // Bind the ggv_spline_data struct
    py::class_<ggv_spline_data>(m, "GGVSplineData")
        .def(py::init<>())
        .def_readwrite("ay", &ggv_spline_data::ay, "Lateral acceleration [m/s²]")
        .def_readwrite("vx", &ggv_spline_data::vx, "Longitudinal velocity [m/s]")
        .def_readwrite("ax_max", &ggv_spline_data::ax_max, "Maximum longitudinal acceleration [m/s²]")
        .def_readwrite("ax_min", &ggv_spline_data::ax_min, "Minimum longitudinal acceleration [m/s²]")
        .def_readwrite("ay_max", &ggv_spline_data::ay_max, "Maximum lateral acceleration [m/s²]")
        .def_readwrite("ay_min", &ggv_spline_data::ay_min, "Minimum lateral acceleration [m/s²]");

    // Bind the FWBW base class (essential methods only)
    py::class_<FWBW>(m, "FWBW")
        .def(py::init<const std::function<real(real, real)>&,
                      const std::function<real(real, real)>&,
                      const gg_range_max_min&>(),
             "Constructor with function pointers",
             py::arg("gg_Upper"), py::arg("gg_Lower"), py::arg("gg_range"))
        .def("compute", &FWBW::compute, 
             "Compute the forward-backward algorithm",
             py::arg("SS"), py::arg("KK"), py::arg("v0"))
        .def("evaluate", &FWBW::evaluate,
             "Evaluate acceleration and velocity at given positions",
             py::arg("SS"), py::arg("AX"), py::arg("AY"), py::arg("V"))
        .def("evalV", &FWBW::evalV,
             "Evaluate velocity at position s",
             py::arg("s"))
        .def("evalAx", &FWBW::evalAx,
             "Evaluate longitudinal acceleration at position s",
             py::arg("s"))
        .def("evalAy", &FWBW::evalAy,
             "Evaluate lateral acceleration at position s", 
             py::arg("s"))
        .def("get_seg_idx", &FWBW::get_seg_idx,
             "Get segment index for position s",
             py::arg("s"))
        .def("evalVmax", &FWBW::evalVmax,
             "Evaluate maximum velocity at position s",
             py::arg("s"));

    // Bind the FB_F1_10 class
    py::class_<FB_F1_10, FWBW>(m, "FB_F1_10")
        .def(py::init<const std::vector<real>&,
                      const std::vector<real>&,
                      const std::vector<real>&,
                      const std::vector<real>&,
                      const std::vector<real>&,
                      const std::vector<real>&>(),
             "Constructor with individual spline vectors",
             py::arg("ay_spline"),
             py::arg("vx_spline"),
             py::arg("ax_bispline_max"),
             py::arg("ax_bispline_min"),
             py::arg("ay_spline_max"),
             py::arg("ay_spline_min"))
        .def(py::init<const ggv_spline_data&>(),
             "Constructor with spline data struct",
             py::arg("spline_data"));

    // Add some useful constants
    m.attr("GRAVITY") = GRAVITY;
    m.attr("PI") = PI;
    m.attr("DEG2RAD") = DEG2RAD;
    m.attr("RAD2DEG") = RAD2DEG;
}
