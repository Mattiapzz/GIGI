cmake_minimum_required(VERSION 3.14)

# Set up the project
project(pygigi)

# Find required packages
find_package(pybind11 REQUIRED)

# Include the main GIGI library directories
include_directories(../src)

# Get the GIGI sources
set(GIGI_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../src")
file(GLOB GIGI_SOURCES ${GIGI_SOURCE_DIR}/*.cc)

# Create the Python module
pybind11_add_module(pygigi 
    pybind_gigi.cpp
    ${GIGI_SOURCES}
)

# Set C++ standard
target_compile_features(pygigi PRIVATE cxx_std_17)

# Link against the third party libraries
set(THIRD_PARTY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../third_party)
target_include_directories(pygigi PRIVATE
    ${THIRD_PARTY_DIR}/cxxopts/include
    ${THIRD_PARTY_DIR}/json/single_include/nlohmann
    ${THIRD_PARTY_DIR}/rapidcsv/src
)

# Set properties for the module
set_target_properties(pygigi PROPERTIES
    CXX_VISIBILITY_PRESET "hidden"
    VISIBILITY_INLINES_HIDDEN YES
)
