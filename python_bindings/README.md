# GIGI Python Bindings

This directory contains Python bindings for the GIGI optimization library using pybind11.

## Overview

The Python bindings provide access to the main GIGI functionality, specifically:

- `FB_F1_10` class for Formula 1-style vehicle dynamics optimization
- Core `FWBW` (Forward-Backward) algorithm methods
- Spline-based G-G diagram representation

## Bound Methods

### FB_F1_10 Class
- **Constructors**:
  - `FB_F1_10(ay_spline, vx_spline, ax_max, ax_min, ay_max, ay_min)` - Individual vectors
  - `FB_F1_10(spline_data)` - Using `GGVSplineData` struct

### FWBW Methods (inherited by FB_F1_10)
- `compute(SS, KK, v0)` - Main forward-backward algorithm
- `evaluate(SS, AX, AY, V)` - Batch evaluation of accelerations and velocities
- `evalV(s)` - Evaluate velocity at position s
- `evalAx(s)` - Evaluate longitudinal acceleration at position s
- `evalAy(s)` - Evaluate lateral acceleration at position s
- `get_seg_idx(s)` - Get segment index for position s
- `evalVmax(s)` - Evaluate maximum velocity at position s

### Helper Structures
- `GGVSplineData` - Container for spline data with fields: `ay`, `vx`, `ax_max`, `ax_min`, `ay_max`, `ay_min`

## Installation

### Prerequisites

1. **Install Python dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

2. **Install pybind11** (if not already installed):
   ```bash
   pip install pybind11
   ```

### Method 1: Using setup.py (Recommended)

```bash
cd python_bindings
python setup.py build_ext --inplace
```

This will create the `pygigi` module in the current directory.

### Method 2: Using CMake

```bash
cd python_bindings
mkdir build
cd build
cmake ..
make -j4
```

### Method 3: pip install (editable)

```bash
cd python_bindings
pip install -e .
```

## Usage Example

```python
import numpy as np
import pygigi

# Create spline data
spline_data = pygigi.GGVSplineData()
spline_data.ay = [0.0, 5.0, 10.0, 15.0, 20.0]  # Lateral acceleration points
spline_data.vx = [10.0, 20.0, 30.0, 40.0, 50.0]  # Velocity points
spline_data.ax_max = [8.0, 7.5, 7.0, ...]  # 2D grid: ay x vx
spline_data.ax_min = [-12.0, -11.5, -11.0, ...]  # Braking limits
spline_data.ay_max = [20.0, 19.0, 18.0, 17.0, 16.0]  # Max lateral g
spline_data.ay_min = [-20.0, -19.0, -18.0, -17.0, -16.0]  # Min lateral g

# Create optimizer
fb_solver = pygigi.FB_F1_10(spline_data)

# Define track
s_points = np.linspace(0, 1000, 100)  # Track positions
curvatures = 0.01 * np.sin(2 * np.pi * s_points / 1000)  # Track curvature
v0 = 20.0  # Initial velocity

# Optimize trajectory
total_time = fb_solver.compute(s_points.tolist(), curvatures.tolist(), v0)

# Evaluate results
velocity_at_500m = fb_solver.evalV(500.0)
ax_at_500m = fb_solver.evalAx(500.0)
ay_at_500m = fb_solver.evalAy(500.0)

print(f"Total time: {total_time:.3f} s")
print(f"Velocity at 500m: {velocity_at_500m:.3f} m/s")
```

## Testing

Run the example script:

```bash
python test_example.py
```

This will run the optimization using the exact same data as the C++ test and display plots showing:
- Velocity profile along the track
- Longitudinal acceleration profile  
- Lateral acceleration profile with track curvature overlay
- G-G diagram showing the acceleration trajectory

## Constants

The following constants are available:
- `pygigi.GRAVITY` - Gravitational acceleration (9.81 m/s²)
- `pygigi.PI` - Pi constant
- `pygigi.DEG2RAD` - Degrees to radians conversion
- `pygigi.RAD2DEG` - Radians to degrees conversion

## Data Types

- `real` (C++) maps to Python `float`
- `integer` (C++) maps to Python `int`
- `std::vector<real>` (C++) maps to Python `list[float]`

## Notes

1. The spline data represents a G-G diagram with:
   - `ay`, `vx`: Grid points for lateral acceleration and velocity
   - `ax_max`, `ax_min`: 2D arrays (flattened) of acceleration limits
   - `ay_max`, `ay_min`: 1D arrays of lateral acceleration limits

2. The `ax_max` and `ax_min` arrays should be flattened 2D arrays where:
   - Size = `len(ay_spline) * len(vx_spline)`
   - Indexing: `ax_max[i*len(vx_spline) + j]` for `ay[i]`, `vx[j]`

3. For trajectory optimization:
   - `SS`: Vector of track positions (arclength)
   - `KK`: Vector of track curvatures at each position
   - `v0`: Initial velocity

## Troubleshooting

1. **Import errors**: Make sure the module is built and in your Python path
2. **Build errors**: Check that all GIGI dependencies are available
3. **Runtime errors**: Verify that spline data dimensions are consistent

## Files

- `pybind_gigi.cpp` - Main binding code
- `CMakeLists.txt` - CMake build configuration
- `setup.py` - Python setuptools configuration
- `requirements.txt` - Python dependencies
- `test_example.py` - Usage example and test
- `README.md` - This file
