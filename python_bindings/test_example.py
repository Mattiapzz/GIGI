#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import pygigi

def example_usage():
    """Demonstrate basic usage of the FB_F1_10 class using exact C++ test data"""
    
    print("FBGA - Forward Backward Generic Acceleration constraints")
    print(" > Running test F1/10 (Python version)")
    
    # 1. Define track data - EXACT SAME AS C++ TEST
    SS_vec = [0.0, 10, 20.0, 30.0, 40.0, 50.0]  # Arc length [m]
    KK_vec = [0.0, 0.0, 0.01, 0.005, -0.01, 0.0]  # Curvature [1/m]
    KK_vec = [k * 10.0 for k in KK_vec]  # Convert to 1/m for consistency with C++ test
    v_initial = 5.0  # Initial velocity [m/s]
    
    print(f" > Track data:")
    print(f" >   SS_vec: {len(SS_vec)} points")
    print(f" >   KK_vec: {len(KK_vec)} points")
    print(f" >   v_initial: {v_initial} m/s")
    for i in range(len(SS_vec)):
        print(f" >   SS[{i}]: {SS_vec[i]}, KK[{i}]: {KK_vec[i]}")
    
    # Interpolate track data to get evaluation points (like C++ does)
    numpt_eval = int(np.ceil(SS_vec[-1]))*10  # 40 points
    ds = SS_vec[-1] / (numpt_eval - 1)
    SS_eval_vec = [i * ds for i in range(numpt_eval)]
    
    # Simple linear interpolation for curvature (matching C++ approach)
    KK_eval_vec = []
    for s in SS_eval_vec:
        # Find interpolation (simple implementation)
        if s <= SS_vec[0]:
            k = KK_vec[0]
        elif s >= SS_vec[-1]:
            k = KK_vec[-1]
        else:
            # Linear interpolation
            for j in range(len(SS_vec)-1):
                if SS_vec[j] <= s <= SS_vec[j+1]:
                    t = (s - SS_vec[j]) / (SS_vec[j+1] - SS_vec[j])
                    k = KK_vec[j] + t * (KK_vec[j+1] - KK_vec[j])
                    break
        KK_eval_vec.append(k)
    
    # 2. Define splines data - EXACT SAME AS C++ TEST
    AY_spline = [-5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0]  # Lateral acceleration [m/s²]
    VX_spline = [0.0, 5.0, 10.0, 15.0, 20.0]  # Longitudinal velocity [m/s]
    
    # Maximum longitudinal acceleration - EXACT SAME AS C++ TEST
    AX_spline_max = [
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
        0,0,0,0,0
    ]
    
    # Minimum longitudinal acceleration - EXACT SAME AS C++ TEST
    AX_spline_min = [
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
    ]
    
    # Lateral acceleration limits - EXACT SAME AS C++ TEST
    AY_spline_max = [5.0, 4.5, 4.5, 4.0, 3.5]
    AY_spline_min = [-5.0, -5, -4, -4.0, -3.0]
    
    # 4. Create FWBW object using spline data struct (matching C++ approach)
    print("Creating FB_F1_10 instance using spline data struct...")
    spline_data = pygigi.GGVSplineData()
    spline_data.ay = AY_spline
    spline_data.vx = VX_spline
    spline_data.ax_max = AX_spline_max
    spline_data.ax_min = AX_spline_min
    spline_data.ay_max = AY_spline_max
    spline_data.ay_min = AY_spline_min
    
    fb110 = pygigi.FB_F1_10(spline_data)
    
    # 5. Compute optimal velocity profile
    print(f"Computing optimal trajectory for {len(SS_eval_vec)} points...")
    print(f"Track length: {SS_vec[-1]} m")
    print(f"Initial velocity: {v_initial} m/s")
    
    total_time = fb110.compute(SS_eval_vec, KK_eval_vec, v_initial)
    
    # 6. Extract results - EXACT SAME AS C++ TEST
    VX_eval_vec = []
    AX_eval_vec = []
    AY_eval_vec = []
    
    for i in range(len(SS_eval_vec)):
        VX_eval_vec.append(fb110.evalV(SS_eval_vec[i]))   # Velocity [m/s]
        AX_eval_vec.append(fb110.evalAx(SS_eval_vec[i]))  # Longitudinal acceleration [m/s²]
        AY_eval_vec.append(fb110.evalAy(SS_eval_vec[i]))  # Lateral acceleration [m/s²]
    
    print(f"Total lap time: {total_time} seconds")
    
    print("Results:")
    print("SS: ", end="")
    for s in SS_eval_vec:
        print(f"{s} ", end="")
    print()
    
    print("KK: ", end="")
    for k in KK_eval_vec:
        print(f"{k} ", end="")
    print()
    
    print("VX: ", end="")
    for v in VX_eval_vec:
        print(f"{v} ", end="")
    print()
    
    print("AX: ", end="")
    for a in AX_eval_vec:
        print(f"{a} ", end="")
    print()
    
    print("AY: ", end="")
    for a in AY_eval_vec:
        print(f"{a} ", end="")
    print()
    
    # Additional comparison info
    print(f"\nComparison with C++ test:")
    print(f"Number of evaluation points: {len(SS_eval_vec)}")
    print(f"Final position: {SS_eval_vec[-1]}")
    print(f"Final velocity: {VX_eval_vec[-1]}")
    print(f"Final curvature: {KK_eval_vec[-1]}")
    
    # Print some constants for verification
    print(f"\nConstants:")
    print(f"Gravity: {pygigi.GRAVITY} m/s²")
    print(f"Pi: {pygigi.PI}")
    print(f"Degrees to radians: {pygigi.DEG2RAD}")


    ## add evaluation at specific time. evalV_t evalS evalAx_t evalAy_t in linsspace 0-total_time
    print("\nEvaluating at specific time points...")
    for t in np.linspace(0, total_time, num=50):
        s = fb110.evalS(t)
        v = fb110.evalV_t(t)
        ax = fb110.evalAx_t(t)
        ay = fb110.evalAy_t(t)
        
        print(f"At t={t:.2f} s: s={s:.2f} m, v={v:.2f} m/s, ax={ax:.2f} m/s², ay={ay:.2f} m/s²")
    
    # Create plots
    print("\nCreating plots...")
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 12))
    
    # Plot 1: Velocity vs Position
    ax1.plot(SS_eval_vec, VX_eval_vec, 'b-', linewidth=2, label='Velocity')
    ax1.set_xlabel('Position s [m]')
    ax1.set_ylabel('Velocity [m/s]')
    ax1.set_title('Optimal Velocity Profile')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # Plot 2: Longitudinal Acceleration vs Position
    ax2.plot(SS_eval_vec, AX_eval_vec, 'r-', linewidth=2, label='Longitudinal Acceleration')
    ax2.set_xlabel('Position s [m]')
    ax2.set_ylabel('Longitudinal Acceleration [m/s²]')
    ax2.set_title('Longitudinal Acceleration Profile')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    ax2.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    
    # Plot 3: Lateral Acceleration vs Position
    ax3.plot(SS_eval_vec, AY_eval_vec, 'g-', linewidth=2, label='Lateral Acceleration')
    ax3.set_xlabel('Position s [m]')
    ax3.set_ylabel('Lateral Acceleration [m/s²]')
    ax3.set_title('Lateral Acceleration Profile')
    ax3.grid(True, alpha=0.3)
    ax3.legend()
    ax3.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    
    # Also plot the track curvature as a secondary axis for reference
    ax3_twin = ax3.twinx()
    ax3_twin.plot(SS_eval_vec, KK_eval_vec, 'orange', linestyle=':', alpha=0.7, label='Curvature')
    ax3_twin.set_ylabel('Curvature [1/m]', color='orange')
    ax3_twin.tick_params(axis='y', labelcolor='orange')
    ax3_twin.legend(loc='upper right')
    
    plt.tight_layout()
    plt.show()
    
    # Also create a 3D acceleration-velocity plot
    fig2 = plt.figure(figsize=(12, 8))
    ax4 = fig2.add_subplot(111, projection='3d')
    
    # Plot the trajectory in 3D space (ay, velocity, ax)
    ax4.plot(AY_eval_vec, VX_eval_vec, AX_eval_vec, 'bo-', markersize=4, linewidth=2, label='Trajectory', alpha=0.8)
    ax4.scatter(AY_eval_vec[0], VX_eval_vec[0], AX_eval_vec[0], color='green', s=100, label='Start', alpha=1.0)
    ax4.scatter(AY_eval_vec[-1], VX_eval_vec[-1], AX_eval_vec[-1], color='red', s=100, label='End', alpha=1.0)
    
    # Create mesh grid for spline data visualization
    ay_grid, vx_grid = np.meshgrid(AY_spline, VX_spline, indexing='ij')
    
    # Reshape the spline data to match the grid (11 ay points x 5 vx points = 55 elements)
    ax_max_grid = np.array(AX_spline_max).reshape(len(AY_spline), len(VX_spline))
    ax_min_grid = np.array(AX_spline_min).reshape(len(AY_spline), len(VX_spline))
    
    # Plot the spline surfaces
    ax4.plot_surface(ay_grid, vx_grid, ax_max_grid, alpha=0.3, color='red', label='Max Acceleration Surface')
    ax4.plot_surface(ay_grid, vx_grid, ax_min_grid, alpha=0.3, color='blue', label='Min Acceleration Surface')
    
    # Plot spline data points for reference
    ay_flat = ay_grid.flatten()
    vx_flat = vx_grid.flatten()
    ax_max_flat = ax_max_grid.flatten()
    ax_min_flat = ax_min_grid.flatten()
    
    # Plot max acceleration points
    ax4.scatter(ay_flat, vx_flat, ax_max_flat, color='red', alpha=0.6, s=20, label='Max Accel Points')
    # Plot min acceleration points  
    ax4.scatter(ay_flat, vx_flat, ax_min_flat, color='blue', alpha=0.6, s=20, label='Min Accel Points')
    
    ax4.set_xlabel('Lateral Acceleration ay [m/s²]')
    ax4.set_ylabel('Velocity v [m/s]')
    ax4.set_zlabel('Longitudinal Acceleration ax [m/s²]')
    ax4.set_title('3D G-G-V Diagram - Trajectory in Acceleration-Velocity Space')
    ax4.legend()
    
    # Set viewing angle for better visualization
    ax4.view_init(elev=20, azim=45)
    
    plt.tight_layout()
    plt.show()
    
    print("Plots displayed. Close the plot windows to continue.")

if __name__ == "__main__":
    try:
        example_usage()
        print("\nExample completed successfully!")
    except Exception as e:
        print(f"Error: {e}")
        print("Make sure the pygigi module is built and installed correctly.")
        import traceback
        traceback.print_exc()
