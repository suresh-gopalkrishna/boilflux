"""Template: invert your own IR recording of a moving-plate quench.

Replace the loading section with your camera export. Frames must already be
converted to temperature with your emissivity calibration.
"""

import numpy as np

import boilflux as bf
from boilflux.plotting import save_all

# ---- 1. load the IR sequence (lab frame) ------------------------------------
frames = np.load("my_ir_frames.npy")            # (n_frames, n_py, n_px) in degC
fps = 200.0
t = np.arange(1, frames.shape[0] + 1) / fps      # time of each frame after quench start [s]
pixel = 0.7e-3                                   # pixel pitch on the plate [m]
X_px = np.arange(frames.shape[2]) * pixel        # lab x of each pixel column [m]
Y_px = np.arange(frames.shape[1]) * pixel        # lab y of each pixel row [m]

# ---- 2. geometry, kinematics, material --------------------------------------
grid = bf.PlateGrid(Lx=0.250, Ly=0.230, Lz=0.005, nx=125, ny=115, nz=6)
motion = bf.PlateMotion(velocity=0.010, x0=-0.050)  # plate edge at lab x = -50 mm at t = 0
material = bf.get_material("AA6082")                # or bf.Material.from_csv("my_props.csv")

# ---- 3. map the frames onto the plate and build the dataset -----------------
Y, mask = bf.lab_frames_to_plate(frames, X_px, Y_px, t, grid, motion)
data = bf.BoilDataset(grid=grid, material=material.name, dt=1 / fps, T0=float(frames[0].max()),
                        Y=Y, mask=mask, motion=motion,
                        meta={"x_jet_lab": 0.075, "experiment": "describe your test here"})
print("estimated IR noise sigma [K]:", bf.estimate_noise_sigma(Y, mask))

# ---- 4. invert --------------------------------------------------------------
loss = bf.SurfaceLoss(h=10.0, emissivity=0.92, T_amb=25.0)
result = bf.AdjointCGM(data, material=material, top_loss=loss,
                       options=bf.CGMOptions(max_iter=300, tau=1.05)).solve()
result.save("my_result.npz")
save_all(data, result, "my_figures")
