# Data and IR cameras

## The dataset format

A `BoilDataset` stores IR-face temperatures on the plate-frame grid:

| Field | Shape / type | Meaning |
|---|---|---|
| `grid` | `PlateGrid` | plate geometry and resolution |
| `material` | str | material name (built-in or label) |
| `dt` | float | frame interval [s]; frame *n* is at \(t_n = n\,\Delta t\), n = 1..nt |
| `T0` | float | initial uniform plate temperature [°C] |
| `Y` | `(nt, ny, nx)` | IR-face temperatures [°C] |
| `mask` | `(nt, ny, nx)` bool | valid observations (field of view, bad pixels) |
| `sigma` | float or None | noise standard deviation [K]; estimated from the data if None |
| `q_true` | `(nt, ny, nx)` or None | ground truth (synthetic data only) [W/m²] |
| `motion` | `PlateMotion` | plate velocity and initial offset |
| `meta` | dict | free-form metadata (e.g. `x_jet_lab`) |

Datasets are saved as compressed `.npz` files with a JSON header (`dataset.save(path)`,
`BoilDataset.load(path)`), so they can be read without BoilFlux.

## From an IR camera

1. **Export temperatures**, not raw counts, using your emissivity calibration (for example
   \(\varepsilon = 0.90\)–\(0.93\) for high-emissivity black paint). Save the frames as a NumPy array
   `(n_frames, rows, cols)`.
2. **Spatial calibration**: give the lab coordinates of the pixel columns (`X_px`) and rows
   (`Y_px`) in metres, from a calibration target or known plate dimensions.
3. **Time origin**: \(t=0\) is the start of quenching (first water contact). Drop frames before it and
   use the frame just before contact for `T0`.
4. **Kinematics**: plate speed and the lab position of the plate edge at \(t=0\) (`PlateMotion`).
5. **Resample**: `Y, mask = bf.lab_frames_to_plate(frames, X_px, Y_px, t, grid, motion)`.
6. **Mask bad pixels**, reflections or the nozzle shadow by setting `mask[...] = False`.

A complete template is in `examples/02_experimental_ir_template.py`.

## Choosing the grid

- `nz`: 5–8 cells through the thickness is usually enough. The unknown sits on one face and the data on the other.
- `nx, ny`: about the IR pixel pitch or coarser. Features smaller than roughly the plate thickness
  cannot be resolved from the back face, so a cell size of \(0.5\)–\(1\times L_z\) is a good default.
- `dt`: the camera frame interval, or an integer multiple of it after temporal averaging.

## Contributing experimental data

Experimental IR quench records with good metadata are the most valuable contribution to this
project, because no public dataset of this kind exists yet. Please include: material and its
property source, plate dimensions, initial temperature, jet diameter, velocity, flow rate and
water temperature, plate speed, camera model, frame rate, pixel pitch, emissivity calibration, and
any independent validation (thermocouples, energy balance). See [Contributing](../community/contributing.md).
