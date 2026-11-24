---
hide:
  - navigation
  - toc
---

<div class="bf-hero" markdown>

<div class="bf-hero-grid" markdown>

<div class="bf-hero-text" markdown>

# BoilFlux

<p class="bf-tagline">Open inverse heat transfer for boiling and quenching. Reconstruct the transient,
spatially resolved boiling heat flux on a quenched, moving hot plate from infrared thermography of its
back face, with a nonlinear 3D model and the adjoint conjugate-gradient method.</p>

[Get started](getting-started.md){ .md-button .md-button--primary }
[Read the theory](theory/index.md){ .md-button }
[View on GitHub](https://github.com/suresh-gopalkrishna/boilflux){ .md-button }

<div class="bf-badges"><span>MIT licence</span><span>Python ≥ 3.10</span><span>3D · nonlinear · moving plate</span><span>adjoint CGM</span><span>benchmarks with ground truth</span></div>

</div>

<div class="bf-hero-art">
<img src="assets/hero-setup.svg" alt="A water jet quenches a hot plate moving to the left. The heat flux on the quenched face is unknown; an IR camera below measures the back-face temperature, from which BoilFlux infers the flux.">
</div>

</div>

</div>

## How it works

<div class="bf-flow">
<div class="bf-step"><span class="bf-step-n">1</span><strong>IR thermography</strong><span>Back-face temperature frames of the moving plate, in the lab frame.</span></div>
<div class="bf-step"><span class="bf-step-n">2</span><strong>Plate frame</strong><span>Resample onto the plate, mask bad pixels, estimate the noise level.</span></div>
<div class="bf-step"><span class="bf-step-n">3</span><strong>Adjoint CGM</strong><span>Nonlinear 3D conduction, Tikhonov regularisation, discrepancy-principle stop.</span></div>
<div class="bf-step bf-step--hot"><span class="bf-step-n">4</span><strong>Heat flux q(x, y, t)</strong><span>Flux and temperature maps on the quenched face.</span></div>
<div class="bf-step bf-step--hot"><span class="bf-step-n">5</span><strong>Boiling regimes</strong><span>Maximum heat flux, DNB and Leidenfrost temperatures, wetting front.</span></div>
</div>

## Why BoilFlux

The heat flux on a quenched surface cannot be measured directly: water jets, vapour films and
temperatures up to 1000 °C destroy any sensor. It has to be inferred from temperatures measured
elsewhere, which is an ill-posed **inverse heat conduction problem (IHCP)**. The methods are well
established, but there has been no maintained, open, boiling-specific software for them. Every
group writes its own code, and results cannot be compared on common data.
BoilFlux is meant to close that gap.

<div class="bf-grid" markdown>

<div class="bf-card" markdown>
<span class="bf-icon">:material-cube-outline:</span>
### Nonlinear 3D model
Finite-volume conduction in a thick plate with temperature-dependent \(k, \rho, c_p\). Use
`ny = 1` for the 2D model of the 2026 paper.
</div>

<div class="bf-card" markdown>
<span class="bf-icon">:material-arrow-right-bold-outline:</span>
### Moving plates
Solved in the plate frame, so there is no advection term. The jet stagnation point moves and the
camera field of view becomes a time-dependent observation mask.
</div>

<div class="bf-card" markdown>
<span class="bf-icon">:material-function-variant:</span>
### Adjoint CGM
Exact discrete adjoint and sensitivity problems, Polak–Ribière directions, Tikhonov
regularisation and stopping by the discrepancy principle.
</div>

<div class="bf-card" markdown>
<span class="bf-icon">:material-camera-outline:</span>
### Built for IR thermography
Resamples lab-frame frames onto the plate, estimates noise robustly, masks bad pixels and models
radiation and convection losses.
</div>

<div class="bf-card" markdown>
<span class="bf-icon">:material-chart-bell-curve:</span>
### Boiling physics
Specht parabolic boiling curve with impingement and free-flow zones. Film, transition and
nucleate boiling and the rewetting front come out of the model rather than being prescribed.
</div>

<div class="bf-card" markdown>
<span class="bf-icon">:material-fire:</span>
### Regime analysis
Maps the maximum heat flux, the DNB and Leidenfrost (rewetting) temperatures, film, transition and
nucleate boiling zones, and wetting-front speeds. See the [physics guide](physics/index.md).
</div>

<div class="bf-card" markdown>
<span class="bf-icon">:material-waves:</span>
### Neural operators for boiling flow
A Fourier neural operator trained on interface-resolved boiling CFD predicts bubbles, temperature and
wall heat flux in milliseconds. See the [prototype](physics/neural-operators.md).
</div>

<div class="bf-card" markdown>
<span class="bf-icon">:material-trophy-outline:</span>
### BF-Bench
Reproducible synthetic quench cases with known ground-truth flux, generated on a refined grid to
avoid the inverse crime, with standard scores.
</div>

</div>

## A moving-plate quench reconstructed

<figure class="bf-figure" markdown>
![Centre-line heat flux, true vs recovered](assets/figures/centerline_xt.png)
<figcaption>Benchmark case <code>M3D-AA6082</code>: AA6082 plate at 500 °C moving at 10 mm/s under a
single round jet. Left: ground-truth heat flux along the centre line in the plate frame.
Right: BoilFlux reconstruction from noisy back-face IR data (σ = 2 K). The dashed line is the
moving stagnation point.</figcaption>
</figure>

## Gallery

<div class="bf-gallery" markdown>

<a class="bf-shot" href="physics/">
<img src="assets/figures/physics/boiling_curve_regimes.png" alt="Boiling curve of a jet-impingement quench with the convection, nucleate, transition and film boiling regimes marked" loading="lazy">
<strong>Boiling curve and regimes</strong>
<span>Where the maximum heat flux (DNB) and the Leidenfrost point sit on a jet-quench boiling curve.</span>
</a>

<a class="bf-shot" href="benchmark/">
<img src="assets/figures/flux_maps.png" alt="True and recovered surface heat-flux maps at four times during a moving-plate quench" loading="lazy">
<strong>True vs recovered flux maps</strong>
<span>Surface heat flux at four instants of the <code>M3D-AA6082</code> benchmark, from noisy IR data.</span>
</a>

<a class="bf-shot" href="physics/regime-analysis/">
<img src="assets/figures/regimes/m3d_feature_maps.png" alt="Surface maps of maximum heat flux, DNB temperature, rewetting time and boiling regime" loading="lazy">
<strong>Regime maps</strong>
<span>Maximum heat flux, DNB temperature, rewetting time and the regime field across the plate.</span>
</a>

<a class="bf-shot" href="physics/neural-operators/">
<img src="assets/figures/operators/rollout_snapshots.png" alt="Bubble and temperature snapshots from boiling CFD compared with a Fourier neural operator rollout" loading="lazy">
<strong>Neural operator for boiling flow</strong>
<span>A Fourier neural operator rolls out bubbles and temperature next to the CFD reference.</span>
</a>

</div>

## Try it in a minute

```bash
pip install boilflux          # or: pip install git+https://github.com/suresh-gopalkrishna/boilflux
boilflux demo                 # synthetic moving-plate quench -> inversion -> figures
boilflux bench                # run all BF-Bench cases and print the score table
```

```python
import boilflux as bf

data = bf.JetQuenchScenario(material="AA6082", T0=500.0, plate_velocity=0.01).generate()
result = bf.AdjointCGM(data).solve()
print(bf.flux_metrics(result.q[result.reliable], data.q_true[result.reliable]))
```

## Built on published science

BoilFlux implements and extends the method of

> S.B. Gopalkrishna, A.K. Nallathambi, E. Specht, *Inverse heat conduction method to estimate heat
> flux during quenching of hot metals with infrared thermography measurements*, **Applied Thermal
> Engineering** 304 (2026) 132481. [doi:10.1016/j.applthermaleng.2026.132481](https://doi.org/10.1016/j.applthermaleng.2026.132481)

The [roadmap](roadmap.md) adds neural-operator priors, a verifier-gated agent for method
selection, uncertainty quantification and an open experimental IR quench dataset.
