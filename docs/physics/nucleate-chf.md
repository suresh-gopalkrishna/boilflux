# Nucleate boiling, CHF and DNB

## Onset of nucleate boiling (ONB)

A vapour embryo trapped in a surface cavity of mouth radius \(r_c\) grows only if the surrounding
liquid is superheated enough to overcome the excess pressure from surface tension (Young–Laplace,
\(\Delta p = 2\sigma/r\)). Combined with the Clausius–Clapeyron relation, the wall superheat required is
approximately

\[
\Delta T_\text{sat} \gtrsim \frac{2\sigma\,T_\text{sat}}{\rho_v h_{fg}\,r_c}.
\]

Hsu's criterion adds the condition that the thermal boundary layer must be thick enough to superheat the
bubble cap. Only cavities in a certain size band are active, and that band widens as superheat increases.
In a jet quench the incoming liquid is **subcooled** and the forced-convection boundary layer is thin,
so incipience is delayed. Subcooled nucleate boiling, where bubbles collapse near the wall, dominates
near the stagnation point.

## Mechanisms of nucleate-boiling heat transfer

High nucleate-boiling heat fluxes come from several mechanisms acting together during each bubble cycle
(nucleation, growth, departure, waiting time):

1. **Microlayer evaporation.** A thin liquid film trapped under a rapidly growing bubble evaporates
   with a very high local flux.
2. **Contact-line evaporation.** Intense evaporation in the nanometre-to-micrometre region where the
   liquid–vapour interface meets the wall.
3. **Transient conduction ("quenching") on departure.** Cold liquid rushes in to fill the dry spot
   left by a departing bubble, and the wall conducts heat into it as into a semi-infinite medium.
4. **Micro-convection and Marangoni flow.** Bubble motion stirs the boundary layer, and surface-tension
   gradients drive flow along the interface in subcooled liquid.

Engineering CFD condenses this into **heat-flux partitioning**, the RPI model of Kurul and Podowski:

\[
q_w = q_\text{conv}\,(1 - A_b) + q_\text{quench}\,A_b + q_\text{evap},
\qquad q_\text{evap} = N_a\, f\, \tfrac{\pi}{6} D_d^3\, \rho_v h_{fg},
\]

with closures for the active nucleation-site density \(N_a\), bubble departure diameter \(D_d\),
departure frequency \(f\), and the fraction \(A_b\) of the wall influenced by bubbles. This is the model
behind OpenFOAM's `wallBoiling` fvModel (see [Modelling](modelling.md)).

The classical pool-boiling correlation of **Rohsenow** (1952) relates heat flux to superheat:

\[
q = \mu_l h_{fg}\left[\frac{g(\rho_l-\rho_v)}{\sigma}\right]^{1/2}
\left(\frac{c_{p,l}\,\Delta T_\text{sat}}{C_{sf}\,h_{fg}\,\Pr_l^{\,n}}\right)^{3}.
\]

The cubic dependence on superheat explains the steep rise of the nucleate branch. The surface–fluid
constant \(C_{sf}\) contains everything the model does not resolve: roughness, wettability, oxide layers.
This is one reason measured quench curves differ from one specimen to another.

## Critical heat flux and departure from nucleate boiling

As the superheat rises, bubbles merge into vapour columns and mushrooms. Beyond a maximum, vapour
removal limits liquid supply to the wall and the flux falls. Several mechanisms have been proposed. They
are not mutually exclusive and dominate under different conditions.

### Hydrodynamic instability (Zuber, Kutateladze)

Counter-flowing vapour jets and liquid become Kelvin–Helmholtz unstable. The jet spacing is set by the
Rayleigh–Taylor wavelength of the vapour–liquid interface. This gives

\[
q_\text{CHF} = C\,\rho_v h_{fg}\left[\frac{\sigma g(\rho_l-\rho_v)}{\rho_v^2}\right]^{1/4},
\qquad C = \pi/24 \approx 0.131 \text{ (Zuber)},\; 0.149 \text{ (Lienhard–Dhir, large plates)}.
\]

For saturated water at 1 atm this gives **1.1–1.3 MW/m²**. This is the pool-boiling CHF.

### Macrolayer dryout and dry-spot theories

Haramura and Katto proposed that a liquid **macrolayer** under the large vapour mushrooms evaporates
before the mushroom departs. Dry-spot (hot-spot) theories, supported by modern high-speed IR
measurements of the wall, attribute the crisis to **irreversible growth of dry patches**: when the
wall under a dry spot heats beyond the rewetting temperature, the liquid cannot return. This links CHF
directly to the Leidenfrost physics of the next page.

### Why jet quenching reaches 5–20 MW/m²

Measured quench MHFs (12–20 MW/m² in the ATE 2026 experiments) are an order of magnitude above pool CHF because:

- **Subcooling.** Part of the energy goes into sensible heating and condensing bubbles collapse.
  Correlations increase CHF roughly linearly with \(\rho_l c_{p,l}\Delta T_\text{sub}/(\rho_v h_{fg})\).
- **Liquid momentum.** The jet continuously replenishes the wall and sweeps vapour away. Impinging-jet
  CHF correlations (for example those of Monde and Katto) scale with jet velocity, diameter and the
  density ratio.
- **Transient wall.** The wall supplies heat from its stored energy at a rate set by its effusivity,
  rather than at a prescribed flux.

## What this means for reconstructions

- \(q_\text{max}\) and \(T_\text{DNB}\) are **local**. They vary with distance from the stagnation point
  (impingement zone vs free-flow zone). BoilFlux maps them over the surface (`regimes.feature_maps`).
- The nucleate peak is sharp in time. Its reconstructed amplitude is limited by the time resolution and the
  through-thickness damping of the plate. Benchmark MHF errors are reported separately for this reason.
