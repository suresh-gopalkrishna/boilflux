# Quenching and rewetting with impinging jets

## Hydrodynamic zones of a round jet on a hot plate

| Zone | Flow | Boiling behaviour |
|---|---|---|
| **Stagnation / impingement zone (IZ)** (≈ 1–2 jet diameters, plus the early wetted core) | flow decelerates, highest pressure and liquid supply, thinnest thermal boundary layer | earliest rewetting, highest \(q_\text{max}\), strongly subcooled nucleate boiling |
| **Radial wall-jet / parallel-flow zone** | thin liquid film spreading radially and slowing down | wetting front propagates outwards; transition boiling just behind the front |
| **Hydraulic jump / splashing** | film thickens abruptly or is ejected as droplets by vapour | weaker, irregular cooling |
| **Free-flow zone (FFZ)** / beyond the front | liquid flows over a vapour film or falls under gravity | film boiling, low flux |

In the QuenchFlux synthetic model these zones are represented by zone-dependent boiling curves (IZ, FFZ)
and a liquid-coverage radius. In the real flow they come from the momentum and mass balance of the
spreading film, which is what [CFD and neural operators](modelling.md) try to predict.

## The wetting front

Rewetting does not occur everywhere at once. A **wetting front** (the boundary between the wetted region
and the film-boiling region) propagates outwards from the stagnation point. Ahead of it, the wall
is still in film boiling. Just behind it, a narrow band of transition and peak nucleate boiling removes
most of the heat.

Classical **conduction-controlled rewetting** models (Yamanouchi; Duffey and Porthouse, developed for
nuclear reflood) treat the front as a moving boundary. Behind it a high heat-transfer coefficient cools
the wall. Ahead of it, axial conduction precools the dry wall to the rewetting temperature. The front velocity
then depends on a Péclet number \(u_f\,s/\alpha\) and a Biot number \(h\,s/k\) of the plate of thickness \(s\):

- thin, highly conducting walls: front speed ∝ \(\sqrt{h/(k s)}\,\alpha / (\text{superheat ratio})\)
- thick walls: two-dimensional conduction, with speed nearly independent of thickness

The superheat ratio \((T_w - T_\text{Leid})/(T_\text{Leid} - T_\text{sat})\) is central: the hotter the wall
above its rewetting temperature, the slower the front. This is consistent with the slow fronts on
800 °C nickel compared with 500 °C AA6082.

## Moving plates

When the plate moves under a stationary jet (run-out tables, strip cooling, moving-plate rigs):

- The cooling pattern is **asymmetric**. On the upstream side the plate arrives hot and meets the front
  head-on. On the downstream side already-cooled metal leaves the jet. The upstream front is pushed against the plate
  motion and can stall when the plate speed equals the front speed in the plate frame.
- A material point experiences a **sequence of regimes in time** that depends on its distance from the
  jet path, not only on its temperature.
- The relevant frame for conduction is the plate frame. For the flow, it is the jet frame. QuenchFlux uses the plate
  frame for the inverse problem (see [moving-plate formulation](../theory/moving-plate.md)).

## Material and surface effects

| Property | Effect |
|---|---|
| thermal diffusivity \(\alpha\) | faster lateral precooling ahead of the front, so faster rewetting (AA6082 compared with Ni) |
| effusivity \(e=\sqrt{k\rho c_p}\) | sets the contact temperature, and through it the Leidenfrost behaviour |
| initial temperature | larger superheat, so a longer film-boiling period and slower front |
| oxide layers, coatings, roughness | can raise \(T_\text{Leid}\) markedly and cause specimen-to-specimen scatter |
| phase transformations (steels) | latent heat appears as an effective \(c_p\) peak; it must be in the property table |

## What a back-face IR inverse solution can and cannot resolve

The plate filters the boundary flux before it reaches the IR face. As a rule of thumb, features
**smaller than about the plate thickness** in space, or **faster than about \(s^2/\alpha\)** in time, are
strongly attenuated:

| Quantity | Resolvable from back-face IR? |
|---|---|
| \(q(x,y,t)\) at the scale of the jet and the front | yes |
| \(q_\text{max}\) and \(T_\text{DNB}\) at the stagnation point | yes, with smoothing of the peak amplitude |
| \(T_\text{Leid}\), rewetting time, front velocity | yes |
| individual bubbles, microlayer, contact-line flux | **no**: needs front-face measurement (for example IR through a transparent heater) or interface-resolving CFD |
| very fast transients (≪ \(s^2/\alpha\)) | only their time integral |

These limits motivate the **combination** pursued in QuenchFlux: inverse reconstruction for the
plate-scale flux, with flow models (CFD, neural operators) for the mechanisms below that scale.
