# RIKsrf
-----------
Ruiz Integral Kinematic (RIK) earthquake source model: Slip rate generation and 1D full wavefield synthetics

Suite of code for earthquake strong ground motion simulations using an advanced kinematic source model based on the model of Ruiz et al. (2011). RIKsrf provides slip rates functions on a finely discretized source that result in synthetics with the desired
omega-squared spectral decay in full (broadband) frequency range.

####Main features of the RIK model:

- The source is composed of randomly distributed overlapping subsources with fractal number-size distribution.
Position of the subsources can be constrained by prior knowledge of major asperities
(stemming, e.g., from slip inversions), or can be completely random. This fractal composition of the source model implies that the slip decays as *k*<sup>-2</sup> at high wavenumbers *k*.

- The rise time is considered to depend on subsource radius, i.e. there is a positive correlation between slip and rise time
as observed in dynamic source modeling.

- The latter two properties ensures omega-squared decay of resulting source time function. 

- Rupture velocity and rise time can follow local S-wave velocity profile, so that the rupture slows down and rise times increase close to the surface, avoiding unrealistically strong ground motions.

- Rupture velocity can be either constant or with have random variations, which results in irregular rupture front while satisfying the
causality principle.

- The generated slip rates can be simply incorporated in any numerical wave propagation code without requiring any cross-over filtering with stochastic Green's functions.

------------

###Content of directories:
 - `src-RIKsrf` - Source codes of the RIK slip rate generator.
 - `src-1Dsynthetics` - Calculation of 1D full wavefield synthetics using Axitra for Green's function calculations.
 - `examples` - Several examples for testing the code including graphics (requires Gnuplot).
 - `docs` - Documentation and papers related to the RIK model.
 - `src-graphics` - Codes for generating graphics such as slip rate snapshots (requires Gnuplot).
