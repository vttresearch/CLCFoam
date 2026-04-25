# CLCFoam
![OpenFOAM dev](https://img.shields.io/badge/OpenFOAM-dev_20250616-brightgreen)

This repository contains materials for the article by Morev et al. [[1]](#1) on the coupling between heterogeneous chemistry and Eulerian-Eulerian multiphase model. The implementation is verified and validated using chemical-looping combustion (CLC) experimental data from literature.

## OpenFOAM codes
The heterogeneous chemistry model, coupled with Eulerian-Eulerian multiphase model, is implemented as [fvModel HeterogeneousShrinkingCoreChemistry](src/fvModels/HeterogeneousShrinkingCoreChemistry). Two types or reactions are implemented, but the number can be easily extended by implementing a new [HeterogeneousShrinkingCoreReaction](src/fvModels/HeterogeneousShrinkingCoreChemistry/HeterogeneousShrinkingCoreReaction.H).

Several [ThermophysicalTransportModels](src/ThermophysicalTransportModels) are implemented, disabling diffusion through patches to enable reactions close to the simulation inlet. 

Simulation case setups are automated using [C++ helper functions](src/common).

## Pre-processing
1. Reaction heat is enforced by correction of the JANAF polynomial of $\mathrm{Fe}_2\mathrm{O}_3$. Correction is done in [janaf_ilmenite.py](src/python/specie_properties/janaf_ilmenite.py) and is verified in [janaf.ipynb](janaf.ipynb).
2. Kinetic parameters are corrected in [scm.py](src/python/scm/scm.py).

## Simulation setups
Three simulation setups, used in the article, are available:
1. [0D](run/0D) - simple single-cell verification case, resembling thermogravimetric analysis conditions;  
2. [labReactor](run/labReactor) - batch reactor setup;  
3. [pilot300Wv2](run/pilot300Wv2) - small 300 W CLC reactor setup.  

These folders contain validation data from the literature and python post-processing scripts. Cases are automated to a degree, so that all required parameters of a simulation can be set in a single place - `caseParamsDict` file. Derived parameters (mass fractions of all species, for example) are calculated in `caseDict` file, using aforementioned helper functions, and then are used in other case dictionaries. `Allrun` script executes `setICBC`, which creates mass fraction field files from `Y.<phase>.template`.

## References

<a id="1">[1]</a>  
  Morev, Ilya, et al. “Coupling of Shrinking Core and Eulerian-Eulerian Models for Chemical Looping Combustion.” Chemical Engineering Science, 2026, p. 123431, [https://doi.org/10.1016/j.ces.2026.123431](https://doi.org/10.1016/j.ces.2026.123431).

<details>
<summary>BibTex</summary>
<p>
 
```
@article{morev2026Coupling,
  title={Coupling of shrinking core and Eulerian-Eulerian models for chemical looping
combustion},
  author = {Ilya Morev and Johanna Kapanen and Juho Peltola and Sirpa Kallio},
  journal = {Chemical Engineering Science},
  pages = {123431},
  year = {2026},
  issn = {0009-2509},
  doi = {https://doi.org/10.1016/j.ces.2026.123431}
}
```
 
</p>
</details>