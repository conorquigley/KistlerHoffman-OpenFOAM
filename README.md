# KistlerHoffman-OpenFOAM

A dynamic contact angle boundary condition for OpenFOAM based on the **Kistler–Hoffman model**, together with a modified `interFoam` solver (`interKistlerFoam`) that exposes the field quantities the boundary condition requires. Compatible with **OpenFOAM v2412** (updated from v7).

Associated publication: *OpenFOAM® Journal Publication — Kistler–Hoffman Validation*, C. Quigley & L. Ó Náraigh, UCD.

---

## Table of Contents

- [Background](#background)
- [Physics](#physics)
  - [Volume of Fluid Method](#volume-of-fluid-method)
  - [The Kistler–Hoffman Contact Angle Model](#the-kistlerhoffman-contact-angle-model)
  - [Contact Line Velocity: The Key Implementation Choice](#contact-line-velocity-the-key-implementation-choice)
- [Repository Structure](#repository-structure)
- [Dependencies](#dependencies)
- [Key Parameters](#key-parameters)
- [Mesh Setup](#mesh-setup)
- [Validation](#validation)
- [Limitations](#limitations)
- [Citation / References](#citation--references)
- [License](#license)

---

## Background

Accurate simulation of droplet impact and spreading depends critically on the dynamic contact angle — the angle formed between the liquid–gas interface and a solid surface at the three-phase contact line. A static or constant contact angle is insufficient for capturing the transient wetting behaviour observed in experiments. As it stands in OpenFOAM, the options for a dynamic contact angle implementation are limited, and many community implementations lack rigorous validation.

This repository provides a validated OpenFOAM implementation of the **Kistler (1993)** dynamic contact angle model. It is designed for droplet impact simulations using the Volume-of-Fluid (VOF) `interFoam` solver (or the bundled `interKistlerFoam` variant). The model has been validated against the experimental and simulation results of Roisman et al. (2008), which were themselves further validated by Esteban et al. (2023).

The key contribution of this work is identifying a practical, easily implemented choice of contact line velocity that produces results in good agreement with published literature, without requiring a custom solver or complex mathematical formulations.

---

## Physics

### Volume of Fluid Method

The VoF method models two-phase flow using a one-fluid formulation of the Navier–Stokes equations (Brackbill et al., 1992). A phase-fraction field α(x) ∈ [0, 1] tracks the interface, where α = 0.5 marks the interface between liquid and gas. The mixture density and viscosity are interpolated as:

```
ρ(x) = ρ_L α(x) + ρ_G [1 − α(x)]
μ(x) = μ_L α(x) + μ_G [1 − α(x)]
```

The full system solved is the incompressible Navier–Stokes equations with a continuum surface tension force F_ST, together with the phase-fraction advection equation:

```
∂α/∂t + u·∇α = 0
```

### The Kistler–Hoffman Contact Angle Model

The dynamic contact angle θ_d is computed from the Kistler–Hoffman function (Kistler, 1993):

```
θ_d = f_H( Ca + f_H⁻¹(θ_e) )
```

where the **Hoffman function** f_H is the empirical function:

```
f_H(x) = arccos{ 1 − 2·tanh[ 5.16·( x / (1 + 1.31·x^0.99) )^0.706 ] }
```

The capillary number Ca encodes the contact-line speed and fluid properties:

```
Ca = μ · U_cl / σ
```

where μ is the dynamic viscosity of the liquid, U_cl is the contact-line velocity, and σ is the surface tension. The equilibrium angle θ_e switches depending on the direction of contact-line motion:

| Contact line state | Reference angle |
|---|---|
| Advancing (U_cl > 0) | Advancing angle θ_A |
| Static (U_cl = 0) | Equilibrium angle θ_E |
| Receding (U_cl < 0) | Receding angle θ_R |

The capillary number acts as a shift factor on the chosen equilibrium angle, allowing the dynamic contact angle to exceed the static advancing or receding limits during rapid impact events.

### Contact Line Velocity: The Key Implementation Choice

The accuracy of the Kistler–Hoffman model depends critically on how U_cl is estimated from the simulation field. A naive approach — sampling the radial velocity from the boundary-layer cells (n = 0) — suffers from two problems:

1. **No-slip condition**: boundary-layer cell velocities are forced toward zero, underestimating the true contact-line speed.
2. **Localised vortices**: during impact, vortices form near the interface and boundary (see Figures 3 and 4 in the paper). These can cause the local velocity vector to point inward even while the droplet as a whole is spreading outward — leading cells in the advancing phase to be misclassified as receding, and vice versa. This produces significant overestimation of the maximum spreading radius.

**This implementation samples the radial velocity n cells above the boundary layer**, where n ≥ 1 (n = 0 being the boundary-layer cells themselves). Sampling above the vortex-dominated region means the radial velocity better reflects the actual outward or inward motion of the contact line. Based on the validation studies in this work, **n = 3 was found to produce the best results**.

One important advantage of this approach is that it can read μ and σ as constants set before the simulation, rather than requiring them to be written dynamically by the solver at each timestep. This means the boundary condition is compatible with the standard, unmodified `interFoam` solver — making it easy to adopt as an additional contact angle option without needing a custom solver build.

---

## Repository Structure

```
KistlerHoffman-OpenFOAM/
│
├── KHContactAngle/
│   └── dynamicKistlerAlphaContactAngle/
│       ├── dynamicKistlerAlphaContactAngleFvPatchScalarField.H   # Class declaration
│       ├── dynamicKistlerAlphaContactAngleFvPatchScalarField.C   # Implementation
│       └── Make/
│           ├── files
│           └── options
│
└── README.md
```

**`KHContactAngle`** — the boundary condition library. Derives from `alphaContactAngleFvPatchScalarField` and overrides the contact angle evaluation to use the Kistler–Hoffman formulation with the n-layer velocity sampling strategy.


## Dependencies

- **OpenFOAM v2412** (ESI / openfoam.com branch). The code has been ported from OpenFOAM v7 and is not guaranteed to compile against other versions without modification.
- A working `wmake` build environment (sourced OpenFOAM shell environment).

---

## Installation

### 1. KHContactAngle boundary condition

Copy the `KHContactAngle` directory into your OpenFOAM user source tree, or any location you prefer:

```bash
cp -r KHContactAngle $FOAM_RUN/../src/
cd $FOAM_RUN/../src/KHContactAngle
./Allmake.sh
```

Or invoke `wmake` directly:

```bash
cd KHContactAngle/dynamicKistlerAlphaContactAngle
wmake libso
```

The compiled shared library (`libdynamicKistlerAlphaContactAngle.so`) will be placed in `$FOAM_USER_LIBBIN`.

---

## Key Parameters

| Parameter | Description | Recommended value |
|---|---|---|
| `thetaA` | Advancing (maximum) contact angle (degrees) | Case-specific |
| `thetaR` | Receding (minimum) contact angle (degrees) | Case-specific |
| `theta0` | Equilibrium contact angle (degrees) | Case-specific |
| `nLayers` | Cell layers above boundary for velocity sampling | Case-specific |
| `mu` | Liquid dynamic viscosity (kg/m/s) | 1e-3 for water |
| `sigma` | Surface tension (N/m) | 0.072 for water/air |

The choice of `nLayers` is the most case-sensitive parameter. A value of n = 3 was found to produce the best results in the validation studies here, but users should test this for their own droplet initial conditions if experimental data is available for comparison.

---

## Mesh Setup

The validation simulations used a **wedge geometry** (axisymmetric) with the following properties:

- Domain: 4 mm × 4 mm, wedge half-angle of 2.8578° / 2
- Base mesh: (100 × 1 × 100) cells
- Refinement: `snappyHexMesh`, giving a total of **328,984 cells**
- Cell distribution: L-shaped refinement zone concentrated along the symmetry axis and the substrate boundary, ensuring the droplet spreading region is well resolved while keeping computational cost low in the far field

Mesh convergence was verified by testing a sequence of progressively finer meshes (Mesh 20 through Mesh 120) until results for the spreading factor β(t) and normalised height h(t)/h(1) converged. The coarsest converged mesh was selected to minimise runtime.

---

## Validation

The model was validated against three droplet impact experiments from Roisman et al. (2008) — water droplets impacting a stainless steel substrate. All three cases used a drop diameter of 2.5 mm, advancing angle θ_A = 120°, and receding angle θ_R = 65°:

| Case | Impact velocity (m/s) | Weber number |
|---|---|---|
| Experiment 1 | 0.16 | 0.88 |
| Experiment 2 | 0.48 | 7.9 |
| Experiment 3 | 0.23 | 1.81 |

Results are compared in terms of the non-dimensional spreading factor β(t) = R(t)/R₀ and the normalised droplet height h(t)/h(1).

Key findings:

- **Advancing phase**: very good agreement with both Roisman's simulations and experimental data across all three cases. This is expected as the contact-line velocity is largest and most clearly resolved during spreading.
- **Receding/pinning phase**: results remain within expected targets, though the true contact-line velocity is more susceptible to contamination from spurious currents during this phase.

---

## Limitations

- **n-layer tuning**: the choice of `nLayers` is practical rather than mathematically derived. While n = 3 produced the best results in the cases studied here, this may not generalise to all droplet initial conditions. Users working with significantly different We or Re should verify this choice against available experimental or reference simulation data.
- **Constant μ and σ**: by default, the boundary condition reads μ and σ as constants, which is appropriate for isothermal, single-component droplets. For cases with variable surface tension (e.g. surfactant-laden drops or thermal effects), modifications to interFOAM will be neeeded to update the mu and sigma registries
- **OpenFOAM version**: tested on v2412 only. The original v7 code required changes to the receding-phase sign convention of the capillary number; other versions may need adaptation.

---

## Citation / References

If you use this boundary condition in published work, please cite the following:

[This work was build from the linked githib repositry.]((https://github.com/franciscotovarmit/interKistlerFoam)) The work done there is foundational to anything that has been achieved by these updates.



**Kistler model:**
> S. F. Kistler, "Hydrodynamics of wetting," in *Wettability*, J. C. Berg, Ed. New York: Marcel Dekker, 1993, p. 311.

**Validation reference:**
> I. V. Roisman, L. Opfer, C. Tropea, M. Raessi, J. Mostaghimi, and S. Chandra, "Drop impact onto a dry surface: Role of the dynamic contact angle," *Colloids and Surfaces A*, vol. 322, no. 1–3, pp. 183–191, 2008.

**Further validation:**
> A. Esteban, P. Gómez, C. Zanzi, J. López, M. Bussmann, and J. Hernández, "A contact line force model for the simulation of drop impacts on solid surfaces using volume of fluid methods," *Computers & Fluids*, vol. 263, 2023.

**VoF surface tension formulation:**
> J. Brackbill, D. Kothe, and C. Zemach, "A continuum method for modeling surface tension," *Journal of Computational Physics*, vol. 100, no. 2, pp. 335–354, 1992.

---

## License

This code is distributed under the **GNU General Public License v2 (or later)**, consistent with the OpenFOAM license under which the base `alphaContactAngleFvPatchScalarField` class is distributed.

See [https://www.gnu.org/licenses/old-licenses/gpl-2.0.html](https://www.gnu.org/licenses/old-licenses/gpl-2.0.html) for full terms.
