# Physical Model

Cosmic rays (CRs) provide important feedback in many astrophysical systems, but the prohibitively wide range of relevant scales makes this feedback difficult to simulate directly.

We derive CR fluid equations from the kinetic description to represent these interactions on macroscopic scales.

This page describes the model's physical scope, equations, and closure assumptions.


## 1. Scope

This module describes CRs as one or more fluid species coupled to an MHD gas. It evolves two fluid moments of the CR distribution: the energy density $`\mathcal{E}_{\rm cr}`$ and energy flux $`\boldsymbol{F}_{\rm cr}`$. The central modeling task is to represent microscopic CR transport and CR–gas coupling at the fluid scale. We consider gyroresonant scattering by Alfvén waves propagating parallel or antiparallel to the magnetic field. Four effective scattering coefficients, $`\sigma_{+}^L,\ \sigma_{+}^R,\ \sigma_{-}^L`$, and $`\sigma_{-}^R`$, distinguish the waves' propagation directions and circular polarizations.

The model is intended for macroscopic simulations in which kinetic-scale wave growth and damping reach local equilibrium on the simulation timescale. A local subgrid closure represents the unresolved saturated wave state.




## 2. Lab-Frame Two-Moment Equations

Define the scalar and vector components of the gas velocity parallel to the magnetic field by

```math
u_{\parallel}=\boldsymbol{u}\cdot\boldsymbol{b},
\qquad
\boldsymbol{u}_{\parallel}=u_{\parallel}\boldsymbol{b},
```

Here $`\boldsymbol{b}=\boldsymbol{B}/|\boldsymbol{B}|`$ is the magnetic-field unit vector. In Athena++ magnetic units, the Alfvén velocity is $`\boldsymbol{v}_A=\boldsymbol{B}/\sqrt{\rho}`$. The two-moment formulation is

```math
\frac{\partial \mathcal{E}_{\rm cr}}{\partial t}
+\boldsymbol{\nabla}\cdot\boldsymbol{F}_{\rm cr}=S_E,
```

```math
\frac{1}{V_m^2}\frac{\partial \boldsymbol{F}_{\rm cr}}{\partial t}
+\boldsymbol{\nabla}\cdot\mathbf{P}_{\rm cr}=\boldsymbol{S}_F.
```


The CR energy source is

```math
\begin{aligned}
S_E={}&-\sigma_+
\left[\boldsymbol{F}_{\rm cr}
-(\mathcal{E}_{\rm cr}+P_{\rm cr})(\boldsymbol{u}+\boldsymbol{v}_A)\right]
\cdot(\boldsymbol{u}_{\parallel}+\boldsymbol{v}_A)\\
&-\sigma_-
\left[\boldsymbol{F}_{\rm cr}
-(\mathcal{E}_{\rm cr}+P_{\rm cr})(\boldsymbol{u}-\boldsymbol{v}_A)\right]
\cdot(\boldsymbol{u}_{\parallel}-\boldsymbol{v}_A)\\
&+c\,\Delta P_{\rm cr}
\left[(\sigma_+^L-\sigma_+^R)(u_{\parallel}+v_A)
+(\sigma_-^L-\sigma_-^R)(u_{\parallel}-v_A)\right]\\
&+A\sigma_{\Omega}(\boldsymbol{F}_{\rm cr}\times\boldsymbol{b})\cdot\boldsymbol{u}.
\end{aligned}
```

The CR flux source is

```math
\begin{aligned}
\boldsymbol{S}_F={}&-\sigma_+
\left[\boldsymbol{F}_{\rm cr}
-(\mathcal{E}_{\rm cr}+P_{\rm cr})(\boldsymbol{u}+\boldsymbol{v}_A)\right]
\cdot\boldsymbol{b}\,\boldsymbol{b}\\
&-\sigma_-
\left[\boldsymbol{F}_{\rm cr}
-(\mathcal{E}_{\rm cr}+P_{\rm cr})(\boldsymbol{u}-\boldsymbol{v}_A)\right]
\cdot\boldsymbol{b}\,\boldsymbol{b}\\
&+c(\sigma_+^L-\sigma_+^R+\sigma_-^L-\sigma_-^R)
\Delta P_{\rm cr}\boldsymbol{b}\\
&+A\sigma_{\Omega}
\left[\boldsymbol{F}_{\rm cr}
-(\mathcal{E}_{\rm cr}+P_{\rm cr})\boldsymbol{u}\right]\times\boldsymbol{b}.
\end{aligned}
```

$`\sigma_+`$ and $`\sigma_-`$ denote the total scattering coefficients for forward- and backward-propagating waves:

```math
\sigma_+=\sigma_+^L+\sigma_+^R,
\qquad
\sigma_-=\sigma_-^L+\sigma_-^R.
```

We consider pressure anisotropy of order $`|\Delta P_{\rm cr}|\sim (v_A/c)P_{\rm cr}`$. At this level, anisotropy-driven self-confinement and energy exchange with the MHD gas can regulate the anisotropy. Although the source terms contain the large physical light speed $`c`$, it multiplies the small anisotropy, giving $`|c\Delta P_{\rm cr}|\sim v_A P_{\rm cr}`$.




## 3. Evolved Variables and Closure

For each CR species, the evolved conserved variables are

```math
\boldsymbol{q}
=\left(\mathcal{E}_{\rm cr},
\frac{F_{{\rm cr},1}}{V_m},
\frac{F_{{\rm cr},2}}{V_m},
\frac{F_{{\rm cr},3}}{V_m}\right)^T.
```

$V_m$ is the reduced numerical speed of light. It replaces the much larger physical light speed $c$ in the hyperbolic timestep constraint, allowing larger timesteps provided that $V_m$ remains well above all relevant physical signal speeds.

In a coordinate frame with $\hat{x}\parallel \boldsymbol{b}$, the CR pressure tensor takes the form

```math
\mathbf{P}_{\rm cr}=\begin{pmatrix}P_{\rm cr,\parallel} & 0 & 0\\ 0 & P_{\rm cr,\perp} & 0\\ 0 & 0 & P_{\rm cr,\perp} \end{pmatrix}
```
where $`P_{\rm cr,\parallel}`$ and $`P_{\rm cr,\perp}`$ are the pressures parallel and perpendicular to the magnetic field.

The CR pressure anisotropy is their difference:
```math
\Delta P_{\rm cr}\equiv P_{\rm cr,\parallel} - P_{\rm cr,\perp}
```

In the lab frame, the tensor components are

```math
P_{{\rm cr},ij}
=\frac{\mathcal{E}_{\rm cr}-\Delta P_{\rm cr}}{3}\delta_{ij}
+\Delta P_{\rm cr} b_i b_j,
```

with

```math
P_{\rm cr, \parallel}+2P_{\perp}={\mathcal{E}_{\rm cr}}.
```

Setting $`\Delta P_{\rm cr}=0`$ recovers the isotropic closure $`\mathbf{P}_{\rm cr}=\mathcal{E}_{\rm cr}\mathbf{I}/3`$. The current recommended `solver_id = 2` transport path uses the isotropic $`\mathcal{E}_{\rm cr}/3`$ interface pressure in its numerical flux. The HLLE solver includes the reconstructed pressure-anisotropy contribution in the CR pressure tensor. By default, however, $`\Delta P_{\rm cr}`$ is set to zero. Users should enable it only when they have a physically justified prescription; see [Pressure Anisotropy](physical_model.md#7-pressure-anisotropy) for details.



## 4. Notation

| Symbol | Meaning | Code representation |
|---|---|---|
| $\mathcal{E}_{\rm cr}$ | CR energy density | `cr_cons(cr_energy_id,...)`, `cr_prim(cr_energy_id,...)` |
| $\boldsymbol{F}_{\rm cr}$ | Lab-frame CR energy flux | Components `cr_flux1_id`, `cr_flux2_id`, and `cr_flux3_id` of `cr_cons` and `cr_prim`; stored as $\boldsymbol{F}_{\rm cr}/V_m$ so that all four state components have the same dimensions |
| $\mathbf{P}_{\rm cr}$ | CR pressure tensor | No direct code representation, calculated through $`\mathcal{E}_{\rm cr}$ and $\Delta P_{\rm cr}`$, see the closure above|
| $V_m$ | Reduced numerical speed of light| Given by `Vmax` under `<cosmic_ray>` in the input file; default $100.0$|
| $\gamma_{\rm cr}$ | CR adiabatic index | Given by `gamma_cr` under `<cosmic_ray>` in the input file; default $4/3$ |
| $\boldsymbol{u}$ | Gas velocity | Hydro primitive components `w(IVX,...)`, `w(IVY,...)`, and `w(IVZ,...)` |
| $v_A$ | Alfvén speed, $\lvert\boldsymbol{B}\rvert/\sqrt{\rho}$ | `vA` in the code |
| $\sigma_\pm^{L,R}$ | Alfvén-wave scattering coefficients | `sigma_pL_array`, `sigma_pR_array`, `sigma_mL_array`, and `sigma_mR_array` |
| $\sigma_{\Omega}$ | Effective Lorentz coefficient associated with the CR gyrofrequency | $\sigma_{\Omega}\equiv A\Omega_{\rm cr}/c^2$; represented by `sigma_Lorentz` in the code |
| $A$ | Dimensionless mean-Lorentz-factor prefactor | absorbed into the selected Lorentz coefficient in the current interface |
| $\Delta P_{\rm cr}$ | CR pressure anisotropy | One scalar per CR species, stored as `deltaPcr(CR_id,k,j,i)`|

## 5. Physical Interpretation

### 5.1 Transport and relaxation in the wave frame

The left-hand sides of our two-moment equations propagate CR energy and flux as a hyperbolic system. The scattering terms then relax the field-aligned flux.

Forward and backward waves can coexist. Each wave population tends to isotropize the CR distribution in its own rest frame, relaxing the CR energy flux toward

```math
\boldsymbol{F}_{\rm cr}
\rightarrow
(\mathcal{E}_{\rm cr}+P_{\rm cr})(\boldsymbol{u}\pm\boldsymbol{v}_A).
```
The corresponding flux-relaxation term on the right-hand side vanishes when $`\boldsymbol{F}_{\rm cr}=(\mathcal{E}_{\rm cr}+P_{\rm cr})(\boldsymbol{u}\pm\boldsymbol{v}_A)`$.

A perpendicular CR flux contributes to a macroscopic current and hence to the Lorentz term. The large CR gyrofrequency rapidly constrains transport relative to the magnetic field, consistent with charged particles being tied to field lines when their gyroradii are small.


### 5.2 Wave branches, force, and work


The four scattering coefficients $\sigma_\pm^{L,R}$ are labeled by:

- `+` or `-`: wave propagation along or opposite to $\boldsymbol{B}$;
- `L` or `R`: left- or right-handed circular polarization.

Separating these branches is necessary when pressure anisotropy selects different resonant polarizations.



The wave contributions to $`\boldsymbol{S}_{F}`$ represent the force exerted on the CR fluid. Let $`\boldsymbol{S}_{F,+}`$ and $`\boldsymbol{S}_{F,-}`$ denote the contributions containing the forward- and backward-wave coefficients, respectively. The work associated with these wave forces can be written as

```math
(\vec{u}+\vec{v}_A)\cdot\boldsymbol{S}_{F,+} + (\vec{u}-\vec{v}_A)\cdot\boldsymbol{S}_{F,-}
```
Each term is the force density dotted with the velocity of the corresponding scatterers. Their sum exactly matches the wave contribution in the $\mathcal{E}_{\rm cr}$ source term $S_E$, as expected.



## 6. Scattering-Coefficient Models

The scattering coefficients encode most of the microphysical transport information.

The module provides two built-in coefficient models and a user-defined callback interface.

### 6.1 Constant coefficients

For controlled tests, set constants for all five coefficients for each CR species:

```ini
scattering_model   = constant
sigma_pL_1         = ...
sigma_pR_1         = ...
sigma_mL_1         = ...
sigma_mR_1         = ...
sigma_Lorentz_1    = ...
```

Equal coefficients in all four Alfvén-wave branches represent direction- and polarization-balanced scattering. This is also the simplest way to approximate extrinsic-turbulence transport, for which one may set

```math
\sigma_+^L=\sigma_+^R=\sigma_-^L=\sigma_-^R=\frac{\sigma_{\rm turb}}{4}.
```

### 6.2 CR streaming-instability coefficients

The built-in streaming model selects the active wave direction from the sign of $`\boldsymbol{b}\cdot\boldsymbol{\nabla}\mathcal{E}_{\rm cr}`$ and divides the scattering coefficient equally between the two polarizations in that direction:

```math
\sigma_{\pm}
=\widetilde{\sigma}_{0,\rm st}
\frac{|\boldsymbol{b}\cdot\boldsymbol{\nabla}{P}_{\rm cr}|}
{{P}_{\rm cr}v_A}\Theta(\mp\boldsymbol{b}\cdot\nabla P_{\rm cr})
```
In input files, $\widetilde{\sigma}_{0,\rm st}$ is `sigma0_streaming_N` for the $N$-th CR species.

This closure captures the gradient dependence expected when CR streaming-instability growth is locally balanced by linear damping. It is a steady-state subgrid prescription; the Alfvén-wave dynamics is not evolved.

### 6.3 User-defined coefficients

More general physics should be implemented through `Mesh::EnrollCRScatteringCoefficients`. A callback receives gas and CR primitive variables and fills

```cpp
sigma_pL, sigma_pR, sigma_mL, sigma_mR, sigma_Lorentz, deltaPcr
```

over the requested cell range. This interface can accommodate user-defined scattering coefficients (and prescriptions for CR pressure anisotropy) that reflect the corresponding microphysics, such as nonlinear Landau, ion-neutral damping, CR pressure-anisotropy instability, spatially varying extrinsic turbulence, or coefficients calibrated from kinetic simulations. See [User-Defined Scattering](parameter_reference.md#4-user-defined-scattering) for the callback interface, and a [test example](examples.md#3-1d-cr-streaming-with-nonlinear-landau-damping) for a concrete implementation.



## 7. Pressure Anisotropy

A two-moment method requires a closure for the CR pressure tensor. The isotropic closure sets $`\mathbf{P}_{\rm cr}=\mathcal{E}_{\rm cr}\mathbf{I}/3`$.

In a magnetized plasma, charged-particle dynamics differ along and across the magnetic field, motivating a more general description with distinct parallel and perpendicular pressures. We characterize their difference by $`\Delta P_{\rm cr}\equiv P_{\rm cr,\parallel}-P_{\rm cr,\perp}`$. This anisotropy affects both the transport flux through $`\nabla\cdot\mathbf{P}_{\rm cr}`$ and the source terms through an additional channel of microphysical CR–gas coupling.

The CR pressure-anisotropy instability (CRPAI) can be excited when $`\frac{|\Delta P_{\rm cr}|}{P_{\rm cr}}\gtrsim\frac{v_A}{c}`$. The resulting waves scatter CRs and regulate their anisotropy. A representative steady-state estimate is $`\frac{\Delta P_{\rm cr}}{P_{\rm cr}}
\sim \pm\frac{v_A}{c}
\frac{\nu_{\rm damp}}{\Omega_{\rm cr}}
\frac{\rho}{\rho_{\rm cr}}`$. For $`\Delta P_{\rm cr}<0`$, a branch selection motivated by quasi-linear theory (QLT) is $`\sigma_+^L=\sigma_-^R>0`$ and $`\sigma_+^R=\sigma_-^L=0`$; the active polarizations reverse for $`\Delta P_{\rm cr}>0`$.

We consider a regime in which this regulation reaches a locally saturated state on timescales shorter than those resolved by the simulation. Accordingly, $`\Delta P_{\rm cr}`$ (kept on the level of $`v_A P_{\rm cr}/c`$) is prescribed through a local subgrid closure rather than evolved as an independent fluid variable. 

Since both $`\Delta P_{\rm cr}$ and $\sigma_{\pm}^{L,R}`$ represent unresolved microphysics and depend on prescriptions, the code assigns them together when setting the scattering properties. This treatment is consistent with our aim of incorporating CRPAI as an additional coupling mechanism at the fluid scale.



## 8. CR Feedback on the Gas

The module provides three feedback modes, corresponding to different levels of coupling between CRs and the MHD gas:

1. No feedback: With `CRFeedback_Flag = false`, the gas receives neither energy nor momentum from the CRs. The CR energy source is also explicitly set to zero, so $`\mathcal{E}_{\rm cr}`$ evolves only through the advection term $`\nabla\cdot\boldsymbol{F}_{\rm cr}`$. This mode is primarily intended for tests, as it simplifies the analysis of the CR subsystem.

2. Momentum feedback: With `CRFeedback_Flag = true` and `CR_EnergyFeedback_Flag = false`, CR momentum transfer updates the gas momentum and the corresponding kinetic energy. The $`\mathcal{E}_{\rm cr}`$ source remains active, but the associated energy exchange is not applied to the gas. This mode approximates situations in which CR heating is balanced by some cooling process or CR energy is inefficiently deposited as gas heat.

3. Full feedback: With both flags set to `true`, CRs exchange both energy and momentum with the gas. The corresponding CR source terms enter the gas energy and momentum equations with opposite signs.

If `CRFeedback_Flag = false` while `CR_EnergyFeedback_Flag = true`, the code issues a warning, ignores the energy-feedback flag, and selects the no-feedback mode.


## 9. Choosing a Physical Model

1. Select $V_m$ well above all relevant acoustic, Alfvén, CR and other characteristic speeds.
2. Choose a built-in scattering model or specify the four $\sigma_\pm^{L,R}$ branches and the Lorentz coefficient through the user callback.
3. Enable $\Delta P_{\rm cr}$ only when a defensible local prescription is available.
4. Select the intended CR–gas coupling: no feedback, momentum feedback, or full energy and momentum feedback, subject to the implementation status above.
5. (Recommended) Verify timestep convergence, especially near extrema where a gradient-dependent streaming coefficient changes sign.

See [Numerical Algorithm](numerical_algorithm.md) for how these quantities enter each stage and [Known Limitations](known_limitations.md) for the current validation boundary.
