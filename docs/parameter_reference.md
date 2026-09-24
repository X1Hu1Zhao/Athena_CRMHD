# Parameter Reference

## 1. Compile-Time Configuration

The CR module adds one Athena++ configure option:

```text
--ncrs=N
```

It sets the number of CR fluid species at compile time. The generated macros are

```cpp
#define NCRS     N      // Number of CR species
#define NCRVARS  (N*4)  // One energy and three flux components per species
#define NSPECIES (N+1)  // Number of fluid species, i.e. CR species plus the gas
```

Use `--ncrs=0` to disable our CR module and `--ncrs=1` for the recommended starting configuration. Include `-b` for the magnetized CR calculations described in this guide; the built-in streaming model explicitly requires magnetic fields.

Example:

```bash
python configure.py --prob=CR_streaming -b -hdf5 --ncrs=1 --hdf5_path=/path/to/hdf5
```

`--ncrs` is the only new configure option introduced by this module. All other options are standard Athena++ options.

## 2. Core `<cosmic_ray>` Parameters

### 2.1 Transport and reconstruction

| Parameter | Type | Default | Meaning |
|---|---:|---:|---|
| `Vmax` | Real | `100.0` | Numerical speed of light $V_m$ in code units. Must exceed all other physical signal speeds. |
| `gamma_cr` | Real | `4/3` | CR adiabatic index used by scattering integrators. Keep `4/3` for the validated path. |
| `cr_xorder` | Integer | Inherits `xorder` | CR reconstruction order, capped at the hydro reconstruction order. Use `1` or the recommended `2`. Higher-order branches are not part of the validated path. |
| `solver_id` | Integer | `2` | `2` selects the recommended adaptive HLLE solver. Currently we only support this solver choice.|
| `crfloor_N` | Real | Athena++ floating-point floor | Minimum CR energy density for $N$-th CR species. |

### 2.2 Scattering model selection

| Parameter | Type | Default | Requirement and meaning |
|---|---:|---:|---|
| `Scattering_Flag` | Boolean | Must be `true` | Controls scattering-coefficient updates in `SetProperties`. It must currently be `true` when the CR module is active. |
| `scattering_model` | String | `constant` | Built-in model: `constant` or `streaming`. Ignored when a user callback is enrolled. |
| `CRFeedback_Flag` | Boolean | Required | Controls CR feedback to the gas. If `false`, the gas receives neither CR energy nor momentum, and the CR-energy scattering source is suppressed. |
| `CR_EnergyFeedback_Flag` | Boolean | `true` | When `CRFeedback_Flag = true`, selects momentum-only feedback (`false`) or full energy-and-momentum feedback (`true`). It is ignored, with a warning, when `CRFeedback_Flag = false`. |
| `scattering_integrator` | String | `2nd-implicit` | Selects the scattering source integrator. Use `2nd-implicit` with `vl2` or `1st-implicit` with `rk1`. The second-order combination is recommended. |



### 2.3 Constant-coefficient model

With

```ini
scattering_model = constant
```

the following values are required for each species `N = 1, ..., NCRS`:

| Parameter | Meaning |
|---|---|
| `sigma_pL_N` | Forward, left-polarized scattering coefficient $\sigma_+^L$ |
| `sigma_pR_N` | Forward, right-polarized scattering coefficient $\sigma_+^R$ |
| `sigma_mL_N` | Backward, left-polarized scattering coefficient $\sigma_-^L$ |
| `sigma_mR_N` | Backward, right-polarized scattering coefficient $\sigma_-^R$ |
| `sigma_Lorentz_N` | Lorentz-coupling coefficient, including any chosen prefactor $A$. A large value is intended to constrain CR transport relative to the magnetic field (by default `sigma_Lorentz_N = 1e8`). |

Example:

```ini
<cosmic_ray>
Scattering_Flag        = true
scattering_model       = constant
sigma_pL_1             = 0.5
sigma_pR_1             = 0.5
sigma_mL_1             = 0.0
sigma_mR_1             = 0.0
sigma_Lorentz_1        = 1e8
CRFeedback_Flag        = true
scattering_integrator  = 2nd-implicit
```

### 2.4 Streaming model

With

```ini
scattering_model = streaming
```

the model reads:

| Parameter | Type | Default | Meaning |
|---|---:|---:|---|
| `sigma0_streaming_N` | Real | `0.0` | Dimensionless streaming prefactor $\widetilde{\sigma}_{0,\rm st}$ for species `N`; must be non-negative. |
| `sigma_Lorentz_N` | Real | `0.0` | Lorentz-coupling coefficient for species `N`. |

At least one `sigma0_streaming_N` must be non-zero, or initialization stops with a fatal error. The largest value among species sets the built-in streaming timestep constraint described in [Timestep Restrictions](known_limitations.md#7-timestep-restrictions).

The model requires magnetic fields and computes

```math
\sigma_{\pm}
=\widetilde{\sigma}_{0,\rm st}
\frac{|\boldsymbol{b}\cdot\boldsymbol{\nabla}P_{\rm cr}|}
{P_{\rm cr}v_A}\Theta(\mp\boldsymbol{b}\cdot\nabla P_{\rm cr}).
```


## 3. Explicit CR Source Parameters

| Parameter | Block | Type | Default | Meaning |
|---|---|---:|---:|---|
| `pressure_anisotropy_flag` | `cosmic_ray` | Boolean | `false` | Master switch for CR pressure anisotropy. If `false`, `deltaPcr` is cleared after `SetProperties`, anisotropic transport is disabled, and the explicit anisotropy source is skipped. |
| `grav_acc1` | `cosmic_ray` | Real | `0.0` | Constant acceleration in coordinate direction 1. |
| `grav_acc2` | `cosmic_ray` | Real | `0.0` | Constant acceleration in coordinate direction 2. |
| `grav_acc3` | `cosmic_ray` | Real | `0.0` | Constant acceleration in coordinate direction 3. |
| `GM` | `problem` | Real | `0.0` | Point-mass parameter used in supported cylindrical or spherical-polar configurations. |


## 4. User-Defined Scattering

A problem generator may replace the built-in coefficient models by enrolling a callback in `Mesh::InitUserMeshData`:

```cpp
EnrollCRScatteringCoefficients(MyCRScatteringCoefficients);
```

The callback type is declared in `src/athena.hpp`:

```cpp
using CRScatteringCoeffFunc = void (*)(
    CosmicRay *cr_pointer, MeshBlock *pmb,
    const AthenaArray<Real> &w,
    const AthenaArray<Real> &cr_prim,
    const AthenaArray<Real> &bcc,
    AthenaArray<Real> &sigma_pL,
    AthenaArray<Real> &sigma_pR,
    AthenaArray<Real> &sigma_mL,
    AthenaArray<Real> &sigma_mR,
    AthenaArray<Real> &sigma_Lorentz,
    AthenaArray<Real> &deltaPcr,
    AthenaArray<Real> &cs_cr,
    int is, int ie, int js, int je, int ks, int ke);
```
Registration overrides the built-in `scattering_model` selection. The callback must assign both the scattering coefficients and the species-dependent pressure anisotropy `deltaPcr(CR_id,k,j,i)` over the requested cells. A non-zero CR pressure anisotropy is retained only when `pressure_anisotropy_flag = true`; otherwise, the array is cleared after the callback returns. See [Set Scattering Properties](numerical_algorithm.md#4-step-1-set-scattering-properties) for the callback contract.

## 5. CR Output

For species `N`, CR datasets use names that describe their physical contents:

| Output name | Stored quantity |
|---|---|
| `cr_N_energy_dens` | CR energy density $\mathcal{E}_{\rm cr}$ |
| `cr_N_flux` | Vector $\boldsymbol{F}_{\rm cr}/V_m$ |
| `cr_N_flux1`, `cr_N_flux2`, `cr_N_flux3` | Individual flux components divided by $V_m$ |
| `cr_N_xyz_flux` | Cartesian-converted flux vector where applicable |
| `cr_N_deltaPcr` | Prescribed pressure anisotropy $\Delta P_{\rm cr}$; included in aggregate `prim` and `cons` output only when `pressure_anisotropy_flag = true` |

The aggregate selector `variable = cons` writes the gas and CR conserved states, whereas `variable = prim` writes the gas and CR primitive states. Because the current CR EOS copies the CR values between these representations after applying the energy floor, their CR datasets are numerically identical. If pressure anisotropy is enabled, either aggregate selector also includes one `cr_N_deltaPcr` scalar per CR species.

Example output block:

```ini
<output2>
file_type = hdf5
variable  = cons
dt        = 0.01
```


Users may define custom output variables through `user_out_var` in the problem generators, allocate them with `AllocateUserOutputVariables(output numbers)`, and enable their output by setting `variable = uov` in the corresponding `<output>` block of the input file. See given generator files for more examples.