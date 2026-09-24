# Limitations

This page describes the current validation scope and the limitations relevant to interpreting or using the code.

## 1. Recommended and Validated Path

The recommended starting configuration is:

```text
NCRS = 1
gamma_cr = 4/3
time integrator = vl2
hydro reconstruction = 2
CR reconstruction = 2
solver_id = 2
scattering_integrator = 2nd-implicit
```

Earlier development tests exercised 1D CR streaming and several coupled CR–MHD eigenmodes, including CR-acoustic, MHD-acoustic, growing, and rapidly damped branches. This does not establish validation for every compile-time option, coordinate system, integrator, or physical closure.

## 2. Single-Species Validation

The infrastructure supports compile-time allocation for `NCRS > 1`, and the primary state, flux, scattering, boundary, and matrix arrays include species indices. The method paper and development regression tests, however, focus on one CR species. Therefore:

- ordinary multi-species transport should be treated as unvalidated;
- multi-species pressure-anisotropy calculations should not be used without code review and dedicated tests;


## 3. Pressure Anisotropy Is an External Closure

The numerical interfaces for a prescribed CR pressure anisotropy are implemented throughout the module. The code stores one `deltaPcr(CR_id,k,j,i)` value per CR species, reconstructs it to cell interfaces for the full pressure-tensor flux, and passes it to the explicit anisotropy source. However, $\Delta P_{\rm cr}$ is not evolved as an independent fluid variable. Users prescribe it alongside the scattering coefficients when the scattering properties are set, because both quantities are local closures for unresolved microphysics.

The principal limitation is therefore the physical prescription for $\Delta P_{\rm cr}$, rather than the availability of numerical interfaces. As a safety control, `pressure_anisotropy_flag = false` clears the array after `SetProperties`, recovers isotropic HLLE transport, and skips the explicit anisotropy source.

Although the numerical interface for $\Delta P_{\rm cr}$ is ready, we recommend using it cautiously because of its uncertain physical prescription.

The current framework is intended for the small-anisotropy regime

```math
\left|\frac{\Delta P_{\rm cr}}{P_{\rm cr}}\right|
\lesssim \frac{v_A}{c}.
```

The underlying physical picture is that, once the relative anisotropy approaches this level, the CR pressure-anisotropy instability (CRPAI) excites gyroresonant waves. These waves scatter the CRs and feedback on their distribution, allowing the anisotropy to self-regulate near this marginal or locally saturated level. The prescribed $\Delta P_{\rm cr}$ should represent such unresolved, self-regulated state.

This assumption need not hold in regimes where a much larger pressure anisotropy is maintained by other physics or imposed for a different modeling purpose. Such problems may require a different closure or a kinetic treatment, and the present model should not be assumed to apply to them.


## 4. Scattering Physics Is Subgrid

The module does not evolve forward/backward Alfvén-wave energies. Constant and streaming coefficients assume that unresolved wave growth and damping can be represented by a local saturated state. This may fail in regimes with:

- Alfvén-wave dark regions;
- non-local wave transport;
- comparable kinetic and fluid timescales;
- simultaneous CRSI and CRPAI without a calibrated closure;
- dominant non-wave scattering, such as transport controlled by magnetic reversals or intermittency.

A user-defined callback can supply a more complete closure, but the user is responsible for its physical validity and numerical stability.

## 5. Reduced-Speed-of-Light Approximation

`Vmax` must remain much larger than all dynamically important gas, Alfvén, CR, and coupled-mode speeds. A value that is merely numerically convenient can alter wave frequencies, source-term relaxation, and feedback.

Repeat at least one representative calculation with a larger `Vmax`. Spatial-resolution convergence alone does not establish convergence with respect to the reduced propagation speed.

The numerical speed of light $V_m$ is constant throughout a run. Reduced-speed-of-light methods have been studied in previous work, but extensions such as spatially varying $V_m$ require further investigation and validation.

## 6. `gamma_cr` Should Remain `4/3`

In the ultrarelativistic limit, the CR adiabatic index is $`4/3`$, giving $`P_{\rm cr}=(\gamma_{\rm cr}-1)\mathcal{E}_{\rm cr}=\mathcal{E}_{\rm cr}/3`$. Parts of the module use this relation explicitly. Changing the input parameter alone therefore does not produce a consistent model with a different adiabatic index.

Keep $`\gamma_{\rm cr}=4/3`$ unless every use of the closure has been audited and the alternative model has been validated. In an input file, write `gamma_cr = 1.3333333333333333`; Athena++ expects a numeric value, not the expression `4/3`.


## 7. Timestep Restrictions

Before the global CFL factor, the current built-in streaming model applies

```math
\Delta t_{\rm scatt}
=\frac{\sqrt{3}\Delta s_{\min}}
{(1+4\widetilde{\sigma}_{0,\rm st})V_m}.
```

This guards the sign-changing streaming closure near CR-pressure extrema. It is not a universal timestep for arbitrary callbacks. A user-defined coefficient model, explicit pressure-anisotropy source, or new damping law may require a stronger bound.

For multiple species, the implementation uses the largest $`\widetilde{\sigma}_{0,\rm st}`$.



## 8. Idealized Units and Application Scope

The examples accompanying the method paper use dimensionless code units and are designed to isolate numerical and physical effects. They are not turnkey galaxy, wind, shock, or multiphase-ISM applications.

For a physical application, define consistent units for $`\mathcal{E}_{\rm cr}`$, $`\boldsymbol{F}_{\rm cr}`$, $`V_m`$, and all coupling coefficients.



## 9. Unresolved deviations in the CR flux perturbation.
In the CRMHD growing-wave test, the $`\mathcal{E}_{\rm cr}`$ perturbation agrees well with the linear eigenmode prediction. The CR flux perturbation, however, develops a localized spike, a small drift in its spatial mean. Although the extracted fundamental Fourier component remains close to the predicted mode, the full $`F_{\rm cr}`$ profile is not currently considered quantitatively validated. The origin of these deviations——finite-amplitude effects, the Riemann solver, implicit source integration, or operator splitting——has not yet been isolated. We will investigate this issue in future work.

We suggest that, users requiring accurate CR flux perturbations should perform amplitude, resolution, and timestep convergence tests.
