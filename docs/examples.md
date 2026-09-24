# Quick Start and Examples

This page describes selected examples from the method paper, with build instructions, input-file locations, run commands, and simple plots. Problem generators are in `athena_CRMHD/src/pgen`; the example inputs are in `athena_CRMHD/inputs/CR`.


## 1. General Workflow

Athena++ compiles one problem generator into the executable at a time. For each example:

1. Configure the source tree with the required `--prob` and `--ncrs`.
2. Rebuild Athena++.
3. Create a separate run directory and copy the matching input file into it.
4. Run the executable with that input file.
5. Check the log for fatal errors, NaNs, and unexpectedly small timesteps.
6. Inspect the output and compare it with the expected behavior.

### 1.1 Set the paths

In your terminal, replace the two paths below with your local paths. Run this block once and keep using the same terminal for the examples:

```bash
cd /path/to/athena_CRMHD
CR_SOURCE="$PWD"
HDF5_PREFIX=/path/to/hdf5
```

`HDF5_PREFIX` is the installation prefix containing the HDF5 headers and libraries. The examples require a C++ compiler, `make`, Python, and an HDF5 installation. Compilation commands return to `CR_SOURCE` so that the working directory is unambiguous.

Each run uses its own directory under `run/`. For another configuration or a repeat run, choose a new directory or archive the old outputs first. Running again in the same directory can overwrite output files.

### 1.2 Analytical comparisons from the command line

The supplied `vis/python/plot_cr_examples.py` plots Athena++ output as colored solid lines and the analytical predictions as black dashed lines. It reads each snapshot's actual time and cell-center coordinates and rejects an initial CR-energy profile that does not match the selected reference; you do not need to edit the plotting code.

Install the plotting dependencies in your Python environment:

```bash
python -m pip install numpy matplotlib h5py
python "$CR_SOURCE/vis/python/plot_cr_examples.py" --help
```

Run each plot command from the corresponding run directory. The mode name selects the expected HDF5 filename pattern, output PNG name, and comparison times automatically. Optional filename and time arguments remain available through `--help`. The reference formulas assume the benchmark parameters stated in each section; if you change those parameters, the dashed curves are no longer the corresponding prediction.

All three example inputs use Athena++'s built-in CR output fields:

| Field | Quantity |
|---|---|
| `x1v` | Cell-center coordinate along the one-dimensional domain |
| `cr_1_energy_dens` | CR energy density $\mathcal{E}_{\rm cr}$ for species 1 |
| `cr_1_flux1` | Stored CR flux component $F_{{\rm cr},1}/V_m$ for species 1 |
| `rho` | Gas density |

The plotter multiplies `cr_1_flux1` by $V_m$ when displaying the physical CR flux.

Users may define custom output variables through `user_out_var` in the problem generators, allocate them with `AllocateUserOutputVariables(output numbers)`, and enable their output by setting `variable = uov` in the corresponding `<output>` block of the input file.

## 2. 1D CR Streaming

### 2.1 Purpose

`CR_streaming.cpp` initializes a triangular CR-energy profile in a uniform, initially static, magnetized gas. The CR pressure gradient drives the CR streaming instability (CRSI); a scattering-coefficient closure represents the saturated state of the unresolved waves.

The intended CR-only solution develops:

- a central plateau in $\mathcal{E}_{\rm cr}$;
- outward-moving fronts on both sides;
- symmetry about $x=0$ for zero background gas velocity.

This example is intended as a first transport test: it combines the CR state, built-in streaming coefficients, HLLE transport, implicit scattering, output, and the streaming timestep constraint in one dimension. Reproducing the paper's analytic solution also requires disabling gas feedback and suppressing the CR-energy scattering source.



### 2.2 Build

```bash
cd "$CR_SOURCE"
python configure.py --prob=CR_streaming -b -hdf5 --ncrs=1 \
  --hdf5_path="$HDF5_PREFIX"
make clean
make -j
```

### 2.3 Run

Copy the supplied input into a fresh run directory:

```bash
mkdir -p "$CR_SOURCE/run/CR_streaming"
cd "$CR_SOURCE/run/CR_streaming"
cp "$CR_SOURCE/inputs/CR/athinput.CRstreaming" athinput.CRstreaming
```

The distributed input requests frequent HDF5 output. Override those two settings and run to $t=0.30$ without changing the input file:

```bash
"$CR_SOURCE/bin/athena" -i athinput.CRstreaming \
  job/problem_id=CR_streaming time/nlim=-1 time/tlim=0.30 \
  output2/dt=0.15 > streaming.log 2>&1
tail -n 20 streaming.log
```

Check that the log ends at $t=0.30$ and that three `CR_streaming.out2.*.athdf` snapshots were written.

### 2.4 Plot

From the run directory, create a two-panel CR energy density/CR energy density flux comparison with one command:

```bash
python "$CR_SOURCE/vis/python/plot_cr_examples.py" streaming
open streaming_comparison.png
```

For $`\widetilde{\sigma}_{0,\rm st}=1`$ and $`v_A=1`$, the analytical front speed is $`v_{\rm st}=5/3`$. The initial CR energy peak streams away from the central maximum, flattening the central gradient into a plateau that broadens with time while the outer fronts move outward. Accordingly, the flux is negative on the left, positive on the right, and small within the plateau. Because the two-moment method evolves the flux dynamically, the numerical solution shows a smooth central transition rather than the sharp zero-flux region of the idealized reference.

Use CR energy as the quantitative check and flux as a qualitative sign-and-shape check.

## 3. 1D CR Streaming with Nonlinear Landau Damping

### 3.1 Purpose

CR–gas coupling depends on both wave driving and wave damping, which together determine the saturated wave amplitude. In this model, the scattering coefficients represent those microscopic processes.

`CR_streaming_nlldamp.cpp` initializes $`\mathcal{E}_{\rm cr}(x,0)=1/(10+3|x|)`$ and registers `MyCRScatteringCoefficient` as a user-defined nonlinear-Landau-damping closure. The callback evaluates the field-aligned CR-energy gradient and activates the forward- or backward-wave coefficient according to its sign, while setting the opposite branch, the Lorentz coefficient, and `deltaPcr` to zero. Active cells and the first ghost layer use centered gradients; the outer ghost layer copies the nearest computed coefficients to avoid out-of-bounds differencing.

### 3.2 Build

```bash
cd "$CR_SOURCE"
python configure.py --prob=CR_streaming_nlldamp -b -hdf5 --ncrs=1 \
  --hdf5_path="$HDF5_PREFIX"
make clean
make -j
```

### 3.3 Run

Copy the supplied input into a fresh run directory:

```bash
mkdir -p "$CR_SOURCE/run/CR_streaming_nll"
cd "$CR_SOURCE/run/CR_streaming_nll"
cp "$CR_SOURCE/inputs/CR/athinput.CRstreaming_nll" athinput.CRstreaming_nll
```

```bash
"$CR_SOURCE/bin/athena" -i athinput.CRstreaming_nll \
  job/problem_id=CR_streaming_nll > streaming_nll.log 2>&1
tail -n 20 streaming_nll.log
```

### 3.4 Plot

The nonlinear-damping reference solution uses an outward front speed $v_{\rm st}=11/6$. Its central plateau edge is found from an implicit scalar equation, and the outer flux is $\pm v_{\rm st}\mathcal{E}_{\rm cr}$. Create the two-panel comparison at $t=0,0.05,0.10$ with:

```bash
python "$CR_SOURCE/vis/python/plot_cr_examples.py" landau
open landau_comparison.png
```

The plotter first verifies that the numerical initial condition matches the analytical profile. As the CRs stream outward, the central energy gradient is flattened and a widening plateau develops. The analytical CR-energy profile provides the primary quantitative comparison. The flux reference assumes a quasi-steady state and is reliable mainly outside the plateau; near the symmetry point, where this approximation breaks down, the numerical flux is expected to vary smoothly rather than follow the sharp piecewise curve.

## 4. Coupled CR–MHD Waves

### 4.1 Purpose

`CRMHD_wave.cpp` initializes a small-amplitude eigenmode of the coupled CR–MHD system. It tests CR transport, implicit scattering, gas feedback, and gas dynamics against a linear-theory prediction.

For further details on these coupled waves, including their linear theory, dispersion relation, and eigenmode structure, see Appendix A of our [method paper](https://doi.org/10.3847/1538-4365/ae3e7f).

### 4.2 Build

```bash
cd "$CR_SOURCE"
python configure.py --prob=CRMHD_wave -b -hdf5 --ncrs=1 -h5double \
  --hdf5_path="$HDF5_PREFIX"
make clean
make -j
```

The problem parameters are:

| Parameter | Meaning |
|---|---|
| `vA` | Background Alfvén speed |
| `Cs` | Gas sound speed |
| `Ccr` | CR characteristic speed that sets the background CR energy density, $`\mathcal{E}_{\rm cr,0}=9\rho_0 C_{\rm cr}^2/4`$ |
| `u0` | Background gas velocity |
| `Rew` | Real part of the wave angular frequency, normalized by $`kv_A`$, required to initialize the eigenvector |
| `Imw` | Imaginary part of the wave angular frequency, normalized by $`kv_A`$, required to initialize the eigenvector; negative value for damping and positive value for growing |
| `delta_rho_over_rho` | Initial fractional gas-density perturbation; defaults to $`10^{-3}`$ |

The generator uses `Rew` and `Imw` in a dimensionless eigenvector normalization (normalized by $`kv_A`$). Do not interpret them directly as angular frequencies in inverse code-time units without checking that normalization.

The method paper discusses four wave branches in different regimes. The table lists their representative input values:

| Case | δρ/ρ<sub>0</sub> | <i>v</i><sub>A</sub> | <i>C</i><sub>s</sub> | <i>V</i><sub>m</sub> | <i>L</i><sub>box</sub> | σ | `Rew` | `Imw` | Interpretation |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|:--|
| CR-MHD acoustic | 10<sup>−3</sup> | 0.01 | 0.1 | 100 | 99.615 | 0.1 | 100.186 | −10.487 | Acoustic mode for the coupled CR–MHD fluid |
| Damping CR-modified acoustic | 10<sup>−3</sup> | 0.01 | 0.1 | 100 | 1.198 | 0.1 | 9.596 | −2.718 | CR-modified MHD acoustic mode |
| Growing CR-modified acoustic | 10<sup>−3</sup> | 0.1 | 0.01 | 1000 | 1.079 | 0.1 | 0.322 | 0.142 | CR-modified MHD acoustic mode |
| Rapid damping | 10<sup>−8</sup> | 0.1 | 0.01 | 1000 | 0.628 | 0.1 | −1.007 | −99693.278 | Strongly damped branch |


Each row provides the case-specific parameters needed to initialize the corresponding wave. Use these values with the common numerical settings in `athinput.CRMHDwave`. In the input file, set `Vmax` to $`V_m`$, choose the domain length $`x_{1,\max}-x_{1,\min}=L_{\rm box}`$, and split the listed forward-wave scattering coefficient as `sigma_pL_1` $=$ `sigma_pR_1` $=$ $\sigma/2$, with `sigma_mL_1` $=$ `sigma_mR_1` $= 0$.

### 4.3 Run

Use the **growing CR-modified acoustic branch** as the single wave example. The supplied `inputs/CR/athinput.CRMHDwave` is a current-format growing-mode input. Copy it into a fresh run directory:

```bash
mkdir -p "$CR_SOURCE/run/CRMHD_growing"
cd "$CR_SOURCE/run/CRMHD_growing"
cp "$CR_SOURCE/inputs/CR/athinput.CRMHDwave" athinput.CRMHDwave
```

Run a shorter, 128-cell comparison at $t=0,1,2$. The overrides supply the benchmark's full-precision eigenmode values and set $\sigma_{+,1}^L=\sigma_{+,1}^R=0.05$ (total $\sigma_{+,\rm species1}=0.1$):

```bash
"$CR_SOURCE/bin/athena" -i athinput.CRMHDwave \
  job/problem_id=CRMHD_growing mesh/nx1=128 mesh/x1max=1.079007753048546 \
  problem/Rew=0.3218887947075146 problem/Imw=0.14203248351814002 \
  cosmic_ray/sigma_pL_1=0.05 cosmic_ray/sigma_pR_1=0.05 \
  time/tlim=2.0 output2/dt=1.0 > growing.log 2>&1
tail -n 20 growing.log
```

Check that the log ends at $t=2$ and that three `CRMHD_growing.out2.*.athdf` snapshots were written. This shortened run is a quick comparison, not a full paper-resolution reproduction.

### 4.4 Plot

Plot $`\mathcal{E}_{\rm cr}-\langle\mathcal{E}_{\rm cr}\rangle`$ against the growing eigenmode:

```bash
python "$CR_SOURCE/vis/python/plot_cr_examples.py" growing
open growing_comparison.png
```

The black dashed curves show the linear eigenmode prediction $`\mathrm{Re}[{\delta\mathcal{E}}_{\rm cr}e^{i(kx-\omega t)}]`$, with $`k=2\pi/L_{\rm box}`$ and $`\omega=kv_A`$(`Rew`+$`i`$`Imw`). The real part of $\omega$ determines the phase propagation, while a positive imaginary part produces exponential growth. At each output time, the plotter subtracts the numerical spatial mean $\langle\mathcal{E}_{\rm cr}\rangle$, isolating the wave perturbation from the uniform background shift. Agreement in phase and amplitude therefore tests both the eigenmode initialization and its subsequent evolution in the coupled CR–MHD system.

## 5. Problem Generators from the Method Paper

We list below the CR problem generators in `athena_CRMHD/src/pgen`. Inclusion in this catalog does not imply that a complete current-format input and reproduction workflow are available for every generator.

| Problem generator | Geometry | Main test | Expected comparison |
|---|---|---|---|
| `CR_streaming.cpp` | 1D Cartesian | Linear-damping streaming closure | Triangular-profile analytic solution |
| `CR_streaming_nlldamp.cpp` | 1D Cartesian | Streaming with nonlinear Landau damping | Nonlinear-damping reference solution |
| `CR_diffusion.cpp` | 1D Cartesian | diffusion-like homogenious scattering | Gaussian diffusion solution |
| `CR_circB.cpp` | 2D Cartesian | Transport along circular field lines | Field-aligned transport and symmetry |
| `CR_streaming_sph.cpp` | Spherical polar | Radial streaming in curvilinear geometry | Spherical analytic scaling |
| `CRaniso_cylin.cpp` | Cylindrical | CRPAI with prescribed anisotropy | Bessel-function solution |
| `CRMHD_wave.cpp` | 1D Cartesian | Coupled CR–MHD waves | Analytic eigenmode predictions |

The curvilinear tests in the paper use monopolar magnetic-field configurations and manually freeze MHD evolution. Those setups require problem-specific code changes beyond the standard input interface, so command-line reproduction instructions are not included here.

## 6. General Athena++ Functionality

This guide focuses on the CR module and its accompanying examples. For general Athena++ functionality, including boundary-condition configuration, MPI and OpenMP parallelism, mesh refinement, coordinate systems, output and restart files, and other framework features, please consult the [original Athena++ project](https://github.com/PrincetonUniversity/athena). The CR module uses the standard Athena++ infrastructure, so the corresponding Athena++ procedures apply unless this documentation explicitly states otherwise.
