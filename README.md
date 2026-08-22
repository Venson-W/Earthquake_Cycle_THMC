# Earthquake Cycle THMC

A numerical earthquake-cycle simulator for investigating the development, evolution, and interaction of **shear localization and delocalization** in fault zones.

The simulator is intended for studying earthquake-cycle processes involving coupled thermo-hydro-mechanical-chemical (THMC) behavior.

## Overview

Fault deformation may transition between distributed and localized shear during an earthquake cycle. This project provides a computational framework for exploring these transitions and their effects on fault-zone mechanics.

The simulator can be used to investigate:

- Formation and evolution of localized shear zones
- Shear delocalization during interseismic or postseismic periods
- Earthquake nucleation and cyclic slip
- Stress and deformation evolution within fault zones
- Coupled thermal, hydraulic, mechanical, and chemical processes
- Sensitivity to material properties, loading conditions, and model parameters

## Physical Model

The model considers the evolution of fault-zone deformation under applied loading. A representative shear-stress relation may be written as

$$
\tau = \sigma_{\mathrm{eff}} f,
$$

where

- $\tau$ is the shear stress,
- $\sigma_{\mathrm{eff}}$ is the effective normal stress, and
- $f$ is the friction coefficient.

The effective normal stress is

$$
\sigma_{\mathrm{eff}} = \sigma_n - p,
$$

where $\sigma_n$ is the total normal stress and $p$ is the pore-fluid pressure.

A characteristic measure of shear strain is

$$
\gamma = \frac{\partial u}{\partial y},
$$

where $u$ is the displacement parallel to the fault and $y$ is the coordinate normal to the fault.

> The precise governing equations, constitutive laws, and numerical methods implemented by this repository should be documented in the corresponding source files and publications.

## Repository Structure

```text
Earthquake_Cycle_THMC/
├── README.md
├── LICENSE
├── src/            # Main source code
├── input/          # Model parameters and input files
├── examples/       # Example simulations
├── scripts/        # Preprocessing and postprocessing tools
├── tests/          # Verification and test cases
└── output/         # Simulation results (normally ignored by Git)
```

The actual directory structure may differ from this example.

## Requirements

Document the required compiler, interpreter, numerical libraries, and visualization tools here. For example:

- Supported operating system: Windows, Linux, or macOS
- Compiler or interpreter: `[specify name and version]`
- Required numerical libraries: `[specify dependencies]`
- Optional visualization tools: `[specify tools]`

## Installation

Clone the repository:

```bash
git clone git@github.com:YOUR_GITHUB_USERNAME/Earthquake_Cycle_THMC.git
cd Earthquake_Cycle_THMC
```

Then build or install the simulator using the appropriate command:

```bash
# Replace this example with the actual build command
make
```

## Running a Simulation

Prepare or select an input configuration and run the simulator:

```bash
# Replace this example with the actual execution command
./earthquake_cycle_thmc input/example.in
```

Simulation results should be written to a designated output directory.

## Input Parameters

Typical inputs may include:

| Parameter | Description | Unit |
|---|---|---|
| $\sigma_n$ | Applied normal stress | Pa |
| $p$ | Pore-fluid pressure | Pa |
| $\tau$ | Shear stress | Pa |
| $V$ | Slip or loading velocity | m/s |
| $G$ | Shear modulus | Pa |
| $\eta$ | Effective viscosity | Pa s |
| $T$ | Temperature | K |
| $t_{\mathrm{end}}$ | Final simulation time | s |

Update this table to match the parameters used by the code.

## Output

Depending on the selected model configuration, output may include:

- Displacement and velocity
- Shear strain and shear-strain rate
- Stress evolution
- Pore-fluid pressure
- Temperature
- State or internal variables
- Localization metrics
- Earthquake timing and recurrence intervals

Large output files should not normally be committed to Git. Add their directories or filename patterns to `.gitignore`.

## Examples

Example problems should be placed in the `examples/` directory. Each example should document:

1. Its physical purpose
2. Required input files
3. Execution command
4. Expected output
5. Relevant plotting or postprocessing commands

## Citation

The paper is in preparation.

## Contributing

Contributions, bug reports, and suggestions are welcome. For substantial changes, please open an issue first to discuss the proposed modification.

When contributing:

1. Create a new branch.
2. Make and test your changes.
3. Document any new parameters or functionality.
4. Submit a pull request with a clear description.

## License

This project is licensed under the MIT License. See the [`LICENSE`](LICENSE) file for details.

## Disclaimer

This software is provided for research purposes without warranty. Users are responsible for verifying that the numerical methods, parameter choices, and model assumptions are appropriate for their applications.

## Contact

For questions or collaboration inquiries, contact:

**Yuhan Wang**  
Swiss Seismological Service (SED) and ETH Zurich
yuhan.wang@sed.ethz.ch
