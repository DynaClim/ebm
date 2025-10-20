<img src="images/SNSF.png" align="left"> <img src="images/AGPLv3.png" align="right">

# **E**nergy **B**alance **M**odel

**ebm** is a simple numerical simulator written in Rust that models the balance between incoming energy from an external source (e.g. the Sun) and outgoing energy from a reservoir (e.g. a planet) over a specified time period. It is used to study long-term evolution of planetary climates.


Features:

 - 1-D model that resolves the evolution of surface temperature as a function of planetary latitude and time.

 - Meridional heat transport: Energy redistribution is parameterised as temperature-driven diffusion, capturing large-scale atmospheric transport processes. The diffusion coefficient is linked to the planetary rotation rate (see [Williams & Kasting, 1997](https://doi.org/10.1006/icar.1997.5759) and [Dressing et al. 2010](https://doi.org/10.1088/0004-637X/721/2/1295)).

 - Variable resolution: The number of latitudes to calculate climate changes for is specified by the grid size.

 - Variable surface planet conditions: planet composition (aqua planet, land planet, hybrid) and surface albedo (temperature and spectrum-dependent, following [Spiegel et al., 2008](https://doi.org/10.1086/588089), Dressing et al. 2010 and [Gilmore, 2014](https://doi.org/10.1093/mnras/stu302)).

 - Variable atmospheric conditions: atmospheric composition can be specified through the temperature-dependent outgoing longwave radiation profile, following [Turbet et al., 2022](https://doi.org/10.1038/s41586-021-03873-w) and [Chaverot et al., 2023](https://doi.org/10.1051/0004-6361/202346936).

 - Variable orbital parameters: eccentricity and obliquity can be constant, interpolated from a profile, or derived from sinusoidal excitation.

The physics has been validated by the following contributors:

 - [Siddharth Bhatnagar, UniGE](https://www.unige.ch/sciences/astro/exoplanets/en/team/scientific-collaborators/bhatnagar-siddharth)
 - [Marine Leyvraz, UniBE](https://www.climate.unibe.ch/about_us/team/leyvraz_marine/index_eng.html)
 - [Emeline Bolmont, UniGE](https://www.emelinebolmont.com/)
 - [Maura Brunetti, UniGE](https://www.unige.ch/gap/nonlinear/people/brunetti-maura)
 - [Jérôme Kasparian, UniGE](https://www.unige.ch/gap/nonlinear/people/jerome-kasparian)

## Installation
### Requirements

- Rust: [see rustup](https://www.rustup.rs/)

Clone the repository, build, and install:

```bash
cargo install --path .
```
The executable will be copied into `$HOME/.cargo/bin/`. You may need to add this directory to your `$PATH`.

## Quickstart Example
```
# Create simulation case(s) into the "my_cases" directory.
python3 scripts/setup.py my_cases

# Launch all simulations from "my_cases" directory, putting results into "my_output" directory.
ebm -b --output-format jsonl my_cases my_output
```

## Usage

Initial conditions are specified into configuration files (`.conf`) as JSON, which is read by `ebm` to start a simulation. Output data from the simulation is in JSONL format. A simple python script is included to help in generating the JSON initial conditions.

### Create a JSON case
Set the planet properties into the configuration file(s) to the desired values and create the JSON initial conditions file(s) into the `config_files` directory.

`python3 scripts/setup.py config_files`

This generates a JSON case (.conf) for each set of input conditions into the designated directory. The total number of simulations will be the combinatorial product of the variables of the planet.

### Start a simulation
Simulations can be launched individually, or in batch mode. The user must specify the location of the input file(s) (`.conf`) and the desired output location (`output_directory`).

#### Start a single simulation
`ebm input_file.conf output_directory --output-format jsonl`

With input arguments:

- `input_file.conf`: Path to the desired JSON input configuration file (e.g. `./relative/path/to/config.conf` or `/home/$USER/simulations/example.conf`).
- `output_directory`: Path to the desired output destination. Will be automatically created if it does not already exist.

#### Start multiple simulations (batch mode)
Batch mode launches simulations in parallel.

`ebm -b input_directory output_directory --output-format jsonl`

With input arguments:

- `input_directory`: Path to the directory containing one or more input configuration files (`.conf`).
- `output_directory`: Path to the desired output destination. Will be automatically created if it does not already exist.

### Output
Launching a simulation produces three output files:

- `input_file.conf`: A copy of the input configuration file used to launch this simulation (JSON).
- `simulation.log`: A log file containing status information for each simulation. Will be appended to if it already exists.
- `input_file.jsonl`: Output for the simulation (JSONL).

For example:
```bash
ebm example.conf simulations --output-format jsonl
```
will produce the following structure:
(A `run_n` sub directory will be created for each simulation, with `n` as the smallest non pre-existing numerical suffix.)
```
simulations/
└── simulation.log
└── run_0
    ├── example.conf
    └── example.jsonl
```

### Analyse a simulation
The output data is available in `/output/run_n/temperature.jsonl`

## Design and Structure

- constants.rs: Physical constants
- lib.rs: Specifies the structure of the `Universe` for the config file and implements the `Integrator::System` trait required by the integrator.
- main.rs: Robust example of using the ebm `lib` crate.
- physics.rs: Public interface is through the `Universe` struct, which provides a `derive()` method used at every iteration of the integrator (`time` and `temperatures`), to solve the ODEs.

## Benchmarking

Benchmarking is handled by divan:

```bash
cargo bench --features divan
```

## Notes

This work has been carried out within the framework of the NCCR PlanetS supported by the Swiss National Science Foundation under grants 51NF40_182901 and 51NF40_205606.
