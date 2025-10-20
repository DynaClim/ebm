"""
Creates an initial conditions file for input into EBM.

For the quantities of eccentricity and obliquity,
the source must be specified as one of the following types:

Constant: Value does not vary with time.
Example:
"obliquity": {
    "Constant": 0.0
}


Damping: Value decreases with time.
Example:
"obliquity": {
    "Damping": {
        # time it takes for the value to linearly decrease from initial to 0 (year)
        "damping_timescale": 1.,
        # starting value
        "initial": 1.
    }
}

Sin: Value excited following a sine time evolution.
Example:
"obliquity": {
    "Sin": {
        # Amplitude of the oscillation of the value, 0 < amplitude < 1, no unit
        "amplitude": 1.,
        # Period of the oscillation of the value (Earth days)
        "period": 1.,
        # starting value
        "initial": 1.
    }
}

Interpolated: Value interpolated by time from a precomputed grid. Arrays can be any length, but must be of equal length.
Example:
"obliquity": {
    "Interpolate": {
        # time points
        "x_vals": [1., 2., 3.],
        # value points
        "y_vals": [4., 5., 6.],
    }
}

"""

import json
import sys
import itertools
from pathlib import Path

TEMPERATURE_GRID_DIMENSION = 146

"""
Earth like planet
A simulation input config file is generated for all combinations of the values of the planet.
"""


def planet_base():
    return {
        # creates one simulation for each temperature value, with a uniform distribution.
        # K
        "initial_temperatures": uniform_temperatures([300, 100]),
        # alternatively, construct your own custom temperature profiles:
        # "initial_temperatures": [[100., 200., 300.], [100., 100., 100.]],
        "albedo_low_temperature": [0.77],
        "albedo_high_temperature": [0.28],
        "temperature_centre_transition": [268.0],
        "temperature_width_transition": [5.0],
        "insolation_factor": [1.0, 2.0],
        "fiducial_diffusion_coefficient": [0.5394],
        "eccentricity": [
            *constant([0.0167]),
            # Can add interpolate, damping, sin, etc structures here
        ],
        "obliquity": [
            *constant([23.44]),
            # Can add interpolate, damping, sin, etc structures here
        ],
        "orbital_period": [365.0],
        # (land_fraction, water_faction): ensure land_fraction + water_fraction == 1.
        "land_water_fraction": [(0.3, 0.7)],
    }


def config_base(system):
    return {
        "resume": False,
        # years
        "initial_time": 0.0,
        "final_time": 10.0,
        "system": system,
        "integrator": {
            "Dopri5": {
                # Controller for selecting adaptive step size.
                "step_size_controller": {
                    # Relative tolerance for computation of adaptive step size.
                    "relative_tolerance": 1.0e-14,
                    # Absolute tolerance for computation of adaptive step size.
                    "absolute_tolerance": 1.0e-14,
                    # Initial step size. Will be computed automatically if set to 0 when starting the integrator.
                    "step_size": 0.0,
                    # Minimum factor between successive steps.
                    "step_size_factor_min": 1.0 / 3.0,
                    # Maximum factor between successive steps.
                    "step_size_factor_max": 6.0,
                    # Safety factor for computation of adaptive step size.
                    "step_size_error_factor": 0.9,
                    # Maximum step size. Will be computed automatically if set to 0 when starting the integrator.
                    "step_size_max": 0.0,
                    # Coefficient of step size calculation.
                    "alpha": 1.0 / 5.0,
                    # Coefficient of step size calculation.
                    # Step size control stability increases with 0 < beta <= 0.04
                    "beta": 0.0,
                },
                # Determines whether to test for stiffness as a criteria for aborting.
                "stiffness_test": "Disabled",
                # Determines whether to test for a maximum number of integration steps as a criteria for aborting.
                "max_integration_steps": None,
                # Dense or Sparse solution output.
                # Determines whether to output at every accepted timestep, or only at specific interval (interpolated from actual steps).
                "solution_output": {
                    # "Sparse"
                    "Dense": {
                        # Interval to output at, in seconds (simulation time)
                        # E.g. once per day would be 86400.
                        "increment": 60 * 60 * 24 * 30
                    }
                },
            }
        },
    }


#### No need to edit below here #####

OLR_INTERPOLATOR = {
    "x_vals": [
        267.63502454991800,
        274.5090016366610,
        280.2373158756140,
        286.53846153846200,
        292.26677577741400,
        296.8494271685760,
        300.5728314238950,
        306.0147299509000,
        311.7430441898530,
        315.75286415711900,
        320.0490998363340,
        324.3453355155480,
        327.7823240589200,
        334.0834697217680,
        342.10310965630100,
        350.69558101473000,
        358.14238952536800,
        363.58428805237300,
        368.4533551554830,
        373.03600654664500,
        379.3371522094930,
        387.3567921440260,
        398.24058919803600,
        408.5515548281510,
        524.8363338788870,
        539.1571194762690,
        557.4877250409170,
        581.5466448445170,
        602.7414075286420,
        629.0916530278230,
        653.1505728314240,
        670.9083469721770,
        696.1129296235680,
        719.0261865793780,
        743.6579378068740,
        762.5613747954170,
        776.8821603927990,
        787.1931260229130,
        893.1669394435350,
        898.8952536824880,
        904.0507364975450,
        938.9934533551560,
        988.2569558101480,
        1055.8510638297900,
    ],
    "y_vals": [
        218.53317102860600,
        233.44491783323200,
        246.65246500304300,
        261.56421180766900,
        275.1978088861840,
        285.84905660377400,
        293.0919050517350,
        298.63055386488100,
        304.169202678028,
        308.42970176506400,
        312.47717589774800,
        308.42970176506400,
        304.169202678028,
        294.79610468654900,
        281.1625076080340,
        264.9726110772980,
        249.63481436396800,
        238.98356664637900,
        231.31466828971400,
        226.20206938527100,
        220.66342057212400,
        216.40292148508800,
        212.35544735240400,
        208.7340231284240,
        208.7340231284240,
        208.7340231284240,
        208.7340231284240,
        208.30797321972000,
        208.30797321972000,
        207.8819233110160,
        207.8819233110160,
        207.45587340231300,
        207.45587340231300,
        207.24284844796100,
        207.45587340231300,
        207.45587340231300,
        207.45587340231300,
        207.45587340231300,
        242.3919659160070,
        250.91296409007900,
        258.5818624467440,
        293.51795496043800,
        406.42118076689000,
        559.7991479001830,
    ],
}


def make_config_files(output_path):
    planet = planet_base()
    # All planets created from the combinatorial expansion of all planet values.
    combis = [x for x in itertools.product(*planet.values())]
    for i, planet_values in enumerate(combis):
        # Create the system for integration from the planet values.
        system = make_system(*planet_values)
        # If the temperature profile is uniform, include the starting value in the initial conditions filename.
        temperature = (
            int(system["temperatures"][0])
            if len(set(system["temperatures"])) == 1
            else "mixed"
        )
        # Create the initial conditions structure.
        config = config_base(system)
        duration = config["final_time"]
        # include the initial temperature and simulation duration in the initial conditions filename.
        config_name = f"{output_path}/{temperature}K_{int(duration)}year_{i}.conf"
        print(f"Making config: {config_name}")
        with open(config_name, "x") as f:
            f.write(json.dumps(config, indent=4))


def make_system(
    initial_temperatures,
    albedo_low_temperature,
    albedo_high_temperature,
    temperature_centre_transition,
    temperature_width_transition,
    insolation_factor,
    fiducial_diffusion_coefficient,
    eccentricity,
    obliquity,
    orbital_period,
    land_water_fraction,
):
    """Creates the system dictionary from provided values."""
    return {
        "temperatures": initial_temperatures,
        "planet": {
            "albedo_low_temperature": albedo_low_temperature,
            "albedo_high_temperature": albedo_high_temperature,
            "temperature_centre_transition": temperature_centre_transition,
            "temperature_width_transition": temperature_width_transition,
            "insolation_factor": insolation_factor,
            "fiducial_diffusion_coefficient": fiducial_diffusion_coefficient,
            "eccentricity": eccentricity,
            "obliquity": obliquity,
            "outgoing_longwave_radiation": OLR_INTERPOLATOR,
            "orbital_period": orbital_period,
            "land_fraction": land_water_fraction[0],
            "water_fraction": land_water_fraction[1],
        },
    }


def main():
    if len(sys.argv) == 2:
        output_path = sys.argv[1]
        # Create the output directory path if it doesn't already exist.
        Path(output_path).mkdir(parents=True, exist_ok=True)

        make_config_files(output_path)
    else:
        print("usage: python3 setup.py output_dir")


def uniform_temperatures(temperatures):
    """Creates a 1D grid of temperatures for each provided value."""
    return [[temp for x in range(TEMPERATURE_GRID_DIMENSION)] for temp in temperatures]


def constant(values):
    """Creates an instance of a constant planet value (obliquity, eccentricity, etc) for each provided value."""
    out = []
    for val in values:
        out.append({"Constant": val})
    return out


if __name__ == "__main__":
    main()
