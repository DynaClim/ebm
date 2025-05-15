use anyhow::Result;
use ebm_lib::{SECONDS_IN_YEAR, Simulation, Universe};

use rayon::prelude::*;

fn main() -> Result<()> {
    // Launch each simuluation (input config file) in parallel.
    let simulations = Simulation::<Universe>::new()?;
    simulations
        .into_par_iter()
        .map(|mut simulation| {
            let mut initial_time = simulation.initial_time;
            let mut final_time = simulation.final_time;
            let initial_temperatures = simulation.system.initial_temperatures();

            if !simulation.resume {
                simulation.system.initialise()?;
                // Convert the time units from years to seconds for a new simulation.
                initial_time *= SECONDS_IN_YEAR;
                final_time *= SECONDS_IN_YEAR;
            }

            simulation.launch(initial_time, final_time, &initial_temperatures)?;

            Ok(())
        })
        .collect::<Result<()>>()?;

    Ok(())
}
