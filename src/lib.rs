#![feature(slice_as_chunks)]
#![feature(portable_simd)]
// Macros for simple math operations: {min, max, etc}
#[macro_use]
extern crate math_macros;
pub use math_macros::*;

mod constants;
mod physics;

pub use constants::*;
pub use physics::{Planet, Universe};
pub use sci_file::OutputFile;
pub use simulation::{InputConfig, Integrator, Simulation, System};

use anyhow::Result;
use serde_json::json;

/// Struct that contains data used by physics module and sink for the integrator solout.
#[derive(Debug)]
pub struct EnergyBalanceModel {
    /// [`Universe`] contains parameters and workspace used by physics module.
    pub data: Universe,
    /// Data output from the integrator is written into [`OutputFile`].
    pub output: OutputFile,
}

impl System for EnergyBalanceModel {
    type Data = Universe;
    type Output = OutputFile;
    fn new(output: OutputFile, data: Universe) -> Self {
        Self { data, output }
    }
    // Function that calculates the derivates of temperatures with respect to time
    fn derive(
        &mut self,
        time: f64,
        temperatures: &[f64],
        output: &mut [f64],
    ) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
        Ok(self.data.derive(time, temperatures, output)?)
    }
    // Function that outputs each step of the solution
    fn solout(
        &mut self,
        time: f64,
        temperatures: &[f64],
    ) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
        // Convert to JSON
        let j = json!({"time": time / SECONDS_IN_YEAR,
        "temperatures": &temperatures,
        });
        Ok(self.output.write_json_line(&j)?)
    }
}
