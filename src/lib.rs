// Macros for simple math operations: {min, max, etc}
#[macro_use]
extern crate math_macros;
pub use math_macros::*;

mod constants;
mod physics;

pub use constants::*;
pub use physics::{Planet, Universe};
pub use simulation::{Integrator, Simulation, System};

use anyhow::Result;

impl System for Universe {
    // Function that calculates the derivates of temperatures with respect to time
    fn derive(
        &mut self,
        time: f64,
        temperatures: &[f64],
        output: &mut [f64],
    ) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
        Ok(self.derive(time, temperatures, output)?)
    }

    // Function that calculates the derivates of temperatures with respect to time
    fn update(
        &mut self,
        time: f64,
        temperatures: &[f64],
    ) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
        self.update(time, temperatures);
        Ok(())
    }
}
