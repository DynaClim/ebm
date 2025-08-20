#![feature(portable_simd)]
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

impl System<f64> for Universe {
    // Function that calculates the derivates of temperatures with respect to time
    fn derive(
        &mut self,
        time: f64,
        temperatures: &[f64],
        output: &mut [f64],
    ) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
        self.derive(time, temperatures, output)?;
        Ok(())
    }

    fn update(
        &mut self,
        time: f64,
        temperatures: &[f64],
    ) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
        self.update_for_output(time, temperatures);
        Ok(())
    }
}

#[cfg(any(test, feature = "divan"))]
mod tests;
