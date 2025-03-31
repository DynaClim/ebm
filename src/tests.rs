use super::*;
//use pretty_assertions::assert_eq;
use sci_file::{OutputFile, deserialize_json_from_path};
use serde::{Deserialize, Serialize};
use simulation::{InputConfig, Integrator, System};
use std::path::Path;
use std::path::PathBuf;

#[derive(Clone)]
struct Test {
    pub data: Universe,
}

// Mock System implementation without writing any output.
impl System for Test {
    type Data = Universe;
    type Output = OutputFile;

    fn new(_: OutputFile, data: Universe) -> Self {
        Self { data }
    }

    fn derive(
        &mut self,
        time: f64,
        y: &[f64],
        dy: &mut [f64],
    ) -> Result<(), Box<(dyn std::error::Error + Send + Sync + 'static)>> {
        self.data.derive(time, y, dy)?;
        Ok(())
    }

    fn solout(
        &mut self,
        _time: f64,
        _y: &[f64],
    ) -> Result<(), Box<(dyn std::error::Error + Send + Sync + 'static)>> {
        Ok(())
    }
}

fn build_simulation(config_path: PathBuf) -> InputConfig<Universe> {
    // Parse the config file.
    let mut config: InputConfig<Universe> = deserialize_json_from_path(&config_path).unwrap();
    config.initial_time *= SECONDS_IN_YEAR;
    config.final_time *= SECONDS_IN_YEAR;

    config.universe.initialise().unwrap();

    // Initial values for the integrator.
    let y = config.universe.initial_temperatures();
    config
        .integrator
        .initialise(config.initial_time, config.final_time, &y);

    config
}

#[allow(dead_code)]
fn test_simulation(config: PathBuf) -> Vec<f64> {
    let mut config = build_simulation(config);

    let mut system = Test {
        data: config.universe,
    };

    // Run the full integration.
    let stats = config.integrator.integrate(&mut system).unwrap();
    dbg!(stats);

    // Collect the final y values.
    config.integrator.y().to_vec()
}

// Float roundtrip in only supported in JSON (not CSV)
// So store the expected output as JSON with one field.
#[derive(Deserialize, Serialize)]
struct Compare {
    expected: Vec<f64>,
}

#[allow(dead_code)]
fn compare_or_create(path: impl AsRef<Path> + std::fmt::Display, result: &[f64]) {
    if let Ok(saved_data) = deserialize_json_from_path::<Compare>(&path) {
        // Saved file exists, compare the results.
        assert_eq!(saved_data.expected, result);
    } else {
        // Saved file does not exist, save the results.
        let mut writer = OutputFile::new(&path).unwrap();
        let result = Compare {
            expected: result.to_vec(),
        };
        writer.write_json(&result).unwrap();
        panic!("comparison file `{path}` did not exist, so it was created");
    }
}

#[test]
fn example_300k_1year() {
    let result = test_simulation("examples/300K_1year.json.conf".into());
    compare_or_create("examples/300K_1year.expected", &result);
}

#[cfg(feature = "divan")]
#[divan::bench]
// This is a test case for performance benchmarking.
fn bench_300k_1year_145d(bencher: divan::Bencher) {
    // Load the config file and create the simulation/integrator once.
    let config = build_simulation("examples/300K_1year.json.conf".into());

    let system = Test {
        data: config.universe,
    };

    // The simulation mutates the integrator, so clone it for each iteration of the bench test.
    bencher.bench(|| config.integrator.clone().integrate(&mut system.clone()));
}
