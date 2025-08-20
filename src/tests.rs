use super::*;
//use pretty_assertions::assert_eq;
use sci_file::{OutputWriter, read_json_from_file};
use simulation::{Integrator, simulation::InputConfig};
use std::path::Path;
use std::path::PathBuf;

fn build_simulation(config_path: PathBuf) -> InputConfig<Universe> {
    // Parse the config file.
    let mut config: InputConfig<Universe> = read_json_from_file(&config_path).unwrap();
    config.initial_time *= SECONDS_IN_YEAR;
    config.final_time *= SECONDS_IN_YEAR;

    config.system.initialise().unwrap();

    // Initial values for the integrator.
    let y = config.system.temperatures();
    config
        .integrator
        .initialise(config.initial_time, config.final_time, &y);

    config
}

fn test_simulation(config: PathBuf) -> Universe {
    let mut config = build_simulation(config);

    // Run the full integration.
    let (x, y) = config.integrator.integrate(&mut config.system).unwrap();
    config.system.update(x, y).unwrap();
    config.system
}

fn compare_or_create(path: impl AsRef<Path> + std::fmt::Display, result: &Universe) {
    match read_json_from_file::<Universe>(&path) {
        Ok(expected) => {
            // Saved file exists, compare the results.
            // We roundtrip our `Universe` through serde before comparison
            // to reset fields that are not serialized (serde skip_serializing)
            // (i.e. interpolation data read from file, internal buffers).
            let tmp = serde_json::to_string(&result).unwrap();
            let result: Universe = serde_json::from_str(&tmp).unwrap();
            assert_eq!(expected, result);
        }
        Err(err) => {
            match err {
                sci_file::Error::FileIo(_) => {
                    // Saved file does not exist save the results.
                    let mut writer = OutputWriter::new(&path).unwrap();
                    writer.write(&result).unwrap();
                    panic!("comparison file `{path}` did not exist, so it was created");
                }
                _ => {
                    dbg!(&err);
                    panic!(
                        "the comparison file `{path}` is corrupt or has invalid structure. if it contains 'null' values, the value was probably NaN or inifinity"
                    );
                }
            }
        }
    }
}

#[test]
fn example_300k_1year() {
    let result = test_simulation("examples/300K_1year.json.conf".into());
    compare_or_create("examples/300K_1year_expected.json", &result);
}

#[cfg(feature = "divan")]
#[divan::bench]
// This is a test case for performance benchmarking.
fn bench_300k_1year_145d(bencher: divan::Bencher) {
    // Load the config file and create the simulation/integrator once.
    let config = build_simulation("examples/300K_1year.json.conf".into());

    let system = config.system.clone();
    let integrator = config.integrator.clone();
    // The simulation mutates the integrator, so clone it for each iteration of the bench test.
    bencher.bench(|| integrator.clone().integrate(&mut system.clone()));
}
