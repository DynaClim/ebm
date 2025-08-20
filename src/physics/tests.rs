use super::*;
//use pretty_assertions::assert_eq;

fn test_universe(system_size: usize) -> Universe {
    let planet = test_planet();
    let mut universe = Universe::new();
    universe.planet = planet;
    universe.temperatures = vec![300.; system_size];
    universe.initialise().unwrap();
    universe
}

fn test_planet() -> Planet {
    let x_vals = vec![
        267.635024549918,
        274.509001636661,
        280.237315875614,
        286.538461538462,
        292.266775777414,
        296.849427168576,
        300.572831423895,
        306.0147299509,
        311.743044189853,
        315.752864157119,
        320.049099836334,
        324.345335515548,
        327.78232405892,
        334.083469721768,
        342.103109656301,
        350.69558101473,
        358.142389525368,
        363.584288052373,
        368.453355155483,
        373.036006546645,
        379.337152209493,
        387.356792144026,
        398.240589198036,
        408.551554828151,
        524.836333878887,
        539.157119476269,
        557.487725040917,
        581.546644844517,
        602.741407528642,
        629.091653027823,
        653.150572831424,
        670.908346972177,
        696.112929623568,
        719.026186579378,
        743.657937806874,
        762.561374795417,
        776.882160392799,
        787.193126022913,
        893.166939443535,
        898.895253682488,
        904.050736497545,
        938.993453355156,
        988.256955810148,
        1055.85106382979,
    ];
    let y_vals = vec![
        218.533171028606,
        233.444917833232,
        246.652465003043,
        261.564211807669,
        275.197808886184,
        285.849056603774,
        293.091905051735,
        298.630553864881,
        304.169202678028,
        308.429701765064,
        312.477175897748,
        308.429701765064,
        304.169202678028,
        294.796104686549,
        281.162507608034,
        264.972611077298,
        249.634814363968,
        238.983566646379,
        231.314668289714,
        226.202069385271,
        220.663420572124,
        216.402921485088,
        212.355447352404,
        208.734023128424,
        208.734023128424,
        208.734023128424,
        208.734023128424,
        208.30797321972,
        208.30797321972,
        207.881923311016,
        207.881923311016,
        207.455873402313,
        207.455873402313,
        207.242848447961,
        207.455873402313,
        207.455873402313,
        207.455873402313,
        207.455873402313,
        242.391965916007,
        250.912964090079,
        258.581862446744,
        293.517954960438,
        406.42118076689,
        559.799147900183,
    ];

    let mut interpolator = Interpolator::<f64>::new();
    interpolator.init(&x_vals, &y_vals);
    Planet {
        // will be computed during initialisation
        heat_capacity_reciprocal: 0.0,

        albedo_low_temperature: 0.77,
        albedo_high_temperature: 0.28,
        temperature_centre_transition: 268.0,
        temperature_width_transition: 5.0,
        insolation_factor: 1.0,
        fiducial_diffusion_coefficient: 0.5394,
        eccentricity: Source::Constant(0.0167),
        obliquity: Source::Constant(23.44),
        outgoing_longwave_radiation: interpolator,
        orbital_period: 365.0,
        land_fraction: 0.3,
        water_fraction: 0.7,
        ..Default::default()
    }
}

#[test]
fn test_derive() {
    let system_size = 5;
    let mut universe = test_universe(system_size);
    let temperatures = vec![300.0; system_size];
    let expected = vec![300.; system_size];
    let mut result = vec![0.0; system_size];
    let time = 1234567.;
    universe.derive(time, &temperatures, &mut result).unwrap();
    assert_eq!(expected, result);
}

#[test]
fn test_insolation() {
    panic!();
}

#[test]
fn test_temperature_derivatives() {
    let system_size = 5;
    let temperatures = vec![300.; system_size];
    let mut universe = test_universe(system_size);

    // first and second derivatives
    let expected_first = vec![1., 2., 3., 4., 5.];
    let expected_second = vec![1., 2., 3., 4., 5.];
    universe.temperature_derivative(&temperatures);
    assert_eq!(expected_first, universe.dx);
    assert_eq!(expected_second, universe.dx2);
}

#[test]
fn test_sinusoidal_excitation() {
    panic!();
}

#[test]
fn test_damping() {
    panic!();
}

#[test]
fn test_albedo() {
    panic!();
}

#[test]
fn test_eccentricity() {
    panic!();
}

#[test]
fn test_obliquity() {
    panic!();
}

#[test]
fn test_declination() {
    panic!();
}

#[test]
fn test_stellar_longitude() {
    panic!();
}

// Planet tests
#[test]
fn test_merged_ir_cooling() {
    let temperature = 300.;
    let planet = test_planet();
    let result = planet.merged_ir_cooling(temperature).unwrap();

    let expected = 5.;
    assert_eq!(expected, result);
}

#[test]
fn test_heat_capacity() {
    let temperature = 300.;
    let planet = test_planet();
    let result = planet.heat_capacity(temperature);

    let expected = 5.;
    assert_eq!(expected, result);
}

#[test]
fn test_energy_transport() {
    let dx = 0.;
    let dx2 = 0.;
    let lat_tan = 0.;
    let planet = test_planet();
    let result = planet.energy_transport(dx, dx2, lat_tan);
    let expected = 5.;
    assert_eq!(expected, result);
}

// Benches
#[cfg(feature = "divan")]
#[divan::bench]
fn bench_stellar_longitude(bencher: divan::Bencher) {
    let planet = test_planet();
    let temperature = 300.;
    bencher.bench(|| planet.stellar_longitude(temperature));
}

#[cfg(feature = "divan")]
#[divan::bench]
fn bench_albedo(bencher: divan::Bencher) {
    let planet = test_planet();
    let temperature = 300.;
    // Warmup
    for _ in 0..256_000 {
        let _ = planet.albedo(temperature);
    }
    bencher.bench(|| planet.albedo(temperature));
}
