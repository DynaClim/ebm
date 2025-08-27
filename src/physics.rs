use crate::constants::{
    COLDICE_HEAT_CAPACITY, CRITICAL_TEMPERATURE, ICE_HEAT_CAPACITY, LAND_HEAT_CAPACITY,
    OCEAN_HEAT_CAPACITY, PERIHELION_LONGTITUDE_DEG, PI, SECONDS_IN_DAY, SECONDS_IN_YEAR,
    SOLAR_CONSTANT_AT_1_AU, STEFAN_BOLTZMANN,
};
use anyhow::{Context, Result, bail};
use itertools::izip;
use sci_file::{Interpolator, LinSpace};
use serde::{Deserialize, Serialize};

const ICE_HEAT_CAPACITY_RECIPROCAL: f64 = ICE_HEAT_CAPACITY.recip();
const COLDICE_HEAT_CAPACITY_RECIPROCAL: f64 = COLDICE_HEAT_CAPACITY.recip();
const PERIHELION_LONGTITUDE_RAD: f64 = PERIHELION_LONGTITUDE_DEG.to_radians();

/// Contains all static data required for calculation of physics in this module.
/// Provides the derivation function [derive] as the public interface.
#[derive(Deserialize, Serialize, Clone, Debug, Default, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct Universe {
    #[serde(default)]
    time: f64,
    /// User specified planetary parameters
    /// NOTE: this field is marked only as default to allow easier testing. Values _MUST_ be provided by the user.
    #[serde(skip_serializing, default)]
    planet: Planet,
    temperatures: Vec<f64>,

    // Precomputed latitude arrays
    #[serde(skip)]
    latitude: Vec<Trig>,

    // Preallocated arrays, used as scratch for calculations.
    #[serde(skip)]
    albedos: Vec<f64>,
    #[serde(skip)]
    source: Vec<f64>,
    #[serde(skip)]
    transport: Vec<f64>,
    #[serde(skip)]
    sink: Vec<f64>,
    #[serde(skip)]
    dx: Vec<f64>,
    #[serde(skip)]
    dx2: Vec<f64>,
}

impl Universe {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn new_from(planet: Planet, temperatures: Vec<f64>) -> Self {
        Self {
            planet,
            temperatures,
            ..Self::default()
        }
    }

    pub fn initialise(&mut self) -> Result<()> {
        let system_size = self.temperatures.len();
        self.latitude = Trig::new_vec(system_size);
        self.albedos = vec![0.0; system_size];
        self.source = vec![0.0; system_size];
        self.transport = vec![0.0; system_size];
        self.sink = vec![0.0; system_size];
        self.dx = vec![0.0; system_size];
        self.dx2 = vec![0.0; system_size];
        self.planet.initialise()
    }

    #[must_use]
    pub fn temperatures(&self) -> Vec<f64> {
        self.temperatures.clone()
    }

    pub fn update_for_output(&mut self, time: f64, temperatures: &[f64]) {
        self.time = time / SECONDS_IN_YEAR;
        self.temperatures.copy_from_slice(temperatures);
    }

    //TODO Sid reference this function
    pub fn derive(&mut self, time: f64, temperatures: &[f64], d_temp: &mut [f64]) -> Result<()> {
        // First (K rad^-1) and second (K rad^-2) derivative of the temperature (K rad^-1)
        self.temperature_derivative(temperatures);

        // Compute the sink
        for (olr, temp) in izip!(&mut self.sink, temperatures) {
            *olr = self.planet.merged_ir_cooling(*temp)?;
        }

        // Compute the current albedo for each latitude
        for (albedo, temp) in izip!(&mut self.albedos, temperatures) {
            *albedo = self.planet.albedo(*temp);
        }

        // Compute the source
        self.insolation(time)?;

        for (transp, dx_i, dx2_i, lat) in
            izip!(&mut self.transport, &self.dx, &self.dx2, &self.latitude)
        {
            *transp = self.planet.energy_transport(*dx_i, *dx2_i, lat.tan);
        }

        for (outval, source_i, transp_i, sink_i, temp_i) in izip!(
            d_temp,
            &self.source,
            &self.transport,
            &self.sink,
            temperatures
        ) {
            *outval = self.planet.heat_capacity(*temp_i) * (source_i + transp_i - sink_i);
        }

        Ok(())
    }

    // Absorbed Stellar Radiation (ASR) function, Spiegel et al. 2008
    // Computes insolation received at the top of the atmosphere TOA as a function of the global orbital parameters defined.
    // Daily averaged insolation (W m^-2), Hartmann 2016.
    fn insolation(&mut self, time: f64) -> Result<()> {
        let decl = self.planet.declination(time)?;
        let decl_tmp = abs!(decl.rad) - PI / 2.;
        let eccentricity = self.planet.eccentricity(time)?;

        // Square of the ratio of the semi-major axis over the instantneous distance given by the Kepler's third law for eccentric orbits.
        // see expression of the daily averaged insolation, Hartmann 2016.
        let distance_ratio_squared = (1.
            + eccentricity * cos!(self.planet.stellar_longitude(time) - PERIHELION_LONGTITUDE_RAD))
        .powi(2)
            / (1. - eccentricity.powi(2)).powi(2);

        let factor = self.planet.insolation_factor * distance_ratio_squared;

        for (outval, lat, albedo_i) in izip!(&mut self.source, &self.latitude, &self.albedos) {
            // Conditions determining the value of the hour angle (rad) as defined in Berger 1978 and Hartmann 2016.
            // Cosine of the zenith angle for any latitude, season and time of the day, Hartmann 2016.
            let zenith_angle_cos = {
                if (abs!(lat.rad) + decl_tmp) < 0.0 {
                    let hour_angle = acos!(-lat.tan * decl.tan);
                    hour_angle * lat.sin * decl.sin + sin!(hour_angle) * lat.cos * decl.cos
                } else if lat.rad * decl.rad > 0.0 {
                    PI * lat.sin * decl.sin
                } else {
                    0.0
                }
            };

            // Absorbed stellar radiation (W m^-2), Spiegel et al.
            *outval = (factor * zenith_angle_cos) * (1.0 - *albedo_i);
        }
        Ok(())
    }

    fn first_temperature_derivative(
        temp_prev: f64,
        temp_next: f64,
        lat_prev: f64,
        lat_next: f64,
    ) -> f64 {
        (temp_next - temp_prev) / (lat_next - lat_prev)
    }

    fn second_temperature_derivative(
        temp_prev: f64,
        temp: f64,
        temp_next: f64,
        lat_prev: f64,
        lat_next: f64,
    ) -> f64 {
        4. * (temp_next - (2. * temp) + temp_prev) / (lat_next - lat_prev).powi(2)
    }

    // Computation of the first and second derivative of the temperature as a function of the sine of latitude (cf 1D time-dependent diffusion equation).
    // Spiegel et al. 2008
    fn temperature_derivative(&mut self, temperatures: &[f64]) {
        let temp_len = temperatures.len();
        assert!(
            temp_len == self.latitude.len()
                && temp_len == self.dx.len()
                && temp_len == self.dx2.len()
        );

        for (dx, dx2, temperature_window, latitude_window) in izip!(
            &mut self.dx[1..temp_len],
            &mut self.dx2[1..temp_len],
            temperatures.windows(3),
            self.latitude.windows(3),
        ) {
            let [temp_prev, temp, temp_next] = temperature_window else {
                unreachable!();
            };
            let [lat_prev, _, lat_next] = latitude_window else {
                unreachable!();
            };

            *dx = Self::first_temperature_derivative(
                *temp_prev,
                *temp_next,
                lat_prev.rad,
                lat_next.rad,
            );
            *dx2 = Self::second_temperature_derivative(
                *temp_prev,
                *temp,
                *temp_next,
                lat_prev.rad,
                lat_next.rad,
            );
        }

        self.dx2[0] = 4. * (2. * temperatures[1] - (2. * temperatures[0]))
            / (2. * (self.latitude[1].rad - self.latitude[0].rad)).powi(2);

        self.dx2[temp_len - 1] = 4.
            * (2. * temperatures[temp_len - 2] - (2. * temperatures[temp_len - 1]))
            / (2. * (self.latitude[temp_len - 1].rad - self.latitude[temp_len - 2].rad)).powi(2);
    }
}

// Defines the sources of calculation used for each quantity (e.g. eccentricity, obliquity).
#[derive(Deserialize, Serialize, Clone, Debug, PartialEq)]
pub enum Source {
    // Value does not vary with time.
    Constant(f64),
    // Value decreases with time.
    Damping {
        // time it takes for the value to linearly decrease from initial to 0 (year)
        damping_timescale: f64,
        initial: f64,
    },
    // Value excited following a sine time evolution.
    Sin {
        // Amplitude of the oscillation of the value, 0 < amplitude < 1, no unit
        amplitude: f64,
        // Period of the oscillation of the value (Earth days)
        period: f64,
        initial: f64,
    },
    // Value interpolated from precomputed grid.
    Interpolate(Interpolator<f64>),
}

impl Default for Source {
    fn default() -> Self {
        Self::Constant(0.0)
    }
}

impl Source {
    fn calculate(&self, time: f64) -> Result<f64> {
        match &self {
            Source::Constant(initial) => Ok(*initial),
            Source::Damping {
                damping_timescale,
                initial,
            } => Ok(Self::damping(time, *initial, *damping_timescale)),
            Source::Sin {
                amplitude,
                period,
                initial,
            } => Ok(Self::sinusoidal_excitation(
                time, *amplitude, *period, *initial,
            )),
            Source::Interpolate(interpolator) => {
                let (_x, y) = interpolator.interpolate(time)?;
                Ok(y)
            }
        }
    }

    fn sinusoidal_excitation(time: f64, amplitude: f64, period: f64, initial: f64) -> f64 {
        amplitude * sin!((2. * PI) / (period * SECONDS_IN_DAY) * time) + initial
    }

    fn damping(time: f64, initial: f64, damping_timescale: f64) -> f64 {
        if time < damping_timescale {
            -(initial / damping_timescale) * time + initial
        } else {
            0.0
        }
    }
}

#[derive(Deserialize, Serialize, Clone, Debug, Default, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct Planet {
    /// Planet albedo at low temperature, no unit
    /// Albedo of the surface of the planet below the temperature of water phase transition
    /// (combined effect of ice/snow and clouds at low temperature)
    albedo_low_temperature: f64,
    /// Planet albedo at high temperature, no unit
    /// Albedo of the surface of the planet above the temperature of water phase transition
    /// (effective albedo of clouds and surface at high temperature)
    albedo_high_temperature: f64,
    /// Temperarure of the center of the albedo transition (between high a low values) (K)
    temperature_centre_transition: f64,
    /// Temperature range of the transition of albedo around phase transition (K)
    temperature_width_transition: f64,
    // Internal cache
    #[serde(skip)]
    gilmore_albedo_term1: f64,
    #[serde(skip)]
    gilmore_albedo_term2: f64,

    /// Diurnally average the insolation factor of the planet.
    /// Hartmann 2016
    insolation_factor: f64,
    /// Diffusion coefficient (J mˆ-2 sˆ-1 Kˆ-1 x (omega_p/omega_earth)ˆ-2)
    /// This value is usually tuned to reproduce Earth’s present latitudinal temperature gradient.
    /// Spiegel et al. 2008, Williams & Kasting 1997
    /// This value depends on the planetary angular spin frequency, needs to be adapted depending on the planet considered.
    /// Spiegel et al. 2008
    fiducial_diffusion_coefficient: f64,

    /// Time evolution of the eccentrity of the planetary orbit, no unit
    eccentricity: Source,
    /// Time evolution of the obliquity of the planteary orbit (deg)
    obliquity: Source,
    /// Interpolator for user specified Outgoing Longwave Radiation (OLR) values (W m^-2)
    outgoing_longwave_radiation: Interpolator<f64>,

    /// Orbital period of the planet (d)
    orbital_period: f64,
    /// Proportion of the planet's surface that is land, no unit. Note [`Planet.land_fraction`] + [`Planet.water_fraction`] must equal 1.
    land_fraction: f64,
    /// Proportion of the planet's surface that is aqua, no unit. Note [`Planet.land_fraction`] + [`Planet.water_fraction`] must equal 1.
    water_fraction: f64,

    // Calculated internally
    #[serde(skip)]
    heat_capacity_reciprocal: f64,
}

impl Planet {
    pub fn new() -> Self {
        Self::default()
    }

    fn initialise(&mut self) -> Result<()> {
        #[allow(clippy::float_cmp)]
        // User supplied fractions should equal exactly 1.0.
        if abs!(1.0 - (self.land_fraction + self.water_fraction)) > f64::EPSILON {
            bail!(
                "Land fraction + Water fraction does not equal 1: {} + {} = {}",
                self.land_fraction,
                self.water_fraction,
                self.land_fraction + self.water_fraction
            );
        }

        // Correct the diurnally averaged insolation (W m^-2) for the Earth at 1 AU to account for the insolation received
        // by a planet located at TRAPPIST-1e semi-major axis and orbiting a Sun with weaker irradiation.
        // TRAPPIST-1e irradiation: Gillon et al. 2017
        // Spiegel et al. 2008
        self.insolation_factor *= SOLAR_CONSTANT_AT_1_AU / PI;

        // Convert orbital period from days to seconds.
        self.orbital_period *= SECONDS_IN_DAY;

        self.heat_capacity_reciprocal = 1.
            / (self.water_fraction * OCEAN_HEAT_CAPACITY + self.land_fraction * LAND_HEAT_CAPACITY);

        self.gilmore_albedo_term1 =
            f64::midpoint(self.albedo_low_temperature, self.albedo_high_temperature);
        self.gilmore_albedo_term2 =
            (self.albedo_low_temperature - self.albedo_high_temperature) / 2.;

        // Convert obliquity from degrees to radians.
        match &mut self.obliquity {
            Source::Constant(initial) => {
                *initial = initial.to_radians();
            }
            // Convert obliquity damping timescale from years to seconds.
            Source::Damping {
                initial,
                damping_timescale,
            } => {
                *initial = initial.to_radians();
                *damping_timescale *= SECONDS_IN_YEAR;
            }
            _ => (),
        }

        // Convert eccentricity damping timescale from years to seconds.
        if let Source::Damping {
            damping_timescale, ..
        } = &mut self.eccentricity
        {
            *damping_timescale *= SECONDS_IN_YEAR;
        }

        Ok(())
    }

    // Gilmore 2014, equation 5, modified with a pade approximation of tanh(x)
    // for performance (~1.6x speedup)
    // https://mathr.co.uk/blog/2017-09-06_approximating_hyperbolic_tangent.html
    //
    // In general, tanh⁡(x) asymptotically tends towards ± 1.
    // x ≈ 0 => tanh(x) ≈ x
    // and the further x is from 0, the closer tanh(x) is to 1.
    // x ≈ ±3 => tanh⁡(x) ≈ ± 0.9950 (very close to the asymptodes).
    //
    // So a tanh approximation only needs to be accurate between (-3, 3)
    // and we clamp x to ± 1 outside of this range.
    //
    // In practice, we simplify the output of the albdeo function
    // and clamp directly on to `albedo_low_temperature` and `albedo_high_temperature`
    // outside of `temperature_centre_transition` ± 3 * `temperature_width_transition` range
    // and use the pade tanh approximation inside the (-3, 3) range
    // (which makes the effective transition width ≈ 6ΔT).
    fn albedo(&self, temperature: f64) -> f64 {
        let x =
            (temperature - self.temperature_centre_transition) / self.temperature_width_transition;
        if x < -3. {
            // low clamp
            self.albedo_low_temperature
        } else if x > 3. {
            // high clamp
            self.albedo_high_temperature
        } else {
            // tanh approximation
            let y = x * ((27. + x.powi(2)) / (27. + 9. * x.powi(2)));

            self.gilmore_albedo_term1 - y * self.gilmore_albedo_term2
        }
    }

    fn eccentricity(&self, time: f64) -> Result<f64> {
        self.eccentricity.calculate(time).context("Eccentricity")
    }

    fn obliquity(&self, time: f64) -> Result<f64> {
        self.obliquity.calculate(time).context("Obliquity")
    }

    // Declination angle latitude of the substellar point as a function of time.
    // Berger 1978
    fn declination(&self, time: f64) -> Result<Trig> {
        let declination_radians =
            asin!(sin!(self.obliquity(time)?) * sin!(self.stellar_longitude(time)));
        Ok(Trig::new(declination_radians))
    }

    // Stellar longitude as a function of time
    // Williams & Kasting 1997
    fn stellar_longitude(&self, time: f64) -> f64 {
        (time % self.orbital_period) * (2.0 * PI / self.orbital_period)
    }

    // Calculates the outgoing longwave radiation of the planet as a function of its temperature.
    fn merged_ir_cooling(&self, temperature: f64) -> Result<f64> {
        // For T < 322.5 K:
        //      Modified Stefan-Boltzmann law to account for Earth's atmospheric composition, Dressing et al. 2010
        // For 322.5 < T < 787.2 K:
        //      Interpolation of GCM simulations, Chaverot et al. 2023
        // For 893.2 < T < 1055.9:
        //      Interpolation of GCM simulations, Turbet et al. 2022
        // The Stefan-Boltzmann relationship is not valid above 322.5 K
        let olr = if temperature < 322.5 {
            // 0.5925 (no unit) is the coefficient accounting for atmospheric composition
            (STEFAN_BOLTZMANN * temperature.powi(4))
                / (1. + (0.5925 * (temperature / CRITICAL_TEMPERATURE).powi(3)))
        } else {
            let (_x, y) = self
                .outgoing_longwave_radiation
                .interpolate(temperature)
                .context("Outgoing Longwave Radiation")?;
            y
        };

        Ok(olr)
    }

    // Calculates effective surface heat capacity depending on the temperature and the surface type (ocean, land, ice, cold ice)
    // Spiegel et al. 2008; Williams & Kasting 1997
    // Can be tweaked to include other materials. If this is done, a corresponding tweak has to be made to the albedo function.
    fn heat_capacity(&self, temperature: f64) -> f64 {
        // All warm surfaces have the same ocean to land ratio.
        if temperature > CRITICAL_TEMPERATURE {
            self.heat_capacity_reciprocal
        // Earth like planet.
        } else if temperature >= CRITICAL_TEMPERATURE - 10. {
            ICE_HEAT_CAPACITY_RECIPROCAL
        } else {
            COLDICE_HEAT_CAPACITY_RECIPROCAL
        }
    }

    // 1D time-dependent diffusion equation, Spiegel et al. 2008
    // Calculates variation of the temperature as a function of time (K s^-1)
    fn energy_transport(&self, dx: f64, dx2: f64, lat_tan: f64) -> f64 {
        self.fiducial_diffusion_coefficient * (dx2 - (lat_tan * dx))
    }
}

// Radian and associated trigonometric values.
#[derive(Deserialize, Serialize, Clone, Debug, Default, PartialEq)]
struct Trig {
    rad: f64,
    cos: f64,
    sin: f64,
    tan: f64,
}

impl Trig {
    fn new(x: f64) -> Self {
        Self {
            cos: cos!(x),
            sin: sin!(x),
            tan: tan!(x),
            rad: x,
        }
    }

    // Create the array and precompute all the trigonometric values
    fn new_vec(system_size: usize) -> Vec<Self> {
        LinSpace::new(90.0_f64, -90.0_f64, system_size)
            .map(|x| Self::new(x.to_radians()))
            .collect()
    }
}

// References:
// Berger 1978, https://doi.org/10.1175/1520-0469(1978)035<2362:LTVODI>2.0.CO;2
// Chaverot et al. 2023, https://doi.org/10.1051/0004-6361/202346936
// Dressing et al. 2010, https://doi.org/10.1088/0004-637X/721/2/1295
// Gillon et al. 2017, https://doi.org/10.1038/nature21360
// Gilmore 2014, https://doi.org/10.1093/mnras/stu302
// Hartmann 2016, https://doi.org/10.1016/B978-0-12-328531-7.00002-5
// Kasting 1997, https://doi.org/10.1006/icar.1997.5759
// Spiegel et al. 2008, https://doi.org/10.1086/588089
// Turbet et al. 2022, https://doi.org/10.1038/s41586-021-03873-w
// Williams & Kasting 1997, https://doi.org/10.1006/icar.1997.5759
