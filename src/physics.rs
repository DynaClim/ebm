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
#[derive(Deserialize, Serialize, Clone, Debug)]
#[serde(deny_unknown_fields)]
pub struct Universe {
    /// User specified planetary parameters
    planet: Planet,
    initial_temperatures: Vec<f64>,

    // Precomputed latitude arrays
    #[serde(default)]
    latitude: Vec<Trig>,

    // Preallocated arrays, used as scratch for calculations.
    #[serde(default)]
    albedos: Vec<f64>,
    #[serde(default)]
    source: Vec<f64>,
    #[serde(default)]
    transport: Vec<f64>,
    #[serde(default)]
    sink: Vec<f64>,
    #[serde(default)]
    dx: Vec<f64>,
    #[serde(default)]
    dx2: Vec<f64>,
}

impl Universe {
    pub fn initialise(&mut self) -> Result<()> {
        let system_size = self.initial_temperatures.len();
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
    pub fn initial_temperatures(&self) -> Vec<f64> {
        self.initial_temperatures.clone()
    }

    //TODO Sid reference this function
    pub fn derive(&mut self, time: f64, temperatures: &[f64], d_temp: &mut [f64]) -> Result<()> {
        // First (K rad^-1) and second (K rad^-2) derivative of the temperature (K rad^-1)
        self.d2t_dx2(temperatures);
        // Compute the sink
        self.planet
            .merged_ir_cooling(temperatures, &mut self.sink)?;

        // Compute the current albedo for each latitude
        self.planet.albedo(temperatures, &mut self.albedos);

        // Compute the source
        self.insolation(time)?;
        self.energy_transport();

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
                } else if lat.rad.is_normal() && (lat.rad.signum() == decl.rad.signum()) {
                    // check if lat.rad * decl.rad > 0.0, but since we don't care about
                    // the result, we can infer if the product would be positive.
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

    // 1D time-dependent diffusion equation, Spiegel et al. 2008
    // Calculates variation of the temperature as a function of time (K s^-1)
    fn energy_transport(&mut self) {
        for (transp, dx_i, dx2_i, lat) in
            izip!(&mut self.transport, &self.dx, &self.dx2, &self.latitude)
        {
            *transp = self.planet.fiducial_diffusion_coefficient * (dx2_i - (lat.tan * dx_i));
        }
    }

    // Computation of the second derivative of the temperature as a function of the sine of latitude (cf 1D time-dependent diffusion equation).
    // Spiegel et al. 2008
    fn d2t_dx2(&mut self, temperatures: &[f64]) {
        let temp_len = temperatures.len();
        assert!(temp_len == self.latitude.len() && temp_len == self.dx2.len());
        let temperatures = &temperatures[..temp_len];
        let latitude = &self.latitude[..temp_len];
        let dx = &mut self.dx[..temp_len];
        let dx2 = &mut self.dx2[..temp_len];

        dx2[0] = 4. * (2. * temperatures[1] - (2. * temperatures[0]))
            / (2. * (latitude[1].rad - latitude[0].rad)).powi(2);

        for i in 1..temp_len - 1 {
            dx[i] = (temperatures[i + 1] - temperatures[i - 1])
                / (latitude[i + 1].rad - latitude[i - 1].rad);

            dx2[i] = 4. * (temperatures[i + 1] - (2. * temperatures[i]) + temperatures[i - 1])
                / (latitude[i + 1].rad - latitude[i - 1].rad).powi(2);
        }

        dx2[temp_len - 1] = 4.
            * (2. * temperatures[temp_len - 2] - (2. * temperatures[temp_len - 1]))
            / (2. * (latitude[temp_len - 1].rad - latitude[temp_len - 2].rad)).powi(2);
    }
}

// Defines the sources of calculation used for each quantity (e.g. eccentricity, obliquity).
#[derive(Deserialize, Serialize, Clone, Debug)]
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

// TODO complete documentation of these variables
#[derive(Deserialize, Serialize, Clone, Debug)]
#[serde(deny_unknown_fields)]
pub struct Planet {
    /// Planet albedo at low temperature, no unit
    /// Albedo of the surface of the planet below the temperature of water phase transition
    /// Gilmore 2014
    albedo_low_temperature: f64,
    /// Planet albedo at high temperature, no unit
    /// Albedo of the surface of the planet above the temperature of water phase transition
    /// Gilmore 2014
    albedo_high_temperature: f64,
    /// Temperarure of the center of the albedo transition (between high a low values) (K)
    /// Gilmore 2014
    /// TODO Write this comment elsewhere: We suppose 1 bar atm and take 268 K according to the considerations of Spiegel et al. 2008.
    temperature_centre_transition: f64,
    /// Temperature range of the transition of albedo around water phase transition (K)
    /// Gilmore 2014
    /// TODO Write this comment elsewhere: We take 25 K as for the Earth in Gilmore 2014
    temperature_width_transition: f64,

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
    // Calculated internally as 1. / (water mass heat capacity + land mass heat capacity)
    #[serde(skip)]
    heat_capacity_reciprocal: f64,
}

impl Planet {
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

    // Calculates planet albedo according to Gilmore 2014, Equation 5.
    fn albedo(&self, temperatures: &[f64], albedos: &mut [f64]) {
        // Coefficients can be tuned to account for the spectral type of the star.
        // albedo_low_temperature:  combined effect of ice/snow and clouds at low temperature, no unit
        // albedo_high_temperature: effective albedo of clouds and surface at high temperature, no unit
        // planet.temperature_centre_transition: center of the smooth transition between the two constant values (K)
        // planet.temperature_width_transition:  width of the transition between the two constant values (K)
        for (albedo, temp) in izip!(albedos, temperatures) {
            *albedo = (f64::midpoint(self.albedo_low_temperature, self.albedo_high_temperature))
                - ((self.albedo_low_temperature - self.albedo_high_temperature) / 2.)
                    * tanh!(
                        (temp - self.temperature_centre_transition)
                            / self.temperature_width_transition
                    );
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
    fn merged_ir_cooling(&self, temperatures: &[f64], sink: &mut [f64]) -> Result<()> {
        // For T < 322.5 K:
        //      Modified Stefan-Boltzmann law to account for Earth's atmospheric composition, Dressing et al. 2010
        // For 322.5 < T < 787.2 K:
        //      Interpolation of GCM simulations, Chaverot et al. 2023
        // For 893.2 < T < 1055.9:
        //      Interpolation of GCM simulations, Turbet et al. 2022
        for (olr, temp) in izip!(sink, temperatures) {
            *olr = {
                // The Stefan-Boltzmann relationship is not valid above 322.5 K
                if *temp < 322.5 {
                    // 0.5925 (no unit) is the coefficient accounting for atmospheric composition
                    (STEFAN_BOLTZMANN * temp.powi(4))
                        / (1. + (0.5925 * (temp / CRITICAL_TEMPERATURE).powi(3)))
                } else {
                    let (_x, y) = self
                        .outgoing_longwave_radiation
                        .interpolate(*temp)
                        .context("Outgoing Longwave Radiation")?;
                    y
                }
            };
        }
        Ok(())
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
}

// Trigs and associated trigonometric values.
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

    // Create the latitude array and precompute all the trigonometric values
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
