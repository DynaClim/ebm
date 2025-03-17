pub use astro_const::constants::*;
pub use core::f64::consts::PI;

/// Stefan-Boltzmann constant (W mˆ-2 Kˆ-4)
/// Spiegel et al.
/// Dressing et al. 2010
/// TODO  Reference to be completed, neither of them give the value explicitly
pub const STEFAN_BOLTZMANN: f64 = 5.67037e-8;

/// Critical temperature (K), Phase transition from water to ice at pressure P = 1 (atm)
/// Spiegel et al. 2008
/// Dressing et al. 2010
pub const CRITICAL_TEMPERATURE: f64 = 273.0;

/// Solar constant (W m^-2)
/// Spiegel et al. 2008
pub const SOLAR_CONSTANT_AT_1_AU: f64 = 1360.0;

/// Effective heat capacity of land (J mˆ-2 Kˆ-1)
/// Williams and Kasting 1997
pub const LAND_HEAT_CAPACITY: f64 = 5.25e6;
/// Effective heat capacity of ocean: 40.0 * [`LAND_HEAT_CAPACITY`] (J mˆ-2 Kˆ-1)
/// Williams and Kasting 1997
pub const OCEAN_HEAT_CAPACITY: f64 = 40. * LAND_HEAT_CAPACITY;
/// Effective heat capacity of ice when 263 K < T < 273 K: 9.2 * [`LAND_HEAT_CAPACITY`] (J mˆ-2 Kˆ-1)
/// Williams and Kasting 1997
pub const ICE_HEAT_CAPACITY: f64 = 9.2 * LAND_HEAT_CAPACITY;
/// Effective heat capacity of cold ice when T < 263 K: 2.0 * [`LAND_HEAT_CAPACITY`] (J mˆ-2 Kˆ-1)
/// Williams and Kasting 1997
pub const COLDICE_HEAT_CAPACITY: f64 = 2.0 * LAND_HEAT_CAPACITY;

/// Albedo (IR) value for ice surface for TRAPPIST-1 planets, no units
/// Rushby et al. 2020 (Table 1)
pub const M_DWARF_TR1_ICE: f64 = 0.134;
/// Albedo (IR) value for snow surface for TRAPPIST-1 planets, no units
/// Rushby et al. 2020 (Table 1)
pub const M_DWARF_TR1_SNOW: f64 = 0.418;
/// Representative albedo of open ocean surface in the East Antarctic sea ice zone in spring and summer, for NIR (> 700 nm) band, no units
/// Brandt et al. 2005
pub const M_DWARF_TR1_OCEAN: f64 = 0.06;

/// Orbital period of TRAPPIST-1e (d)
/// Agol et al. 2021 (Table 2)
pub const ORBITAL_PERIOD_TR1_E: f64 = 6.101_013;
/// Orbital period of TRAPPIST-1e (s)
/// Agol et al. 2021 [`ORBITAL_PERIOD_TR1_E`] * [`SECONDS_IN_DAY`]
pub const ORBITAL_PERIOD_TR1_E_SEC: f64 = ORBITAL_PERIOD_TR1_E * SECONDS_IN_DAY;

/// Longitude of perihelion of the Earth
/// E. J. Rose, CLIMLAB Python toolkit
/// E. J. Rose, Climate Laboratory Book
/// Berger and Loutre 1991
pub const PERIHELION_LONGTITUDE_DEG: f64 = 281.37;

// Agol et al. 2021, https://doi.org/10.3847/PSJ/abd022
// Berger and Loutre 1991 https://doi.org/10.1016/0277-3791(91)90033-Q
// Brandt et al. 2005, https://doi.org/10.1175/JCLI3489.1
// Dressing et al. 2010, https://doi.org/10.1088/0004-637X/721/2/1295
// E. J. Rose, CLIMLAB Python toolkit, https://doi.org/10.21105/joss.00659
// E. J. Rose, Climate Laboratory, https://brian-rose.github.io/ClimateLaboratoryBook/courseware/insolation.html
// Rushby et al. 2020, https://doi.org/10.3847/1538-4357/abbe04
// Spiegel et al. 2008, https://doi.org/10.1086/588089
// Williams and Kasting 1997, https://doi.org/10.1006/icar.1997.5759
