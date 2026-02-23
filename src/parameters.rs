//! Parameters for the simulation, including physical constants, potential parameters, and algorithm parameters.

/////////////////////////////////////////////////////
////////// ISM 2 : Lennard-Jones potential //////////
/////////////////////////////////////////////////////
pub const R_STAR: f64 = 3.0;
pub const EPSILON_STAR: f64 = 0.2;

/////////////////////////////////////////////////
////////// ISM 3 : Periodic conditions //////////
/////////////////////////////////////////////////

/// Cutoff radius for interactions
pub const R_CUT: f64 = 10.0;

/// Smoothing width around cutoff, to ensure continuity of the potential and forces near the cutoff radius. The smoothing is applied between R_MIN and R_MAX, which are defined as R_CUT +/- CUTOFF_SMOOTHING_WIDTH.
pub const CUTOFF_SMOOTHING_WIDTH: f64 = 0.1 * R_CUT;

/// Start of smoothing region around cutoff
pub const R_MIN: f64 = R_CUT - CUTOFF_SMOOTHING_WIDTH;

/// End of smoothing region around cutoff
pub const R_MAX: f64 = R_CUT + CUTOFF_SMOOTHING_WIDTH;

/// Pre-computed inverse of the smoothing width, for efficiency in the smoothing function
pub const INVERSE_DIFF_MAX_MIN: f64 = 1.0 / (R_MAX - R_MIN);

/// A distance large enough to be considered "infinite" and come back to non-periodic conditions in tests
pub const FAR_AWAY: f64 = 99999999.0;

/// Box side length for periodic conditions
pub const BOX_SIDE: f64 = 42.0;

/////////////////////////////////////////////////////////////////
////////// ISM 4 : Movement, Velocity-Verlet algorithm //////////
/////////////////////////////////////////////////////////////////

/// Time step for the Velocity-Verlet algorithm, in femtoseconds
pub const DELTA_TIME: f64 = 1.0;

pub const CONVERSION_FORCE: f64 = 0.0001 * 4.186;
pub const PARTICLE_MASS: f64 = 18.0;
pub const R_CONSTANT: f64 = 0.00199;

/// Initial/Target temperature for the system
pub const T_0: f64 = 300.0;

/// Berendsen Thermostat correction factor
pub const GAMMA: f64 = 0.01;

/// A number of steps large enough to never apply the thermostat, if we want to run a simulation without thermostat
pub const NEVER_APPLY_THERMOSTAT: usize = usize::MAX;

/////////////////////////////////////////////////////
////////// ISM 8 : Faster neighbor queries //////////
/////////////////////////////////////////////////////

/// Maximum velocity anticipated for a particle in the system, used to determine the skin size for Verlet lists and the buffer size for domain decomposition
pub const MAX_PARTICLE_VELOCITY: f64 = R_STAR;

/// Frequency of reconstruction of Verlet lists, in number of iterations.
pub const REBUILD_VERLET_LISTS_FREQUENCY: usize = 5;

/// Number of steps to give the simulation to spread particles evenly until applying the Verlet list optimization
pub const PRE_VERLET_LISTS_STEPS: usize = 10;
