use std::path::Path;

use mlom::parameters::FAR_AWAY;
use mlom::periodic_conditions::neighboring_3d_symmetries;
use mlom::{algebra::Vector3, system::System};
use mlom::{assert_approx_eq, assert_vector_approx_eq};

#[test]
fn sum_of_forces_is_null() {
	let system = System::from_file(Path::new("dataset/particles.xyz"), 0);

	// Non periodic conditions
	system.compute_forces();
	assert_vector_approx_eq!(System::sum_of_forces(&system.compute_forces()), Vector3::zero());

	// Periodic conditions
	assert_vector_approx_eq!(System::sum_of_forces_periodic(&system.compute_forces_periodic()), Vector3::zero());
}

/// NOTE: Only works with the clear cut potential, not the smooth cut potential, because I hard-coded it like an idiot.
#[test]
fn if_nb_sym_1_then_equivalent_to_non_periodic() {
	let system = System::from_file(Path::new("dataset/particles.xyz"), 0);
	let u_lj_non_periodic = system.microscopic_energy();
	let u_lj_periodic = system.microscopic_energy_periodic(&[Vector3::zero()], FAR_AWAY);
	assert_approx_eq!(u_lj_non_periodic, u_lj_periodic);
}

/// NOTE: Only works with the clear cut potential, not the smooth cut potential, because I hard-coded it like an idiot.
#[test]
fn if_far_away_then_equivalent_to_non_periodic() {
	let system = System::from_file(Path::new("dataset/particles.xyz"), 0);
	let u_lj_non_periodic = system.microscopic_energy();
	let u_lj_periodic = system.microscopic_energy_periodic(&neighboring_3d_symmetries(FAR_AWAY), FAR_AWAY);
	assert_approx_eq!(u_lj_non_periodic, u_lj_periodic);
}

#[test]
fn right_initial_temperature() {
	let system = System::from_file(Path::new("dataset/particles.xyz"), 0);
	let (_, t_init) = system.kinetic_energy_and_temperature();
	assert_approx_eq!(t_init, 300.0);
}
