use std::path::Path;

use mlom::parameters::{BOX_SIDE, FAR_AWAY, R_CUT};
use mlom::periodic_conditions::neighboring_3d_translations;
use mlom::{algebra::Vector3, system::System};
use mlom::{assert_approx_eq, assert_vector_approx_eq};

#[test]
fn sum_of_forces_is_null() {
	let system = System::from_file(Path::new("dataset/particles.xyz"), 0);

	// Non periodic conditions
	system.compute_forces();
	assert_vector_approx_eq!(System::sum_of_forces(&system.compute_forces()), Vector3::zero());

	// Periodic conditions
	system.compute_forces_periodic(&neighboring_3d_translations(5.0), R_CUT);
	assert_vector_approx_eq!(
		System::sum_of_forces_periodic(&system.compute_forces_periodic(&neighboring_3d_translations(BOX_SIDE), R_CUT)),
		Vector3::zero()
	);
}

#[test]
fn if_nb_sym_1_then_equivalent_to_non_periodic() {
	let system = System::from_file(Path::new("dataset/particles.xyz"), 0);
	let u_lj_non_periodic = system.microscopic_energy();
	let u_lj_periodic = system.microscopic_energy_periodic(&[Vector3::zero()], FAR_AWAY);
	assert_approx_eq!(u_lj_non_periodic, u_lj_periodic);
}

#[test]
fn if_far_away_then_equivalent_to_non_periodic() {
	let system = System::from_file(Path::new("dataset/particles.xyz"), 0);
	let u_lj_non_periodic = system.microscopic_energy();
	let u_lj_periodic = system.microscopic_energy_periodic(&neighboring_3d_translations(FAR_AWAY), FAR_AWAY);
	assert_approx_eq!(u_lj_non_periodic, u_lj_periodic);
}
