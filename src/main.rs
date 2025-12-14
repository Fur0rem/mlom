use mlom::{
	algebra::Vector3,
	parameters::{BOX_SIDE, FAR_AWAY, R_CUT},
	periodic_conditions::neighboring_3d_translations,
	system::System,
};
use std::path::Path;

fn main() {
	let system = System::from_file(Path::new("dataset/particles.xyz"), 0);
	println!("{}", system.microscopic_energy());

	// Periodic conditions with 1 symmetry
	println!("{}", system.microscopic_energy_periodic(&[Vector3::zero()], BOX_SIDE));

	// Periodic conditions without cutoff
	println!(
		"{}",
		system.microscopic_energy_periodic(&neighboring_3d_translations(FAR_AWAY), BOX_SIDE)
	);

	// Periodic conditions with cutoff
	println!(
		"{}",
		system.microscopic_energy_periodic(&neighboring_3d_translations(BOX_SIDE), R_CUT)
	);
}
