use mlom::{parameters::NEVER_APPLY_THERMOSTAT, system::System};
use std::path::Path;

fn main() {
	let mut system = System::from_file(Path::new("dataset/particles.xyz"), 0);

	let (ke_init, t_init) = system.kinetic_energy_and_temperature();
	println!("INIT: K = {}, T = {}", ke_init, t_init);

	// system.simulate(1000, NEVER_APPLY_THERMOSTAT, "plots/v18.png");
	// system.simulate(2000, 5, "plots/v15_with_thermostat.png");

	system.simulate_with_verlet_lists(3000, NEVER_APPLY_THERMOSTAT, "plots/smooth.png");
}
