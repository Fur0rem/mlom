use mlom::system::System;
use std::path::Path;

fn main() {
	let mut system = System::from_file(Path::new("dataset/particles.xyz"), 0);
	system.energy_evolution(1000, "plots/jvaismetuer.png");
}
