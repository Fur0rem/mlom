use mlom::system::System;
use std::path::Path;

fn main() {
	let system = System::from_file(Path::new("dataset/particles.xyz"), 0);
	println!("{}", system.microscopic_energy());
}
