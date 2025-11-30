use std::path::Path;

use mlom::{algebra::Vector3, system::System};

#[test]
fn sum_of_forces_is_null() {
	let system = System::from_file(Path::new("dataset/particles.xyz"), 0);
	assert_eq!(system.sum_of_forces(), Vector3::zero());
}
