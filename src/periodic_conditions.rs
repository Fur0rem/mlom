use crate::{
	algebra::Vector3,
	energy::energy_between_particles,
	parameters::*,
	system::{Particle, System},
};

/// Computes the neighbors in 3D of the simulation box.
pub fn neighboring_3d_translations(box_side: f64) -> Vec<Vector3> {
	let mut symmetries = Vec::with_capacity(27);
	for x in -1..=1 {
		for y in -1..=1 {
			for z in -1..=1 {
				let x = x as f64 * box_side;
				let y = y as f64 * box_side;
				let z = z as f64 * box_side;

				let translation = Vector3::from(x, y, z);
				symmetries.push(translation);
			}
		}
	}

	return symmetries;
}

impl System {
	/// Compute the microscopic energy in the system, according to the Lennard-Jones potential, with periodic conditions.
	pub fn microscopic_energy_periodic(&self, translations: &[Vector3], radius_cut: f64) -> f64 {
		let mut total = 0.0;
		for sym in translations {
			for i in 0..self.nb_particles_total() {
				for j in 0..self.nb_particles_total() {
					if i == j && *sym == Vector3::zero() {
						continue;
					}

					// Compute translated particle j
					let particle_j_with_symmetry = Particle {
						coordinates: (self.particles[j].coordinates + *sym).as_point(),
						momentum: self.particles[j].momentum,
					};

					total += energy_between_particles(&self.particles[i], &particle_j_with_symmetry, radius_cut);
				}
			}
		}

		return (4.0 * EPSILON_STAR * total) / 2.0;
	}

	/// Compute the sum of all the forces applied to particles in the system, with periodic conditions.
	pub fn sum_of_forces_periodic(forces: &Vec<Vector3>) -> Vector3 {
		let mut sx = 0.0;
		let mut sy = 0.0;
		let mut sz = 0.0;

		for f in forces {
			sx += f.x();
			sy += f.y();
			sz += f.z();
		}

		return Vector3::from(sx, sy, sz);
	}
}
