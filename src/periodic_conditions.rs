use crate::{
	algebra::Vector3,
	parameters::{EPSILON_STAR, R_STAR},
	system::System,
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
				for j in (i + 1)..self.nb_particles_total() {
					// Compute translated particle j
					let particle_j_with_symmetry = self.particles[j].coordinates + sym;

					// Apply cut above given radius
					let dist_ij_squared = self.particles[i].coordinates.distance_to_squared(&particle_j_with_symmetry);
					assert!(dist_ij_squared != 0.0);
					if dist_ij_squared > radius_cut.powi(2) {
						continue;
					}

					// The usual energy term
					let r_star_over_r_ij_pow2 = R_STAR.powi(2) / dist_ij_squared;
					let r_star_over_r_ij_pow6 = r_star_over_r_ij_pow2.powi(3);
					let r_star_over_r_ij_pow12 = r_star_over_r_ij_pow6.powi(2);
					let u_ij = EPSILON_STAR * (r_star_over_r_ij_pow12 - (2.0 * r_star_over_r_ij_pow6));
					total += u_ij;
				}
			}
		}

		return 4.0 * total;
	}

	/// Compute the forces between pairs of particles, with periodic conditions.
	/// The new force that a particle j applies on particle i is the sum of forces of all its symmetries.
	pub fn compute_forces_periodic(&mut self, translations: &[Vector3], radius_cut: f64) {
		for i in 0..self.nb_particles_total() {
			for j in 0..self.nb_particles_total() {
				self.forces[i][j] = Vector3::zero();
			}
		}

		for i in 0..self.nb_particles_total() {
			for j in (i + 1)..self.nb_particles_total() {
				for sym in translations {
					// Compute translated particle j
					let particle_j_with_symmetry = self.particles[j].coordinates + sym;

					if self.particles[i].coordinates == particle_j_with_symmetry {
						continue;
					}

					// Apply cut above given radius
					let dist_ij_squared = self.particles[i].coordinates.distance_to_squared(&particle_j_with_symmetry);
					if dist_ij_squared > radius_cut.powi(2) {
						continue;
					}

					// Gradient function for any coordinate
					let r_star_over_r_ij_pow2 = R_STAR.powi(2) / dist_ij_squared;
					let r_star_over_r_ij_pow8 = r_star_over_r_ij_pow2.powi(3);
					let r_star_over_r_ij_pow14 = r_star_over_r_ij_pow8.powi(7);
					let gradient = |c_i: f64, c_j: f64| {
						-48.0 * EPSILON_STAR * (r_star_over_r_ij_pow14 - r_star_over_r_ij_pow8) * (c_i - c_j)
					};

					// Apply gradient in the x, y, and z directions
					let (x_i, y_i, z_i) = self.particles[i].xyz();
					let (x_j, y_j, z_j) = self.particles[j].xyz();
					self.forces[i][j] += Vector3::from(gradient(x_i, x_j), gradient(y_i, y_j), gradient(z_i, z_j));
				}
			}
		}
	}
}
