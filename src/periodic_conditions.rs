use crate::{
	algebra::Vector3,
	parameters::{EPSILON_STAR, R_STAR},
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
					let particle_j_with_symmetry = (self.particles[j].coordinates + sym).as_point();

					// FIXME: Sometimes they're at the same place, no idea why
					if self.particles[i].coordinates == particle_j_with_symmetry {
						dbg!(&i, &j, &sym);
						dbg!(particle_j_with_symmetry);
						dbg!(self.particles[i]);
						dbg!(sym);
						dbg!(i);
						dbg!(j);
						continue;
					}

					// Apply cut above given radius
					let dist_ij_squared = self.particles[i].coordinates.distance_to_squared(&particle_j_with_symmetry);
					assert!(dist_ij_squared != 0.0);
					if dist_ij_squared < 0.0001 {
						println!("dist = {dist_ij_squared}");
						assert!(dist_ij_squared > 0.01);
					}
					if dist_ij_squared > radius_cut.powi(2) {
						continue;
					}

					// The usual energy term (ISM 3)
					let r_star_over_r_ij = R_STAR / dist_ij_squared.sqrt();
					let u_ij = r_star_over_r_ij.powi(12) - 2.0 * r_star_over_r_ij.powi(6);
					total += u_ij;
				}
			}
		}

		return (4.0 * EPSILON_STAR * total) / 2.0;
	}

	/// Compute the forces between pairs of particles, with periodic conditions.
	/// The new force that a particle j applies on particle i is the sum of forces of all its symmetries.
	pub fn compute_forces_periodic(&self, translations: &[Vector3], radius_cut: f64) -> Vec<Vec<Vec<Vector3>>> {
		let mut forces =
			vec![vec![vec![Vector3::zero(); self.nb_particles_total()]; self.nb_particles_total()]; translations.len()];

		for (sym_idx, sym) in translations.iter().enumerate() {
			for i in 0..self.nb_particles_total() {
				for j in 0..self.nb_particles_total() {
					if i == j && *sym == Vector3::zero() {
						continue;
					}

					// Compute translated particle j
					let particle_j_with_symmetry = (self.particles[j].coordinates + sym).as_point();

					if self.particles[i].coordinates == particle_j_with_symmetry {
						continue;
					}

					// Apply cut above given radius
					let dist_ij_squared = self.particles[i].coordinates.distance_to_squared(&particle_j_with_symmetry);
					if dist_ij_squared > radius_cut.powi(2) {
						continue;
					}

					if dist_ij_squared < 0.0001 {
						println!("dist = {dist_ij_squared}");
						continue;
						// assert!(dist_ij_squared > 0.01);
					}

					// Gradient function for any coordinate
					forces[sym_idx][i][j] += self.energy_gradient(
						&self.particles[i],
						&Particle {
							coordinates: particle_j_with_symmetry,
							momentum: self.particles[j].momentum,
						},
					);
				}
			}
		}

		return forces;
	}

	/// Compute the sum of all the forces between pairs of particles in the system with periodic conditions
	pub fn sum_of_forces_periodic(forces: &Vec<Vec<Vec<Vector3>>>) -> Vector3 {
		let mut sx = 0.0;
		let mut sy = 0.0;
		let mut sz = 0.0;

		for sym_idx in 0..forces.len() {
			for i in 0..forces[sym_idx].len() {
				for j in 0..forces[sym_idx][i].len() {
					let f = forces[sym_idx][i][j];
					sx += f.x();
					sy += f.y();
					sz += f.z();
				}
			}
		}

		return Vector3::from(sx, sy, sz);
	}
}
