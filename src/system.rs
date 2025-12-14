//! A system of particles

use std::{fs::File, io::Read, path::Path};

use crate::{
	algebra::{Point3, Vector3},
	parameters::{EPSILON_STAR, R_STAR},
};

/// A particle in the system
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Particle {
	/// The coordinates of the particle
	pub(crate) coordinates: Point3,
}

impl Particle {
	/// Parse a particle from a string
	/// The format should be an unsigned integer, followed by 3 floating point numbers, all separated by whitespace
	///
	/// # Arguments
	///
	/// * `s` - The string to parse
	pub fn parse(s: &str) -> Self {
		let mut parts = s.split_whitespace();
		let _type: usize = parts.next().unwrap().parse().unwrap();
		let x: f64 = parts.next().unwrap().parse().unwrap();
		let y: f64 = parts.next().unwrap().parse().unwrap();
		let z: f64 = parts.next().unwrap().parse().unwrap();

		return Self {
			coordinates: Point3::from(x, y, z),
		};
	}

	/// The x coordinate of the particle
	pub fn x(&self) -> f64 {
		self.coordinates.x()
	}

	/// The y coordinate of the particle
	pub fn y(&self) -> f64 {
		self.coordinates.y()
	}

	/// The z coordinate of the particle
	pub fn z(&self) -> f64 {
		self.coordinates.z()
	}

	/// The coordinates of the particle
	pub fn xyz(&self) -> (f64, f64, f64) {
		(self.x(), self.y(), self.z())
	}

	/// Compute the distance to another [particle](Self), squared
	/// More optimized than calling `distance_to` then squaring the result.
	///
	/// # Arguments
	///
	/// * `rhs` - The [particle](Self) to compute the squared distance to
	pub fn distance_to_squared(&self, rhs: &Self) -> f64 {
		self.coordinates.distance_to_squared(&rhs.coordinates)
	}

	/// Compute the distance to another [particle](Self)
	///
	/// # Arguments
	///
	/// * `rhs` - The [particle](Self) to compute the distance to
	pub fn distance_to(&self, rhs: &Self) -> f64 {
		return self.distance_to_squared(rhs).sqrt();
	}
}

/// A system of [particles](Particle)
#[derive(Debug, Clone, PartialEq)]
pub struct System {
	/// The [particles](Particle) in the system
	pub(crate) particles:          Vec<Particle>,
	/// The number of local particles (unused for now)
	pub(crate) nb_particles_local: usize,
}

impl System {
	/// Parse a system from a file
	pub fn from_file(path: &Path, nb_particles_local: usize) -> Self {
		// Read the file
		let mut file = File::open(path).unwrap();
		let mut contents = String::new();
		file.read_to_string(&mut contents).unwrap();

		// Ignore first line
		let lines = contents.lines().skip(1);

		// Parse the rest of the lines
		let mut particles = Vec::new();
		for line in lines {
			particles.push(Particle::parse(line));
		}

		assert!(nb_particles_local < particles.len());
		return Self {
			particles,
			nb_particles_local,
		};
	}

	/// Get the total number of particles in the [system](Self)
	pub fn nb_particles_total(&self) -> usize {
		self.particles.len()
	}

	/// Get the local number of particles in the [system](Self)
	pub fn nb_particles_local(&self) -> usize {
		self.nb_particles_local
	}

	/// Compute the distance between 2 [particles](Particle) of the [system](Self), squared.
	///
	/// # Arguments
	///
	/// * `i` - The index of the first particle
	/// * `j` - The index of the second particle
	pub fn distance_between_squared(&self, i: usize, j: usize) -> f64 {
		self.particles[i].distance_to_squared(&self.particles[j])
	}

	/// Reference (unoptimized) method for computing the microscopic energy in the system, according to the Lennard-Jones potential.
	#[allow(unused)]
	fn microscopic_energy_reference(&self) -> f64 {
		let mut total = 0.0;
		for i in 0..self.nb_particles_total() {
			for j in 0..self.nb_particles_total() {
				if i == j {
					continue;
				}
				let r_ij = self.distance_between_squared(i, j).sqrt();
				let u_ij = EPSILON_STAR * ((R_STAR / r_ij).powi(12) - 2.0 * (R_STAR / r_ij).powi(6));
				total += u_ij;
			}
		}

		return 2.0 * total;
	}

	/// Compute the microscopic energy in the system, according to the Lennard-Jones potential.
	pub fn microscopic_energy(&self) -> f64 {
		let mut total = 0.0;
		for i in 0..self.nb_particles_total() {
			for j in 0..self.nb_particles_total() {
				if i == j {
					continue;
				}
				let r_star_over_r_ij_pow2 = R_STAR.powi(2) / self.distance_between_squared(i, j);
				let r_star_over_r_ij_pow6 = r_star_over_r_ij_pow2.powi(3);
				let r_star_over_r_ij_pow12 = r_star_over_r_ij_pow6.powi(2);
				let u_ij = EPSILON_STAR * (r_star_over_r_ij_pow12 - (2.0 * r_star_over_r_ij_pow6));
				total += u_ij;
			}
		}

		return 2.0 * total;
	}

	/// Compute the forces between pairs of particles
	pub fn compute_forces(&self) -> Vec<Vec<Vector3>> {
		let mut forces = vec![vec![Vector3::zero(); self.nb_particles_total()]; self.nb_particles_total()];
		for i in 0..self.nb_particles_total() {
			for j in 0..self.nb_particles_total() {
				if i == j {
					// Force between a particle and itself is 0
					continue;
				}

				// Gradient function for any coordinate
				let r_star_over_r_ij_pow2 = R_STAR.powi(2) / self.distance_between_squared(i, j);
				let r_star_over_r_ij_pow8 = r_star_over_r_ij_pow2.powi(3);
				let r_star_over_r_ij_pow14 = r_star_over_r_ij_pow8.powi(7);
				let gradient = |c_i: f64, c_j: f64| {
					-48.0 * EPSILON_STAR * (r_star_over_r_ij_pow14 - r_star_over_r_ij_pow8) * (c_i - c_j)
				};

				// Apply gradient in the x, y, and z directions
				let (x_i, y_i, z_i) = self.particles[i].xyz();
				let (x_j, y_j, z_j) = self.particles[j].xyz();
				forces[i][j] = Vector3::from(gradient(x_i, x_j), gradient(y_i, y_j), gradient(z_i, z_j));
			}
		}

		return forces;
	}

	/// Compute the sum of all the forces between pairs of particles in the system
	pub fn sum_of_forces(forces: &Vec<Vec<Vector3>>) -> Vector3 {
		let mut sx = 0.0;
		let mut sy = 0.0;
		let mut sz = 0.0;

		for i in 0..forces.len() {
			for j in 0..forces.len() {
				let f = forces[i][j];
				sx += f.x();
				sy += f.y();
				sz += f.z();
			}
		}

		return Vector3::from(sx, sy, sz);
	}

	/// Get a reference to the particles in the system
	pub fn particles(&self) -> &[Particle] {
		&self.particles
	}
}

#[cfg(test)]
mod tests {
	use crate::assert_approx_eq;

	use super::*;

	#[test]
	fn check_energy_optimizations() {
		let system = System::from_file(Path::new("dataset/particles.xyz"), 0);
		assert_approx_eq!(system.microscopic_energy_reference(), system.microscopic_energy());
	}
}
