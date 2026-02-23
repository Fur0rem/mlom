//! The system of particles

use std::{fs::File, io::Read, path::Path};

use crate::{
	algebra::{Point3, Vector3},
	parameters::*,
};

/// A particle in the system
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Particle {
	/// The coordinates of the particle
	pub(crate) coordinates: Point3,
	/// The momentum of the particle
	pub(crate) momentum: Vector3,
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

		let momentum = Vector3::zero(); // 0 for now

		return Self {
			coordinates: Point3::from(x, y, z),
			momentum,
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

	/// Compute the kinetic moment of the [particle](Self)
	pub fn kinetic_moment(&self) -> Vector3 {
		return self.momentum;
	}

	/// Put the particle back in the box
	pub fn put_back_in_box(&mut self) {
		self.coordinates.x = (self.coordinates.x + BOX_SIDE / 2.0).rem_euclid(BOX_SIDE) - BOX_SIDE / 2.0;
		self.coordinates.y = (self.coordinates.y + BOX_SIDE / 2.0).rem_euclid(BOX_SIDE) - BOX_SIDE / 2.0;
		self.coordinates.z = (self.coordinates.z + BOX_SIDE / 2.0).rem_euclid(BOX_SIDE) - BOX_SIDE / 2.0;

		assert!(
			self.coordinates.x <= BOX_SIDE / 2.0,
			"Particle x coordinate is out of bounds: {}",
			self.coordinates.x
		);
		assert!(
			self.coordinates.x >= -BOX_SIDE / 2.0,
			"Particle x coordinate is out of bounds: {}",
			self.coordinates.x
		);
		assert!(
			self.coordinates.y <= BOX_SIDE / 2.0,
			"Particle y coordinate is out of bounds: {}",
			self.coordinates.y
		);
		assert!(
			self.coordinates.y >= -BOX_SIDE / 2.0,
			"Particle y coordinate is out of bounds: {}",
			self.coordinates.y
		);
		assert!(
			self.coordinates.z <= BOX_SIDE / 2.0,
			"Particle z coordinate is out of bounds: {}",
			self.coordinates.z
		);
		assert!(
			self.coordinates.z >= -BOX_SIDE / 2.0,
			"Particle z coordinate is out of bounds: {}",
			self.coordinates.z
		);
	}
}

/// A system of [particles](Particle)
#[derive(Debug, Clone, PartialEq)]
pub struct System {
	/// The [particles](Particle) in the system
	pub(crate) particles: Vec<Particle>,
	/// The number of local particles (unused for now)
	pub(crate) nb_particles_local: usize,
}

impl System {
	/// Parse a system from a string
	///
	/// # Arguments
	/// * `s` - The string to parse
	/// * `nb_particles_local` - Unused.
	pub fn from_str(s: &str, nb_particles_local: usize) -> Self {
		// Ignore first line
		let lines = s.lines().skip(1);

		// Parse the rest of the lines
		let mut particles = Vec::new();
		for line in lines {
			particles.push(Particle::parse(line));
		}

		// Create the system
		assert!(nb_particles_local < particles.len());
		let mut system = Self {
			particles,
			nb_particles_local,
		};

		// Initialize the particle momentums
		system.init_particles_momentums();

		return system;
	}

	/// Parse a system from a file
	///
	/// # Arguments
	/// * `path` - The path to the file to parse
	/// * `nb_particles_local` - Unused.
	pub fn from_file(path: &Path, nb_particles_local: usize) -> Self {
		// Read the file
		let mut file = File::open(path).unwrap();
		let mut contents = String::new();
		file.read_to_string(&mut contents).unwrap();

		return Self::from_str(&contents, nb_particles_local);
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

	/// Compute the distance between 2 [particles](Particle) of the [system](Self).
	///
	/// # Arguments
	///
	/// * `i` - The index of the first particle
	/// * `j` - The index of the second particle
	pub fn distance_between(&self, i: usize, j: usize) -> f64 {
		self.particles[i].distance_to_squared(&self.particles[j]).sqrt()
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
