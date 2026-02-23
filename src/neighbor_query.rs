use crate::{
	algebra::Vector3,
	energy::force_between_particles,
	parameters::*,
	system::{Particle, System},
};

/// Estimate the maximum number of neighbors per particle based on the density and cut radius, to dimension neighbor lists.
pub fn max_number_of_neighbors(nb_particles: usize, box_side: f64, cut_radius: f64) -> usize {
	// Compute a conservative "skin" distance that particles can travel and therefore still be relevant for the next rebuild window.
	// Two particles can approach by up to 2 * vmax * n_steps * dt per rebuild period,
	let skin = 2.0 * MAX_PARTICLE_VELOCITY * REBUILD_VERLET_LISTS_FREQUENCY as f64 * DELTA_TIME;
	let effective_cut = cut_radius + skin;

	let density = nb_particles as f64 / box_side.powi(3);
	let volume_cut_sphere = (4.0 / 3.0) * std::f64::consts::PI * effective_cut.powi(3);
	let safe_factor = 2.0;

	// Estimate cannot exceed total particles-1 (a particle cannot have itself as neighbor)
	let estimate = (density * volume_cut_sphere * safe_factor).ceil() as usize;
	println!(
		"Estimated max number of neighbors per particle: {}, with effective cut radius: {:.3}",
		estimate, effective_cut
	);
	let capped = std::cmp::min(estimate, nb_particles - 1);
	return std::cmp::max(1, capped);
}

#[derive(Debug, Clone, PartialEq)]
pub struct Neighbor {
	pub index: usize,
	pub symmetry: Vector3,
}

#[derive(Debug, Clone, PartialEq)]
pub struct VerletList {
	/// Flattened list of neighbors for all particles, assumed to be of size `nb_particles * max_number_of_neighbors`
	pub neighbors: Vec<Neighbor>,

	/// Number of neighbors per particle, used to index into the flattened `neighbors` vector
	pub number_of_neighbors_per_particle: Vec<usize>,
}

impl VerletList {
	pub fn neighbors_of_particle(&self, particle_index: usize, max_number_of_neighbors: usize) -> &[Neighbor] {
		let start = particle_index * max_number_of_neighbors;
		let end = start + self.number_of_neighbors_per_particle[particle_index];
		return &self.neighbors[start..end];
	}

	pub fn build(system: &System, translations: &[Vector3], cut_radius: f64, max_number_of_neighbors: usize) -> Self {
		let mut neighbors = vec![
			Neighbor {
				index: 0,
				symmetry: Vector3::zero()
			};
			system.nb_particles_total() * max_number_of_neighbors
		];
		let mut number_of_neighbors_per_particle = vec![0; system.nb_particles_total()];
		for i in 0..system.nb_particles_total() {
			let mut count_neighbors = 0;
			for j in i + 1..system.nb_particles_total() {
				for sym in translations {
					let particle_j_with_symmetry = (system.particles[j].coordinates + *sym).as_point();
					let dist_ij_squared =
						system.particles[i].coordinates.distance_to_squared(&particle_j_with_symmetry);

					// Put in the neighbor list if within cut radius
					if dist_ij_squared < cut_radius.powi(2) {
						if count_neighbors >= max_number_of_neighbors {
							panic!(
								"Number of neighbors for particle {} ({}) exceeds the maximum anticipated number of neighbors ({}).",
								i, count_neighbors, max_number_of_neighbors
							);
						}
						neighbors[i * max_number_of_neighbors + count_neighbors] =
							Neighbor { index: j, symmetry: *sym };
						count_neighbors += 1;
					}
				}
			}
			number_of_neighbors_per_particle[i] = count_neighbors;
		}
		return Self {
			neighbors,
			number_of_neighbors_per_particle,
		};
	}
}

impl System {
	pub fn compute_forces_periodic_with_neighbor_lists(
		&self, verlet_list: &VerletList, max_number_of_neighbors: usize,
	) -> Vec<Vector3> {
		let mut forces = vec![Vector3::zero(); self.nb_particles_total()];

		// Iterate particle pairs and apply equal-and-opposite forces.
		for i in 0..self.nb_particles_total() {
			for neighbor in verlet_list.neighbors_of_particle(i, max_number_of_neighbors) {
				let j = neighbor.index;
				let sym = neighbor.symmetry;

				// Skip self-interaction for the original particle
				if i == j && sym == Vector3::zero() {
					continue;
				}

				// Process each physical pair once: only (i, j) with sym but not (j, i) with -sym
				if !(j > i || (j == i && sym != Vector3::zero())) {
					continue;
				}

				// Compute translated particle j
				let j_with_sym = Particle {
					coordinates: (self.particles[j].coordinates + sym).as_point(),
					momentum: self.particles[j].momentum,
				};

				// Apply equal and opposite forces to particles i and j
				let force_ij = force_between_particles(&self.particles[i], &j_with_sym);
				forces[i] += force_ij;
				forces[j] -= force_ij;
			}
		}

		return forces;
	}
}
