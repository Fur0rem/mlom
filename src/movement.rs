use crate::smoothing::quintic_smoothing_and_derivative;
use crate::system::Particle;
use crate::{algebra::Vector3, parameters::*, periodic_conditions::neighboring_3d_translations, system::System};
use plotters::prelude::*;
use plotters::prelude::{RED, WHITE};

impl System {
	/// Compute the degrees of liberty of the system. Noted N_dl
	pub fn degrees_of_liberty(&self) -> f64 {
		(3 * self.nb_particles_total() - 3) as f64
	}

	/// Calibrate the kinetic momentum of the particles to have the desired temperature
	pub fn recalibrate_according_to_temperature(&mut self) {
		let initial_ke = self.kinetic_energy_and_temperature().0;
		let desired_ke = self.degrees_of_liberty() * R_CONSTANT * T_0;

		let scale = (desired_ke / initial_ke).sqrt();

		for p in self.particles.iter_mut() {
			p.momentum *= scale;
		}
	}

	/// Calibrate the kinetic momentum of the particles to conserve the kinetic momentum of the center of mass
	pub fn recalibrate_according_to_center_of_mass(&mut self) {
		// Compute the average momentum of the particles
		let mut sum_momentum = Vector3::zero();
		for particle in self.particles() {
			sum_momentum += particle.momentum;
		}
		let avg_momentum = sum_momentum / self.nb_particles_total() as f64;

		// Recalibrate the momentum of each particle by removing the average momentum
		for particle in self.particles.iter_mut() {
			particle.momentum -= avg_momentum;
		}
	}

	pub fn init_particles_momentums(&mut self) {
		// Step 1: Set momentums to random vectors in unit cube
		for particle in self.particles.iter_mut() {
			particle.momentum = Vector3::random_in_unit_cube();
		}

		// Step 2: Recalibrate the momentums to have the right initial temperature
		self.recalibrate_according_to_temperature();
		self.recalibrate_according_to_center_of_mass();
		self.recalibrate_according_to_temperature();
	}

	pub fn kinetic_energy_and_temperature(&self) -> (f64, f64) {
		// Compute kinetic energy: K = sum_i (p_i^2 / m)
		let mut sum_kinetic_energy = 0.0;
		for particle in self.particles() {
			let p = particle.momentum;
			sum_kinetic_energy += p.x().powi(2) + p.y().powi(2) + p.z().powi(2);
		}
		sum_kinetic_energy /= PARTICLE_MASS;

		// Convert kinetic energy to real units
		let kinetic_energy = sum_kinetic_energy / (2.0 * CONVERSION_FORCE);

		// Compute temperature: T = 1 / (N_dl * R) * K
		let temperature = kinetic_energy / (self.degrees_of_liberty() * R_CONSTANT);

		return (kinetic_energy, temperature);
	}

	pub fn compute_forces_periodic(&self) -> Vec<Vector3> {
		let mut forces = vec![Vector3::zero(); self.nb_particles_total()];

		// Iterate particle pairs and apply equal-and-opposite forces.
		for i in 0..self.nb_particles_total() {
			for j in 0..self.nb_particles_total() {
				for sym in neighboring_3d_translations(BOX_SIDE) {
					// Skip self-interaction for the original particle
					if i == j && sym == Vector3::zero() {
						continue;
					}

					// Process each physical pair once: only (i, j) with sym but not (j, i) with -sym
					if !(j > i || (j == i && sym != Vector3::zero())) {
						continue;
					}

					let j_with_sym = Particle {
						coordinates: (self.particles[j].coordinates + sym).as_point(),
						momentum: self.particles[j].momentum,
					};
					let dist_ij_squared = self.particles[i].distance_to_squared(&j_with_sym);

					// Cutoff radius
					if dist_ij_squared > R_MAX * R_MAX {
						continue;
					}

					let r_r2 = R_STAR.powi(2) / dist_ij_squared;
					let r_r6 = r_r2 * r_r2 * r_r2;
					let r_r12 = r_r6 * r_r6;

					// P5 smoothing factors (keeps U and F continuous near R_CUT)
					let radius = dist_ij_squared.sqrt();
					assert!(radius > 0.0);
					let (p5, p5_derivative) = quintic_smoothing_and_derivative(radius);

					// Main gradient term
					let inv_r2 = 1.0 / dist_ij_squared;
					let gradient = -48.0 * EPSILON_STAR * (r_r12 - r_r6) * inv_r2 * p5;

					// Distance components
					let dx = self.particles[i].x() - j_with_sym.x();
					let dy = self.particles[i].y() - j_with_sym.y();
					let dz = self.particles[i].z() - j_with_sym.z();

					// P5 derivative term
					let factor = p5_derivative * (4.0 * EPSILON_STAR * (r_r12 - 2.0 * r_r6)) / radius;

					// Total gradient for each coordinate
					let grad_x = gradient * dx + factor * dx;
					let grad_y = gradient * dy + factor * dy;
					let grad_z = gradient * dz + factor * dz;

					// Apply equal and opposite forces
					forces[i].x += grad_x;
					forces[i].y += grad_y;
					forces[i].z += grad_z;

					forces[j].x -= grad_x;
					forces[j].y -= grad_y;
					forces[j].z -= grad_z;
				}
			}
		}

		return forces;
	}

	/// Update the momentums of the particles according to the forces applied to them, for the velocity Verlet algorithm
	#[inline(always)]
	pub fn velocity_verlet_momentums_update(&mut self, forces: &Vec<Vector3>) {
		// 1st equation: half time step update of the kinetic momentum
		// d(r_i) / d_t is momentum
		// F_i = -Nabla U => momentum(t) + 1/2 Nabla U(t) dt = momentum(t) - 1/2 F_i dt
		for i in 0..self.nb_particles_total() {
			self.particles[i].momentum = self.particles[i].momentum - forces[i] * DELTA_TIME * 0.5 * CONVERSION_FORCE;
		}
	}

	/// Update the coordinates of the particles according to their momentums, for the velocity Verlet algorithm
	#[inline(always)]
	pub fn velocity_verlet_coordinates_update(&mut self) {
		// According to Newton's equations: m_i * momentum = F_i
		for i in 0..self.nb_particles_total() {
			let offset = (self.particles[i].momentum * DELTA_TIME) / PARTICLE_MASS;
			self.particles[i].coordinates = (self.particles[i].coordinates + offset).as_point();
			self.particles[i].put_back_in_box();
		}
	}

	/// A step in the simulation
	/// Applies the velocity Verlet algorithm to update the coordinates and momentums of the particles, with periodic conditions.
	/// If `correct_with_temperature` is Some, it applies the Berendsen thermostat correction with the given target temperature after the full step update.
	pub fn step(&mut self, correct_with_temperature: Option<f64>) {
		// // INFO: max force magnitude and max particle momentum before update
		// let forces = self.compute_forces_periodic();
		// let max_force = forces.iter().map(|f| f.norm()).fold(0.0, f64::max);
		// let max_momentum_before = self.particles.iter().map(|p| p.momentum.norm()).fold(0.0, f64::max);
		// println!("INFO: max_force = {}, max_momentum_before = {}", max_force, max_momentum_before);

		// // INFO: minimal pair distance (considering periodic images)
		// let mut min_pair_dist2 = std::f64::INFINITY;
		// let mut min_pair = (0usize, 0usize);
		// for sym in neighboring_3d_translations(BOX_SIDE) {
		// 	for i in 0..self.nb_particles_total() {
		// 		for j in 0..self.nb_particles_total() {
		// 			if i == j {
		// 				continue;
		// 			}
		// 			let particle_j_with_symmetry = (self.particles[j].coordinates + sym).as_point();
		// 			let dist2 = self.particles[i].coordinates.distance_to_squared(&particle_j_with_symmetry);
		// 			if dist2 < min_pair_dist2 {
		// 				min_pair_dist2 = dist2;
		// 				min_pair = (i, j);
		// 			}
		// 		}
		// 	}
		// }
		// let min_pair_distance = min_pair_dist2.sqrt();
		// println!("INFO: min_pair_distance = {}, min_pair = {:?}", min_pair_distance, min_pair);

		let forces = self.compute_forces_periodic();
		self.velocity_verlet_momentums_update(&forces);

		self.velocity_verlet_coordinates_update();

		let forces = self.compute_forces_periodic();
		self.velocity_verlet_momentums_update(&forces);

		// INFO: max force magnitude and max particle momentum after update
		// let max_force_after = forces.iter().map(|f| f.norm()).fold(0.0, f64::max);
		// let max_momentum_after = self.particles.iter().map(|p| p.momentum.norm()).fold(0.0, f64::max);
		// println!("INFO: max_force_after = {}, max_momentum_after = {}", max_force_after, max_momentum_after);

		// Correct the momentums with the Berendsen thermostat if requested
		if let Some(target_temperature) = correct_with_temperature {
			let current_temperature = self.kinetic_energy_and_temperature().1;
			let factor = GAMMA * ((target_temperature / current_temperature) - 1.0);
			for p in self.particles.iter_mut() {
				p.momentum.x += factor * p.momentum.x;
				p.momentum.y += factor * p.momentum.y;
				p.momentum.z += factor * p.momentum.z;
			}
		}
	}

	/// Compute the total energy of the system (kinetic + potential) and its temperature
	pub fn total_energy_and_temperature(&self) -> (f64, f64) {
		// Calculate kinetic energy
		let (kinetic_energy, temperature) = self.kinetic_energy_and_temperature();

		// Calculate potential energy using the periodic conditions
		let potential_energy = self.microscopic_energy_periodic(&neighboring_3d_translations(BOX_SIDE), R_CUT);

		return (kinetic_energy + potential_energy, temperature);
	}

	/// Simulate the system for a given number of steps, applying the velocity Verlet algorithm with periodic conditions.
	/// If `correct_each` is not 0, it applies the Berendsen thermostat correction every `correct_each` steps with the target temperature T_0.
	/// It saves the evolution of the total energy and temperature during the simulation in a plot at the given path `save_to`.
	pub fn simulate(&mut self, nb_steps: usize, correct_each: usize, save_to: &str) {
		let mut energies = vec![];
		let mut temperatures = vec![];
		for step in 0..nb_steps {
			// Correct only every `correct_each` steps
			if step % correct_each == 0 && step != 0 {
				self.step(Some(T_0));
			} else {
				self.step(None);
			}
			let (total_energy, temperature) = self.total_energy_and_temperature();
			println!(
				"Step {:<5}: Total energy = {:<14.8}  |  Temperature = {:<15.6}",
				step, total_energy, temperature
			);
			energies.push(total_energy);
			temperatures.push(temperature);
		}

		let root = BitMapBackend::new(save_to, (1600, 600)).into_drawing_area();
		let (left, right) = root.split_horizontally(800);

		// Energy plot
		left.fill(&WHITE).unwrap();
		let mut chart_energy = ChartBuilder::on(&left)
			.caption("Energy Evolution", ("sans-serif", 40).into_font())
			.margin(20)
			.x_label_area_size(30)
			.y_label_area_size(40)
			.build_cartesian_2d(
				0..nb_steps,
				*energies.iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap()
					..*energies.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap(),
			)
			.unwrap();

		chart_energy.configure_mesh().draw().unwrap();

		chart_energy
			.draw_series(plotters::series::LineSeries::new(
				energies.iter().enumerate().map(|(x, y)| (x, *y)),
				&RED,
			))
			.unwrap()
			.label("Total Energy")
			.legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));
		chart_energy
			.configure_series_labels()
			.background_style(&WHITE.mix(0.8))
			.draw()
			.unwrap();

		// Temperature plot
		right.fill(&WHITE).unwrap();
		let mut chart_temp = ChartBuilder::on(&right)
			.caption("Temperature Evolution", ("sans-serif", 40).into_font())
			.margin(20)
			.x_label_area_size(30)
			.y_label_area_size(40)
			.build_cartesian_2d(
				0..nb_steps,
				*temperatures.iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap()
					..*temperatures.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap(),
			)
			.unwrap();

		chart_temp.configure_mesh().draw().unwrap();

		chart_temp
			.draw_series(plotters::series::LineSeries::new(
				temperatures.iter().enumerate().map(|(x, y)| (x, *y)),
				&BLUE,
			))
			.unwrap()
			.label("Temperature")
			.legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));
		chart_temp
			.configure_series_labels()
			.background_style(&WHITE.mix(0.8))
			.draw()
			.unwrap();
	}
}
