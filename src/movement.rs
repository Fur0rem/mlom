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
		if initial_ke <= 0.0 {
			return;
		}

		// Desired kinetic energy: K_desired = (N_dl / 2) * k_B * T0
		let desired_ke = 0.5 * self.degrees_of_liberty() * R_CONSTANT * T_0;
		let scale = (desired_ke / initial_ke).sqrt();

		for p in self.particles.iter_mut() {
			p.momentum *= scale;
		}
	}

	/// Calibrate the kinetic momentum of the particles to conserve the kinetic momentum of the center of mass
	pub fn recalibrate_according_to_center_of_mass(&mut self) {
		// Compute center of mass
		let mut center = Vector3::zero();
		for p in self.particles() {
			center += p.momentum;
		}
		center /= self.nb_particles_total() as f64;

		// Modify particles momentums
		for p in self.particles.iter_mut() {
			p.momentum -= center;
		}
	}

	pub fn init_particles_momentums(&mut self) {
		// Step 1: Random vectors in unit cube
		for particle in self.particles.iter_mut() {
			particle.momentum = Vector3::random_in_unit_cube();
		}

		// Step 2: Recalibrate
		self.recalibrate_according_to_temperature();
		self.recalibrate_according_to_center_of_mass();
		self.recalibrate_according_to_temperature();
	}

	pub fn kinetic_energy_and_temperature(&self) -> (f64, f64) {
		let mut sum_p2 = 0.0;
		for particle in self.particles() {
			let p = particle.kinetic_moment();
			sum_p2 += p.x().powi(2) + p.y().powi(2) + p.z().powi(2);
		}

		// Physical kinetic energy: K = sum(p^2) / (2 m)
		let kinetic_energy = sum_p2 / (2.0 * PARTICLE_MASS);

		// Temperature: K = (N_dl / 2) * k_B * T  =>  T = 2K / (N_dl * k_B)
		let temperature = 2.0 * kinetic_energy / (self.degrees_of_liberty() * R_CONSTANT);

		return (kinetic_energy, temperature);
	}

	pub fn forces_applied_to_particles(forces: &Vec<Vec<Vec<Vector3>>>) -> Vec<Vector3> {
		let nb_particles = forces[0].len();
		let mut flattened_forces = vec![Vector3::zero(); nb_particles];
		for sym_idx in 0..forces.len() {
			for i in 0..nb_particles {
				for j in 0..nb_particles {
					flattened_forces[i] += forces[sym_idx][i][j];
				}
			}
		}
		return flattened_forces;
	}

	pub fn step(&mut self) {
		// Compute forces applied to each particle
		let forces = self.compute_forces_periodic(&neighboring_3d_translations(BOX_SIDE), R_CUT);
		let forces = Self::forces_applied_to_particles(&forces);

		// INFO: max force magnitude and max particle momentum before update
		let max_force = forces.iter().map(|f| f.norm()).fold(0.0, f64::max);
		let max_momentum_before = self.particles.iter().map(|p| p.momentum.norm()).fold(0.0, f64::max);
		println!("INFO: max_force = {}, max_momentum_before = {}", max_force, max_momentum_before);

		// INFO: minimal pair distance (considering periodic images)
		let mut min_pair_dist2 = std::f64::INFINITY;
		let mut min_pair = (0usize, 0usize);
		for sym in neighboring_3d_translations(BOX_SIDE) {
			for i in 0..self.nb_particles_total() {
				for j in 0..self.nb_particles_total() {
					if i == j {
						continue;
					}
					let particle_j_with_symmetry = (self.particles[j].coordinates + sym).as_point();
					let dist2 = self.particles[i].coordinates.distance_to_squared(&particle_j_with_symmetry);
					if dist2 < min_pair_dist2 {
						min_pair_dist2 = dist2;
						min_pair = (i, j);
					}
				}
			}
		}
		let min_pair_distance = min_pair_dist2.sqrt();
		println!("INFO: min_pair_distance = {}, min_pair = {:?}", min_pair_distance, min_pair);

		// 1st equation: half time step update of the kinetic momentum
		// d(r_i) / d_t is momentum
		// F_i = -Nabla U => momentum(t) - 1/2 Nabla U(t) dt = momentum(t) + 1/2 F_i dt
		// TODO: Verify why it explodes when it's += but stays stable with -=
		for (i, p) in self.particles.iter_mut().enumerate() {
			p.momentum = p.momentum - 0.5 * forces[i] * DELTA_TIME * CONVERSION_FORCE;
		}

		// 2nd equation: full time step update of the positions
		// According to Newton's equations: m_i * momentum = F_i
		for p in self.particles.iter_mut() {
			let velocity = p.momentum / PARTICLE_MASS;
			p.coordinates = (p.coordinates + velocity * DELTA_TIME).as_point();
		}

		// 3rd equation: full time step update of the kinetic momentum
		// Before, compute the energy at the next time step and forces applied to each particle
		// TODO: Same as 1st equation
		let forces = self.compute_forces_periodic(&neighboring_3d_translations(BOX_SIDE), R_CUT);
		let forces = Self::forces_applied_to_particles(&forces);

		// INFO: max force (after position update)
		let max_force_after = forces.iter().map(|f| f.norm()).fold(0.0, f64::max);
		println!("INFO: max_force_after = {}", max_force_after);

		for (i, p) in self.particles.iter_mut().enumerate() {
			p.momentum = p.momentum - 0.5 * forces[i] * DELTA_TIME * CONVERSION_FORCE;
		}

		// INFO: max particle momentum after full update
		let max_momentum_after = self.particles.iter().map(|p| p.momentum.norm()).fold(0.0, f64::max);
		println!("INFO: max_momentum_after = {}", max_momentum_after);

		// Periodic conditions: put the particles in the box
		for p in self.particles.iter_mut() {
			p.put_back_in_box();
		}
	}

	pub fn total_energy(&self) -> f64 {
		let (kinetic_energy, _temp) = self.kinetic_energy_and_temperature();

		// Calculate potential energy using the periodic conditions
		let potential_energy = self.microscopic_energy_periodic(&neighboring_3d_translations(BOX_SIDE), R_CUT);

		println!("k: {kinetic_energy}, t: {_temp}, p: {potential_energy}");

		return kinetic_energy + potential_energy;
	}

	pub fn energy_evolution(&mut self, nb_steps: usize, save_to: &str) {
		let mut energies = vec![];
		for step in 0..nb_steps {
			self.step();
			println!("Step {}: Total energy = {}", step, self.total_energy());
			let total_energy = self.total_energy();
			energies.push(total_energy);
		}

		let root = BitMapBackend::new(save_to, (800, 600)).into_drawing_area();
		root.fill(&WHITE).unwrap();
		let mut chart = ChartBuilder::on(&root)
			.caption("Energy Evolution", ("sans-serif", 50).into_font())
			.margin(20)
			.x_label_area_size(30)
			.y_label_area_size(40)
			.build_cartesian_2d(
				0..nb_steps,
				*energies.iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap()
					..*energies.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap(),
			)
			.unwrap();

		chart.configure_mesh().draw().unwrap();

		chart.draw_series(plotters::series::LineSeries::new(
			energies.iter().enumerate().map(|(x, y)| (x, *y)),
			&RED,
		))
		.unwrap()
		.label("Total Energy")
		.legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));
		chart.configure_series_labels().background_style(&WHITE.mix(0.8)).draw().unwrap();
	}
}
