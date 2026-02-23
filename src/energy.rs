use crate::{algebra::Vector3, parameters::*, system::Particle};

/// P5 smoothing function for the potential and its derivative, to ensure continuity of the potential and forces near the cutoff radius.
pub fn quintic_smoothing_and_derivative(radius: f64) -> (f64, f64) {
	let mut p5 = 1.0;
	let mut p5_derivative = 0.0;
	if radius > R_MIN && radius < R_MAX {
		let r1 = (radius - R_MIN) * (INVERSE_DIFF_MAX_MIN);
		let r2 = r1 * r1;
		let r3 = r1 * r2;
		let r4 = r2 * r2;
		let r5 = r2 * r3;
		p5 = 1.0 - 10.0 * r3 + 15.0 * r4 - 6.0 * r5;
		p5_derivative = -30.0 * INVERSE_DIFF_MAX_MIN * r2 + 60.0 * INVERSE_DIFF_MAX_MIN * r3 - 30.0 * INVERSE_DIFF_MAX_MIN * r4;
	}

	return (p5, p5_derivative);
}

/// P5 smoothing function
pub fn quintic_smoothing(radius: f64) -> f64 {
	let mut p5 = 1.0;
	if radius > R_MIN && radius < R_MAX {
		let r1 = (radius - R_MIN) * (INVERSE_DIFF_MAX_MIN);
		let r2 = r1 * r1;
		let r3 = r1 * r2;
		let r4 = r2 * r2;
		let r5 = r2 * r3;
		p5 = 1.0 - 10.0 * r3 + 15.0 * r4 - 6.0 * r5;
	}

	return p5;
}

/// Derivative of the P5 smoothing function
pub fn quintic_smoothing_derivative(radius: f64) -> f64 {
	let mut p5_derivative = 0.0;
	if radius > R_MIN && radius < R_MAX {
		let r1 = (radius - R_MIN) * (INVERSE_DIFF_MAX_MIN);
		let r2 = r1 * r1;
		let r3 = r1 * r2;
		let r4 = r2 * r2;
		p5_derivative = -30.0 * INVERSE_DIFF_MAX_MIN * r2 + 60.0 * INVERSE_DIFF_MAX_MIN * r3 - 30.0 * INVERSE_DIFF_MAX_MIN * r4;
	}

	return p5_derivative;
}

#[allow(dead_code)]
fn energy_between_particles_smooth(particle_i: &Particle, particle_j: &Particle, radius_cut: f64) -> f64 {
	// Apply cut above given radius
	let dist_ij_squared = particle_i.distance_to_squared(particle_j);
	assert!(dist_ij_squared > 0.0001);

	// Apply cut above given radius
	if dist_ij_squared > radius_cut.powi(2) {
		return 0.0;
	}

	// The usual energy term
	let r_star_over_r_ij_pow_2 = (R_STAR * R_STAR) / dist_ij_squared;
	let r_star_over_r_ij_pow6 = r_star_over_r_ij_pow_2 * r_star_over_r_ij_pow_2 * r_star_over_r_ij_pow_2;
	let r_star_over_r_ij_pow12 = r_star_over_r_ij_pow6 * r_star_over_r_ij_pow6;
	let u_ij = r_star_over_r_ij_pow12 - (2.0 * r_star_over_r_ij_pow6);

	// Apply P5 smoothing near cutoff (keeps U and F continuous near R_CUT)
	let radius = dist_ij_squared.sqrt();
	return u_ij * quintic_smoothing(radius);
}

#[allow(dead_code)]
fn force_between_particles_smooth(particle_i: &Particle, particle_j: &Particle) -> Vector3 {
	let dist_ij_squared = particle_i.distance_to_squared(particle_j);

	// Cutoff radius
	if dist_ij_squared > R_MAX * R_MAX {
		return Vector3::zero();
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
	let dx = particle_i.x() - particle_j.x();
	let dy = particle_i.y() - particle_j.y();
	let dz = particle_i.z() - particle_j.z();

	// P5 derivative term
	let factor = p5_derivative * (4.0 * EPSILON_STAR * (r_r12 - 2.0 * r_r6)) / radius;

	// Total gradient for each coordinate
	let grad_x = gradient * dx + factor * dx;
	let grad_y = gradient * dy + factor * dy;
	let grad_z = gradient * dz + factor * dz;

	return Vector3::from(grad_x, grad_y, grad_z);
}

#[allow(dead_code)]
fn energy_between_particles_clear_cut(particle_i: &Particle, particle_j: &Particle, radius_cut: f64) -> f64 {
	// Apply cut above given radius
	let dist_ij_squared = particle_i.distance_to_squared(particle_j);
	assert!(dist_ij_squared > 0.0001);

	// Apply cut above given radius
	if dist_ij_squared > radius_cut.powi(2) {
		return 0.0;
	}

	// The usual energy term
	let r_star_over_r_ij_pow_2 = (R_STAR * R_STAR) / dist_ij_squared;
	let r_star_over_r_ij_pow6 = r_star_over_r_ij_pow_2 * r_star_over_r_ij_pow_2 * r_star_over_r_ij_pow_2;
	let r_star_over_r_ij_pow12 = r_star_over_r_ij_pow6 * r_star_over_r_ij_pow6;
	return r_star_over_r_ij_pow12 - (2.0 * r_star_over_r_ij_pow6);
}

#[allow(dead_code)]
fn force_between_particles_clear_cut(particle_i: &Particle, particle_j: &Particle) -> Vector3 {
	let dist_ij_squared = particle_i.distance_to_squared(particle_j);

	// Cutoff radius
	if dist_ij_squared > R_MAX * R_MAX {
		return Vector3::zero();
	}

	let r_r2 = R_STAR.powi(2) / dist_ij_squared;
	let r_r6 = r_r2 * r_r2 * r_r2;
	let r_r12 = r_r6 * r_r6;

	// Main gradient term
	let inv_r2 = 1.0 / dist_ij_squared;
	let gradient = -48.0 * EPSILON_STAR * (r_r12 - r_r6) * inv_r2;

	// Distance components
	let dx = particle_i.x() - particle_j.x();
	let dy = particle_i.y() - particle_j.y();
	let dz = particle_i.z() - particle_j.z();

	// Total gradient for each coordinate
	let grad_x = gradient * dx;
	let grad_y = gradient * dy;
	let grad_z = gradient * dz;

	return Vector3::from(grad_x, grad_y, grad_z);
}

pub fn energy_between_particles(particle_i: &Particle, particle_j: &Particle, radius_cut: f64) -> f64 {
	return energy_between_particles_smooth(particle_i, particle_j, radius_cut);
}

pub fn force_between_particles(particle_i: &Particle, particle_j: &Particle) -> Vector3 {
	return force_between_particles_smooth(particle_i, particle_j);
}
