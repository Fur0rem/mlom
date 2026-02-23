use crate::parameters::{INVERSE_DIFF_MAX_MIN, R_MAX, R_MIN};

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
