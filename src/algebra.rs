//! Linear algebra operations needed for the simulation

use rand::{Rng, distr::Uniform, rng};

/// Checks if 2 floats are approximately equal.
#[macro_export]
macro_rules! assert_approx_eq {
	($L: expr, $R: expr) => {
		let difference = $L - $R;
		assert!(difference.abs() < 1e-5, "{}", format!("{} =/= {}", $L, $R));
	};
}

/// Checks if 2 vectors are approximately equal.
#[macro_export]
macro_rules! assert_vector_approx_eq {
	($L: expr, $R: expr) => {
		let difference = $L - $R;
		assert!(difference.x() < 1e-5, "{}", format!("{:?} =/= {:?}", $L, $R));
		assert!(difference.y() < 1e-5, "{}", format!("{:?} =/= {:?}", $L, $R));
		assert!(difference.z() < 1e-5, "{}", format!("{:?} =/= {:?}", $L, $R));
	};
}

/// A point in 3D space
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point3 {
	/// x coordinate
	pub(crate) x: f64,
	/// y coordinate
	pub(crate) y: f64,
	/// z coordinate
	pub(crate) z: f64,
}

impl Point3 {
	/// Create a point from the given x y z coordinates
	pub fn from(x: f64, y: f64, z: f64) -> Self {
		Self { x, y, z }
	}

	/// The x coordinate of the point
	pub fn x(&self) -> f64 {
		self.x
	}

	/// The y coordinate of the point
	pub fn y(&self) -> f64 {
		self.y
	}

	/// The z coordinate of the point
	pub fn z(&self) -> f64 {
		self.z
	}

	/// Create a point at the origin (0, 0, 0)
	pub fn origin() -> Self {
		Self { x: 0.0, y: 0.0, z: 0.0 }
	}

	/// Compute the distance to another [point](Self), squared
	pub fn distance_to_squared(&self, rhs: &Self) -> f64 {
		let dist_squared = (self.x - rhs.x).powi(2) + (self.y - rhs.y).powi(2) + (self.z - rhs.z).powi(2);

		// In the context of a simulation, having a distance of 0 is a bad sign.
		// So even if it doesn't make sense to check if in here specifically, I am likely to forget to check in the rest of the code.
		assert!(dist_squared != 0.0);

		return dist_squared;
	}

	/// Compute the distance to another [point](Self)
	pub fn distance_to(&self, rhs: &Self) -> f64 {
		self.distance_to_squared(rhs).sqrt()
	}
}

/// A vector in 3D space
#[derive(Clone, Copy, PartialEq, PartialOrd)]
pub struct Vector3 {
	/// x component
	x: f64,
	/// y component
	y: f64,
	/// z component
	z: f64,
}

impl std::fmt::Debug for Vector3 {
	fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
		write!(f, "(x: {}, y: {}, z: {})", self.x, self.y, self.z)
	}
}

impl Vector3 {
	/// Create a vector from the given x y z components
	pub fn from(x: f64, y: f64, z: f64) -> Self {
		Self { x, y, z }
	}

	/// Create a random vector in the unit cube, i.e. each component is in the range [-1, 1]
	pub fn random_in_unit_cube() -> Self {
		let distribution = Uniform::new_inclusive(-1.0f64, 1.0f64).unwrap();
		let mut rng = rng();

		return Vector3 {
			x: rng.sample(distribution),
			y: rng.sample(distribution),
			z: rng.sample(distribution),
		};
	}

	/// The x component of the vector
	pub fn x(&self) -> f64 {
		self.x
	}

	/// The y component of the vector
	pub fn y(&self) -> f64 {
		self.y
	}

	/// The z component of the vector
	pub fn z(&self) -> f64 {
		self.z
	}

	/// Create a vector of norm 0
	pub fn zero() -> Self {
		Self { x: 0.0, y: 0.0, z: 0.0 }
	}

	/// Compute the squared norm of the vector
	pub fn norm_squared(&self) -> f64 {
		self.x.powi(2) + self.y.powi(2) + self.z.powi(2)
	}

	/// Compute the norm of the vector
	pub fn norm(&self) -> f64 {
		(self.x.powi(2) + self.y.powi(2) + self.z.powi(2)).sqrt()
	}

	pub fn as_point(self) -> Point3 {
		Point3 {
			x: self.x,
			y: self.y,
			z: self.z,
		}
	}
}

macro_rules! add_vec_and_point {
	($L: ty, $R: ty) => {
		impl std::ops::Add<$R> for $L {
			type Output = Vector3;
			fn add(self, rhs: $R) -> Self::Output {
				Vector3::from(self.x() + rhs.x(), self.y() + rhs.y(), self.z() + rhs.z())
			}
		}
	};
}

macro_rules! impl_add_point_vector {
	($L: ident, $R: ident) => {
		add_vec_and_point!($R, $R);
		add_vec_and_point!($R, $L);
		add_vec_and_point!(&$R, $L);
		add_vec_and_point!($R, &$L);
		add_vec_and_point!(&$R, &$L);
	};
}

impl_add_point_vector!(Point3, Vector3);
impl_add_point_vector!(Vector3, Point3);

macro_rules! add_assign_vec_and_point {
	($L: ty, $R: ty) => {
		impl std::ops::AddAssign<$R> for $L {
			fn add_assign(&mut self, rhs: $R) {
				self.x += rhs.x();
				self.y += rhs.y();
				self.z += rhs.z();
			}
		}
	};
}

macro_rules! impl_add_assign_point_vector {
	($L: ident, $R: ident) => {
		add_assign_vec_and_point!($R, $R);
		add_assign_vec_and_point!($R, &$R);
		add_assign_vec_and_point!($R, $L);
		add_assign_vec_and_point!($R, &$L);
	};
}

impl_add_assign_point_vector!(Point3, Vector3);
impl_add_assign_point_vector!(Vector3, Point3);

macro_rules! sub_vec_and_point {
	($L: ty, $R: ty) => {
		impl std::ops::Sub<$R> for $L {
			type Output = Vector3;
			fn sub(self, rhs: $R) -> Self::Output {
				Vector3::from(self.x() - rhs.x(), self.y() - rhs.y(), self.z() - rhs.z())
			}
		}
	};
}

macro_rules! impl_sub_point_vector {
	($L: ident, $R: ident) => {
		sub_vec_and_point!($R, $R);
		sub_vec_and_point!($R, $L);
		sub_vec_and_point!(&$R, $L);
		sub_vec_and_point!($R, &$L);
		sub_vec_and_point!(&$R, &$L);
	};
}

impl_sub_point_vector!(Point3, Vector3);
impl_sub_point_vector!(Vector3, Point3);

macro_rules! sub_assign_vec_and_point {
	($L: ty, $R: ty) => {
		impl std::ops::SubAssign<$R> for $L {
			fn sub_assign(&mut self, rhs: $R) {
				self.x -= rhs.x();
				self.y -= rhs.y();
				self.z -= rhs.z();
			}
		}
	};
}

macro_rules! impl_sub_assign_point_vector {
	($L: ident, $R: ident) => {
		sub_assign_vec_and_point!($R, $R);
		sub_assign_vec_and_point!($R, &$R);
		sub_assign_vec_and_point!($R, $L);
		sub_assign_vec_and_point!($R, &$L);
	};
}

impl_sub_assign_point_vector!(Point3, Vector3);
impl_sub_assign_point_vector!(Vector3, Point3);

impl std::ops::Mul<f64> for Vector3 {
	type Output = Vector3;
	fn mul(self, rhs: f64) -> Self::Output {
		Self::Output {
			x: self.x * rhs,
			y: self.y * rhs,
			z: self.z * rhs,
		}
	}
}

impl std::ops::Mul<Vector3> for f64 {
	type Output = Vector3;
	fn mul(self, rhs: Vector3) -> Self::Output {
		Self::Output {
			x: rhs.x * self,
			y: rhs.y * self,
			z: rhs.z * self,
		}
	}
}

impl std::ops::MulAssign<f64> for Vector3 {
	fn mul_assign(&mut self, rhs: f64) {
		self.x *= rhs;
		self.y *= rhs;
		self.z *= rhs;
	}
}

impl std::ops::Div<f64> for Vector3 {
	type Output = Vector3;
	fn div(self, rhs: f64) -> Self::Output {
		Self::Output {
			x: self.x * rhs,
			y: self.y * rhs,
			z: self.z * rhs,
		}
	}
}

impl std::ops::DivAssign<f64> for Vector3 {
	fn div_assign(&mut self, rhs: f64) {
		self.x *= rhs;
		self.y *= rhs;
		self.z *= rhs;
	}
}

#[cfg(test)]
mod tests {
	use super::*;

	#[test]
	fn point_distance_to_origin() {
		let o = Point3::origin();
		assert_eq!(o.x(), 0.0);
		assert_eq!(o.y(), 0.0);
		assert_eq!(o.z(), 0.0);

		let a = Point3::from(1.0, 2.0, 2.0);
		// difference = (1,1,2) -> squared = 1+4+4 = 9 -> distance = 3
		let dist_sq = a.distance_to_squared(&o);
		let dist = a.distance_to(&o);
		let eps = 1e-12;
		assert!((dist_sq - 9.0).abs() < eps);
		assert!((dist - 3.0).abs() < eps);
	}

	#[test]
	fn point_distance() {
		let p = Point3::from(1.0, 2.0, 3.0);
		let q = Point3::from(4.0, 6.0, 3.0);
		// difference = (3,4,0) -> squared = 9+16+0 = 25 -> distance = 5
		assert!((p.distance_to_squared(&q) - 25.0).abs() < 1e-12);
		assert!((p.distance_to(&q) - 5.0).abs() < 1e-12);
	}

	#[test]
	fn vector_zero_and_norms() {
		let z = Vector3::zero();
		assert_eq!(z.x(), 0.0);
		assert_eq!(z.y(), 0.0);
		assert_eq!(z.z(), 0.0);
		assert!((z.norm_squared() - 0.0).abs() < 1e-12);
		assert!((z.norm() - 0.0).abs() < 1e-12);
	}
}
