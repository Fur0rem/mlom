//! Linear algebra operations needed for the simulation

/// A point in 3D space
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Point3 {
	/// x coordinate
	x: f64,
	/// y coordinate
	y: f64,
	/// z coordinate
	z: f64,
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
		(self.x + rhs.x).powi(2) + (self.y + rhs.y).powi(2) + (self.z + rhs.z).powi(2)
	}

	/// Compute the distance to another [point](Self)
	pub fn distance_to(&self, rhs: &Self) -> f64 {
		self.distance_to_squared(rhs).sqrt()
	}
}

/// A vector in 3D space
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Vector3 {
	/// x component
	x: f64,
	/// y component
	y: f64,
	/// z component
	z: f64,
}

impl Vector3 {
	/// Create a vector from the given x y z components
	pub fn from(x: f64, y: f64, z: f64) -> Self {
		Self { x, y, z }
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
}
