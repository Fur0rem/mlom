use bevy::input::mouse::{MouseMotion, MouseWheel};
use bevy::prelude::*;
use mlom::system::System;
use std::path::Path;

#[derive(Component)]
struct OrbitCamera {
	radius:  f32,
	angle_x: f32,
	angle_y: f32,
}

#[derive(Component)]
struct Billboard;

fn main() {
	App::new()
		.add_plugins(DefaultPlugins)
		.add_systems(Startup, setup)
		.add_systems(Update, orbit_camera)
		.add_systems(PostUpdate, update_billboards)
		.run();
}

fn setup(
	mut commands: Commands, mut meshes: ResMut<Assets<Mesh>>, mut materials: ResMut<Assets<StandardMaterial>>,
	asset_server: Res<AssetServer>,
) {
	let system = System::from_file(Path::new("../dataset/particles.xyz"), 0);

	let texture_handle = asset_server.load("particle.png");

	// Create shared mesh and material for particles
	let mesh_handle = meshes.add(Rectangle::new(0.8, 0.8));
	let material_handle = materials.add(StandardMaterial {
		base_color: Color::WHITE,
		base_color_texture: Some(texture_handle),
		alpha_mode: AlphaMode::Mask(0.5),
		unlit: true,
		..default()
	});

	// Spawn a circle for each particle
	for particle in system.particles() {
		let (x, y, z) = particle.xyz();
		commands.spawn((
			Mesh3d(mesh_handle.clone()),
			MeshMaterial3d(material_handle.clone()),
			Transform::from_xyz(x as f32, y as f32, z as f32),
			Billboard,
		));
	}

	// Add a camera
	let radius = 40.0;
	commands.spawn((
		Camera3d::default(),
		OrbitCamera {
			radius,
			angle_x: 0.7,
			angle_y: 0.5,
		},
		Transform::from_xyz(radius * 0.5, radius * 0.5, radius * 0.5).looking_at(Vec3::ZERO, Vec3::Y),
	));
}

fn update_billboards(
	camera_query: Query<&Transform, With<Camera>>, mut billboard_query: Query<&mut Transform, (With<Billboard>, Without<Camera>)>,
) {
	let Ok(camera_transform) = camera_query.single()
	else {
		return;
	};

	for mut transform in billboard_query.iter_mut() {
		transform.rotation = camera_transform.rotation;
	}
}

fn orbit_camera(
	mut mouse_motion: MessageReader<MouseMotion>, mut scroll_event: MessageReader<MouseWheel>,
	input_mouse: Res<ButtonInput<MouseButton>>, mut query: Query<(&mut Transform, &mut OrbitCamera)>,
) {
	let mut rotation_move = Vec2::ZERO;
	let mut scroll = 0.0;

	if input_mouse.pressed(MouseButton::Left) {
		for ev in mouse_motion.read() {
			rotation_move += ev.delta;
		}
	}

	for ev in scroll_event.read() {
		scroll += ev.y;
	}

	for (mut transform, mut orbit) in query.iter_mut() {
		if rotation_move.length_squared() > 0.0 {
			orbit.angle_x -= rotation_move.x * 0.01;
			orbit.angle_y -= rotation_move.y * 0.01;
			// Clamp vertical angle to avoid flipping
			orbit.angle_y = orbit.angle_y.clamp(-1.5, 1.5);
		}

		if scroll.abs() > 0.0 {
			orbit.radius -= scroll * 2.0;
			orbit.radius = orbit.radius.max(2.0);
		}

		// Update transform based on angles and radius
		let rot = Quat::from_euler(EulerRot::YXZ, orbit.angle_x, orbit.angle_y, 0.0);
		transform.translation = rot * Vec3::Z * orbit.radius;
		transform.look_at(Vec3::ZERO, Vec3::Y);
	}
}
