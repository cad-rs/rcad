struct CameraUniform {
    view_proj: mat4x4<f32>,
};

@group(0) @binding(0)
var<uniform> camera: CameraUniform;

struct VertexInput {
    @location(0) position: vec3<f32>,
};

struct VertexOutput {
    @builtin(position) clip_position: vec4<f32>,
    @location(0) world_position: vec3<f32>,
};

@vertex
fn vs_main(in: VertexInput) -> VertexOutput {
    var out: VertexOutput;
    out.world_position = in.position;
    out.clip_position = camera.view_proj * vec4<f32>(in.position, 1.0);
    return out;
}

@fragment
fn fs_main(in: VertexOutput) -> @location(0) vec4<f32> {
    // Use a more visible diagnostic color
    let normal = normalize(cross(dpdx(in.world_position), dpdy(in.world_position)));
    let light_dir = normalize(vec3<f32>(1.0, 1.0, 1.0));
    let diff = max(dot(normal, light_dir), 0.2);
    let base_color = vec3<f32>(0.0, 0.5, 1.0); // Bright blue
    return vec4<f32>(base_color * diff, 1.0);
}
