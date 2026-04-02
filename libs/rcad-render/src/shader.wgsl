struct CameraUniform {
    view_proj: mat4x4<f32>,
    eye_pos: vec4<f32>,
};

struct MaterialUniform {
    color: vec4<f32>,
    // flags.x > 0.5 means unlit flat color output
    flags: vec4<f32>,
};

@group(0) @binding(0)
var<uniform> camera: CameraUniform;

@group(1) @binding(0)
var<uniform> material: MaterialUniform;

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
fn fs_main(in: VertexOutput, @builtin(front_facing) front_facing: bool) -> @location(0) vec4<f32> {
    if (material.flags.x > 0.5) {
        return vec4<f32>(material.color.xyz, material.color.w);
    }

    var normal = normalize(cross(dpdx(in.world_position), dpdy(in.world_position)));
    if (!front_facing) {
        normal = -normal;
    }

    let light_dir = normalize(vec3<f32>(0.45, 0.85, 0.35));
    let view_dir = normalize(camera.eye_pos.xyz - in.world_position);
    let half_dir = normalize(light_dir + view_dir);

    let base_color = material.color.xyz;
    let ambient = 0.18;
    let diffuse = max(dot(normal, light_dir), 0.0);
    let specular = pow(max(dot(normal, half_dir), 0.0), 24.0) * 0.20;

    let lit = base_color * (ambient + diffuse * 0.82) + vec3<f32>(specular);
    return vec4<f32>(lit, 1.0);
}
