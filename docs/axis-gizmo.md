# Role
You are an expert Rust Graphics Engineer specializing in low-level GPU programming using the `wgpu` crate.

# Objective
Create a self-contained Rust application using `winit` and `wgpu` that renders a "3D Viewport Axis Gizmo" similar to engineering CAD software (like ANSYS or Blender).

# Visual Description (High Fidelity)
Please implement the following visual elements exactly as described:

1. **Coordinate System**:
   - Right-handed 3D coordinate system.
   - **X-Axis (Red)**: Pointing Right.
   - **Y-Axis (Green)**: Pointing Up.
   - **Z-Axis (Blue)**: Pointing Out of the screen (towards camera).

2. **Geometry**:
   - **Axes Shafts**: Three cylinders originating from (0,0,0). They should be thick enough to be clearly visible, not thin lines.
   - **Axis Tips**: Cones placed at the end of each cylinder, matching the cylinder's radius.
   - **Origin Node**: A small, semi-transparent or opaque gray metallic cube (or sphere) located at (0,0,0) where the axes intersect. It should slightly obscure the back of the axes to give a 3D feel.
   - **Background**: A clean white background with a faint gray perspective grid on the XZ plane (Ground Plane).

3. **Material & Lighting**:
   - **Shading Model**: Use a standard Phong shading model or simple PBR in the WGSL shader.
   - **Lighting**: A directional light source coming from the top-right-front to create specular highlights on the glossy cones (visible in the reference image).
   - **Colors**:
     - X: Bright Red (`#FF0000` or similar)
     - Y: Bright Green (`#00FF00` or similar)
     - Z: Bright Blue (`#0000FF` or similar)
     - Center: Metallic Gray (`#808080`)

4. **Camera Control**:
   - Implement a basic "Orbit Camera" (Arcball) using mouse input.
   - Left-click drag to rotate the view around the gizmo.
   - Scroll to zoom in/out.
   - The Gizmo should always remain centered in the screen (0,0 in NDC), or the camera orbits around the world origin.

