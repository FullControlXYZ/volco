import numpy as np
import trimesh
import plotly.graph_objects as go
from skimage import measure


def color_mesh(mesh, color_scheme='cyan_blue'):
    """
    Apply colors to the mesh based on height (z-value).
    Returns the colored mesh.
    
    Parameters:
    -----------
    mesh : trimesh.Trimesh or trimesh.Scene
        The mesh to color
    color_scheme : str
        The color scheme to use ('cyan_blue' or 'viridis')
        
    Returns:
    --------
    trimesh.Trimesh or trimesh.Scene
        The colored mesh
    """
    # Check if mesh is a Scene object (from box representation)
    if isinstance(mesh, trimesh.Scene):
        # For box representation (Scene object), we can't easily color individual vertices
        # Return the original scene without coloring
        return mesh
        
    # For Trimesh objects (from marching cubes)
    if not hasattr(mesh, 'vertices') or len(mesh.vertices) == 0:
        return mesh

    z_values = mesh.vertices[:, 2]
    # Normalize z values to 0-1 range for coloring
    z_normalized = (z_values - z_values.min()) / (z_values.max() - z_values.min() if z_values.max() > z_values.min() else 1)
    
    # Create a color array
    colors = np.zeros((len(mesh.vertices), 4))
    
    if color_scheme == 'cyan_blue':
        colors[:, 0] = 0
        colors[:, 1] = z_normalized                   # Green decreases but stays bright
        colors[:, 2] = 1.0                            # Blue stays at maximum brightness
        colors[:, 3] = 1.0                            # Alpha channel (fully opaque)
    elif color_scheme == 'viridis':
        # Approximate viridis colormap
        colors[:, 0] = 0.267 + 0.733 * (1 - z_normalized)  # Red
        colors[:, 1] = 0.004 + 0.996 * z_normalized        # Green
        colors[:, 2] = 0.329 + 0.413 * z_normalized        # Blue
        colors[:, 3] = 1.0                                 # Alpha
    
    # Apply colors to the mesh
    mesh.visual.vertex_colors = (colors * 255).astype(np.uint8)
    
    return mesh


def visualize_with_trimesh(mesh):
    """
    Create a 3D visualization of the mesh using trimesh.
    
    Parameters:
    -----------
    mesh : trimesh.Trimesh or trimesh.Scene
        The mesh to visualize
        
    Returns:
    --------
    trimesh.Scene
        A trimesh scene that can be displayed
    """
    # If mesh is already a Scene (from box representation), use it directly
    if isinstance(mesh, trimesh.Scene):
        scene = mesh
    else:
        # For Trimesh objects (from marching cubes), create a scene
        scene = trimesh.Scene(mesh)
    
    # Calculate a good distance for camera position
    try:
        # For Trimesh objects
        if not isinstance(mesh, trimesh.Scene):
            distance = mesh.bounding_box.volume ** (1/3) * 4
        else:
            # For Scene objects, try to get a reasonable distance
            # Get the extents of all geometries in the scene
            bounds = scene.bounds
            if bounds is not None:
                extents = bounds[1] - bounds[0]
                distance = max(extents) * 2
            else:
                distance = 10  # Default distance if bounds can't be determined
    except Exception:
        # Fallback distance if calculation fails
        distance = 10
    
    # Set camera to a position similar to the Plotly view
    # Looking from approximately -x, -y direction with positive z up
    scene.set_camera(angles=(np.pi/4, 0, 0), distance=distance)
    
    # For trimesh 3.21.7, we can only adjust basic rendering properties
    # Increase the ambient light to make the mesh more visible
    try:
        # Attempt to make the mesh brighter by increasing ambient light
        scene.graph[scene.camera.name][0] = 4.0  # Increase ambient light
    except (AttributeError, KeyError, IndexError, TypeError):
        # If that doesn't work, we'll continue without modifying the lighting
        pass
    
    return scene


def visualize_with_plotly(mesh):
    """
    Create a 3D visualization of the mesh using Plotly.
    
    Parameters:
    -----------
    mesh : trimesh.Trimesh or trimesh.Scene
        The mesh to visualize
        
    Returns:
    --------
    plotly.graph_objects.Figure
        A plotly figure that can be displayed
    """
    # Check if mesh is a Scene object (from box representation)
    if isinstance(mesh, trimesh.Scene):
        # For Scene objects, we need to extract all meshes and combine them
        # This is a simplified approach that may not work for all Scene objects
        # Create a warning message
        fig = go.Figure()
        fig.add_annotation(
            text="Plotly visualization is not supported for box representation.<br>Use trimesh visualizer instead.",
            xref="paper", yref="paper",
            x=0.5, y=0.5,
            showarrow=False,
            font=dict(size=20, color="red")
        )
        
        # Set basic layout
        fig.update_layout(
            width=800,
            height=800
        )
        
        return fig
    
    # For Trimesh objects (from marching cubes)
    if mesh is not None:
        verts = mesh.vertices
        faces = mesh.faces
        
        # Extract coordinates
        x_mesh = verts[:, 0]
        y_mesh = verts[:, 1]
        z_mesh = verts[:, 2]
        
        # Use z values for coloring
        intensity = z_mesh
        
        # Create mesh3d visualization
        fig = go.Figure([go.Mesh3d(
            x=verts[:, 0], y=verts[:, 1], z=verts[:, 2],
            i=faces[:, 0], j=faces[:, 1], k=faces[:, 2],
            opacity=1,
            vertexcolor=(mesh.visual.vertex_colors[:, :3]
                         if hasattr(mesh.visual, 'vertex_colors') else None),
            colorscale=('viridis'
                        if not hasattr(mesh.visual, 'vertex_colors') else None),
            intensity=(verts[:, 2]
                       if not hasattr(mesh.visual, 'vertex_colors') else None),
            showscale=False,
            flatshading=True,
            lighting=dict(ambient=1, diffuse=0.8, specular=0.2,
                          roughness=0, fresnel=0.1),
            lightposition=dict(x=100, y=200, z=0)
        )])

        common_axes_dict = dict(
            showbackground=True,
            backgroundcolor='black',
            gridcolor='gray',
            color='white'
        )
        # Improve layout
        fig.update_layout(
            paper_bgcolor='black',
            plot_bgcolor='black',
            scene=dict(                
                xaxis_title='X (mm)',
                yaxis_title='Y (mm)',
                zaxis_title='Z (mm)',
                xaxis=common_axes_dict,
                yaxis=common_axes_dict,
                zaxis=common_axes_dict,
                bgcolor='black',
                aspectmode='data',  # Maintain aspect ratio based on data
                camera=dict(
                    eye=dict(x=-1.8, y=-1.8, z=1.0),  # Position camera further away for a more zoomed out view
                    up=dict(x=0, y=0, z=1)      # Keep z-axis pointing up
                )
            ),
            width=800,
            height=400,
            margin=dict(l=10, r=10, t=10, b=10),
        )
        
        return fig
    
    # If mesh is None, return an empty figure
    return go.Figure()