"""
Example script demonstrating how to load FEA results from a file.

This script shows how to:
1. Load FEA results from a previously saved file
2. Extract and display key information from the results
3. Visualize the loaded results
"""

import os
import numpy as np
from volco_fea import load_fea_results, visualize_fea, export_visualization


def main():
    """Load and visualize previously saved FEA results."""
    print("Loading FEA results from file...")
    
    # Create output directory if it doesn't exist
    os.makedirs("Results_volco/fea", exist_ok=True)
    
    # Try to load results from different formats
    formats = ['pickle', 'json', 'hdf5']
    loaded_results = None
    
    for fmt in formats:
        file_path = f'Results_volco/fea/results_{fmt}'
        if fmt == 'pickle':
            file_path += '.pkl'
        elif fmt == 'json':
            file_path += '.json'
        elif fmt == 'hdf5':
            file_path += '.h5'
        
        if os.path.exists(file_path):
            print(f"Loading results from {file_path}...")
            try:
                loaded_results = load_fea_results(file_path)
                print(f"Successfully loaded results from {file_path}")
                break
            except Exception as e:
                print(f"Error loading {file_path}: {e}")
    
    if loaded_results is None:
        print("No saved results found. Please run fea_save_results.py first.")
        return
    
    # Display key information from the loaded results
    print("\nFEA Results Summary:")
    print(f"Number of nodes: {loaded_results['nodes'].shape[0]}")
    print(f"Number of elements: {loaded_results['elements'].shape[0]}")
    print(f"Maximum displacement: {loaded_results['max_displacement']:.6f} mm")
    print(f"Maximum von Mises stress: {loaded_results['max_von_mises']:.2f} MPa")
    
    # Create visualization from the loaded results
    print("\nCreating visualization from loaded results...")
    viz = visualize_fea(
        nodes=loaded_results['nodes'],
        elements=loaded_results['elements'],
        displacements=loaded_results['displacements'],
        von_mises=loaded_results['von_mises'],
        result_type='von_mises',
        scale_factor=10.0,
        show_undeformed=False
    )
    
    # Save the visualization
    print("Saving visualization...")
    export_visualization(
        viz,
        "Results_volco/fea/loaded_von_mises.html"
    )
    
    # Create a displacement visualization
    print("Creating displacement visualization...")
    disp_viz = visualize_fea(
        nodes=loaded_results['nodes'],
        elements=loaded_results['elements'],
        displacements=loaded_results['displacements'],
        von_mises=loaded_results['von_mises'],
        result_type='displacement',
        scale_factor=10.0,
        show_undeformed=False
    )
    
    # Create a visualization with both original and deformed meshes
    print("Creating visualization with original and deformed meshes...")
    dual_viz = visualize_fea(
        nodes=loaded_results['nodes'],
        elements=loaded_results['elements'],
        displacements=loaded_results['displacements'],
        von_mises=loaded_results['von_mises'],
        result_type='von_mises',
        scale_factor=10.0,
        show_undeformed=True,
        original_opacity=0.3
    )
    
    # Save the dual visualization
    export_visualization(
        dual_viz,
        "Results_volco/fea/loaded_von_mises_with_undeformed.html"
    )
    
    # Save the displacement visualization
    export_visualization(
        disp_viz,
        "Results_volco/fea/loaded_displacement.html"
    )
    
    # Create a comparison plot showing stress distribution using Plotly
    print("Creating stress distribution histogram...")
    import plotly.graph_objects as go
    
    fig = go.Figure()
    fig.add_trace(go.Histogram(
        x=loaded_results['von_mises'],
        nbinsx=50,
        marker_color='blue',
        opacity=0.7
    ))
    
    fig.update_layout(
        title='Distribution of von Mises Stress',
        xaxis_title='von Mises Stress (MPa)',
        yaxis_title='Frequency',
        width=800,
        height=500
    )
    
    fig.write_html("Results_volco/fea/stress_distribution.html")
    
    print("\nExample complete. Visualizations saved to Results_volco/fea/")
    print("Files generated:")
    print("  - Results_volco/fea/loaded_von_mises.html")
    print("  - Results_volco/fea/loaded_displacement.html")
    print("  - Results_volco/fea/loaded_von_mises_with_undeformed.html")
    print("  - Results_volco/fea/stress_distribution.html")


if __name__ == "__main__":
    main()