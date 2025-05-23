"""
Input/Output module for FEA results.

This module provides functionality to save and load FEA results to/from files.
"""

import os
import json
import pickle
import datetime
import numpy as np
from typing import Dict, Any, Optional, Union, Tuple

# Try to import h5py for HDF5 support
try:
    import h5py
    HDF5_AVAILABLE = True
except ImportError:
    HDF5_AVAILABLE = False


def save_results(
    results: Dict[str, Any],
    filename: str,
    format: str = 'pickle',
    include_visualization: bool = False,
    metadata: Optional[Dict[str, Any]] = None
) -> str:
    """
    Save FEA results to a file.
    
    Args:
        results: Dictionary containing FEA results
        filename: Output filename
        format: File format ('pickle', 'json', or 'hdf5')
        include_visualization: Whether to include visualization in the saved file
        metadata: Additional metadata to include in the saved file
        
    Returns:
        Path to the saved file
    
    Raises:
        ValueError: If an unsupported format is specified
        ImportError: If HDF5 format is requested but h5py is not available
    """
    # Create directory if it doesn't exist
    os.makedirs(os.path.dirname(os.path.abspath(filename)), exist_ok=True)
    
    # Prepare data to save
    data = {
        'results': {},
        'metadata': {
            'timestamp': datetime.datetime.now().isoformat(),
            'version': '1.0.0'
        }
    }
    
    # Add user metadata if provided
    if metadata is not None:
        data['metadata'].update(metadata)
    
    # Copy results, excluding visualization by default
    for key, value in results.items():
        if key != 'visualization' or include_visualization:
            data['results'][key] = value
    
    # Save in the specified format
    if format.lower() == 'pickle':
        with open(filename, 'wb') as f:
            pickle.dump(data, f)
    elif format.lower() == 'json':
        # Convert numpy arrays to lists for JSON serialization
        json_data = _prepare_for_json(data)
        with open(filename, 'w') as f:
            json.dump(json_data, f, indent=2)
    elif format.lower() == 'hdf5':
        if not HDF5_AVAILABLE:
            raise ImportError("HDF5 format requires h5py. Install it with 'pip install h5py'.")
        
        with h5py.File(filename, 'w') as f:
            # Save metadata as attributes
            for key, value in data['metadata'].items():
                f.attrs[key] = value
            
            # Create results group
            results_group = f.create_group('results')
            _save_to_hdf5(results_group, data['results'])
    else:
        raise ValueError(f"Unsupported format: {format}. Use 'pickle', 'json', or 'hdf5'.")
    
    return filename


def load_results(
    filename: str
) -> Dict[str, Any]:
    """
    Load FEA results from a file.
    
    Args:
        filename: Path to the file containing saved FEA results
        
    Returns:
        Dictionary containing the loaded FEA results
        
    Raises:
        ValueError: If the file format is not recognized
        ImportError: If HDF5 format is detected but h5py is not available
    """
    # Determine file format based on extension
    _, ext = os.path.splitext(filename)
    
    if ext.lower() == '.pkl' or ext.lower() == '.pickle':
        with open(filename, 'rb') as f:
            data = pickle.load(f)
    elif ext.lower() == '.json':
        with open(filename, 'r') as f:
            data = json.load(f)
        # Convert lists back to numpy arrays
        data = _convert_from_json(data)
    elif ext.lower() == '.h5' or ext.lower() == '.hdf5':
        if not HDF5_AVAILABLE:
            raise ImportError("HDF5 format requires h5py. Install it with 'pip install h5py'.")
        
        with h5py.File(filename, 'r') as f:
            # Load metadata
            metadata = dict(f.attrs)
            
            # Load results
            results = {}
            _load_from_hdf5(f['results'], results)
            
            data = {
                'metadata': metadata,
                'results': results
            }
    else:
        raise ValueError(f"Unrecognized file format: {ext}. Supported formats: .pkl, .pickle, .json, .h5, .hdf5")
    
    # Return the results with metadata
    return data


def _prepare_for_json(data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Prepare data for JSON serialization by converting numpy arrays to lists.
    
    Args:
        data: Dictionary containing data with numpy arrays
        
    Returns:
        Dictionary with numpy arrays converted to lists
    """
    result = {}
    
    for key, value in data.items():
        if isinstance(value, dict):
            result[key] = _prepare_for_json(value)
        elif isinstance(value, np.ndarray):
            result[key] = {
                '__type__': 'ndarray',
                'dtype': str(value.dtype),
                'shape': value.shape,
                'data': value.tolist()
            }
        elif isinstance(value, (np.int32, np.int64)):
            result[key] = int(value)
        elif isinstance(value, (np.float32, np.float64)):
            result[key] = float(value)
        else:
            result[key] = value
    
    return result


def _convert_from_json(data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Convert JSON data back to original format, recreating numpy arrays.
    
    Args:
        data: Dictionary containing data loaded from JSON
        
    Returns:
        Dictionary with lists converted back to numpy arrays
    """
    result = {}
    
    for key, value in data.items():
        if isinstance(value, dict):
            if '__type__' in value and value['__type__'] == 'ndarray':
                # Reconstruct numpy array
                result[key] = np.array(value['data'], dtype=np.dtype(value['dtype']))
            else:
                result[key] = _convert_from_json(value)
        else:
            result[key] = value
    
    return result


def _save_to_hdf5(group, data: Dict[str, Any]) -> None:
    """
    Recursively save data to an HDF5 group.
    
    Args:
        group: HDF5 group to save to
        data: Dictionary containing data to save
    """
    for key, value in data.items():
        if isinstance(value, dict):
            subgroup = group.create_group(key)
            _save_to_hdf5(subgroup, value)
        elif isinstance(value, np.ndarray):
            group.create_dataset(key, data=value)
        else:
            # Store as attribute if it's a simple type
            try:
                group.attrs[key] = value
            except TypeError:
                # If it can't be stored as an attribute, pickle it
                group.create_dataset(key, data=np.void(pickle.dumps(value)))


def _load_from_hdf5(group, result: Dict[str, Any]) -> None:
    """
    Recursively load data from an HDF5 group.
    
    Args:
        group: HDF5 group to load from
        result: Dictionary to store loaded data
    """
    # Load datasets
    for key in group.keys():
        item = group[key]
        if isinstance(item, h5py.Group):
            result[key] = {}
            _load_from_hdf5(item, result[key])
        elif isinstance(item, h5py.Dataset):
            if item.dtype == np.dtype('object'):
                # This might be a pickled object
                result[key] = pickle.loads(item[()])
            else:
                result[key] = item[()]
    
    # Load attributes
    for key, value in group.attrs.items():
        result[key] = value