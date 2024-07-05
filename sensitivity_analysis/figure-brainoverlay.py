# script to derive figure with brain overlay of the four different estimators
# average, pool.se, pool.gls, pool.gls1

# import relevant packages
import pandas as pd
import numpy as np
import nibabel as nib
from nibabel import freesurfer
from glob import glob

# extract current pwd
pwd = '/staff/mganz/Ganz/gitrepos/multiverse_tools/sensitivity_analysis/'
# get the fsaverage subject of freesurfer on our cluster
fsaveragePath = '/usr/local/nru/freesurfer/fs71/subjects/fsaverage/'

def extract_estimator_results(resultsFile,estimator,rename_region):
    """
    Extract the values from the resultsFile for a given estimator

    Parameters:
    - resultsFile: File with estimator results from Brice's statistical modelling
    - estimator: string with one of the four estimators either 'average', 'pool.se', 'pool.gls', 'pool.gls1'
    - rename_region: dataframe wih renaming that brice originally created in R

    Returns:
    - region_value_pairs: 
    """
    # Load data into pandas dataframe
    results = pd.read_csv(resultsFile)
    # Extract only the relevant rows for the four estimators
    rowLabelsEst = np.where(results['type']==estimator)
    # Extract the estimate for those rows from the results
    estimatorVal = results['estimate'][rowLabelsEst[0]]
    # Extract the unique regions from the data
    regions = results['region'].unique()
    # Pair the region with each value in the results file
    region_value_pairs = {name: value for name, value in zip(regions, estimatorVal)}
    # Create the mapping dictionary from combined to full names
    combined_to_full_mapping = dict(zip(rename_region['combined'], rename_region['full']))
    # Create a new dictionary with renamed keys
    renamed_region_value_pairs = {combined_to_full_mapping.get(key, key): value for key, value in region_value_pairs.items()}

    return renamed_region_value_pairs 


def find_region_index(names, region_name):
    """
    Find the index of the specified region name or substring in the names list.

    Parameters:
    - names (list): List of region names.
    - region_name (str): The region name or substring to find in the list.

    Returns:
    - index (int): The index of the region name in the list, or None if not found.
    """
    if not isinstance(region_name, str):
        raise TypeError("The substring should be a string.")
    
    decoded_names = []
    for i, name in enumerate(names):
        if isinstance(name, bytes):
            decoded_names.append(name.decode('utf-8'))
        elif isinstance(name, str):
            decoded_names.append(name)
        else:
            raise TypeError(f"All elements in names_list should be strings or bytes. Found {type(name)} at index {i}.")

    
    matching_indices = [i for i, name in enumerate(decoded_names) if region_name in name]
    if matching_indices:
        return matching_indices[0]  # Return the first match
    return None


def create_freesurfer_overlay(annot_file, surf_file, region_value_pairs, output_file):
    """
    Create a FreeSurfer overlay using an annotation file and a dictionary of region-value pairs.

    Parameters:
    - annot_file (str): Path to the FreeSurfer annotation file (.annot).
    - region_value_pairs (dict): Dictionary of region names and their corresponding estimator values.
    - output_file (str): Path to save the output NIfTI file.

    Returns:
    - None
    """
    # Load the annotation file
    labels, ctab, names = nib.freesurfer.read_annot(annot_file)

    # Load the surface file to create a correct overlay
    surface = nib.freesurfer.read_geometry(surf_file)
    vertices, faces = surface
    num_vertices = vertices.shape[0]

    # Initialize overlay array with zeros (or any background value)
    overlay = np.zeros(num_vertices,dtype=np.float32)

    # Find the correct region indices and apply the values
    for region_name, value in region_value_pairs.items():
        index = find_region_index(names, region_name)
        if index is not None:
            overlay[labels == index] = value
        else:
            print(f"Warning: The region '{region_name}' was not found in the annotation file.")

    
    # For testing - load existing overlay file
    # overlay = nib.freesurfer.read_morph_data('/usr/local/nru/freesurfer/fs71/subjects/fsaverage/surf/lh.curv')
    
    # Write an overlay file 
    nib.freesurfer.write_morph_data(output_file, overlay)


# Actual script starting

# Utilize the same name matching as Brice has done in application.R line 42 to get the correct labels for the overlay
data = {
    'combined': ["amygdala", "thalamus", "putamen", "caudate", "ACC", "hippocampus", "OFC", "SFC", "OC", "STG", "insula", "ITG", "PC", "EC"],
        'full': ["amygdala", "thalamus", "putamen", "caudate", "anteriorcingulate", "hippocampus", "orbitofrontal", "superiorfrontal", "occipital", "superiortemporal", "insula", "inferiortemporal", "parietal", "entorhinal"]
    #'position':["subcortical","subcortical","subcortical","subcortical","cortical","subcortical","cortical","cortical","cortical","cortical","cortical","cortical","cortical","cortical"]
}
# Create the DataFrame
rename_region = pd.DataFrame(data)

# Extract pipeline results first
resultsFile = pwd+'results/data-gg-forest.csv'

estimators = ['average','pool.se','pool.gls','pool.gls1']
for estimator in estimators:
    region_value_pairs = extract_estimator_results(resultsFile,estimator,rename_region)

    # Create the overlay on the surfaces
    annot_file = fsaveragePath + '/label/lh.aparc.annot'  # Path to your annotation file
    surf_file = fsaveragePath + '/surf/lh.pial'  # Path to pial surface file
    output_file = pwd+'figures/lh.overlay_'+estimator  # Path to save the output overlay file

    create_freesurfer_overlay(annot_file, surf_file, region_value_pairs, output_file)


 





