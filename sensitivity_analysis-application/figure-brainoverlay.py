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
    - renamed_region_value_pairs: dictionary with full names as in FreeSurfer and the region value 
    """
    # Load data into pandas dataframe
    results = pd.read_csv(resultsFile, sep=',')
    # Extract only the relevant rows for the four estimators
    rowLabelsEst = np.where(results['type']==estimator)
    # Extract the estimate for those rows from the results
    estimatorVal = results['estimate'][rowLabelsEst[0]]
    # Extract the unique regions from the data
    regions = results['region'].unique()
    # Pair the region with each value in the results file
    region_value_pairs = {name: value for name, value in zip(regions, estimatorVal)}
    # Create the mapping dictionary from resultsName to overlayName
    resultsName_to_overlayName_mapping = dict(zip(rename_region['resultsName'], rename_region['overlayName']))
    # Create a new dictionary with renamed keys
    renamed_region_value_pairs = {resultsName_to_overlayName_mapping.get(key, key): value for key, value in region_value_pairs.items()}

    return renamed_region_value_pairs 

def extract_atlas_results(atlasFile,tracer,rename_region):
    """
    Extract the values from the resultsFile for a given estimator

    Parameters:
    - atlasFile: File with atlas results per region
    - tracer: string with traceer name used in the atlas
    - rename_region: dataframe wih renaming that brice originally created in R

    Returns:
    - renamed_region_value_pairs: dictionary with full names as in FreeSurfer and the region value 
    """
    # Load data into pandas dataframe
    atlas = pd.read_csv(atlasFile, sep=',')
    # reduce it to only the DASB column
    # Extract the estimate for those rows from the results
    atlasVal = atlas[tracer]
    # Extract the unique regions from the data
    regions = atlas['Region'].unique()
    # Pair the region with each value in the results file
    region_atlas_pairs = {name: value for name, value in zip(regions, atlasVal)}
    # Create the mapping dictionary from atlasName to overlayName
    atlasName_to_overlayName_mapping = dict(zip(rename_region['atlasName'], rename_region['overlayName']))
    # Create a new dictionary with renamed keys
    renamed_region_atlas_pairs = {atlasName_to_overlayName_mapping.get(key, key): value for key, value in region_atlas_pairs.items()}

    return renamed_region_atlas_pairs 


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


def create_freesurfer_overlay(annot_file, surf_file, region_value_pairs, atlas_value_pairs, output_file):
    """
    Create a FreeSurfer overlay using an annotation file and a dictionary of region-value pairs.

    Parameters:
    - annot_file (str): Path to the FreeSurfer annotation file (.annot).
    - region_value_pairs (dict): Dictionary of region names and their corresponding estimator values.
    - atlas_value_pairs (dict): Dictionary of region names and their corresponding atlas values for a specific tracer.
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
        atlasValue = eval(atlas_value_pairs[region_name])
        if index is not None:
            overlay[labels == index] = value - atlasValue
        else:
            print(f"Warning: The region '{region_name}' was not found in the annotation file.")

    
    # For testing - load existing overlay file
    # overlay = nib.freesurfer.read_morph_data('/usr/local/nru/freesurfer/fs71/subjects/fsaverage/surf/lh.curv')
    
    # Write an overlay file 
    nib.freesurfer.write_morph_data(output_file, overlay)


# Actual script starting

# Utilize the same name matching as Brice has done in application.R line 42 to get the correct labels for the overlay
data = {
    'resultsName': ["amygdala", "thalamus", "putamen", "caudate", "ACC", "hippocampus", "OFC", "SFC", "OC", "STG", "insula", "ITG", "PC", "EC"],
    'overlayName': ["amygdala", "thalamus", "putamen", "caudate", "rostralanteriorcingulate", "hippocampus", "medialorbitofrontal", "superiorfrontal", "lateraloccipital", "superiortemporal", "insula", "inferiortemporal", "superiorparietal", "entorhinal"]
    #'position':["subcortical","subcortical","subcortical","subcortical","cortical","subcortical","cortical","cortical","cortical","cortical","cortical","cortical","cortical","cortical"]
}
# Create the DataFrame
rename_region = pd.DataFrame(data)

# Utilize the same name matching as Brice has done in application.R line 120 to get the correct labels for the overlay
data2 = {
    'atlasName': ["Amygdala", "Thalamus", "Putamen", "Caudate", "Rostral Anterior", "Hippocampus", "Medial Orbito Frontal", "Superior Frontal", "Lateral Occipital", "Superior Temporal", "Insula", "Inferior Temporal", "Superior Parietal", "Entorhinal"],
    'overlayName': ["amygdala", "thalamus", "putamen", "caudate", "rostralanteriorcingulate", "hippocampus", "medialorbitofrontal", "superiorfrontal", "lateraloccipital", "superiortemporal", "insula", "inferiortemporal", "superiorparietal", "entorhinal"]
    #'position':["subcortical","subcortical","subcortical","subcortical","cortical","subcortical","cortical","cortical","cortical","cortical","cortical","cortical","cortical","cortical"]
}
# Create the DataFrame
rename_region2 = pd.DataFrame(data2)

# Extract pipeline results first
resultsFile = pwd+'results/data-gg-forest.csv'

# Extract atlas results first
atlasFile = pwd+'data/bp_table.csv'

estimators = ['average','pool.se','pool.gls','pool.gls1']
for estimator in estimators:
    region_value_pairs = extract_estimator_results(resultsFile,estimator,rename_region)
    atlas_value_pairs = extract_atlas_results(atlasFile,'[11C]DASB',rename_region2)

    # Create the overlay on the surfaces
    annot_file = fsaveragePath + '/label/lh.aparc.annot'  # Path to your annotation file
    surf_file = fsaveragePath + '/surf/lh.pial'  # Path to pial surface file
    output_file = pwd+'figures/lh.overlay_atlasDiff_'+estimator  # Path to save the output overlay file

    create_freesurfer_overlay(annot_file, surf_file, region_value_pairs,atlas_value_pairs, output_file)


 





