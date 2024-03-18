import prody
from scipy.spatial import ConvexHull
import reader
from ripser import ripser
from persim import plot_diagrams
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import glob
import plotly.express as px
def get_ca_atoms(prody_structure):
    """
    Get the alpha carbon atoms from a pdb file.
    """
    # Select the alpha carbon atoms
    ca_atoms = prody_structure.select('calpha')
    return ca_atoms



def get_ca_coordinates(ca_atoms):
    """
    Get the coordinates of the alpha carbon atoms.
    """
    # Get the coordinates of the alpha carbon atoms
    ca_coordinates = ca_atoms.getCoords()
    return ca_coordinates

def get_k_grouping(all_coordinates, k):
    """
    Get the k grouping of the coordinates.
    """
    return
def read_pdbs_from_file(file_path):
    pdbs = []
    pdb_sidechain = {}
    pdb_conforms = {}
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if parts[9] == "NA" or parts[9] == "10_PDB_validation":
                continue
            pdb = parts[9]
            pdbs.append(pdb)
            pdb_sidechain[pdb[0:4]] = pdb[4]
            pdb_id = pdb[0:4]
            pdb_conforms[pdb_id] = parts[10]
    return pdbs, pdb_sidechain, pdb_conforms




if __name__ == "__main__":
    # Read the pdb file
    files = glob.glob("*.pdb")
        
    structures = []
    list_of_files = files

    with open('vector_labels.csv','w') as f:
        for file in files:
            f.write(file + '\n')
    for file in files:
        structures.append(reader.read_pdb(file))
    atoms = [get_ca_atoms(structure) for structure in structures]
    coords = [get_ca_coordinates(atom) for atom in atoms] 
    vectors = [ripser(coord)['dgms'] for coord in coords]
    
    data = []
    # Create vector list
    for vector in vectors:
        newfeature = []
        for feature in vector[1]:
            newfeature.append(feature[1]-feature[0])
        data.append(newfeature)

    # find Max shape
    min_val = len(data[0])
    for vec in data:
        if len(vec) < min_val:
            min_val = len(vec)

    # Fix vectors to be sorted and same lenght
    for i in range(len(data)):
        data[i] = sorted(data[i], reverse=True)[:5]
         

    data = np.array(data)
    np.save("data.csv", data)

    # Assuming you have your data stored in a numpy array 'data'
    # Each row represents a sample and each column represents a feature
    # Make sure your data is properly preprocessed (e.g., mean-centered, scaled)
    
    pdb_conforms = read_pdbs_from_file('paper_data.txt')[2] 
    ''' 
    print(pdb_conforms)
    for label in labels:
        print(label)
        print(pdb_conforms[label])
        print(state[pdb_conforms[label]])
    colors = [state[pdb_conforms[label]] for label in labels]
    '''
    # Initialize PCA with the number of components you want to keep
    num_components = 4  # You can change this value based on your requirement
    pca = PCA(n_components=num_components)
    # Fit PCA to your data
    pca.fit(data)
    
    # Transform the data to its principal components
    transformed_data = pca.transform(data)
    
    # Access the principal components (eigenvectors)
    principal_components = pca.components_
    
    # Access the explained variance ratio of each principal component
    explained_variance_ratio = pca.explained_variance_ratio_
    
    # Print the explained variance ratio
    print("Explained Variance Ratio:", explained_variance_ratio)
    
    # Print the principal components
    print("Principal Components:")
    for i, component in enumerate(principal_components):
        print(f"Principal Component {i + 1}: {component}")
    
    # Print the transformed data
    print("Transformed Data:")
    print(transformed_data)

    components = pca.fit_transform(data)
    labels = {
        str(i): f"PC {i+1} ({var:.1f}%)"
        for i, var in enumerate(pca.explained_variance_ratio_ * 100)
    }
    state = {"DFGin": "red", "DFGout": "blue", "DFGinter":"green", "NA":"black"}    
    





    fig = px.scatter_matrix(
        components,
        labels=labels,
        dimensions=range(4)

    )
    fig.update_traces(diagonal_visible=False)
    fig.show()
