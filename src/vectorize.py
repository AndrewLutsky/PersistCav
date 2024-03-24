from sklearn.preprocessing import StandardScaler
import reader
from ripser import ripser
from persim import plot_diagrams
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import glob
import plotly.express as px
import prody 
import numpy as np
from persim import PersImage
from persim import PersistenceImager
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
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

def vectorize(file):
    """
    Vectorizes a particular structure
    """
    struct = reader.read_pdb(file)
    atoms = get_ca_atoms(struct)
    coords = get_ca_coordinates(atoms)
    vectors = ripser(coords)['dgms']
    return vectors


def generate_descriptor(file, max_vector_size):
    """
    Generates a descriptor for a particular structure.
    """
    vector = vectorize(file)
    descriptor = vector[1]
    second_desc = vector[0]
    descriptor = descriptor[0:max_vector_size]
    second_desc = second_desc[0:max_vector_size]

    new_desc = descriptor + second_desc
    #for i in range(len(new_desc)):
    #    new_desc[i][0] = w

    descriptor = [coord[1] for coord in new_desc]

    return sorted(descriptor)

def get_persistance_image(file):
    """
    Get the persistance image for a particular structure.
    """
    vector = vectorize(file)
    vector[0] = vector[0][0:-1]
    vector = vector[1]
    pimgr = PersistenceImager(pixel_size=0.5)
    pimgr.fit(np.asarray(vector))
    img = pimgr.transform(vector, skew=True)
    print(img)
    pimgr.plot_image(img)

def get_dist_bd_line(coords):
    """
    Get the distance of the coordinates from the birth death line
    """
    x = coords[0]
    y = coords[1]
    d = (y-x) / np.sqrt(2)
    return d

def create_pca_plot(vectors, n_components=3, labels=None):
    """
    Create a PCA plot for the given list of vectors.

    Parameters:
    vectors (list of lists): A list of vectors (each vector is a list of numbers).
    n_components (int): Number of principal components to keep.

    Returns:
    PCA plot displayed using matplotlib.
    """

    # Step 1: Standardize the dataset (important for PCA)
    scaler = StandardScaler()
    standardized_vectors = scaler.fit_transform(vectors)

    # Step 2: Perform PCA
    pca = PCA(n_components=n_components)
    pca_results = pca.fit_transform(standardized_vectors)

    # Step 3: Plot the results
    plt.figure(figsize=(8, 6))
    plt.scatter(pca_results[:, 0], pca_results[:, 1])
    for i, txt in enumerate(labels):
        plt.annotate(txt, (pca_results[i, 0], pca_results[i, 1]))
    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.title('PCA Plot')
    plt.grid(True)
    plt.show()




files = glob.glob('../tests/comparison/*.pdb')
vectors = []
file = files[0]
get_persistance_image(file)
plt.tight_layout()
plt.show()
#    vectors.append(generate_descriptor(file, 5))


#create_pca_plot(vectors, labels=files)
