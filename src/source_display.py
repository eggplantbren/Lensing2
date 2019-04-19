import dnest4.classic as dn4
import matplotlib.pyplot as plt
import numba
import numpy as np

@numba.jit
def evaluate_blobs(x, y, blobs):
    """
    Evaluate a surface brightness profile
    """

    img = np.zeros(x.shape)
    for i in range(blobs.shape[0]):

        # Unpack blob parameters
        xc, yc, M, w = blobs[i, :]
        rsq = (x - xc)**2 + (y - yc)**2
        inside = rsq <= w**2
        img[inside] += 2.0*M/(np.pi*w**2)*(1.0 - rsq[inside]/w**2)

    return img

# Load the posterior samples
print("Loading posterior samples...", end="", flush=True)
posterior_sample = dn4.my_loadtxt("posterior_sample.txt")
indices = dn4.load_column_names("posterior_sample.txt")["indices"]
print("done.")

# Create a grid based on the properties of the run
import yaml
f = open("run_files.yaml")
metadata_filename = yaml.load(f, Loader=yaml.SafeLoader)["metadata_file"]
f.close()
f = open(metadata_filename)
the_dict = yaml.load(f, Loader=yaml.SafeLoader)["dimensions"]
x_min, x_max, y_min, y_max = the_dict["x_min"], the_dict["x_max"],\
                             the_dict["y_min"], the_dict["y_max"]
nj, ni = the_dict["nj"], the_dict["ni"]
f.close()

# Higher resolution than the inference used
nj *= 10
ni *= 10
if nj*ni > 100000000:
    exit()

# Create new grid
x_range = x_max - x_min
x_min = x_min + 0.25*x_range
x_max = x_max - 0.25*x_range
y_range = y_max - y_min
y_min = y_min + 0.25*y_range
y_max = y_max - 0.25*y_range

dx = (x_max - x_min)/(nj - 1)
dy = (y_max - y_min)/(ni - 1)
x = np.linspace(x_min + 0.5*dx, x_max - 0.5*dx, nj)
y = np.linspace(y_min + 0.5*dy, y_max - 0.5*dy, ni)
x, y = np.meshgrid(x, y)
y = y[::-1, :]


tot = np.zeros((ni, nj))

# Loop over the posterior samples
for i in range(posterior_sample.shape[0]):

    # Extract the blob parameters
    num_blobs = int(posterior_sample[i, indices["num_source_blobs"]])
    blobs = []
    for j in range(num_blobs):
        xc = posterior_sample[i, indices["source_blob_x[{j}]".format(j=j)]]
        yc = posterior_sample[i, indices["source_blob_y[{j}]".format(j=j)]]
        M  = posterior_sample[i, indices["source_blob_mass[{j}]".format(j=j)]]
        w  = posterior_sample[i, indices["source_blob_width[{j}]".format(j=j)]]
        blobs.append([xc, yc, M, w])
    blobs = np.array(blobs)

    # Approximate center of mass (treats blob widths as zero)
    x_com = np.sum(blobs[:,0]*blobs[:,2]) / np.sum(blobs[:,2])
    y_com = np.sum(blobs[:,1]*blobs[:,2]) / np.sum(blobs[:,2])

    # Re-center blobs
    blobs[:,0] -= x_com
    blobs[:,1] -= y_com

    # Evaluate and display the source model
    img = evaluate_blobs(x, y, blobs)
    tot += img

    print(i+1, "/", posterior_sample.shape[0])
#    plt.imshow(img, origin="lower", extent=[x_min, x_max, y_min, y_max])
#    plt.show()

plt.imshow(tot/posterior_sample.shape[0],
           extent=[x_min, x_max, y_min, y_max])
plt.title("Posterior mean source")
plt.show()

