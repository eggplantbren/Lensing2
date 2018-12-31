"""
A little script to load the output and plot the model images
and the residuals. The Python here is hacky and old, please
forgive the inelegance.
"""
from pylab import *
import os
import dnest4.deprecated as dn4

print("WARNING! This will delete\
 movie.mp4 and the Frames/ directory, if these exist.")
ch = input("Continue? y/n: ")
if ch != "y" and ch != "Y":
    exit()

def blob_density(x, y, params):
    xc, yc, mass, width = params
    rsq = (x - params[0])**2 + (y - params[1])**2
    widthsq = width**2
    f = zeros(x.shape)
    f[rsq < widthsq] = 2*mass/pi*(1 - rsq[rsq < widthsq]/widthsq)/widthsq
    return f

def stretch(img):
    """
    A generic stretch function that just saturates at the 99th percentile
    of nonzero values
    """
    if img.max() == 0:
        return img

    values = np.sort(img.flatten())
    values = values[values > 0.0]
    peak = values[int(0.99*len(values))]
    img[img >= peak] = peak
    return img

os.system("rm -rf Frames/ movie.mp4")
os.mkdir("Frames")

mass_units = 1.0

output = dn4.my_loadtxt('posterior_sample.txt')
indices = dn4.load_column_names("posterior_sample.txt")["indices"]

# Open run_files.yaml to get data filenames used for the run
import yaml
f = open("run_files.yaml")
run_files = yaml.load(f)
f.close()
a, b, c, d = run_files["metadata_file"],\
             run_files["image_file"],\
             run_files["sigma_file"],\
             run_files["psf_file"]

# Load images etc
data = loadtxt(b)
sig = loadtxt(c)
not_masked = (sig < 1E100)

# Load metadata
f = open(a)
metadata_dict = yaml.load(f)
f.close()
# Put it into a list because silly old me did it that way
metadata = [0 for i in range(0, 9)]
metadata[0] = metadata_dict["dimensions"]["ni"]
metadata[1] = metadata_dict["dimensions"]["nj"]
metadata[2] = metadata_dict["dimensions"]["x_min"]
metadata[3] = metadata_dict["dimensions"]["x_max"]
metadata[4] = metadata_dict["dimensions"]["y_min"]
metadata[5] = metadata_dict["dimensions"]["y_max"]
metadata[6] = metadata_dict["psf"]["num_pixels"]
metadata[7] = metadata_dict["computation"]["nrays"]
metadata[8] = int(metadata_dict["psf"]["is_highres"])

# A grid
dx = (metadata[3] - metadata[2])/(metadata[7]*metadata[1])
dy = (metadata[5] - metadata[4])/(metadata[7]*metadata[0])
xgrid = linspace(metadata[2] + 0.5*dx, metadata[3] - 0.5*dx,\
            metadata[7]*metadata[1])
ygrid = linspace(metadata[4] + 0.5*dy, metadata[5] - 0.5*dy,\
            metadata[7]*metadata[0])

[xgrid, ygrid] = meshgrid(xgrid, ygrid)
ygrid = ygrid[::-1, :]

total = zeros((metadata[0]*metadata[7], metadata[1]*metadata[7]))
magnification = zeros(output.shape[0])
all_substructures_x = array([])
all_substructures_y = array([])
substructure_num_in_image = []
substructure_mass_in_image = array([])

max_num_blobs = int(output[0, indices["max_num_lens_blobs"]])

figure(figsize=(14, 9))
for i in range(0, output.shape[0]):
    clf()
    x = output[i, :]

    # Extract substructure information
    n_substructures = x[indices["num_lens_blobs"]]
    x_substructures = x[indices["lens_blob_x[0]"]:indices["lens_blob_x[0]"] + max_num_blobs]
    y_substructures = x[indices["lens_blob_y[0]"]:indices["lens_blob_y[0]"] + max_num_blobs]
    m_substructures = x[indices["lens_blob_mass[0]"]:indices["lens_blob_mass[0]"] + max_num_blobs]
    w_substructures = x[indices["lens_blob_width[0]"]:indices["lens_blob_width[0]"] + max_num_blobs]

    # Remove substructures out of image boundaries (don't plot these)
#    good = logical_and(x_substructures > metadata[2],
#            x_substructures < metadata[3])
#    good = logical_and(good, y_substructures > metadata[4])
#    good = logical_and(good, y_substructures < metadata[5])
    good = (m_substructures > 0.)
    x_substructures = x_substructures[good]
    y_substructures = y_substructures[good]
    m_substructures = m_substructures[good]
    w_substructures = w_substructures[good]
    x_nie, y_nie = x[indices["xc"]], x[indices["yc"]]

    # Extract images
    # For MyModel2 (sersic source only), replace 468 with 66
    # Sersic source model parameters are columns 59-65
    src = x[indices["source[0][0]"]:indices["source[0][0]"] + metadata[0]*metadata[1]*metadata[7]**2]
    src = src.reshape((metadata[0]*metadata[7], metadata[1]*metadata[7]))

    img1 = x[indices["source[0][0]"] + metadata[0]*metadata[1]*metadata[7]**2:indices["source[0][0]"] + 2*metadata[0]*metadata[1]*metadata[7]**2]
    img1 = img1.reshape((metadata[0]*metadata[7], metadata[1]*metadata[7]))

    img2 = x[indices["source[0][0]"] + 2*metadata[0]*metadata[1]*metadata[7]**2:]
    img2 = img2.reshape((metadata[0], metadata[1]))

    subplot(2,3,1)
    imshow(src, extent=metadata[2:6], interpolation='nearest', cmap='viridis')
    title('Model Source ' + str(i+1))
    axis(metadata[2:6])

    subplot(2,3,2)
    imshow(log(img1),
           extent=metadata[2:6], interpolation='nearest', cmap='viridis')
    title('Lens profile ' + str(i+1) + ' (log stretch)')
    axis(metadata[2:6])

    subplot(2,3,3)
    imshow(img2, extent=metadata[2:6], interpolation='nearest', cmap='viridis')
    # Plot center of NIE
    plot(x_nie, y_nie, 'wo', markersize=15, alpha=0.5)
    # Substructures
    plot(x_substructures, y_substructures, 'w*', markersize=15, alpha=0.5)
    title('Model Image ' + str(i+1))
    axis(metadata[2:6])

    subplot(2,3,4)
    substructure_density = zeros(xgrid.shape)
    for j in range(0, int(n_substructures)):
        substructure_density += blob_density(xgrid, ygrid,\
                [x_substructures[j], y_substructures[j],\
                    m_substructures[j], w_substructures[j]])
    imshow(substructure_density, cmap="viridis", interpolation="nearest",\
            extent=metadata[2:6])
    title("Substructure Map")

    subplot(2,3,5)
    temp = data*not_masked
    imshow(data*not_masked, extent=metadata[2:6],\
            vmin=-0.1*temp.max(), vmax=temp.max(),\
            interpolation='nearest', cmap='viridis')
    title('Data')
    axis(metadata[2:6])

    subplot(2,3,6)
    sigma = np.ones(sig.shape)
    sigma[not_masked] = sqrt(sig[not_masked]**2 + x[0]**2 + x[1]*(img2[not_masked] - img2[not_masked].min()))
    imshow((img2 - data)*not_masked, extent=metadata[2:6], interpolation='nearest', cmap='coolwarm')
    title('Residuals')
    axis(metadata[2:6])
    draw()

    savefig('Frames/' + '%0.6d'%(i+1) + '.png', bbox_inches='tight')
    print('Frames/' + '%0.6d'%(i+1) + '.png')

    total += src
    magnification[i] = 2.5*log10(metadata[7]**2*img2.sum()/src.sum())

    # Accumulate these before converting to pixels
    all_substructures_x = hstack([all_substructures_x, x_substructures])
    all_substructures_y = hstack([all_substructures_y, y_substructures])

    inside = (x_substructures > metadata[2]) &\
             (x_substructures < metadata[3]) &\
             (y_substructures > metadata[4]) &\
             (y_substructures < metadata[5])

    substructure_num_in_image.append(sum(inside))
    substructure_mass_in_image = hstack([substructure_mass_in_image,\
      m_substructures[inside].sum()])

show()

rcParams["font.family"] = "serif"
rcParams["font.size"] = 16
rc("text", usetex=True)

figure(1)
mean_source = total/output.shape[0]
imshow(mean_source, interpolation='nearest', cmap='viridis')
title('Posterior Mean Source')

try:
    figure(2)
    hist(magnification, 50, alpha=0.5, color="k")
    xlabel('Magnification (magnitudes)')
    title('Magnification = {a:.3f} +- {b:.3f}'.format(a=magnification.mean(), b=magnification.std()))
    show()
except:
    pass

figure(3)
plot(output[:,indices["b"]], mass_units*output[:,indices["lens_blob_mass[0]"]:indices["lens_blob_mass[0]"]+max_num_blobs].sum(axis=1),\
                    'k.', alpha=0.2, label="Total")
plot(output[:,indices["b"]], mass_units*array(substructure_mass_in_image),\
                    'g.', alpha=0.2, label="Center within image")
xlabel('SPEMD Einstein Radius')
ylabel('Total substructure mass')
legend(loc="upper right", numpoints=1)

# Plot the true masses
#truth = loadtxt('Data/mock_truth.txt')
#plot(pi*truth[2]**2, truth[39:49].sum(), 'r*', markersize=20)
savefig('masses.pdf', bbox_inches='tight')
show()

try:
    width=0.6
    hist(output[:,indices["num_lens_blobs"]], bins=arange(0, 51) - 0.5*width, width=width, alpha=0.2, color="k",
            label="Total")
    hist(substructure_num_in_image, bins=arange(0, 51) - 0.5*width, width=width, alpha=0.2,
            label="Center within image", color="g")
    xlim([-0.5, max_num_blobs+0.5])
    xlabel('$N_{\\rm lens}$')
    ylabel('Number of Posterior Samples')
    legend(loc="upper right")
    savefig('N_lens.pdf', bbox_inches='tight')
    show()
except:
    pass

figure(4, figsize=(8, 8))
temp = data*not_masked
imshow(data*not_masked, extent=metadata[2:6],\
        vmin=-0.1*temp.max(), vmax=temp.max(),\
        interpolation='nearest', cmap='viridis')
plot(all_substructures_x, all_substructures_y, 'w.', markersize=3, alpha=0.1)
gca().set_xticks([])
gca().set_yticks([])
axis(metadata[2:6])
title('Substructure Positions')
savefig('substructures.pdf', bbox_inches='tight')
show()

# Make movie
os.system('ffmpeg -r 10 -i Frames/%06d.png -c:v h264 movie.mp4')

