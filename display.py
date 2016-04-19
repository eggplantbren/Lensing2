"""
A little script to load the output and plot the model images
and the residuals.
"""
from pylab import *
import os
import dnest4.deprecated as dn4

mass_units = 1.0
os.system("rm Frames/*.png")
os.system("rm movie.mkv")

output = atleast_2d(dn4.my_loadtxt('posterior_sample.txt'))
data = loadtxt('Data/mock_image.txt')
sig = loadtxt('Data/mock_sigma.txt')
not_masked = (sig < 1E100)
metadata = loadtxt('Data/mock_metadata.txt')

total = zeros((metadata[0]*metadata[7], metadata[1]*metadata[7]))
magnification = zeros(output.shape[0])
all_substructures_x = array([])
all_substructures_y = array([])

figure(figsize=(12, 8))
hold(False)
for i in range(0, output.shape[0]):
	x = output[i, :]

	# Extract substructure information
	n_substructures = x[19]
	x_substructures = x[20:30]
	y_substructures = x[30:40]
	m_substructures = x[40:50]

	# Remove substructures out of image boundaries (don't plot these)
	good = logical_and(x_substructures > metadata[2],
			x_substructures < metadata[3])
	good = logical_and(good, y_substructures > metadata[4])
	good = logical_and(good, y_substructures < metadata[5])
	good = logical_and(good, m_substructures > 0.)
	x_substructures = x_substructures[good]
	y_substructures = y_substructures[good]
	m_substructures = m_substructures[good]
	x_nie, y_nie = x[6], x[7]

	# Extract images
	# For MyModel2 (sersic source only), replace 468 with 66
	# Sersic source model parameters are columns 59-65
	src = x[469:469 + metadata[0]*metadata[1]*metadata[7]**2]
	src = src.reshape((metadata[0]*metadata[7], metadata[1]*metadata[7]))

	img1 = x[469 + metadata[0]*metadata[1]*metadata[7]**2:469 + 2*metadata[0]*metadata[1]*metadata[7]**2]
	img1 = img1.reshape((metadata[0]*metadata[7], metadata[1]*metadata[7]))

	img2 = x[469 + 2*metadata[0]*metadata[1]*metadata[7]**2:]
	img2 = img2.reshape((metadata[0], metadata[1]))

	subplot(2,3,1)
	imshow(src, extent=metadata[2:6], interpolation='nearest', cmap='viridis')
	title('Model Source ' + str(i+1))
	axis(metadata[2:6])

	subplot(2,3,2)
	imshow(img1, extent=metadata[2:6], interpolation='nearest', cmap='viridis')
	title('Unblurred Image ' + str(i+1))
	axis(metadata[2:6])

	subplot(2,3,3)
	imshow(img2, extent=metadata[2:6], interpolation='nearest', cmap='viridis')
	hold(True)
	# Plot center of NIE
	plot(x_nie, y_nie, 'wo', markersize=15, alpha=0.5)
	# Substructures
	plot(x_substructures, y_substructures, 'g*', markersize=15, alpha=0.5)
	title('Model Image ' + str(i+1))
	hold(False)
	axis(metadata[2:6])

	subplot(2,2,3)
	imshow(data*not_masked, extent=metadata[2:6], interpolation='nearest', cmap='viridis')
	title('Data')
	axis(metadata[2:6])

	subplot(2,2,4)
	sigma = np.ones(sig.shape)
	sigma[not_masked] = sqrt(sig[not_masked]**2 + x[0]**2 + x[1]*img2[not_masked])
	imshow(((img2 - data)/sigma)*not_masked, extent=metadata[2:6], interpolation='nearest', cmap='viridis')
	title('Standardised Residuals')
	axis(metadata[2:6])
	draw()

	savefig('Frames/' + '%0.6d'%(i+1) + '.png', bbox_inches='tight')
	print('Frames/' + '%0.6d'%(i+1) + '.png')

	total += src
	magnification[i] = 2.5*log10(metadata[7]**2*img2.sum()/src.sum())

	# Accumulate these before converting to pixels
	all_substructures_x = hstack([all_substructures_x, x_substructures])
	all_substructures_y = hstack([all_substructures_y, y_substructures])

show()

os.system('ffmpeg -r 10 -i Frames/%06d.png -c:v h264 -b:v 4192k movie.mkv')

figure(1)
mean_source = total/output.shape[0]
imshow(mean_source, interpolation='nearest', cmap='viridis')
title('Posterior Mean Source')

figure(2)
hist(magnification, 50, alpha=0.5)
xlabel('Magnification (magnitudes)')
title('Magnification = {a:.3f} +- {b:.3f}'.format(a=magnification.mean(), b=magnification.std()))
show()

figure(3)
rc("font", size=16, family="serif", serif="Computer Sans")
rc("text", usetex=True)
plot(mass_units*pi*output[:,2]**2, mass_units*output[:,39:49].sum(axis=1),\
					'bo', markersize=5, alpha=0.2)
xlabel('SIE mass (within critical ellipse)')
ylabel('Total substructure mass')

# Plot the true masses
#truth = loadtxt('Data/mock_truth.txt')
#hold(True)
#plot(pi*truth[2]**2, truth[39:49].sum(), 'r*', markersize=20)
savefig('masses.pdf', bbox_inches='tight')
show()

width=0.6
hist(output[:,19], bins=arange(0, 10) - 0.5*width, width=width, alpha=0.5)
xlim([-0.2, 10.2])
gca().set_xticks(arange(0, 11))
xlabel('$N_{\\rm lens}$')
ylabel('Number of Posterior Samples')
savefig('N_lens.pdf', bbox_inches='tight')
show()

figure(4, figsize=(8, 8))
imshow(data*not_masked, extent=metadata[2:6], interpolation="nearest", cmap='viridis')
hold(True)
plot(all_substructures_x, all_substructures_y, 'w.', markersize=3)
gca().set_xticks([])
gca().set_yticks([])
axis(metadata[2:6])
title('Substructure Positions')
savefig('substructures.pdf', bbox_inches='tight')
show()

