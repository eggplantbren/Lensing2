"""
A little script to load the output and plot the model images
and the residuals.
"""
from pylab import *
import os

saveFrames = False # For making movies
if saveFrames:
	os.system('rm Frames/*.png')

output = atleast_2d(loadtxt('posterior_sample.txt'))
data = loadtxt('Data/mock_image.txt')
sig = loadtxt('Data/mock_sigma.txt')
metadata = loadtxt('Data/mock_metadata.txt')

total = zeros((metadata[0]*metadata[7], metadata[1]*metadata[7]))
magnification = zeros(output.shape[0])
all_substructures_x = array([])
all_substructures_y = array([])

figure(figsize=(12, 8))
ion()
hold(False)
for i in xrange(0, output.shape[0]):
	x = output[i, :]

	# Extract substructure information
	n_substructures = x[18]
	x_substructures = x[19:29]
	y_substructures = x[29:39]
	m_substructures = x[39:49]

	# Remove substructures out of image boundaries (don't plot these)
	good = logical_and(x_substructures > metadata[2],
			x_substructures < metadata[3])
	good = logical_and(good, y_substructures > metadata[4])
	good = logical_and(good, y_substructures < metadata[5])
	good = logical_and(good, m_substructures > 0.)
	x_substructures = x_substructures[good]
	y_substructures = y_substructures[good]
	m_substructures = m_substructures[good]

	# Convert x and y to pixel coordinates for overplotting
	dx = (metadata[3] - metadata[2])/metadata[1]
	dy = (metadata[5] - metadata[4])/metadata[0]
	x_substructures = (x_substructures - metadata[2])/dx - 0.5
	y_substructures = (metadata[5] - y_substructures)/dy - 0.5
	x_nie = (x[5] - metadata[2])/dx - 0.5
	y_nie = (metadata[5] - x[6])/dy - 0.5

	# Extract images
	# For MyModel2 (sersic source only), replace 468 with 66
	# Sersic source model parameters are columns 59-65
	src = x[468:468 + metadata[0]*metadata[1]*metadata[7]**2]
	src = src.reshape((metadata[0]*metadata[7], metadata[1]*metadata[7]))

	img1 = x[468 + metadata[0]*metadata[1]*metadata[7]**2:468 + 2*metadata[0]*metadata[1]*metadata[7]**2]
	img1 = img1.reshape((metadata[0]*metadata[7], metadata[1]*metadata[7]))

	img2 = x[468 + 2*metadata[0]*metadata[1]*metadata[7]**2:-2]
	img2 = img2.reshape((metadata[0], metadata[1]))

	subplot(2,3,1)
	imshow(src, interpolation='nearest')
	title('Model Source ' + str(i+1))

	subplot(2,3,2)
	imshow(img1, interpolation='nearest')
	title('Unblurred Image ' + str(i+1))

	subplot(2,3,3)
	imshow(img2, interpolation='nearest')
	hold(True)
	# Plot center of NIE
	plot(x_nie, y_nie, 'yo', markersize=10)
	# Substructures
	plot(x_substructures, y_substructures, 'wo', alpha=0.5)
	xlim([-0.5, metadata[1] - 0.5])
	ylim([metadata[0] - 0.5, -0.5])
	title('Model Image ' + str(i+1))
	hold(False)

	subplot(2,2,3)
	imshow(data, interpolation='nearest')
	title('Data')

	subplot(2,2,4)
	sigma = sqrt(sig**2 + x[0]**2 + x[1]*img2)
	imshow((img2 - data)/sigma, interpolation='nearest')
	title('Standardised Residuals')
	draw()

	if saveFrames:
		savefig('Frames/' + '%0.4d'%(i+1) + '.png', bbox_inches='tight')
		print('Frames/' + '%0.4d'%(i+1) + '.png')

	total += src
	magnification[i] = 2.5*log10(metadata[7]**2*img2.sum()/src.sum())

	# Accumulate these before converting to pixels
	all_substructures_x = hstack([all_substructures_x, x_substructures])
	all_substructures_y = hstack([all_substructures_y, y_substructures])

ioff()
show()

figure(1)
mean_source = total/output.shape[0]
imshow(mean_source, interpolation='nearest')
title('Posterior Mean Source')

figure(2)
hist(magnification, 50, alpha=0.5)
xlabel('Magnification (magnitudes)')
title('Magnification = {a:.3f} +- {b:.3f}'.format(a=magnification.mean(), b=magnification.std()))
show()

figure(3)
rc("font", size=16, family="serif", serif="Computer Sans")
rc("text", usetex=True)
plot(output[:,2]**2, output[:,39:48].sum(axis=1), 'bo', markersize=5,\
                       alpha=0.2)
xlabel('SIE mass (within critical ellipse)')
ylabel('Total substructure mass')
savefig('masses.pdf', bbox_inches='tight')
show()

hist(output[:,18], bins=arange(0, 11) - 0.5*0.3, width=0.3, alpha=0.3)
xlim([-0.2, 10.2])
gca().set_xticks(arange(0, 11))
xlabel('$N_{\\rm lens}$')
ylabel('Number of Posterior Samples')
savefig('N_lens.pdf', bbox_inches='tight')
show()

figure(4, figsize=(8, 8))
imshow(data, interpolation="nearest", cmap="Oranges")
hold(True)
plot(all_substructures_x, all_substructures_y, 'k.', markersize=3)
gca().set_xticks([])
gca().set_yticks([])
xlim([-0.5, metadata[1] - 0.5])
ylim([metadata[0] - 0.5, -0.5])
title('Substructure Positions')
savefig('substructures.pdf', bbox_inches='tight')
show()

