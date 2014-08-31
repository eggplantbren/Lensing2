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
data = loadtxt('Data/test_image.txt')
sig = loadtxt('Data/test_sigma.txt')
metadata = loadtxt('Data/test_metadata.txt')

figure(figsize=(12, 8))
ion()
hold(False)
for i in xrange(0, output.shape[0]):
	x = output[i, :]

	# Extract substructure information
	n_substructures = x[17]
	x_substructures = x[18:28]
	y_substructures = x[28:38]
	m_substructures = x[38:48]

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
	src = x[466:466 + metadata[0]*metadata[1]*metadata[7]**2]
	src = src.reshape((metadata[0]*metadata[7], metadata[1]*metadata[7]))
	img = x[466 + (metadata[0]*metadata[1]*metadata[7]**2):-2]
	img = img.reshape((metadata[0], metadata[1]))

	subplot(2,2,1)
	imshow(src, interpolation='nearest')
	title('Model Source ' + str(i+1))

	subplot(2,2,2)
	imshow(img, interpolation='nearest')
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
	sigma = sqrt(sig**2 + x[0]**2 + x[1]*img)
	imshow((img - data)/sigma, interpolation='nearest')
	title('Standardised Residuals')
	draw()

	if saveFrames:
		savefig('Frames/' + '%0.4d'%(i+1) + '.png', bbox_inches='tight')
		print('Frames/' + '%0.4d'%(i+1) + '.png')

ioff()
show()

