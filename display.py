"""
A little script to load the output and plot the model images
and the residuals.
"""
from pylab import *

output = atleast_2d(loadtxt('posterior_sample.txt'))
data = loadtxt('Data/harder_image.txt')
metadata = loadtxt('Data/test_metadata.txt')

figure(figsize=(12, 8))
ion()
hold(False)
for i in xrange(0, output.shape[0]):
	x = output[i, :]

	# Extract substructure information
	n_sumstructures = x[17]
	x_substructures = x[17:27]
	y_substructures = x[27:37]
	m_substructures = x[37:47]

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

	# Extract images
	src = x[466:466 + metadata[0]*metadata[1]*metadata[7]**2]
	src = src.reshape((metadata[0]*metadata[7], metadata[1]*metadata[7]))
	img = x[466 + (metadata[0]*metadata[1]*metadata[7]**2):]
	img = img.reshape((metadata[0], metadata[1]))

	subplot(2,2,1)
	imshow(src)
	title('Model Source ' + str(i+1))

	subplot(2,2,2)
	imshow(img)
	hold(True)
	plot(x_substructures, y_substructures, 'wo', alpha=0.5)
	xlim([-0.5, metadata[1] - 0.5])
	ylim([metadata[0] - 0.5, -0.5])
	title('Model Image ' + str(i+1))
	hold(False)

	subplot(2,2,3)
	imshow(data)
	title('Data')

	subplot(2,2,4)
	sigma = sqrt(x[0]**2 + x[1]*img)
	imshow((img - data)/sigma)
	title('Standardised Residuals')
	draw()

ioff()
show()

