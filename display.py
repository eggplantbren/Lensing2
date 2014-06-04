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
	src = x[466:466 + metadata[0]*metadata[1]*metadata[7]**2]
	src = src.reshape((metadata[0]*metadata[7], metadata[1]*metadata[7]))
	img = x[466 + (metadata[0]*metadata[1]*metadata[7]**2):]
	img = img.reshape((metadata[0], metadata[1]))

	subplot(2,2,1)
	imshow(src)
	title('Model Source ' + str(i+1))

	subplot(2,2,2)
	imshow(img)
	title('Model Image ' + str(i+1))

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

