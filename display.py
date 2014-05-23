"""
A little script to load the output and plot the model images
and the residuals.
"""
from pylab import *

output = atleast_2d(loadtxt('posterior_sample.txt'))
data = loadtxt('Data/test_image.txt')

figure(figsize=(10, 6))
ion()
hold(False)
for i in xrange(0, output.shape[0]):
	x = output[i, 8:]
	img = x.reshape((100, 100))

	subplot(1,3,1)
	imshow(img)
	title('Model Image')

	subplot(1,3,2)
	imshow(data)
	title('Data')

	subplot(1,3,3)
	imshow(img - data)
	title('Residuals')
	draw()

ioff()
show()

