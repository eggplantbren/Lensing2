"""
A little script to load the output and plot the model images
and the residuals.
"""
from pylab import *

output = atleast_2d(loadtxt('posterior_sample.txt'))
data = loadtxt('Data/test_image.txt')

figure(figsize=(12, 8))
ion()
hold(False)
for i in xrange(0, output.shape[0]):
	x = output[i, :]
	src = x[58:58 + 200**2]
	src = src.reshape((200, 200))
	img = x[58 + 200**2:]
	img = img.reshape((100, 100))

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
	imshow(img - data)
	title('Residuals')
	draw()

ioff()
show()

