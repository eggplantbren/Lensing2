from pylab import *

output = atleast_2d(loadtxt('posterior_sample.txt'))
data = loadtxt('Data/test_image.txt')

ion()
hold(False)
for i in xrange(0, output.shape[0]):
	x = output[i, 5:]
	imshow(x.reshape((100, 100)) - data)
	draw()

ioff()
show()

