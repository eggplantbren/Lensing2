from pylab import *

output = atleast_2d(loadtxt('posterior_sample.txt'))

ion()
hold(False)
for i in xrange(0, output.shape[0]):
	x = output[i, :]
	imshow(x.reshape((100, 100)))
	draw()

ioff()
show()

