from pylab import *

output = loadtxt('output.txt')

ion()
hold(False)
for i in xrange(0, output.shape[0]):
	x = output[i, :]
	imshow(x.reshape((198, 198)))
	draw()

ioff()
show()

