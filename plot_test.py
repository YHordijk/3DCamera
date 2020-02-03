import numpy as np
import modules.plot as plot


p = plot.Plot()






for i in range(10):
	x = np.arange(-1,1.05,0.05)
	y = i*x**4
	p.plot(x,y, style='line')


p.show((400,400))