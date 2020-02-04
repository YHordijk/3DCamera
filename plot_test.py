import numpy as np
import modules.plot as plot
import modules.colour_maps as cmap


p = plot.Plot()





f = lambda x,y: x-y

x = np.array([np.arange(-5,7.5, 2.5)],)
y = np.array([np.arange(-5,7.5, 2.5)],)
z = f(x.T, y)





colour = cmap.Viridis()
p.plot_3d(x,y,z, colour_map=colour, style='contour')


# for i in range(5):
# 	x = np.arange(-1,1.1,0.1)
# 	y = 2*i*x**2
# 	p.plot(x,y, style='scatter')

p.title = f'Voorbeeld heatmap met kleurenschema {colour}'
p.show()