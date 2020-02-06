import numpy as np
import modules.plot as plot
import modules.colour_maps as cmap


p = plot.Plot()





f = lambda x,y: x**3 + y**2

x = np.array([np.arange(-1,1.04, 0.04)],)
y = np.array([np.arange(-1,1.04, 0.04)],)
z = f(x.T, y)





colour = cmap.CoolWarm()
p.plot_3d(x,y,z, colour_map=colour, style='contour')


# for i in range(7):
# 	x = np.arange(-i+1,i+1+.1,0.1)
# 	y = x**3
# 	p.plot(x,y, style='line', width=2)

p.x_label = 'x'
p.y_label = 'y'
p.title = f'Voorbeeld heatmap met kleurenschema {colour}'
p.show()