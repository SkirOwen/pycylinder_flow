import numpy as np
import scipy.io
import matplotlib.pyplot as plt



maxT   = 5000;  # total number of iterations
tPlot  = 1 * 50;  # cycles

# load('TestSymVars.mat')
data = scipy.io.loadmat('TestSymVars.mat')

x = data['x'] - data['obst_x']
y = data['y'] - data['obst_y']
d = 2 * data['obst_r']


## visualize data

myContours = np.linspace(0, 1.4 * data['uMax'], 11)

for cycle in range(2000, maxT + 1, tPlot): #maxT
	# fileName = sprintf('Re100Data/velocity_%d.mat', cycle);
	fileName = f"Testvelocity_{cycle}.mat"
	test_data = scipy.io.loadmat(fileName)
	u = np.sqrt(test_data['vx'] ** 2 + test_data['vy'] ** 2)

	plt.figure(2)
	plt.contourf(x / d, y / d, u, myContours, edgecolor='none')
	plt.gca().set_aspect('equal')
	plt.xlabel('x/d')
	plt.ylabel('y/d')
	plt.title(f'Time step {cycle}')
	plt.colorbar()
	plt.clim(min(myContours), max(myContours))
	plt.show()

	# figure(2)
	# contourf(x/d,y/d,u,myContours,'edgecolor','none'), hold on
	# daspect([1 1 1])
	# xlabel('x/d')
	# ylabel('y/d')
	# title(sprintf('Time step %d', cycle))
	# #plot(x(300,50)/d,y(300,50)/d,'ko')
	# colormap(viridisCMap)
	# colorbar
	# caxis([min(myContours),max(myContours)])
# end

## load data at a given point and plot
# k = 1;
# clear up
# for cycle = 5000:tPlot:maxT
#     fileName = sprintf('Re100Data/velocity_%d.mat', cycle);
#     load(fileName)
#     u = sqrt(vx.^2 + vy.^2);
#     up(k) = u(300,50);
#     k = k+1;
# end

##
load_data = scipy.io.loadmat('u_vs_t_downstream_v1.mat')

ti = np.arange(1, len(load_data['up']) + 1)
plt.figure(10)
plt.plot(ti, load_data['up'], '-r')
plt.xlabel('Timestep')
plt.ylabel('Flow speed, m/s')
plt.grid(True)
plt.show()

# ti = (1:length(up));
# figure(10), clf
# plot(ti,up,'-r')
# xlabel('Timestep')
# ylabel('Flow speed, m/s')
# grid on

