# Py Cylinder Flow

`cylinder.py` runs the flow field over a grid of size `(lx, ly)` 
(preset as quite course in the file and you probably want to refine).  
To reach a stable regime `max_t > 2000`.

It spits out a bunch of data files containing the horizontal `vz`
and vertical `vy` velocity on the grid for each time step.

You can then use `plot_cylinder.m` to visualise the flow in increasing time steps.
Note, here you only want to take the flow from timesteps 
2000+ since the early time is literally 'turning on' the upstream flow
(from zero) so there is some funky behaviour we don't want to consider there.
Seems like you have to manually copy over maxT and tPlot from whatever 
you selected in the previous case.
I'd check the visualisation to make sure you're in the 'uniform flow' 
regime and not the setup time region.

You can play with the obstacle size and position (obst_{x,y,r}),
and the Reynolds number to generate more data.

## Original code
The original (unmodified) code in MATLAB is located in `matlab_original`.


## License