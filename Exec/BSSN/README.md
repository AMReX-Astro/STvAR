# BSSN setup

## Visualization

The `parallel_vis.py` script can be used in an MPI environment to make plots
and/or animations from a time series of multiple plotfiles.

To get a list of command line options, use the `-h` flag.

For example, the following uses 4 MPI tasks to process all the plotfiles in the
current directory with the following options:

- `--axis z`: generate slices along `z=0`. The `-a` flag can take `x`, `y`, or `z` options.
- `--field alpha`: make slice plots using the `alpha` field in the plotfiles.
- `--resolution 1200 1200`: sets the image resolution for the slices to 1200 x 1200.
- `--framerate 1`: generate the animation at 1 frame per second.
- `--save_slices`: save the images generated from each plotfile in the current directory.

```
mpiexec -n 4 python3 parallel_vis.py plt[0-9][0-9][0-9][0-9][0-9][0-9][0-9] --axis z --field alpha --framerate 1 --save_slices --resolution 1200 1200
```
