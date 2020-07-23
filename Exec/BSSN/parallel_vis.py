import yt

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import matplotlib.animation as animation

import numpy as np

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("infiles", type=str, nargs="+", help="Name of input plotfiles to process.")
parser.add_argument("-a", "--axis", type=str, default="z", help="Name of axis to slice across (default is 'z')")
parser.add_argument("-f", "--field", type=str, default="alpha", help="Name of field to plot (default is 'alpha')")
parser.add_argument("-fps", "--framerate", type=int, default=20, help="Framerate for saving the animation (default is 20 [frames/second])")
parser.add_argument("-r", "--resolution", type=int, nargs=2, default=[1024, 1024], help="Resolution of the image frames (default is 1024 x 1024 images)")
parser.add_argument("-s", "--save_slices", action="store_true", help="If supplied, save all the image slices used for the animation.")
args = parser.parse_args()

if __name__ == "__main__":
    # turn on MPI parallelism in yt
    yt.enable_parallelism()

    # this time series holds all the plotfiles we want to process
    time_series = yt.load(args.infiles)

    # shared storage across MPI
    storage = {}

    # use MPI processes to loop through the plotfiles
    for sto, ds in time_series.piter(storage=storage):
        slc = ds.slice(args.axis, 0.0)

        # slice customizations could go here ...


        # convert slice to raw image
        frb = slc.to_frb(ds.domain_width[0], args.resolution, center=ds.domain_center)
        image = np.array(frb[args.field])

        # save the raw image in the shared storage
        sto.result = image
        sto.result_id = str(ds)
        ds.index.clear_all_data()

    # only generate output from the root MPI rank
    if yt.is_root():
        # create a figure to work with
        fig = plt.figure()

        # make a list of all the artists for animation later
        image_list = []

        # loop through the raw images from each plotfile
        # and make matplotlib artists
        for result_id, result in sorted(storage.items()):
            # generate the matplotlib artists from the raw image
            image = plt.imshow(result, animated=True)
            image_list.append([image])
            # save each of the figures
            if args.save_slices:
                plt.savefig("{}_slice_z_{}.png".format(result_id, args.field))

        # animate the figures into a movie
        frame_interval_millisecond = 1.0e3 / float(args.framerate)
        ani = animation.ArtistAnimation(fig, image_list, 
                                        interval=frame_interval_millisecond,
                                        blit=True, repeat_delay=0)

        # save the movie
        plt.show() # this is necessary to set up the figure canvas
        ani.save("{}.mp4".format(args.field), fps=args.framerate)


