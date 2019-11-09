#!/usr/bin/env python
import os
import shutil
import subprocess
import numpy as np

def get_abs_error(fcompare_output):
    field_errors = {}
    fcompare_lines = fcompare_output.splitlines()
    fcompare_lines = fcompare_lines[5:]
    for line in fcompare_lines:
        lss = line.strip().split()
        field_errors[lss[0]] = lss[1]
    return field_errors

def get_error(directory, cfl_list, output_plotfile):
    error_data = {}
    error_data["cfl"] = []

    sorted_cfl = np.sort(np.array(cfl_list))
    sorted_cfl = sorted_cfl[::-1]
    cfl_ref = sorted_cfl[-1]
    cfl_compare = sorted_cfl[:-1]

    plt_ref = os.path.join(directory, "cfl-{}".format(cfl_ref), output_plotfile)

    for cfl in cfl_compare:
        plt_compare = os.path.join(directory, "cfl-{}".format(cfl), output_plotfile)
        comparison = subprocess.run("fcompare {} {}".format(plt_ref, plt_compare), stdout=subprocess.PIPE, shell=True)
        comparison = comparison.stdout.decode("UTF-8")
        error = get_abs_error(comparison)
        error_data["cfl"].append(cfl)
        for k in error.keys():
            if not k in error_data.keys():
                error_data[k] = []
            error_data[k].append(error[k])
    return error_data

def write_error(directory, error_data):
    error_file = os.path.join(directory, "errors.txt")
    f = open(error_file, 'w')
    f.write(directory + "\n")
    f.write("CFL")
    fields = []
    for k in error_data.keys():
        if k != "cfl":
            f.write("   {}".format(k))
            fields.append(k)
    f.write("\n")
    nlines = len(error_data["cfl"])
    for i in range(nlines):
        to_write = [error_data["cfl"][i]]
        for field in fields:
            to_write.append(error_data[field][i])
        to_write = [str(x) for x in to_write]
        line_to_write = "   ".join(to_write)
        f.write("{}\n".format(line_to_write))

def run_order(directory, cfl_list, run_command, args, output_plotfile, initial_plotfile=None):
    os.mkdir(directory)

    for cfl in cfl_list:
        plt_dir = os.path.join(directory, "cfl-{}".format(cfl))
        plt_dest = os.path.join(plt_dir, output_plotfile)
        os.mkdir(plt_dir)
        subprocess.call("{} {} {}".format(run_command, args, "cfl={}".format(cfl)), shell=True)
        shutil.move(output_plotfile, plt_dest)
        if initial_plotfile:
            shutil.rmtree(initial_plotfile)

    error = get_error(directory, cfl_list, output_plotfile)
    write_error(directory, error)

if __name__ == "__main__":
    run = "./main2d.gnu.ex inputs_waveeq nsteps=10000000 plot_int=100000000"
    output_plt = "plt_End_Simulation"
    initial_plt = "plt0000000"

    # FE Run parameters
    directory_fe = "fe"
    cfl_fe = np.array([1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125])
    args_fe = "end_time=0.001"

    # RK FE Run parameters
    directory_rk_fe = "rk_fe"
    cfl_rk_fe = np.array([1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125])
    args_rk_fe = "end_time=0.001 integration.type=1 integration.rk.type=1"

    # RK4 Run parameters
    directory_rk4 = "rk4"
    cfl_rk4 = np.array([1000, 500, 250, 125, 62.5, 31.25])
    args_rk4 = "end_time=1.0 integration.type=1 integration.rk.type=2"

    run_order(directory_fe, cfl_fe, run, args_fe, output_plt, initial_plt)
    run_order(directory_rk_fe, cfl_rk_fe, run, args_rk_fe, output_plt, initial_plt)
    run_order(directory_rk4, cfl_rk4, run, args_rk4, output_plt, initial_plt)