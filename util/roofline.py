# @file roofline.py
# @author Jonas Kühne (jonas.kuehne@proton.me)
# @brief Generates roofline plots based on measurement results

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D
import numpy as np

mpl.rcParams['font.family'] = 'Nimbus Sans'

# set paths
res_dir = "results"
measure_dir = res_dir + "/measurements"
roofline_dir = res_dir + "/rooflines"
fontsize = 45
min_dot_size = 250
max_dot_size = 1000

# load peak perf and bandwidths
bandwidth_specs = open(measure_dir + "/bandwidth.txt", "r")
peakperf_specs = open(measure_dir + "/peakperf.txt", "r")
bandwidth = []
perf = []

line = bandwidth_specs.readline()
while line:
    bandwidth.append(float(line.strip()))
    line = bandwidth_specs.readline()

line = peakperf_specs.readline()
while line:
    perf.append(float(line.strip()))
    line = peakperf_specs.readline()

# specifications for plot
bits = ["8_32", "16_32", "64_64_IFMA", "64_64"]
libs = ["APIB", "APIB_Sliced", "Eigen", "MKL"]
shapes = ["square", "thin_1", "thin_2"]
colors = ["crimson", "seagreen", "peru", "navy"]

# creates roofline plots
def create_rooflines(threading, bits, libs, shapes, colors, scale_font=1.0):
    for shape in shapes:
        height = len(bits)
        if (len(bits) > 1):
            height += 1
        if (shape == "slices"):
            height *= 3
        _, axs = plt.subplots(nrows=len(bits), ncols=len(threading), figsize=(15*len(threading), 10*(height)), squeeze=False, constrained_layout=True)
        for bit_ind, bit in enumerate(bits):
            for thread_ind, thread in enumerate(threading):
                run_name = []
                run_perf = []
                run_op_int = []
    
                eff_libs = []
                dim = []
    
                for lib in libs:
                    try:
                        with open(measure_dir + "/" + shape + "/" + lib + "_" + bit + "_" + str(thread) + "_threads" +  ".txt", "r") as runs:
                            run_name_lib = []
                            run_perf_lib = []
                            run_op_int_lib = []
                            dim = []
                            line = runs.readline()
                            while line:
                                dims = line.strip()
                                
                                dim.append(dims)
                                run_name_lib.append(dims)
                                line = runs.readline()
                                run_perf_lib.append(float(line.strip()))
                                line = runs.readline()
                                run_op_int_lib.append(float(line.strip()) / thread)
                                line = runs.readline()
                            run_perf.append(run_perf_lib)
                            run_op_int.append(run_op_int_lib)
                            eff_libs.append(lib)
                            if (shape == "slices"):
                                if (lib == "APIB_Sliced"):
                                    run_name = run_name_lib
                            else:
                                run_name = run_name_lib
    
                    except FileNotFoundError:
                        continue
        
                max_op_int = 0
                for a in run_op_int:
                    max_op_int = max(max_op_int, max(a))

                max_perf = 0
                for a in run_perf:
                    max_perf = max(max_perf, max(a))
        
                thread_perf = perf[bit_ind] * thread
        
                intersect = thread_perf / bandwidth[thread - 1]
    
                axs[bit_ind][thread_ind].tick_params(labelsize=scale_font*fontsize)
                # Move x-tick labels slightly down
                for tick in axs[bit_ind][thread_ind].get_xticklabels():
                    tick.set_y(-0.025)  # Negative moves it down; adjust value as needed
        
                xlim = 2*(max(max_op_int, intersect))
                if (shape == "slices"):
                    axs[bit_ind][thread_ind].set_xlim(left=0.05, right = 2*(max(max_op_int, intersect)))
                else:
                    axs[bit_ind][thread_ind].set_xlim(left=0.1, right = 2*(max(max_op_int, intersect)))
                if (shape=="thin_2"):
                    axs[bit_ind][thread_ind].set_ylim(bottom=0.05, top=1.2*thread_perf)
                else:
                    axs[bit_ind][thread_ind].set_ylim(bottom=0.5, top=1.2*thread_perf)
        
                axs[bit_ind][thread_ind].set_xscale("log")
        
                axs[bit_ind][thread_ind].set_yscale("log")
        
                axs[bit_ind][thread_ind].tick_params(axis="both", which="minor", labelleft=True, labelbottom=True)
                axs[bit_ind][thread_ind].set_facecolor('gainsboro')
                
                if (thread == 1):
                    title_threads = "thread"
                else:
                    title_threads = "threads"

    
                if (bit_ind == 0):
                    axs[bit_ind][thread_ind].set_title(str(thread) + " " + title_threads, fontsize=1.75*scale_font*fontsize, pad=85)
                    if (shape == "slices"):
                        y=1.035
                    else:
                        y=1.075
                    axs[bit_ind][thread_ind].text(
                        x=0,
                        y=y,
                        s="Performance (ops/cycle)",
                        transform=axs[bit_ind][thread_ind].transAxes,
                        ha="left", 
                        va="top",
                        fontsize=1.25*scale_font*fontsize
                    )
                    markerlines = []
                    markertext = []
                    if (thread_ind == 0):
                        for ind, el in enumerate(libs):
                            markertext.append(el.replace("_", " "))
                            markerlines.append(Line2D([0], [0], color=colors[ind], marker="o", markersize=25, linestyle='None'))
                    if (shape=="slices"):
                        bit_min = run_name[0].split(")")[0] + ")" 
                        bit_max = run_name[len(run_name) - 1].split(")")[0] + ")" 
                        markertext.append("b$_A$=" + bit_min)
                        markertext.append("b$_A$=" + bit_max)
                    else:
                        if (shape=="thin_2"):
                            n_min = run_name[0].split(",")[1].split(",")[0]
                            n_max = run_name[len(run_name) - 1].split(",")[1].split(",")[0]
                            markertext.append("m=" + n_min)
                            markertext.append("m=" + n_max)
                        else:
                            n_min = run_name[0].split(",")[0]
                            n_max = run_name[len(run_name) - 1].split(",")[0]
                            markertext.append("n=" + n_min)
                            markertext.append("n=" + n_max)
                    markerlines.append(Line2D([0], [0], color="black", marker="o", markersize=15, linestyle='None'))
                    markerlines.append(Line2D([0], [0], color="black", marker="o", markersize=35, linestyle='None'))
                    axs[bit_ind][thread_ind].legend(markerlines, markertext, fontsize=1*scale_font*fontsize, loc='best')

                if (bit_ind == len(bits) - 1):
                    axs[bit_ind][thread_ind].set_xlabel("Operational Intensity (ops/byte)", fontsize=1.25*scale_font*fontsize)

                if (thread_ind == 0):
                    axs[bit_ind][thread_ind].text(
                        x=-0.15,
                        y=0.5,
                        s=bit.replace("_", "->", 1).replace("_", " "),
                        transform=axs[bit_ind][thread_ind].transAxes,
                        ha="center", 
                        va="center",
                        rotation=90,
                        fontsize=1.75*scale_font*fontsize
                )
        
                axs[bit_ind][thread_ind].grid(which="both", linestyle='solid', color="white", lw=2)
                axs[bit_ind][thread_ind].plot([intersect, xlim], [thread_perf, thread_perf], color="black", label="_nolegend_", lw=8)
                axs[bit_ind][thread_ind].plot([0, intersect], [0, thread_perf], color="black", label="_nolegend_", lw=8)
                
                for lib_ind, lib_el in enumerate(eff_libs):
                    if (lib_el == "APIB" or lib_el == "APIB_Sliced"):
                        scalar = 1.5
                    else:
                        scalar = 1

                    axs[bit_ind][thread_ind].plot(run_op_int[lib_ind], run_perf[lib_ind], color=colors[lib_ind], lw=scalar*5)
                    if (len(run_perf[lib_ind]) > 1):
                        n_sizes = np.array([i for i in range(len(run_perf[lib_ind]))], dtype=float)
                        min_val = np.min(n_sizes)
                        max_val = np.max(n_sizes)
                        scaled_sizes = (n_sizes - min_val) / (max_val - min_val)
                        marker_sizes = [(min_dot_size + el*(max_dot_size - min_dot_size)) for el in scaled_sizes]
                    else:
                        marker_sizes = [max_dot_size]

                    axs[bit_ind][thread_ind].scatter(run_op_int[lib_ind], run_perf[lib_ind], color=colors[lib_ind], s=marker_sizes, zorder=5)

        plt.savefig(roofline_dir + "/" + str(min(threading)) + "_" + str(max(threading)) + "_threads_" + shape + ".svg")
        plt.close()

# primary results
create_rooflines([1, 12], bits, libs, shapes, colors)
create_rooflines([1, 12], ["8_32"], ["APIB_Sliced", "Eigen", "MKL"], ["slices"], colors[1:])
# additional results for appendix
create_rooflines([4, 8], bits, libs, shapes, colors)
create_rooflines([4, 8], ["8_32"], ["APIB_Sliced", "Eigen", "MKL"], ["slices"], colors[1:])
# LUT benchmark
create_rooflines([1], ["8_32"], ["APIB", "LUT"], ["lut"], colors, scale_font=0.5)
