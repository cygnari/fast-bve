import os
import math

file = open('namelist.txt', 'r')

use_fast = False
use_amr = False
use_remesh = False
force_conservative = False
vor_fix = False

for line in file:
    string_break = line.rfind('=')
    part1 = line[:string_break]
    part2 = line[string_break+1:]
    if (part1 == 'use_amr'):
        if (part2.strip() == '1'):
            use_amr = True
            use_remesh = True
    if (part1 == 'use_remesh'):
        if (part2.strip() == '1'):
            use_remesh = True
    if (part1 == 'use_fast'):
        if (part2.strip() == '1'):
            use_fast = True
    if (part1 == 'force_conservative'):
        if (part2.strip() == '1'):
            force_conservative = True
    if (part1 == 'vor_fix'):
        if (int(part2) >= 1):
            vor_fix = True
        if (int(part2) >= 2):
            vor_limiter = True
    if (part1 == 'out_path'):
        outpath = part2.strip()
    if (part1 == 'end_time'):
        end_time = float(part2)
    if (part1 == 'dynamics_levels_min'):
        dynamics_levels_min = int(part2)
    if (part1 == 'initial_vor_condition'):
        initial_vor_condition = part2.strip()
    if (part1 == 'vor_force'):
        vor_force = part2.strip()
    if (part1 == 'icp1'):
        icp1 = int(part2)
    if (part1 == 'icp2'):
        icp2 = float(part2)
    if (part1 == 'frp1'):
        frp1 = int(part2)
    if (part1 == 'frp2'):
        frp2 = float(part2)
    if (part1 == 'fast_sum_theta'):
        fast_sum_theta = float(part2)
    if (part1 == 'interp_degree'):
        interp_degree = int(part2)
    if (part1 == 'fast_sum_tree_levels'):
        # fast_sum_tree_levels = int(part2)
        if (int(part2) == -1):
            fast_sum_tree_levels = dynamics_levels_min - 3
        else:
            fast_sum_tree_levels = int(part2)
    if (part1 == 'fast_sum_theta'):
        fast_sum_theta = float(part2)
    if (part1 == 'amr_levels_max'):
        amr_levels_max = int(part2)

initial_particles = pow(4, dynamics_levels_min - 1) * 10 + 2
output_filename = str(initial_particles) + "_" + initial_vor_condition + "_"
# output_filename += "_"
if (icp1 > 0):
    output_filename += str(icp1) + "_"
if (icp2 > 0):
    precision = max(int(math.ceil(-math.log10(icp2))), 0);
    output_filename += str('{:.{n}f}'.format(icp2, n=precision)) + '_'
if (vor_force != "none"):
    output_filename += vor_force + "_"
    if (frp1 > 0):
        output_filename += str(frp1) + "_"
    if (frp2 > 0):
        precision = max(int(math.ceil(-math.log10(frp2))), 0);
        output_filename += str('{:.{n}f}'.format(frp2, n=precision)) + '_'
if (use_fast):
    output_filename += "fast_" + str(fast_sum_tree_levels) + "_" + '{:.{n}f}'.format(fast_sum_theta, n=1) + "_"
else:
    output_filename += "direct_"
if (use_amr):
    output_filename += "amr_" + str(amr_levels_max) + "_"
if (use_remesh):
    output_filename += "remesh_"
if (force_conservative):
    output_filename += "fixer_"
precision = max(int(math.ceil(-math.log10(end_time))), 0);
output_filename += '{:.{n}f}'.format(end_time, n=precision)
output_path = outpath + output_filename
isExist = os.path.exists(output_path)
if (not isExist):
    os.makedirs(output_path)
# print(output_filename)
