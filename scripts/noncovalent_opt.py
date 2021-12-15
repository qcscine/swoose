import os
import sys
import time
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.optimize import differential_evolution, OptimizeResult
import scine_utilities as utils

"""
This script optimizes the global non-covalent parameters of SFAM, namely
a1, s8, a2, beta, and the atomic charges scaling factor.
The order of these parameters will be kept constant throughout this script.
Some settings to set (usually only the first one needs to be modified):
"""
# The bounds for the noncovalent parameters in order that they appear in the SFAM parameter file
limits_for_noncovalent_parameters = ((0.0, 1.0), (2.0, 6.0), (4.0, 8.0), (5.0, 12.0), (0.5, 1.5))
# The tolerance of the optimization algorithm for each parameter
tolerance = 0.1
# Directory where the reference data is stored
data_dir = "data_for_noncovalent_parametrization"
# Directory to which the figures that are created by this script are stored
figures_dir = "figures"
# Filename for all of the SFAM parameter files
param_filename = "Parameters.dat"
# Indices of the structures to ignore from the trajectories to get less (well-distanced) data points
structs_to_ignore = {9, 11, 13, 15, 17, 18, 20, 21, 23, 24, 25, 27, 28, 29}
# Number of structures per trajectory
n_structures_per_system = 30
# Settings applied to the SFAM molecular mechanics calculation
mm_settings = {"mm_parameter_file": "parameters_tmp.dat",
               "mm_connectivity_file": "Connectivity.dat",
               "non_covalent_cutoff": 100.0}
"""
End of settings. Code starts.
"""


def get_reference_data() -> dict:
    ref_energy_files = glob(os.path.join(data_dir, "*", "*act.dat"))
    ref_data = {}
    for filename in ref_energy_files:
        energies = []
        indices = []
        with open(filename) as file:
            for line in file:
                idx = int(float(line.split()[0]))
                if idx in structs_to_ignore:
                    continue
                energy = float(line.split()[1])
                energies.append(energy)
                indices.append(idx)
        last_energy = energies[-1]
        energies = [utils.KCALPERMOL_PER_HARTREE * (energy - last_energy) for energy in energies]
        system_name = os.path.normpath(filename).split(os.path.sep)[1]
        ref_data[system_name] = {"energies": energies}
        ref_data[system_name].update({"indices": indices})
    return ref_data


def plot_ref_curve(ref_data: dict, system: str, spline: InterpolatedUnivariateSpline) -> None:
    system_data = ref_data[system]
    fig = plt.figure()
    x = np.linspace(1, 30, 1000)
    plt.plot(system_data["indices"], system_data["energies"], 'ro', ms=5)
    plt.plot(x, spline(x))
    plt.xlabel("structure number")
    plt.ylabel("energy (kcal/mol)")
    try:
        os.makedirs(os.path.join(figures_dir, "ref_figs"))
    except FileExistsError:
        pass
    fig.savefig(os.path.join(figures_dir, "ref_figs", system + "_ref.pdf"))
    plt.close(fig)


def plot_comparison(system: str, ref_spline: InterpolatedUnivariateSpline, mm_spline: InterpolatedUnivariateSpline) -> None:
    x = np.linspace(1, n_structures_per_system, 1000)
    fig = plt.figure()
    plt.plot(x, ref_spline(x), color="orange", label="reference")
    plt.plot(x, mm_spline(x), color="blue", label="SFAM")
    plt.legend()
    plt.xlabel("structure number")
    plt.ylabel("energy (kcal/mol)")
    try:
        os.makedirs(os.path.join(figures_dir, "comparisons"))
    except FileExistsError:
        pass
    fig.savefig(os.path.join(figures_dir, "comparisons", system + "_comparison.pdf"))
    plt.close(fig)


def update_ref_with_min_info(ref_data: dict) -> None:
    for system in ref_data:
        system_data = ref_data[system]
        spline = InterpolatedUnivariateSpline(system_data["indices"], system_data["energies"], k=4)
        plot_ref_curve(ref_data, system, spline)
        roots = spline.derivative().roots()
        min_x = roots[0]
        min_y = spline(min_x).item()
        system_data.update({"min_x": min_x, "min_y": min_y, "spline": spline})


def write_tmp_param_file(original_file: str, params: tuple) -> None:
    with open(original_file) as original, open(mm_settings["mm_parameter_file"], "w") as tmp:
        for line in original:
            if line.startswith("a1"):
                tmp.write("a1  " + str(params[0]) + "\n")
            elif line.startswith("s8"):
                tmp.write("s8  " + str(params[1]) + "\n")
            elif line.startswith("a2"):
                tmp.write("a2  " + str(params[2]) + "\n")
            elif line.startswith("beta"):
                tmp.write("beta  " + str(params[3]) + "\n")
            elif line.startswith("scaling_factor_for_atomic_charges"):
                tmp.write("scaling_factor_for_atomic_charges  " + str(params[4]) + "\n")
            else:
                tmp.write(line)


def cost_function(*args) -> float:
    start_time = time.time()
    global global_counter
    global_counter += 1
    params = args[0]
    print("In cost function, params:", params)
    print("Iteration:", global_counter)
    data = args[1]
    plotting = args[2]
    system_dirs = glob(os.path.join(data_dir, "*"))
    total_cost = 0
    for system_dir in system_dirs:
        system_cost = get_cost_for_system(system_dir, data, params, plotting)
        total_cost += system_cost
    print("Total cost (kcal^2/mol^2):", total_cost)
    print("Time (seconds):", time.time() - start_time)
    print("-------------------------------------")
    return total_cost


def get_cost_for_system(system_dir: str, data: dict, params: tuple, plotting: bool) -> float:
    system_name = os.path.normpath(system_dir).split(os.path.sep)[1]
    base_dir = os.getcwd()
    os.chdir(system_dir)
    write_tmp_param_file(param_filename, params)
    energies = []
    calculator = utils.core.get_calculator("SFAM", "Swoose")
    log = utils.core.Log.silent()
    calculator.log = log
    calculator.settings.update(mm_settings)
    xyz_file = glob("*_1.xyz")[0]
    calculator.structure = utils.io.read(xyz_file)[0]
    for i in data[system_name]["indices"]:
        energy = calculate_mm_energy_for_system(i, calculator)
        energies.append(energy)
    os.remove(mm_settings["mm_parameter_file"])
    os.chdir(base_dir)

    last_energy = energies[-1]
    energies = [energy - last_energy for energy in energies]

    spline = InterpolatedUnivariateSpline(data[system_name]["indices"], energies, k=4)
    if plotting:
        plot_comparison(system_name, data[system_name]["spline"], spline)
    roots = spline.derivative().roots()
    try:
        min_x = roots[0]
        min_y = spline(min_x).item()
    except IndexError:
        min_x, min_y = n_structures_per_system, 0.0

    diff_min_x = min_x - data[system_name]["min_x"]
    diff_min_y = min_y - data[system_name]["min_y"]
    mean_squared_error = 0.0
    for e1, e2 in zip(energies, data[system_name]["energies"]):
        mean_squared_error += (e1 - e2)**2
    mean_squared_error /= len(energies)

    return diff_min_x**2 + diff_min_y**2 + mean_squared_error


def calculate_mm_energy_for_system(structure_idx: int, calculator: utils.core.Calculator) -> float:
    allxyz_file = glob("*.allxyz")[0]
    trj = utils.io.read_trajectory(utils.io.TrajectoryFormat.Xyz, allxyz_file)
    calculator.positions = trj[structure_idx - 1]
    results = calculator.calculate()
    return results.energy * utils.KCALPERMOL_PER_HARTREE


def optimize(ref_data: dict) -> OptimizeResult:
    opt_result = differential_evolution(cost_function, limits_for_noncovalent_parameters,
                                        args=(ref_data, False), tol=tolerance)
    return opt_result


if __name__ == "__main__":
    global_counter = 0
    ref = get_reference_data()
    update_ref_with_min_info(ref)
    t1 = time.time()
    result = optimize(ref)
    t2 = time.time()
    print("\n\nOPTIMIZATION DONE!\n")
    print("Optimized variables:", result.x)
    print("RMSE:", np.sqrt(result.fun / (len(ref) * 3)), "kcal/mol")
    print("Time for optimization:", round(t2 - t1, 2), "seconds.")
    print("Creating more figures...")
    emptiness = open(os.devnull, "w")
    sys.stdout = emptiness
    cost_function(result.x, ref, True)
    sys.stdout = sys.__stdout__
    print("Done.")
