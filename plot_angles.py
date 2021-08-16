#!/usr/bin/python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import sys
from pathlib import Path
import uproot
import awkward as ak
import argparse
from tqdm import tqdm
from scipy.special import gammaln
from simple_term_menu import TerminalMenu
import vector

def plot_menu(root):
    plot_menu_title = "Select a fit_results.txt file:"
    fits = [p for p in root.glob("*::fit_results.txt")]
    plot_menu_items = ["Cancel"] + [fit.name[:-17] for fit in fits]
    plot_menu_cursor = "> "
    plot_menu_cursor_style = ("fg_red", "bold")
    plot_menu_style = ("bg_black", "fg_green")
    plot_menu = TerminalMenu(
        menu_entries=plot_menu_items,
        title=plot_menu_title,
        menu_cursor=plot_menu_cursor,
        menu_cursor_style=plot_menu_cursor_style,
        menu_highlight_style=plot_menu_style
    )
    selection_index = plot_menu.show()
    if selection_index == 0:
        print("No plot file chosen")
        sys.exit(1)
    else:
        return fits[selection_index - 1]


parser = argparse.ArgumentParser(description="Plotting tools for PartialWaveAnalysis")
parser.add_argument("path", help="path to fit directory")
parser.add_argument("-u", "--uncorrected", action='store_true', help="plot intensities without acceptance correction")
parser.add_argument("-d", "--data", action='store_true', help="plot true data histogram behind intensities")
parser.add_argument("-b", "--background", action='store_true', help="plot background histogram behind intensities")
parser.add_argument("-s", "--subtracted", action='store_true', help="plot background-subtracted data behind intensities")
parser.add_argument("-l", "--label", default="K_SK_S", help="LaTeX formated string of particles for x-axis label (default: \"K_SK_S\")")
parser_group = parser.add_mutually_exclusive_group()
parser_group.add_argument("--bins", type=int, help="set number of bins for data")
parser_group.add_argument("--optimize", action='store_true', help="optimize data binning (experimental, currently doesn't work)")

if len(sys.argv) == 1: # if the user doesn't supply any arguments, print the help string and exit
    parser.print_help(sys.stderr)
    sys.exit(1)
args = parser.parse_args()
data_colors = ['blue', 'red', 'green']
colors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown", "tab:pink", "tab:olive", "tab:cyan"]

fit_results = plot_menu(Path(args.path)).resolve()
xlabel = rf"m$({args.label})$ GeV/$c^2$"

pdf = matplotlib.backends.backend_pdf.PdfPages(f"angles_{fit_results.stem.split('::')[0]}.pdf")

df = pd.read_csv(fit_results, delimiter='\t', index_col=False)
df.sort_values(['Bin', 'likelihood', 'total_intensity_err_AC'], ascending=[True, False, True], inplace=True)

def mask_first(x):
    result = np.zeros_like(x)
    result[0] = 1
    return result

mask = df.groupby(['Bin'])['Bin'].transform(mask_first).astype(bool)
df_filtered = df.loc[mask]
print(df_filtered.columns)
print(df_filtered.head())
bin_df = pd.read_csv(fit_results.parent / 'bin_info.txt', delimiter='\t')


amplitudes = [column[:-3] for column in df.columns.to_list() if column.endswith('_re')]
wave_set = set([amp[:-1] for amp in amplitudes])
wave_dict = {'S': 0, 'P': 1, 'D': 2, 'F': 3, 'G': 4}
waves_sorted = sorted(list(wave_set), key=lambda wave: 100 * wave_dict[wave[0]] + (-1 if wave[-1] == '-' else 1) * int(wave[1]))
# kind of sneaky wave of sorting waves by L-letter and M-number without creating a whole new class
waves_l = []
waves_m = []
waves_r = []
waves_s = []
for wave in waves_sorted:
    if wave_dict[wave[0]] == 0:
        waves_l += [0, 0]
        waves_m += [0, 0]
    else:
        waves_l += [wave_dict[wave[0]], wave_dict[wave[0]]]
        waves_m += [(1 if wave[2] == '+' else -1) * int(wave[1]), (1 if wave[2] == '+' else -1) * int(wave[1])]
    waves_r += [1, -1]
    reflectivity = 1 if wave[-1] == '+' else -1
    waves_s += [reflectivity * 1, reflectivity * -1]

print(waves_l)
print(waves_m)
print(waves_r)
print(waves_s)

# Get invariant mass from data files
bin_dirs = [bin_dir for bin_dir in fit_results.parent.iterdir() if bin_dir.is_dir()]
phis = []
costhetas = []
weights = []
print("Gathering real data for histograms")
for i, bin_dir in tqdm(enumerate(bin_dirs)):
    if args.data or args.background or args.subtracted:
        bin_phis = np.array([])
        bin_costhetas = np.array([])
        bin_weights = np.array([])
        data_files = [data_file for data_file in bin_dir.iterdir() if data_file.name.endswith(".root") and "DATA" in data_file.name]
        for data_file in data_files:
            with uproot.open(data_file) as df:
                branches = df['kin'].arrays()
                recoil = vector.array({"E": branches['E_FinalState'][:,0],
                                       "px": branches['Px_FinalState'][:,0],
                                       "py": branches['Py_FinalState'][:,0],
                                       "pz": branches['Pz_FinalState'][:,0]})
                p1 = vector.array({"E": branches['E_FinalState'][:,1],
                                   "px": branches['Px_FinalState'][:,1],
                                   "py": branches['Py_FinalState'][:,1],
                                   "pz": branches['Pz_FinalState'][:,1]})
                p2 = vector.array({"E": branches['E_FinalState'][:,2],
                                   "px": branches['Px_FinalState'][:,2],
                                   "py": branches['Py_FinalState'][:,2],
                                   "pz": branches['Pz_FinalState'][:,2]})
                beam = vector.array({"E": branches['E_Beam'],
                                     "px": branches['Px_Beam'],
                                     "py": branches['Py_Beam'],
                                     "pz": branches['Pz_Beam']})
                resonance = p1 + p2
                resonance_boost_vector = resonance.to_beta3()
                beam_res = beam.boost(-resonance_boost_vector)
                p1_res = p1.boost(-resonance_boost_vector)
                recoil_res = recoil.boost(-resonance_boost_vector)

                beam_vect = vector.array({"px": beam.x, "py": beam.y, "pz": beam.z})
                p1_res_vect = vector.array({"px": p1_res.x, "py": p1_res.y, "pz": p1_res.z})
                recoil_vect = vector.array({"px": recoil.x, "py": recoil.y, "pz": recoil.z})
                recoil_res_vect = vector.array({"px": recoil_res.x, "py": recoil_res.y, "pz": recoil_res.z})

                z = -1 *  recoil_res_vect.unit() # Helicity frame
                y = beam_vect.unit().cross(-recoil_vect.unit()).unit()
                x = y.cross(z)
                
                angles = vector.array({"x": p1_res_vect.dot(x), "y": p1_res_vect.dot(y), "z": p1_res_vect.dot(z)})
                costheta = angles.costheta
                phi = angles.phi

                bin_costhetas = np.append(bin_costhetas, costheta)
                bin_phis = np.append(bin_phis, phi)
                bin_weights = np.append(bin_weights, list(branches['Weight']))
        costhetas.append(bin_costhetas)
        phis.append(bin_phis)
        weights.append(bin_weights)
if args.data:
    print(f"Data has {len(phis)} events")

for bin_n in range(len(bin_df)):
    fig = plt.figure()
    plt.hist(costhetas[bin_n], range=(-1, 1), bins=20)
    plt.title(f"Bin {bin_n}")
    plt.xlabel(r"$\cos(\theta_{GJ})$")
    pdf.savefig(dpi=300)
    plt.close()
pdf.close()
