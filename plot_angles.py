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
from scipy.special import gammaln, sph_harm
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
parser.add_argument("-c", "--corrected", action='store_true', help="plot acceptance-corrected data")
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
amps_sorted = []
waves_reflectivity = []
waves_l = []
waves_m = []
waves_r = []
waves_s = []
for wave in waves_sorted:
    wave_pos = wave + "+"
    if wave_pos in amplitudes:
        if wave_dict[wave[0]] == 0:
            waves_l += [0, 0]
            waves_m += [0, 0]
        else:
            waves_l += [wave_dict[wave[0]], wave_dict[wave[0]]]
            waves_m += [(1 if wave[2] == '+' else -1) * int(wave[1]), (1 if wave[2] == '+' else -1) * int(wave[1])]
        waves_r += [1, -1]
        waves_s += [1, -1]
        amps_sorted += [wave_pos, wave_pos]
    wave_neg = wave + "-"
    if wave_neg in amplitudes:
        if wave_dict[wave[0]] == 0:
            waves_l += [0, 0]
            waves_m += [0, 0]
        else:
            waves_l += [wave_dict[wave[0]], wave_dict[wave[0]]]
            waves_m += [(1 if wave[2] == '+' else -1) * int(wave[1]), (1 if wave[2] == '+' else -1) * int(wave[1])]
        waves_r += [1, -1]
        waves_s += [-1, 1]
        amps_sorted += [wave_neg, wave_neg]

print(amps_sorted)
print(waves_l) # list of tuples, each one has two elements which describe the two
print(waves_m) # places where the amplitude shows up in two sums
print(waves_r)
print(waves_s)

def get_angles(tag="DATA", pol=""):
    # Get invariant mass from data files
    bin_dirs = [bin_dir for bin_dir in fit_results.parent.iterdir() if bin_dir.is_dir()]
    phis = []
    costhetas = []
    weights = []
    print("Gathering data for histograms")
    for i, bin_dir in tqdm(enumerate(bin_dirs)):
        bin_phis = np.array([])
        bin_costhetas = np.array([])
        bin_weights = np.array([])
        data_files = [data_file for data_file in bin_dir.iterdir() if data_file.name.endswith(".root") and tag in data_file.name and pol in data_file.name]
        for data_file in data_files:
            with uproot.open(data_file) as df:
                branches = df['kin'].arrays()
                recoil_lab = vector.array({"E": branches['E_FinalState'][:,0],
                                           "px": branches['Px_FinalState'][:,0],
                                           "py": branches['Py_FinalState'][:,0],
                                           "pz": branches['Pz_FinalState'][:,0]})
                p1_lab = vector.array({"E": branches['E_FinalState'][:,1],
                                       "px": branches['Px_FinalState'][:,1],
                                       "py": branches['Py_FinalState'][:,1],
                                       "pz": branches['Pz_FinalState'][:,1]})
                p2_lab = vector.array({"E": branches['E_FinalState'][:,2],
                                       "px": branches['Px_FinalState'][:,2],
                                       "py": branches['Py_FinalState'][:,2],
                                       "pz": branches['Pz_FinalState'][:,2]})
                beam_lab = vector.array({"E": branches['E_Beam'],
                                         "px": branches['Px_Beam'],
                                         "py": branches['Py_Beam'],
                                         "pz": branches['Pz_Beam']})
                # boost to COM
                com = recoil_lab + p1_lab + p2_lab
                com_boost_vector = com.to_beta3()
                recoil = recoil_lab.boost(-com_boost_vector)
                p1 = p1_lab.boost(-com_boost_vector)
                p2 = p2_lab.boost(-com_boost_vector)
                beam = beam_lab.boost(-com_boost_vector)

                # boost to resonance
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
                if 'Weight' in branches.keys():
                    bin_weights = np.append(bin_weights, list(branches['Weight']))
                else:
                    bin_weights = np.append(bin_weights, np.ones_like(costheta))
        costhetas.append(bin_costhetas)
        phis.append(bin_phis)
        weights.append(bin_weights)
    return costhetas, phis, weights

ct_data, phi_data, weights_data = get_angles(pol="PARA_0")
ct_gen, phi_gen, weights_gen = get_angles("GEN", pol="PARA_0")
ct_acc, phi_acc, weights_acc = get_angles("ACCEPT", pol="PARA_0")


plt.rcParams["figure.figsize"] = (40, 10)
plt.rcParams["font.size"] = 24
for bin_n in range(len(bin_df)):
    fig, (ax, ax2) = plt.subplots(ncols=2)
    nbins = 20
    weight_total = np.zeros_like(ct_acc[bin_n])
    for event in range(len(ct_acc[bin_n])):
        weight_positivere = 0
        weight_positiveim = 0
        weight_negativere = 0
        weight_negativeim = 0
        for i, amp in enumerate(amps_sorted):
            costheta = ct_acc[bin_n]
            phi = phi_acc[bin_n]
            fit_amplitude = df_filtered[amp + '_re'].iloc[bin_n] + 1j * df_filtered[amp + '_im'].iloc[bin_n]
            if waves_r[i] == 1:
                if waves_s[i] == 1:
                    # scipy does spherical harmonics as ylm(m, l, phi, theta) for some silly reason
                    zlm = np.real(sph_harm(waves_m[i], waves_l[i], phi[event], costheta[event]))
                    factor = np.sqrt(1 + waves_s[i] * 0.33) # Pgamma = 0.3 for now
                    weight_positivere += factor * zlm * fit_amplitude
                else:
                    zlm = np.real(sph_harm(waves_m[i], waves_l[i], phi[event], costheta[event]))
                    factor = np.sqrt(1 + waves_s[i] * 0.33) # Pgamma = 0.3 for now
                    weight_negativere += factor * zlm * fit_amplitude
            else:
                if waves_s[i] == 1:
                    zlm = np.imag(sph_harm(waves_m[i], waves_l[i], phi[event], costheta[event]))
                    factor = np.sqrt(1 + waves_s[i] * 0.33) # Pgamma = 0.3 for now
                    weight_negativeim += factor * zlm * fit_amplitude
                else:
                    zlm = np.imag(sph_harm(waves_m[i], waves_l[i], phi[event], costheta[event]))
                    factor = np.sqrt(1 + waves_s[i] * 0.33) # Pgamma = 0.3 for now
                    weight_positiveim += factor * zlm * fit_amplitude
        weight_total[event] = np.abs(weight_positivere)**2 + np.abs(weight_negativere)**2 + np.abs(weight_positiveim)**2 + np.abs(weight_negativeim)**2
        weight_total[event] = weight_total[event] / len(ct_gen[bin_n])
    fit_hist, bin_edges = np.histogram(ct_acc[bin_n], bins=nbins, range=(-1, 1), weights=weight_total)
    #ax.bar(bin_edges[:-1], fit_hist, width=np.diff(bin_edges), align="edge") 
    if args.corrected:
        ct_gen_hist, _ = np.histogram(ct_gen[bin_n], bins=nbins, range=(-1, 1), weights=weights_gen[bin_n])
        phi_gen_hist, _ = np.histogram(phi_gen[bin_n], bins=nbins, range=(-np.pi, np.pi), weights=weights_gen[bin_n])
        ct_acc_hist, _ = np.histogram(ct_acc[bin_n], bins=nbins, range=(-1, 1), weights=weights_acc[bin_n])
        phi_acc_hist, _ = np.histogram(phi_acc[bin_n], bins=nbins, range=(-np.pi, np.pi), weights=weights_acc[bin_n])
        ct_acceptance_function = ct_gen_hist / ct_acc_hist
        phi_acceptance_function = phi_gen_hist / phi_acc_hist
        ct_data_hist, bin_edges = np.histogram(ct_data[bin_n], range=(-1, 1), bins=nbins, weights=weights_data[bin_n])
        phi_data_hist, phi_bin_edges = np.histogram(phi_data[bin_n], range=(-np.pi, np.pi), bins=nbins, weights=weights_data[bin_n])
        ct_data_corrected = ct_data_hist * ct_acceptance_function
        phi_data_corrected = phi_data_hist * phi_acceptance_function
        bin_centers = bin_edges[:-1] + np.diff(bin_edges) / 2
        phi_bin_centers = phi_bin_edges[:-1] + np.diff(phi_bin_edges) / 2
        ax.errorbar(bin_centers, ct_data_corrected, yerr=np.sqrt(ct_data_corrected), xerr=np.diff(bin_edges) / 2, color='k', fmt='none')
        ax2.errorbar(phi_bin_centers, phi_data_corrected, yerr=np.sqrt(phi_data_corrected), xerr=np.diff(phi_bin_edges) / 2, color='k', fmt='none')
    else:
        ct_gen_hist, _ = np.histogram(ct_gen[bin_n], bins=nbins, range=(-1, 1), weights=weights_gen[bin_n])
        phi_gen_hist, _ = np.histogram(phi_gen[bin_n], bins=nbins, range=(-np.pi, np.pi), weights=weights_gen[bin_n])
        ct_acc_hist, _ = np.histogram(ct_acc[bin_n], bins=nbins, range=(-1, 1), weights=weights_acc[bin_n])
        phi_acc_hist, _ = np.histogram(phi_acc[bin_n], bins=nbins, range=(-np.pi, np.pi), weights=weights_acc[bin_n])
        ct_data_hist, bin_edges = np.histogram(ct_data[bin_n], range=(-1, 1), bins=nbins, weights=weights_data[bin_n])
        phi_data_hist, phi_bin_edges = np.histogram(phi_data[bin_n], range=(-np.pi, np.pi), bins=nbins, weights=weights_data[bin_n])
        bin_centers = bin_edges[:-1] + np.diff(bin_edges) / 2
        phi_bin_centers = phi_bin_edges[:-1] + np.diff(phi_bin_edges) / 2
        ax.errorbar(bin_centers, ct_data_hist, yerr=np.sqrt(ct_data_hist), xerr=np.diff(bin_edges) / 2, color='k', fmt='none')
        ax2.errorbar(phi_bin_centers, phi_data_hist, yerr=np.sqrt(phi_data_hist), xerr=np.diff(phi_bin_edges) / 2, color='k', fmt='none')
        #ax.errorbar(bin_centers, ct_acc_hist, yerr=np.sqrt(ct_acc_hist), xerr=np.diff(bin_edges) / 2, color='b', fmt='none')
        #ax2.errorbar(phi_bin_centers, phi_acc_hist, yerr=np.sqrt(phi_acc_hist), xerr=np.diff(phi_bin_edges) / 2, color='b', fmt='none')
        #ax.errorbar(bin_centers, ct_gen_hist, yerr=np.sqrt(ct_gen_hist), xerr=np.diff(bin_edges) / 2, color='g', fmt='none')
        #ax2.errorbar(phi_bin_centers, phi_gen_hist, yerr=np.sqrt(phi_gen_hist), xerr=np.diff(phi_bin_edges) / 2, color='g', fmt='none')
    ax.set_ylim(0)
    ax2.set_ylim(0)
    plt.title(f"Bin {bin_n}")
    ax.set_xlabel(r"$\cos(\theta_{HX})$")
    ax2.set_xlabel(r"$\phi_{HX}$")
    pdf.savefig(dpi=300)
    plt.close()
pdf.close()
