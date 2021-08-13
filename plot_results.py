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


def optBINS(data, weights, minM, maxM):
    # Kevin H. Knuth's algorithm for determining optimal binning
    # https://arxiv.org/pdf/physics/0605197.pdf
    print("Computing optimal number of bins")
    N = len(data)
    logp = np.zeros(maxM)
    for M in range(minM, maxM):
        print(f"Testing nbins = {M}", end='')
        n, _ = np.histogram(data, weights=weights)
        part1 = N * np.log(M) + gammaln(M / 2) - gammaln(N + M / 2)
        part2 = - M * gammaln(1/2) + np.sum(gammaln(n + 0.5))
        logp[M] = part1 + part2
        print(f" ... ({logp[M]})")
    optM = np.argmax(logp)
    print(f"Optimal number of bins is {optM}")
    return optM


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

pdf = matplotlib.backends.backend_pdf.PdfPages(f"figs_{fit_results.stem.split('::')[0]}.pdf")

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
bin_width = bin_df['mass'].iloc[1] - bin_df['mass'].iloc[0]
bin_df['Centers'] = np.linspace(bin_df['mass'].iloc[0] + bin_width/2, bin_df['mass'].iloc[-1] - bin_width/2, len(bin_df))
if not args.uncorrected:
    ac_tag = "_AC_INT"
    ac_tag_total = "_AC"
    amplitudes = [column[:-len(ac_tag)] for column in df.columns.to_list() if column.endswith(ac_tag) and not "_err" in column]
else:
    ac_tag = "_INT"
    ac_tag_total = ""
    amplitudes = [column[:-len(ac_tag)] for column in df.columns.to_list() if column.endswith(ac_tag) and not "_err" in column and not "_AC_" in column]

print(f"Plotting Amplitudes: {amplitudes}")

wave_set = set([amp[:-1] for amp in amplitudes])
wave_dict = {'S': 0, 'P': 1, 'D': 2, 'F': 3, 'G': 4}
waves_sorted = sorted(list(wave_set), key=lambda wave: 100 * wave_dict[wave[0]] + (-1 if wave[-1] == '-' else 1) * int(wave[1]))
# kind of sneaky wave of sorting waves by L-letter and M-number without creating a whole new class

# Get invariant mass from data files
bin_dirs = [bin_dir for bin_dir in fit_results.parent.iterdir() if bin_dir.is_dir()]
M2s = []
M2s_bkg = []
M2s_sub = []
weights = []
weights_bkg = []
weights_sub = []
print("Gathering real data for histograms")
for bin_dir in tqdm(bin_dirs):
    if args.data or args.background or args.subtracted:
        data_files = [data_file for data_file in bin_dir.iterdir() if data_file.name.endswith(".root") and "DATA" in data_file.name]
        for data_file in data_files:
            with uproot.open(data_file) as df:
                branches = df['kin'].arrays()
                Final_State_P4 = [np.array([E, Px, Py, Pz]) for E, Px, Py, Pz in zip(np.sum(branches['E_FinalState'][:,1:], axis=1), np.sum(branches['Px_FinalState'][:,1:], axis=1), np.sum(branches['Py_FinalState'][:,1:], axis=1), np.sum(branches['Pz_FinalState'][:,1:], axis=1))] #ignore first particle, which is the proton
                M2s += [P4[0]**2 - np.sum(np.power(P4[1:], 2)) for P4 in Final_State_P4]
                M2s_sub += [P4[0]**2 - np.sum(np.power(P4[1:], 2)) for P4 in Final_State_P4]
            weights += list(branches['Weight'])
            weights_sub += list(branches['Weight'])
    if args.background or args.subtracted:
        bkg_files = [bkg_file for bkg_file in bin_dir.iterdir() if bkg_file.name.endswith(".root") and "BKG" in bkg_file.name]
        for bkg_file in bkg_files:
            with uproot.open(bkg_file) as df:
                branches = df['kin'].arrays()
                Final_State_P4 = [np.array([E, Px, Py, Pz]) for E, Px, Py, Pz in zip(np.sum(branches['E_FinalState'][:,1:], axis=1), np.sum(branches['Px_FinalState'][:,1:], axis=1), np.sum(branches['Py_FinalState'][:,1:], axis=1), np.sum(branches['Pz_FinalState'][:,1:], axis=1))] #ignore first particle, which is the proton
                M2s_bkg += [P4[0]**2 - np.sum(np.power(P4[1:], 2)) for P4 in Final_State_P4]
                M2s_sub += [P4[0]**2 - np.sum(np.power(P4[1:], 2)) for P4 in Final_State_P4]
            weights_bkg += list(branches['Weight'])
            weights_sub += list(np.array(branches['Weight']) * -1)
invariant_mass = np.sqrt(M2s)
invariant_mass_bkg = np.sqrt(M2s_bkg)
invariant_mass_sub = np.sqrt(M2s_sub)
if args.data:
    print(f"Data has {len(invariant_mass)} events")
if args.background or args.subtracted:
    print(f"Background has {len(invariant_mass_bkg)} events")

nbins = len(bin_df)
if args.optimize:
    print("WARNING: Plotting data with an optimized binning algorithm. This is an experimental feature!")
    nbins = optBINS(invariant_mass, weights, nbins, 600)
elif not args.bins == None:
    nbins = args.bins

plt.rcParams["figure.figsize"] = (20, 10)
plt.rcParams["font.size"] = 24

print("Plotting Separate Amplitudes")
for wave in waves_sorted:
    fig = plt.figure()
    print(wave, end='\t')

    # plot real data histograms
    if args.data:
        plt.hist(invariant_mass, bins=nbins, range=(bin_df['mass'].iloc[0], bin_df['mass'].iloc[-1]), weights=weights, color=data_colors[0], label="Data", fill=False, histtype='step', lw=2)
    if args.background:
        plt.hist(invariant_mass_bkg, bins=nbins, range=(bin_df['mass'].iloc[0], bin_df['mass'].iloc[-1]), weights=weights_bkg, color=data_colors[1], label="Background", fill=False, histtype='step', lw=2)
    if args.subtracted:
        plt.hist(invariant_mass_sub, bins=nbins, range=(bin_df['mass'].iloc[0], bin_df['mass'].iloc[-1]), weights=weights_sub, color=data_colors[2], label="Data - Background", fill=False, histtype='step', lw=2)

    # wave = S0, P1+, D2-, etc. -- no reflectivity info
    amp_letter = wave[0]
    amp_m = wave[1]
    if int(amp_m) > 0:
        amp_m_sign = wave[2]
    else:
        amp_m_sign = ""

    wave_pos = wave + "+"
    if wave_pos in amplitudes:
        print("+e", end='\t')
        plt.errorbar(bin_df['Centers'].iloc[df_filtered['Bin']], df_filtered[wave_pos + ac_tag], yerr=df_filtered[wave_pos + "_err" + ac_tag], elinewidth=0.5, fmt='o', color='r', label=r'$+\epsilon$')
    else:
        print("", end='\t')

    wave_neg = wave + "-"
    if wave_neg in amplitudes:
        print("-e", end='\t')
        plt.errorbar(bin_df['Centers'].iloc[df_filtered['Bin']], df_filtered[wave_neg + ac_tag], yerr=df_filtered[wave_neg + "_err" + ac_tag], elinewidth=0.5, fmt='o', color='b', label=r'$-\epsilon$')
    else:
        print("", end='\t')
    # Plot total
    print("tot")
    plt.errorbar(bin_df['Centers'].iloc[df_filtered['Bin']], df_filtered['total_intensity' + ac_tag_total], yerr=df_filtered['total_intensity_err' + ac_tag_total], elinewidth=0.5, fmt='none', color='k')
    plt.hist(bin_df['Centers'].iloc[df_filtered['Bin']], bins=len(bin_df), range=(bin_df['mass'].iloc[0], bin_df['mass'].iloc[-1]), weights=df_filtered['total_intensity' + ac_tag_total], fill=False, histtype='step', color='k', label="Total") # this is a sneaky hack
    plt.title(rf"${amp_letter}_{{{amp_m_sign}{amp_m}}}$")
    plt.ylim(bottom=0)
    plt.xlim(bin_df['Centers'].iloc[0] - 0.1, bin_df['Centers'].iloc[-1] + 0.1)
    plt.ylabel("Intensity")
    plt.xlabel(xlabel)
    plt.legend()
    plt.tight_layout()
    pdf.savefig(fig, dpi=300)


print("Plotting Combined Amplitudes")
print("Positive Reflectivity")
fig = plt.figure()

# plot real data histograms
if args.data:
    plt.hist(invariant_mass, bins=nbins, range=(bin_df['mass'].iloc[0], bin_df['mass'].iloc[-1]), weights=weights, color=data_colors[0], label="Data", fill=False, histtype='step', lw=2)
if args.background:
    plt.hist(invariant_mass_bkg, bins=nbins, range=(bin_df['mass'].iloc[0], bin_df['mass'].iloc[-1]), weights=weights_bkg, color=data_colors[1], label="Background", fill=False, histtype='step', lw=2)
if args.subtracted:
    plt.hist(invariant_mass_sub, bins=nbins, range=(bin_df['mass'].iloc[0], bin_df['mass'].iloc[-1]), weights=weights_sub, color=data_colors[2], label="Data - Background", fill=False, histtype='step', lw=2)

for i, wave in enumerate(waves_sorted):
    amp_letter = wave[0]
    amp_m = wave[1]
    if int(amp_m) > 0:
        amp_m_sign = wave[2]
    else:
        amp_m_sign = ""

    wave_pos = wave + "+"
    if wave_pos in amplitudes:
        print(wave + "\t+e")
        plt.errorbar(bin_df['Centers'].iloc[df_filtered['Bin']], df_filtered[wave_pos + ac_tag], yerr=df_filtered[wave_pos + "_err" + ac_tag], elinewidth=0.5, fmt='o', color=colors[i], label=rf"${amp_letter}_{{{amp_m_sign}{amp_m}}}$")
print("tot")
plt.errorbar(bin_df['Centers'].iloc[df_filtered['Bin']], df_filtered['total_intensity' + ac_tag_total], yerr=df_filtered['total_intensity_err' + ac_tag_total], elinewidth=0.5, fmt='none', color='k')
plt.hist(bin_df['Centers'].iloc[df_filtered['Bin']], bins=len(bin_df), range=(bin_df['mass'].iloc[0], bin_df['mass'].iloc[-1]), weights=df_filtered['total_intensity' + ac_tag_total], fill=False, histtype='step', color='k', label="Total") # this is a sneaky hack
plt.xlim(bin_df['Centers'].iloc[0] - 0.1, bin_df['Centers'].iloc[-1] + 0.1)
plt.ylim(bottom=0)
plt.legend(loc="upper right")
plt.title("Positive Reflectivity")
plt.ylabel("Intensity")
plt.xlabel(xlabel)
plt.tight_layout()
pdf.savefig(fig, dpi=300)

print("Negative Reflectivity")
fig = plt.figure()

# plot real data histograms
if args.data:
    plt.hist(invariant_mass, bins=nbins, range=(bin_df['mass'].iloc[0], bin_df['mass'].iloc[-1]), weights=weights, color=data_colors[0], label="Data", fill=False, histtype='step', lw=2)
if args.background:
    plt.hist(invariant_mass_bkg, bins=nbins, range=(bin_df['mass'].iloc[0], bin_df['mass'].iloc[-1]), weights=weights_bkg, color=data_colors[1], label="Background", fill=False, histtype='step', lw=2)
if args.subtracted:
    plt.hist(invariant_mass_sub, bins=nbins, range=(bin_df['mass'].iloc[0], bin_df['mass'].iloc[-1]), weights=weights_sub, color=data_colors[2], label="Data - Background", fill=False, histtype='step', lw=2)

for i, wave in enumerate(waves_sorted):
    amp_letter = wave[0]
    amp_m = wave[1]
    if int(amp_m) > 0:
        amp_m_sign = wave[2]
    else:
        amp_m_sign = ""

    wave_neg = wave + "-"
    if wave_neg in amplitudes:
        print(wave + "\t\t-e")
        plt.errorbar(bin_df['Centers'].iloc[df_filtered['Bin']], df_filtered[wave_neg + ac_tag], yerr=df_filtered[wave_neg + "_err" + ac_tag], elinewidth=0.5, fmt='o', color=colors[i], label=rf"${amp_letter}_{{{amp_m_sign}{amp_m}}}$")
print("tot")
plt.errorbar(bin_df['Centers'].iloc[df_filtered['Bin']], df_filtered['total_intensity' + ac_tag_total], yerr=df_filtered['total_intensity_err' + ac_tag_total], elinewidth=0.5, fmt='none', color='k')
plt.hist(bin_df['Centers'].iloc[df_filtered['Bin']], bins=len(bin_df), range=(bin_df['mass'].iloc[0], bin_df['mass'].iloc[-1]), weights=df_filtered['total_intensity' + ac_tag_total], fill=False, histtype='step', color='k', label="Total") # this is a sneaky hack
plt.xlim(bin_df['Centers'].iloc[0] - 0.1, bin_df['Centers'].iloc[-1] + 0.1)
plt.ylim(bottom=0)
plt.legend(loc="upper right")
plt.title("Negative Reflectivity")
plt.ylabel("Intensity")
plt.xlabel(xlabel)
plt.tight_layout()
pdf.savefig(fig, dpi=300)

print("Positive and Negative Reflectivity")
fig = plt.figure()

# plot real data histograms
if args.data:
    plt.hist(invariant_mass, bins=nbins, range=(bin_df['mass'].iloc[0], bin_df['mass'].iloc[-1]), weights=weights, color=data_colors[0], label="Data", fill=False, histtype='step', lw=2)
if args.background:
    plt.hist(invariant_mass_bkg, bins=nbins, range=(bin_df['mass'].iloc[0], bin_df['mass'].iloc[-1]), weights=weights_bkg, color=data_colors[1], label="Background", fill=False, histtype='step', lw=2)
if args.subtracted:
    plt.hist(invariant_mass_sub, bins=nbins, range=(bin_df['mass'].iloc[0], bin_df['mass'].iloc[-1]), weights=weights_sub, color=data_colors[2], label="Data - Background", fill=False, histtype='step', lw=2)

for i, wave in enumerate(waves_sorted):
    amp_letter = wave[0]
    amp_m = wave[1]
    if int(amp_m) > 0:
        amp_m_sign = wave[2]
    else:
        amp_m_sign = ""

    wave_pos = wave + "+"
    print(wave, end='\t')
    if wave_pos in amplitudes:
        print("+e", end='\t')
        plt.errorbar(bin_df['Centers'].iloc[df_filtered['Bin']], df_filtered[wave_pos + ac_tag], yerr=df_filtered[wave_pos + "_err" + ac_tag], elinewidth=0.5, fmt='o', color=colors[i], label=rf"${amp_letter}^+_{{{amp_m_sign}{amp_m}}}$")
    else:
        print("", end='\t')

    wave_neg = wave + "-"
    if wave_neg in amplitudes:
        print("-e", end='\t')
        plt.errorbar(bin_df['Centers'].iloc[df_filtered['Bin']], df_filtered[wave_neg + ac_tag], yerr=df_filtered[wave_neg + "_err" + ac_tag], elinewidth=0.5, fmt='s', color=colors[i], label=rf"${amp_letter}^-_{{{amp_m_sign}{amp_m}}}$")
    else:
        print("", end='\t')
    print()
print("tot")
plt.errorbar(bin_df['Centers'].iloc[df_filtered['Bin']], df_filtered['total_intensity' + ac_tag_total], yerr=df_filtered['total_intensity_err' + ac_tag_total], elinewidth=0.5, fmt='none', color='k')
plt.hist(bin_df['Centers'].iloc[df_filtered['Bin']], bins=len(bin_df), range=(bin_df['mass'].iloc[0], bin_df['mass'].iloc[-1]), weights=df_filtered['total_intensity' + ac_tag_total], fill=False, histtype='step', color='k', label="Total") # this is a sneaky hack
plt.xlim(bin_df['Centers'].iloc[0] - 0.1, bin_df['Centers'].iloc[-1] + 0.1)
plt.ylim(bottom=0)
plt.legend(loc="upper right")
plt.title("All Waves")
plt.ylabel("Intensity")
plt.xlabel(xlabel)
plt.tight_layout()
pdf.savefig(fig, dpi=300)

phase_tag = "_PHASE"
phase_diffs = [column[:-len(phase_tag)] for column in df.columns.to_list() if column.endswith(phase_tag) and not "_err" in column]
for phase_diff in phase_diffs:
    fig, ax = plt.subplots()
    ax2 = ax.twinx()
    ax.errorbar(bin_df['Centers'].iloc[df_filtered['Bin']], df_filtered['total_intensity' + ac_tag_total], yerr=df_filtered['total_intensity_err' + ac_tag_total], elinewidth=0.5, fmt='none', color='k')
    ax.hist(bin_df['Centers'].iloc[df_filtered['Bin']], bins=len(bin_df), range=(bin_df['mass'].iloc[0], bin_df['mass'].iloc[-1]), weights=df_filtered['total_intensity' + ac_tag_total], fill=False, histtype='step', color='k', label="Total") # this is a sneaky hack
    ax2.errorbar(bin_df['Centers'].iloc[df_filtered['Bin']], df_filtered[phase_diff + phase_tag], yerr=df_filtered[phase_diff + "_err" + phase_tag], elinewidth=0.5, fmt='o', color='m')
    plt.xlim(bin_df['Centers'].iloc[0] - 0.1, bin_df['Centers'].iloc[-1] + 0.1)
    plt.ylim(bottom=-2*np.pi, top=2*np.pi)
    plt.title(f"Phase Difference {phase_diff}")
    ax.set_ylabel("Intensity")
    ax2.set_ylabel("Phase")
    plt.xlabel(xlabel)
    plt.tight_layout()
    pdf.savefig(fig, dpi=300)
pdf.close()
