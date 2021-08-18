#!/usr/bin/env python3
"""
generate_config.py: This program generates configuration files for AmpTools PWA fits.
    These files are templates for input into divide_data_pol.py and contain tags which
    are later replaced during binning or fitting.

Author: Nathaniel Dene Hoffman - Carnegie Mellon University - GlueX Collaboration
Creation Date: 23 July 2021
"""

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd
from simple_term_menu import TerminalMenu


def clear():
    """Clear terminal screen"""
    if os.name == 'nt':
        os.system('cls')
    else:
        os.system('clear')


class Wave:
    """The Wave object contains the information required to write an AmpTools amplitude
    using the Zlm function basis. This is the recommended method for incorporating
    polarization into your fits.
    """
    REAL = True
    IMAG = False
    POSITIVE = +1
    NEGATIVE = -1
    letter_dict = {0: "S", 1: "P", 2: "D"}

    def __init__(self, reaction: str, l: int, m: int, e: int):
        assert abs(m) <= l, f"-L <= M <= L is required!\tL = {l}\t M = {m}"
        assert abs(e) == 1, f"Reflectivity must be unitary!\t|{e}| != 1"
        assert l >= 0, f"L must be non-negative!\tL = {l}"
        self.reaction = reaction
        self.l = l
        self.m = m
        self.e = e

    def __str__(self):
        if self.m > 0:
            m_sign = "+"
        elif self.m < 0:
            m_sign = "-"
        else:
            m_sign = " "

        if self.e > 0:
            reflectivity_str = "Positive"
        else:
            reflectivity_str = "Negative"

        return (f"{Wave.letter_dict.get(self.l)}-Wave  {self.l} {m_sign}{abs(self.m)}  "
                f"{reflectivity_str} Reflectivity")

    def __hash__(self):
        return hash(repr(self))

    def __repr__(self):
        return str(self)

    def __lt__(self, other):
        if self.l == other.l:
            if self.m == other.m:
                return self.e > other.e
            return self.m < other.m
        return self.l < other.l

    def __gt__(self, other):
        if self.l == other.l:
            if self.m == other.m:
                return self.e < other.e
            return self.m > other.m
        return self.l > other.l

    def __eq__(self, other):
        return (self.l == other.l) and (self.m == other.m) and (self.e == other.e)

    def get_wave_letter(self) -> str:
        """gets the letter associated with an L value"""
        letter = Wave.letter_dict.get(self.l)
        if letter is None:
            print(f"Wave with L = {l} is currently not supported (0 <= L <= 2)")
            return "X"
        return letter

    def get_wave_string(self, real: bool, pol="") -> str:
        """get_wave_string produces a string like KsKs_000::PositiveRe::D2++,
        the details of which are determined by the wave properties.

        :param real: bool: true if the amplitude is paired with Re[Zlm]
        :param pol: The polarization string (Default value = "")

        """
        output = self.reaction + pol + "::"
        if self.e > 0:
            output += "Positive"
        else:
            output += "Negative"

        if real:
            output += "Re::"
        else:
            output += "Im::"

        output += self.get_wave_letter()

        if self.m > 0:
            output += str(self.m) + "+"
        elif self.m < 0:
            output += str(abs(self.m)) + "-"
        else:
            output += str(self.m)

        if self.e > 0:
            output += "+"
        else:
            output += "-"

        return output

    def get_wave(self, init_real=False, cartesian=True):
        """

        :param init_real: fixes Im[amplitude] to zero (Default value = False)
        :param cartesian: uses Cartesian (or Polar) coordinates (Default value = True)

        """
        wave_string_real = self.get_wave_string(Wave.REAL)
        wave_string_real_000 = self.get_wave_string(Wave.REAL, pol="_000")

        wave_string_imag = self.get_wave_string(Wave.IMAG)
        wave_string_imag_000 = self.get_wave_string(Wave.IMAG, pol="_000")

        if self.e > 0:
            amplitude_string = (f"amplitude {wave_string_real} Zlm {self.l} {self.m} "
                                f"+1 +1 LOOPPOLANG LOOPPOLVAL\namplitude "
                                f"{wave_string_imag} Zlm {self.l} {self.m} -1 -1 "
                                f"LOOPPOLANG LOOPPOLVAL\n")
        else:
            amplitude_string = (f"amplitude {wave_string_real} Zlm {self.l} {self.m} "
                                f"+1 -1 LOOPPOLANG LOOPPOLVAL\namplitude "
                                f"{wave_string_imag} Zlm {self.l} {self.m} -1 +1 "
                                f"LOOPPOLANG LOOPPOLVAL\n")

        if cartesian:
            cartesian_str = "cartesian"
        else:
            cartesian_str = "polar"
        if init_real:
            initialize_string = (f"initialize {wave_string_real} {cartesian_str} "
                                 f"@uniform 0.0 real\n")
        else:
            initialize_string = (f"initialize {wave_string_real} {cartesian_str} "
                                 f"@uniform @uniform\n")

        constrain_string = f"constrain {wave_string_real} {wave_string_imag}\n"
        constrain_string += f"constrain {wave_string_real_000} {wave_string_real}\n"
        constrain_string += f"constrain {wave_string_imag_000} {wave_string_imag}\n"

        scale_string = (f"scale {wave_string_real} LOOPSCALE\n"
                        f"scale {wave_string_imag} LOOPSCALE\n")

        return amplitude_string, initialize_string, constrain_string, scale_string


def generate_config(reaction: str, waves: list, n_pos: int, n_neg: int,
                    use_background: bool, use_cartesian: bool) -> None:
    """
    :param reaction: str: reaction name
    :param waves: list: list of Wave objects
    :param n_pos: int: number of positive reflectivities
    :param n_neg: int: number of negative reflectivities
    :param use_background: bool: set True to include the background string
    :param use_cartesian: bool: set True(False) to use Cartesian(Polar) coordinates
    """
    # This text does not change between configs (yet):
    text = f"""define polVal_000 0.3519
define polVal_045 0.3374
define polVal_090 0.3303
define polVal_135 0.3375
define polVal_AMO 0.00001

define polAngle_000 0.0
define polAngle_045 45.0
define polAngle_090 90.0
define polAngle_135 135.0
define polAngle_AMO 0.0

parameter polScale_000 1.0 fixed
parameter polScale_045 1.0
parameter polScale_090 1.0
parameter polScale_135 1.0
parameter polScale_AMO 1.0

fit {reaction}
loop {reaction} {reaction}_000 {reaction}_045 {reaction}_090 {reaction}_135 {reaction}_AMO

loop LOOPDATA @DATAFILE_000 @DATAFILE_045 @DATAFILE_090 @DATAFILE_135 @DATAFILE_AMO
loop LOOPACC @ACCFILE_000 @ACCFILE_045 @ACCFILE_090 @ACCFILE_135 @ACCFILE_AMO
loop LOOPGEN @GENFILE_000 @GENFILE_045 @GENFILE_090 @GENFILE_135 @GENFILE_AMO\n"""

    if use_background: # add background line
        text += "loop LOOPBKG @BKGFILE_000 @BKGFILE_045 @BKGFILE_090 @BKGFILE_135 @BKGFILE_AMO\n"

    text += f"""loop LOOPNIFILE @NIFILE_000 @NIFILE_045 @NIFILE_090 @NIFILE_135 @NIFILE_AMO

loop LOOPPOLANG polAngle_000 polAngle_045 polAngle_090 polAngle_135 polAngle_AMO
loop LOOPPOLVAL polVal_000 polVal_045 polVal_090 polVal_135 polVal_AMO
loop LOOPSCALE [polScale_000] [polScale_045] [polScale_090] [polScale_135] [polScale_AMO]

normintfile {reaction} LOOPNIFILE\n

data {reaction} ROOTDataReader LOOPDATA
genmc {reaction} ROOTDataReader LOOPGEN
accmc {reaction} ROOTDataReader LOOPACC\n"""
    if use_background: # add background line
        text += f"bkgnd {reaction} ROOTDataReader LOOPBKG\n"

    text += f"reaction {reaction} " + input(
        "Enter the particles in the reaction separated by spaces (i.e. gamma Proton Ks1 Ks2): "
    ) + "\n\n" # prompt for particle names

    if n_neg > 0:
        text += f"sum {reaction} NegativeRe\nsum {reaction} NegativeIm\n"
    if n_pos > 0:
        text += f"sum {reaction} PositiveRe\nsum {reaction} PositiveIm\n"

    text += "\n\n# Amplitudes\n\n"
    for wave in waves:
        amp, _, _, _ = wave.get_wave()
        text += amp
    text += "\n\n# Initialize Real Parts\n\n"
    made_neg_real = n_neg == 0
    made_pos_real = n_pos == 0
    for i, wave in enumerate(waves):
        if wave.e > 0:
            if not made_pos_real:
                _, init, _, _ = wave.get_wave(init_real=True,
                                              cartesian=use_cartesian)
                made_pos_real = True
            else:
                _, init, _, _ = wave.get_wave(cartesian=use_cartesian)
        else:
            if not made_neg_real:
                _, init, _, _ = wave.get_wave(init_real=True,
                                              cartesian=use_cartesian)
                made_neg_real = True
            else:
                _, init, _, _ = wave.get_wave(cartesian=use_cartesian)
        text += init
    text += "\n\n# Constrain Real and Imaginary Zlm amplitudes\n\n"
    for wave in waves:
        _, _, const, _ = wave.get_wave()
        text += const
    text += "\n\n# Scale amplitudes\n\n"
    for wave in waves:
        _, _, _, scale = wave.get_wave()
        text += scale

    file_exists = True
    while file_exists:
        filename = input("Enter a name for the new config file: ")
        if not filename.endswith(".cfg"):
            filename += ".cfg"
        file_path = Path(filename).resolve()
        if file_path.is_file():
            print("Please select a different name, that one is already taken!")
        else:
            file_exists = False
    with open(Path(filename).resolve(), 'w') as file:
        file.write(text)
        print(f"Configuration file has been generated and saved to {filename}")


def print_waves(wave_set: set):
    """

    :param wave_set: set: 

    """
    print("╔═════════════════════════════════════════╗")
    wave_list = list(wave_set)
    wave_list.sort()
    s_waves = "\n".join(
        ["║ " + str(wave) + " ║" for wave in wave_list if wave.l == 0])
    p_waves = "\n".join(
        ["║ " + str(wave) + " ║" for wave in wave_list if wave.l == 1])
    d_waves = "\n".join(
        ["║ " + str(wave) + " ║" for wave in wave_list if wave.l == 2])
    if len(s_waves) != 0:
        print(s_waves)
        if len(p_waves) == 0 and len(d_waves) == 0:
            print("╚═════════════════════════════════════════╝")
        else:
            print("╠═════════════════════════════════════════╣")
    if len(p_waves) != 0:
        print(p_waves)
        if len(d_waves) == 0:
            print("╚═════════════════════════════════════════╝")
        else:
            print("╠═════════════════════════════════════════╣")
    if len(d_waves) != 0:
        print(d_waves)
        print("╚═════════════════════════════════════════╝")
    print()


def main():
    """ """
    clear()
    reaction = input("Enter a name for this reaction: ")
    waves = set()

    main_menu_title = f"AmpTools Configuration Generator\n"
    main_menu_items = ["Add Wave", "Remove Wave", "Generate", "Exit"]
    main_menu_cursor = "> "
    main_menu_cursor_style = ("fg_red", "bold")
    main_menu_style = ("bg_black", "fg_green")
    main_menu_exit = False

    main_menu = TerminalMenu(menu_entries=main_menu_items,
                             title=main_menu_title,
                             menu_cursor=main_menu_cursor,
                             menu_cursor_style=main_menu_cursor_style,
                             menu_highlight_style=main_menu_style)

    add_menu_title = "Add Wave\n"
    add_menu_items = ["[s] S Wave", "[p] P Wave", "[d] D Wave"]
    add_menu = TerminalMenu(menu_entries=add_menu_items,
                            title=add_menu_title,
                            menu_cursor=main_menu_cursor,
                            menu_cursor_style=main_menu_cursor_style,
                            menu_highlight_style=main_menu_style,
                            cycle_cursor=True,
                            clear_screen=True)

    p_wave_menu_title = "Select P Orbital M\n"
    p_wave_menu_items = ["+1", "0", "-1"]
    p_wave_menu = TerminalMenu(menu_entries=p_wave_menu_items,
                               title=p_wave_menu_title,
                               menu_cursor=main_menu_cursor,
                               menu_cursor_style=main_menu_cursor_style,
                               menu_highlight_style=main_menu_style,
                               cycle_cursor=True,
                               clear_screen=True)
    d_wave_menu_title = "Select D Orbital M\n"
    d_wave_menu_items = ["+2", "+1", "0", "-1", "-2"]
    d_wave_menu = TerminalMenu(menu_entries=d_wave_menu_items,
                               title=d_wave_menu_title,
                               menu_cursor=main_menu_cursor,
                               menu_cursor_style=main_menu_cursor_style,
                               menu_highlight_style=main_menu_style,
                               cycle_cursor=True,
                               clear_screen=True)

    reflectivity_menu_title = "Select Reflectivity\n"
    reflectivity_menu_items = ["+1", "-1"]
    reflectivity_menu = TerminalMenu(menu_entries=reflectivity_menu_items,
                                     title=reflectivity_menu_title,
                                     menu_cursor=main_menu_cursor,
                                     menu_cursor_style=main_menu_cursor_style,
                                     menu_highlight_style=main_menu_style,
                                     cycle_cursor=True,
                                     clear_screen=True,
                                     multi_select=True,
                                     show_multi_select_hint=True)

    while not main_menu_exit:
        clear()
        if len(waves) != 0:
            print("Current Waves:")
            print_waves(waves)
        main_sel = main_menu.show()

        if main_sel == 0:
            add_sel = add_menu.show()
            if add_sel == 0:
                is_reflectivity_selected = False
                while not is_reflectivity_selected:
                    refl_sels = reflectivity_menu.show()
                    if len(refl_sels) == 0:
                        print("You must select at least one reflectivity!")
                    else:
                        for i in refl_sels:
                            if i == 0:
                                waves.add(Wave(reaction, 0, 0, 1))
                            else:
                                waves.add(Wave(reaction, 0, 0, -1))
                        is_reflectivity_selected = True
            elif add_sel == 1:
                p_wave_sel = p_wave_menu.show()
                wave_m = [1, 0, -1][p_wave_sel]
                is_reflectivity_selected = False
                while not is_reflectivity_selected:
                    refl_sels = reflectivity_menu.show()
                    if len(refl_sels) == 0:
                        print("You must select at least one reflectivity!")
                    else:
                        for i in refl_sels:
                            if i == 0:
                                waves.add(Wave(reaction, 1, wave_m, 1))
                            else:
                                waves.add(Wave(reaction, 1, wave_m, -1))
                        is_reflectivity_selected = True
            elif add_sel == 2:
                d_wave_sel = d_wave_menu.show()
                wave_m = [2, 1, 0, -1, -2][d_wave_sel]
                is_reflectivity_selected = False
                while not is_reflectivity_selected:
                    refl_sels = reflectivity_menu.show()
                    if len(refl_sels) == 0:
                        print("You must select at least one reflectivity!")
                    else:
                        for i in refl_sels:
                            if i == 0:
                                waves.add(Wave(reaction, 2, wave_m, 1))
                            else:
                                waves.add(Wave(reaction, 2, wave_m, -1))
                        is_reflectivity_selected = True
        elif main_sel == 1:
            if len(waves) != 0:
                remove_wave_menu_title = "Remove Waves\n"
                wave_list = list(waves)
                wave_list.sort()
                wave_list.append("Cancel")
                remove_wave_menu_items = [str(wave) for wave in wave_list]
                remove_wave_menu = TerminalMenu(
                    menu_entries=remove_wave_menu_items,
                    title=remove_wave_menu_title,
                    menu_cursor=main_menu_cursor,
                    menu_cursor_style=main_menu_cursor_style,
                    menu_highlight_style=main_menu_style,
                    cycle_cursor=True,
                    clear_screen=True,
                    multi_select=True,
                    show_multi_select_hint=True)
                waves_to_remove = remove_wave_menu.show()
                if waves_to_remove[-1] != len(wave_list) - 1:
                    for to_remove in waves_to_remove:
                        wave_to_remove = wave_list[to_remove]
                        waves.remove(wave_to_remove)
        elif main_sel == 2:
            if len(waves) != 0:
                wave_list = list(waves)
                wave_list.sort()
                n_pos = len([wave for wave in wave_list if wave.e > 0])
                n_neg = len([wave for wave in wave_list if wave.e < 0])
                answered = False
                while not answered:
                    include_bkg = input(
                        "Will you be using background files? [y/n]: ")
                    if include_bkg != 'y' and include_bkg != 'n':
                        clear()
                        print("Please select either 'y' or 'n'!")
                    else:
                        use_background = include_bkg == 'y'
                        answered = True
                answered = False
                while not answered:
                    cartesian = input(
                        "Initialize with polar (p) or Cartesian (c) coordinates? [p/c]: "
                    )
                    if cartesian != 'p' and cartesian != 'c':
                        clear()
                        print("Please select either 'p' or 'c'!")
                    else:
                        use_cartesian = cartesian == 'c'
                        answered = True
                generate_config(reaction, wave_list, n_pos, n_neg,
                                use_background, use_cartesian)
        else:
            main_menu_exit = True
            clear()


if __name__ == "__main__":
    main()
