# PartialWaveAnalysis

Usage: This repo contains all the code I use to perform a partial wave analysis using the AmpTools software provided by IU at GlueX. To utilize this software, you must already have access to the [halld_sim](https://github.com/JeffersonLab/halld_sim) suite of programs, since this code calls upon `FitResults.h` and `fit.C` from the AmpTools library and `split_mass`. ~~You also need access to the [gluex_root_analysis](https://github.com/JeffersonLab/gluex_root_analysis) program [`tree_to_amptools`](https://github.com/JeffersonLab/gluex_root_analysis/tree/master/programs/tree_to_amptools).~~ I no longer advocate the use of `tree_to_amptools` since it seems to not always create the right number of final state particles when running over generated thrown trees.

### Workflow
1. Download and install [halld_sim](https://github.com/JeffersonLab/halld_sim) and [rcdb](https://github.com/JeffersonLab/rcdb). The second one contains the Python library required to run `hadd_rcdb.py`.
2. Run `pip3 install --user -r requirements.txt` or use your favorite Python package installation method to get all of the required Python packages.
3. Use `hadd_rcdb.py <min run number> <max run number> <directory> [basename]` to merge your data into separate `.root` files by polarization.
    * (Optional) Do the same with some "background" data, like out-of-time bins or something.
4. Generate Thrown and Accepted `.root` trees using polarized flux files (use [psflux](https://github.com/JeffersonLab/hd_utilities/tree/master/psflux) on iFarm) and merge by polarization.
5. Turn all `.root` files into AmpTools flat trees (see method below) and all files of each polarization into separate folders for data, thrown, accepted, and (optional) background.
    * Files must include the the strings "AMO", "PARA_0", "PARA_135", "PERP_45", and "PERP_90" somewhere in their names.
6. Run `generate_config.py` to create your first config file with a few selected amplitudes.
7. Run `divide_data_pol.py`, pointing to the config file, the directories for your data, and an output directory (will be made if it doesn't exit). This is also where you specify your maximum and minimum mass bin, as well as the number of bins to use.
8. Run `run_amptools.py`, pointing to the output directory from the previous step. You can also specify the number of times you want to run each fit (iterations) and the seed (for reproducibility). Currently supports SLURM and Python Multiprocessing.
9. Run `plot_results.py` and/or `plot_stats.py` to make nice plots of the waves and the statistical properties of the fits. Both take the resulting fit file from the previous step as their input, and optionally take an x-axis label in LaTeX format as the second argument.
10. Run `run_amptools.py` with the `--bootstrap` flag (you can change the number of iterations here) to bootstrap the fit. This should create more reliable error bars on the amplitudes.
11. Run `plot_bootstrap.py` to see the results of bootstrapping (takes the same arguments as the othr plotting methods).
12. Generate more configs and run `divide_data_pol.py`, only specifying the output directory and the config file to add that wave set to the analysis. Then repeat steps 8-11 with this new config file.

#### Workflow (detailed)
1. Create flat trees from MC Thrown, MC Reconstructed, Data, and (optional) Background ROOT TTrees.
    1. Use [`SetupAmpTools_FlatTree()`](https://github.com/JeffersonLab/gluex_root_analysis/blob/a2c0dddc6e7b3fce28bb1919843ca676c4482975/libraries/DSelector/DSelector.cc#L932) and [`FillAmpTools_FlatTree`](https://github.com/JeffersonLab/gluex_root_analysis/blob/a2c0dddc6e7b3fce28bb1919843ca676c4482975/libraries/DSelector/DSelector.cc#L948) in a DSelector (make sure you also tell the DSelector to fill these trees). Asside from running `SetupAmpTools_FlatTree()`, I also had to do a bit of code to get the proper trees filled:
    ```
    vector<TLorentzVector> locFinalStateP4;
    locFinalStateP4.push_back(locProtonP4);
    locFinalStateP4.push_back(locDecayingKShort1P4);
    locFinalStateP4.push_back(locDecayingKShort2P4);
    dFlatTreeInterface->Fill_Fundamental<Float_t>("Weight", 1.0); // haven't dealt with accidentals yet
    Fill_FlatTree();
    ```
    I insert this right at the end of the combo loop. A bit more work is needed to create these kinematic variables in the first place for the Thrown Trees, but the idea is exactly the same. Note that when this is run, you should avoid generating the default flat trees by running `tree_name->Process("DSelector_name.C+", "DefaultFlatOff")` in `ROOT`. Additionally, I use `hadd_rcdb.py` to divide my files by polarization and generate MC with polarized flux files (no actual polarization dependent physics is included in the MC).
2. Run `generate_config.py` to select which waves you want to include in your config file. The interface should be fairly straightforward. If you want to make yours manually, there are a few tags (preceeded by an `@` symbol) which you can find by looking in the code for `divide_data_pol.py` or `run_amptools.py` scripts.
3. Run `divide_data_pol.py` provided by this repo. Running it without arguments will show a help string with the required arguments and a short description of their usage. This program essentially wraps the `split_mass` program across the AmpTools trees you just generated. `divide_data.py` is essentially deprecated, and it is (probably) more accurate to include polarization data in your fits.
4. Run `run_amptools.py` to actually perform the fits with `fit.C`. Again, running without arguments will display the help message.
    1. Additionally, the `-s <integer>` option allows you to set the seed (default is `seed = 1`). The seed is used to randomly generate starting parameter values to be used for the fit. In theory, if you use the same seed and the same number of bins and iterations, all of those generated values should be consistent across each run, as long as you aren't adding new parameters to your fits.
    2. `run_amptools.py` will also call `get_fit_results`. However, if you just downloaded this repo, you need to build `get_fit_results` by going into that directory and running the command `make`. You might (probably will) get some errors about directories not existing. After creating the directories as suggested, the compilation worked. There is probably a better way to do this, but I just copied the Makefile from Naomi Jarvis, and I don't actually know how it works.
    3. By default, this program will submit fits as SLURM jobs using the `sbatch` command. If you don't use SLURM, run with the option `--parallel Pool` to use Python multiprocessing
        1. The `-p <integer>` option will allow you to specify how many multiprocessing processes to generate. This is limited by the number of cores you have available. Python unfortunately does not support multithreading due to the Python Interpreter Lock, but generally my fits haven't taken very long (31 bins with 20 iterations took about 7 minutes the last time I ran it). Don't worry about not knowing how many cores you have, this option is more to limit the number used if you do know it, otherwise, just don't use this flag and the program will generate up to 60 processes.
    4. The `gather.py` program included just re-gathers the fit data after a fit is performed. It is called at the end of `run_amptools.py`, but can be called again if you accidentally delete a fit file or something.
5. If you've gotten to this step, you will now have an output directory that contains bin folders as well as a file titled `bin_info.txt` and one titled `<name of config>::fit_results.txt`. These are human-readable tab-separated files, so feel free to look inside them and see how they are structured. I wrote a quick script called `plot_results.py` to do some preliminary graphing, but it is very much a work-in-progress. `plot_stats.py` has the same syntax and will generate a set of plots which allow you to see some of the statistical properties of your fits. Both of these plotting methods take an optional second argument, the x-axis label for the plots, which can be written in LaTeX format (you will likely have to use quotes or escape characters if typing this into the shell). The first argument should be the path to the `*fit_results.txt` file.
6. After the initial fit, you can bootstrap to get more accurate error bars. The proceedure is fairly simple. You can run `run_amptools.py` again with the same general settings as before (although feel free to change the number of iterations used in bootstrapping), with the addition of the `--bootstrap` flag. This will generate another file called `<name of config>::bootstrap.txt` which has the same format as `<name of config>::fit_results.txt` but should not be used for anything but error estimation, since it is the result of fits to oversampled data. Finally, `plot_bootstrap.py` will generate some useful plots of this result, and like the other plotting programs, it only needs the fit file (not `*bootstrap.txt`, but the `*fit_results.txt` file) as its main argument.

## Troubleshooting
* > `latest_fit_file = max(fit_files, key=os.path.getctime)
ValueError: max() arg is an empty sequence`

If you see this, check the `/logs/` directory and look at a `.err` file. You'll probably see

`FileNotFoundError: [Errno 2] No such file or directory: '</path/to/fit/>/<bin number>/<iteration number>/<reaction name>.fit' -> '</path/to/fit/>/<bin number>/<iteration number>/<reaction name>::UNKNOWN.fit'`

This typically means your configuration file has an error somewhere, or AmpTools ran into a different error and no fit file was produced. One way to troubleshoot this is to go into that folder, copy in the `.root` files from the parent bin directory, and run the fit manually with `fit -c <config>.cfg`. This will give you the AmpTools error outputs, which might be more helpful for finding a solution.

* > `background_path_AMO = [filepath for filepath in bin_dir.glob(f"*_BKG__{bin_num}.root") if "AMO" in filepath.name][0]
IndexError: list index out of range`

Your config file wants background files, but you didn't include any when you divided the data into mass bins.

* > I ran fit -c <config>.cfg and got `ConfigFileParser ERROR:  Mismatch in loop size:
(31) <config>_<bin>-<iteration>.cfg >>   normintfile <reaction> LOOPNIFILE`
    
This could also point to a different loop, but it generally means that you have some "::" characters in one of your filenames somewhere (this particular one would come from the config file, but if you had that in your background file you'd have a problem too).

### TODO:
1. ~~Need to implement polarization in this analysis~~
2. ~~Need a better way of sending the proper waves in pairs to `get_fit_results` than just hoping the user writes them in a nice order~~
3. ~~`plot_results.py` needs to be more robust and versatile. I want to be able to get separate subplots for each wave, but I want it to look nice also.~~ I also want this program to have some command-line options to specify which waves to plot (or other things to plot)
4. `get_fit_results`
    1. ~~This program needs to be expanded to include more information about the fits, such as phases between waves and the base amplitudes from the fit.~~ Currently none of this code is suited to do any SDMEs, but eventually some separate programs should be added to handle these too.
    2. ~~I should eventually try to implement a way to drop specific waves from the analysis. It might be easier to implement this in `run_amptools.py`.~~ (made it easier to generate more config files for one set of data bins)
5. Bootstrapping methods take standard deviation of all converged runs, ~~but maybe they should fit a gaussian with a mean at the fit value instead?~~ (no, you shouldn't)
6. Create a method of generating fake data for testing purposes.



## Version History
v1.0 (6 August 2021):
* Original set of commits after development
* Generate configuration files with up to J=2 waves
* Divide polarized data into invariant mass bins
* Run multiple fit iterations with a config file for each bin
* Add new config files to existing bin structures
* Plot intensity profiles for each wave
* Plot statistical data for all fits
* Bootstrap fits to get accurate error bars and confidence intervals

v1.1 (13 August 2021):
* SLURM and other parallelism methods can pick up from where they left off on fits
* Added uproot3 support to allow plotting real data behind intensity plots
* Added option to specify tree names other than 'kin' for files in divide_data_pol.py
* Added verbose loading bars with enlighten
* Corrected issue where running a second time with fewer iterations gathered iterations from the previous run
* Added phase difference information to fit results files
* Added plot_amps.py to examine locations of complex amplitudes
* All plotting methods now only require one parameter pointing them to a fit directory
* Fixed incorrect implementation of the Zlm class in generate_config.py
* Added capability to use config files with fewer polarizations
* Added improved logging methods
