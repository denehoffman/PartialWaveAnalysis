/*
 * get_fit_results.cc: Prints out a tab-separated line of fit results from an AmpTools .fit file
 *
 * Author: Nathaniel Dene Hoffman - Carnegie Mellon University - GlueX Collaboration
 * Creation Date: 5 July 2021
 *
 * Makefile was created by Naomi Jarvis (CMU)
 */
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <unistd.h>

#include "IUAmpTools/FitResults.h"

using namespace std;

int main(int argc, char* argv[]) {
    // Parse command line input (fit file, command, param1, param2, ...)
    string fitFile(argv[1]); // 0 is the name of the program
    string commandType(argv[2]); // command
    string intensityStr("intensity");
    string intensityTotalStr("intensityTotal");
    string likelihoodStr("likelihood");
    string realImagStr("realImag");
    string phaseDiffStr("phaseDiff");
    FitResults results(fitFile.c_str()); // create FitResults object
    bool doAcceptanceCorrection = true; // let's hard-code this for now
    if (commandType == intensityStr) {
        vector<string> amplitudes;
        for (int i = 3; i < argc; i++) { 
            string amplitudeString(argv[i]);
            amplitudes.push_back(amplitudeString.c_str());
        }
        // Use AmpTools FitResults to calculate intensity for each individual wave's (or set of waves) contribution
        pair<double, double> intensityResult = results.intensity(amplitudes, doAcceptanceCorrection);
        cout << "#" << intensityResult.first << "\t" << intensityResult.second;
    } else if (commandType == likelihoodStr) {
        // Calculate the likelihood
        cout << "#" << results.likelihood();
    } else if (commandType == intensityTotalStr) {
        // Calculate the total intensity from all waves
        pair<double, double> totalIntensityResult = results.intensity(doAcceptanceCorrection);
        cout << "#" << totalIntensityResult.first << "\t" << totalIntensityResult.second;
    } else if (commandType == realImagStr) {
        // Calculate real and imaginary parts of amplitude
        cout << "#";
        for (int i = 3; i < argc; i++) {
            string amplitudeString(argv[i]);
            cout << results.productionParameter(amplitudeString.c_str()).real() << "\t" << results.productionParameter(amplitudeString.c_str()).imag();
            if (i != argc - 1) {
                cout << "\t";
            }
        }
    } else if (commandType == phaseDiffStr) {
        string amplitudeString1(argv[3]);
        string amplitudeString2(argv[4]);
        pair<double, double> phase = results.phaseDiff(amplitudeString1.c_str(), amplitudeString2.c_str());
        cout << "#" << phase.first << "\t" << phase.second;
    } else {
        cout << "function not found: " << commandType;
    }
}
