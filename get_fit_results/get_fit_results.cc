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
    // Parse command line input (fit file, function name, param1, param2, ...)
    string fitFile(argv[1]); // 0 is the name of the program
    string functionType(argv[2]); // Zlm or TwoPSAngles
    string zlmString("Zlm");
    string twopsanglesString("TwoPSAngles");
    FitResults results(fitFile.c_str()); // create FitResults object
    bool doAcceptanceCorrection = true; // let's hard-code this for now
    if (functionType == zlmString) {
        for (int i = 3; i < argc; i+=2) { // for each parameter (order is ###Im::### ###Re::### so we take two args at once)
            vector<string> amplitudeStringImRe;
            string amplitudeStringIm(argv[i]);
            string amplitudeStringRe(argv[i + 1]);
            amplitudeStringImRe.push_back(amplitudeStringIm.c_str());
            amplitudeStringImRe.push_back(amplitudeStringRe.c_str());
            // Use AmpTools FitResults to calculate intensity for each individual wave's contribution
            pair<double, double> intensityResult = results.intensity(amplitudeStringImRe, doAcceptanceCorrection);
            cout << intensityResult.first << "\t" << intensityResult.second << "\t";
        }
    } else if (functionType == twopsanglesString) {
        for (int i = 3; i < argc; i++) { // for this function, we only have one parameter
            vector<string> amplitudeVect;
            string amplitudeString(argv[i]);
            amplitudeVect.push_back(amplitudeString.c_str());
            // Use AmpTools FitResults to calculate intensity for each individual wave's contribution
            pair<double, double> intensityResult = results.intensity(amplitudeVect, doAcceptanceCorrection);
            cout << intensityResult.first << "\t" << intensityResult.second << "\t";
        }
    } else {
        cout << "function not found: " << functionType << "\t";
    }
    // Calculate the total intensity from all waves
    pair<double, double> totalIntensityResult = results.intensity(doAcceptanceCorrection);
    // Calculate the likelihood
    cout << totalIntensityResult.first << "\t" << totalIntensityResult.second << "\t" << results.likelihood() << "\t";
}
