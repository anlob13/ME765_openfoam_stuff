// Minimal standalone test harness: loads unifiedWallModel from wmLayers_openfoam/
// (the real trained weights for all 4 building blocks + the classifier) and
// checks its predictions against generate_verification_rows.py's reference
// output - one row per regime, same verification methodology used throughout
// this project.
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <string>

#include "unified_wall_model.H"

int main()
{
    Foam::unifiedWallModel model("wmLayers_openfoam");

    std::ifstream f("verification_rows.txt");
    if (!f)
    {
        std::cerr << "Cannot open verification_rows.txt\n";
        return 1;
    }

    std::string line;
    int rowNum = 0;
    double maxDiff = 0.0;

    while (std::getline(f, line))
    {
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        std::vector<Foam::scalar> x(6);
        for (int i = 0; i < 6; ++i) iss >> x[i];

        std::string bar;
        iss >> bar; // "|"

        double refTw, refTwSigma, refAngle, refAngleSigma;
        iss >> refTw >> refTwSigma >> refAngle >> refAngleSigma;

        iss >> bar; // "|"

        std::array<double, 5> refW;
        for (int i = 0; i < 5; ++i) iss >> refW[i];

        std::string trueLabel;
        iss >> trueLabel;

        auto pred = model.predict(x);

        double dTw = std::fabs(pred.twMean - refTw);
        double dTwSigma = std::fabs(pred.twSigma - refTwSigma);
        double dAngle = std::fabs(pred.angleMean - refAngle);
        double dAngleSigma = std::fabs(pred.angleSigma - refAngleSigma);

        double maxWDiff = 0.0;
        for (int i = 0; i < 5; ++i)
        {
            maxWDiff = std::max(maxWDiff, std::fabs(pred.weightsMean[i] - refW[i]));
        }

        double wSum = 0.0;
        for (int i = 0; i < 5; ++i) wSum += pred.weightsMean[i];

        std::cout << "row " << rowNum << " (true=" << trueLabel << "): "
                  << "tw=" << pred.twMean << " (ref " << refTw << ", diff " << dTw << ") "
                  << "twSigma=" << pred.twSigma << " (ref " << refTwSigma << ", diff " << dTwSigma << ") "
                  << "angle=" << pred.angleMean << " (ref " << refAngle << ", diff " << dAngle << ") "
                  << "angleSigma=" << pred.angleSigma << " (ref " << refAngleSigma << ", diff " << dAngleSigma << ") "
                  << "maxWeightDiff=" << maxWDiff << " "
                  << "weightSum=" << wSum << "\n";

        maxDiff = std::max({maxDiff, dTw, dTwSigma, dAngle, dAngleSigma, maxWDiff});
        ++rowNum;
    }

    std::cout << "\nChecked " << rowNum << " rows. Max abs diff vs Python reference: " << maxDiff << "\n";
    return 0;
}
