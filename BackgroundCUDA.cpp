#include "BackgroundCUDA.h"
#include "RCR.h"
#include "Tools.h"
#include "Debugger.h"
#include <iostream>
#include <future>

std::vector<double> clean(Baseline &baseline)
{
    baseline.rejectPoints();
    baseline.returnPoints();
    baseline.setLocalModels();
    std::vector<double> results = baseline.getResults();
    return results;
}

Baseline::Baseline(bool forward, double scatter, double xAnchor, double yAnchor, std::vector<bool> &checks, std::vector<double> &flux, std::vector<double> &dataDumps, std::vector<double> &angDist)
{
    this->forward = forward;
    this->scatter = scatter;
    this->xAnchor = xAnchor;
    this->yAnchor = yAnchor;
    this->BLChecks = checks;
    this->BLFlux = flux;
    this->BLDataDumps = dataDumps;
    this->BLAngDist = angDist;

    BLResults.reserve(BLFlux.size() * 3 * 30);
}
Baseline::Baseline()
{
}
Baseline::~Baseline()
{
}

int Baseline::sufficentPointCheck()
{
    double distinctXHold;
    double distinctYHold;
    int checkCountX, checkCountY;
    int checkCountMin;
    bool testBoolX, testBoolY;
    std::vector<double> distinctXValues, distinctYValues;
    distinctXValues.reserve(4);
    distinctYValues.reserve(4);
    checkCountX = 0, checkCountY = 0;
    // distinctXValues.push_back(xAnchor);
    // COLLECT ORIGINAL FLAGS FOR SUB-SECTION OF POINTS
    for (int i = 0; i < BLChecks.size(); i++)
    {
        distinctXHold = BLAngDist[i];
        distinctYHold = BLFlux[i];
        testBoolX = findDuplicate(distinctXHold, distinctXValues);
        testBoolY = findDuplicate(distinctYHold, distinctYValues);
        if ((BLChecks[i] == true) && (testBoolX == false) && (checkCountX < 4))
        {
            distinctXValues.push_back(distinctXHold);
            checkCountX = distinctXValues.size();
        }
        if ((BLChecks[i] == true) && (testBoolY == false) && (checkCountY < 4))
        {
            distinctYValues.push_back(distinctYHold);
            checkCountY = distinctYValues.size();
        }
    }

    checkCountX = distinctXValues.size();
    checkCountY = distinctYValues.size();

    checkCountMin = std::min(checkCountX, checkCountY);

    return checkCountMin;
}
bool Baseline::findDuplicate(double value, std::vector<double> vectorHold)
{
    bool duplicate = false;
    for (int i = 0; i < vectorHold.size(); i++)
    {
        if (value == vectorHold[i])
        {
            duplicate = true;
        }
    }
    return duplicate;
}

std::vector<double> Baseline::autoFixedWRegression()
{
    // ATTEMPTS QUADRATIC REGRESSION
    // IF A COEFFICENT IS UNSTABLE RESORST TO LINEAR
    // IF LINEAR FALIS NO MODEL
    // LAST INDEX OF toRet CONTAINS MODEL ORDER

    // std::vector<double> quadCoef, linCoef;
    // std::vector<double> toRet;
    std::vector<double> linCoef, quadCoef, nanQuad, nanLin;
    std::vector<double> toRet(5);
    bool linearModel = false;
    bool noModel = false;
    bool sufficentPoints = false;
    int trueCount = 0;
    quadCoef.resize(3, 0);
    linCoef.resize(2, 0);
    nanQuad.resize(3, NAN);
    nanLin.resize(2, NAN);

    // linCoef = fixedRegressionPivot();
    quadCoef = Tools::fixedQuadraticRegressionPivot(BLChecks, BLDataDumps, BLAngDist, BLFlux, xAnchor, yAnchor);

    int i = 0;

    while (i < 3)
    {
        if (quadCoef[i] != quadCoef[i])
        {
            quadCoef = nanQuad;
            linearModel = true;
            linCoef = Tools::fixedRegressionPivot(BLChecks, BLDataDumps, BLAngDist, BLFlux, xAnchor, yAnchor);
            break;
        }
        i++;
    }

    i = 0;
    while (i < 2)
    {
        if (linCoef[i] != linCoef[i])
        {
            noModel = true;
            linCoef = nanLin;
            break;
        }
        i++;
    }

    if (linearModel || noModel)
    {
        trueCount = sufficentPointCheck();

        if (trueCount >= 3 && linearModel)
        {
            Debugger::print("Warn", "Quad Fixed Regression Failed in a Worrysome Way");
        }

        if (trueCount >= 2 && noModel)
        {
            Debugger::print("Warn", "All Fixed Regression Failed in a Worrysome Way");
        }
    }

    toRet[0] = quadCoef[0];
    toRet[1] = quadCoef[1];
    toRet[2] = quadCoef[2];

    toRet[3] = linCoef[0];
    toRet[4] = linCoef[1];

    return toRet;
}
std::vector<double> Baseline::autoWRegression()
{
    // ATTEMPTS QUADRATIC REGRESSION
    // IF A COEFFICENT IS UNSTABLE RESORST TO LINEAR
    // IF LINEAR FALIS NO MODEL
    // LAST INDEX OF toRet CONTAINS MODEL ORDER
    std::vector<double> linCoef, quadCoef, nanQuad, nanLin;
    std::vector<double> toRet(5);
    bool linearModel = false;
    bool noModel = false;
    int trueCount = 0;
    int sufficentPoints = 0;
    quadCoef.resize(3, 0);
    linCoef.resize(2, 0);
    nanQuad.resize(3, NAN);
    nanLin.resize(2, NAN);

    // linCoef = regressionPivot();
    quadCoef = Tools::quadraticRegressionPivot(BLChecks, BLDataDumps, BLAngDist, BLFlux);

    int i = 0;

    while (i < 3)
    {
        if (quadCoef[i] != quadCoef[i])
        {
            quadCoef = nanQuad;
            linearModel = true;
            linCoef = Tools::regressionPivot(BLChecks, BLDataDumps, BLAngDist, BLFlux);
            break;
        }
        i++;
    }

    i = 0;
    while (i < 2)
    {
        if (linCoef[i] != linCoef[i])
        {
            noModel = true;
            linCoef = nanLin;
            break;
        }
        i++;
    }

    if (linearModel || noModel)
    {
        trueCount = sufficentPointCheck();

        if (trueCount >= 3 && linearModel)
        {
            Debugger::print("Warn", "Quad Fixed Regression Failed in a Worrysome Way");
        }

        if (trueCount >= 2 && noModel)
        {
            Debugger::print("Warn", "All Fixed Regression Failed in a Worrysome Way");
        }
    }

    toRet[0] = quadCoef[0];
    toRet[1] = quadCoef[1];
    toRet[2] = quadCoef[2];

    toRet[3] = linCoef[0];
    toRet[4] = linCoef[1];

    return toRet;
}
std::vector<double> Baseline::applyModel(double x, std::vector<double> &coefs)
{
    std::vector<double> modelVal(2);

    // FIRST ATTEMPT QUADRATIC
    // modelVal[0] = coefs[0] + coefs[1] * (x - x0) + coefs[2] * pow(x - x0, 2);
    modelVal[0] = coefs[0] + coefs[1] * (x) + coefs[2] * pow(x, 2);
    // Debugger::print("Info", coefs[0] + coefs[1] * (x - x0) + coefs[2] * pow(x - x0, 2), coefs[0] + coefs[1] * (x) + coefs[2] * pow(x, 2));
    modelVal[1] = 2.0;

    // THEN ATTEMPT LINEAR
    if (modelVal != modelVal)
    {
        // Debugger::print("Info", "Local Model: Failure to use Quadratic Coefficents");
        // modelVal[0] = coefs[3] + coefs[4] * (x - x0);
        modelVal[0] = coefs[3] + coefs[4] * (x);
        modelVal[1] = 1.0;
    }
    else if (modelVal != modelVal)
    {
        // Debugger::print("Info", "Local Model: Failure to use Linear Coefficents");
        modelVal[0] = NAN;
        modelVal[1] = 0.0;
    }

    return modelVal;
}
std::vector<double> Baseline::getResults()
{
    return this->BLResults;
}

double Baseline::rejectPoints()
{
    double largest, delta;
    double splus = 999999, splusTemp;
    int checkCount, counter, counterplus, largeIndex;
    std::vector<double> regressionLine, holderAVec;

    while (splus > scatter)
    {
        checkCount = sufficentPointCheck();
        regressionLine = autoFixedWRegression();

        // BASED ON THE MODEL PARAMETERS, REMOVE THE HIGHEST DEVIATING POINT AND REPEAT
        splus = 0.0;
        largest = 0.0;
        counter = 0;
        counterplus = 0;
        largeIndex = -1;

        for (int k = 0; k < BLChecks.size(); k++)
        {
            if (BLChecks[k])
            {
                counter++;
                holderAVec = applyModel(BLAngDist[k], regressionLine);
                delta = (BLFlux[k] - holderAVec[0]);
                if (BLFlux[k] > holderAVec[0])
                {
                    splus = splus + delta * delta;
                    counterplus++;
                }
                if (delta > largest)
                {
                    largest = delta;
                    largeIndex = k;
                }
            }
        }

        // HALT REJECTION IF THERE ARE ONLY THREE DISTINCT POINTS REMAINING
        if (checkCount >= (int)holderAVec[1])
        {
            splus = sqrt((splus / (counterplus)));
            splusTemp = splus;
            if ((splus > scatter) && (largeIndex > -1))
            {
                BLChecks[largeIndex] = false;
            }
        }
        else
        {
            splusTemp = sqrt((splus) / (counterplus));
            splus = 0;
        }
    }

    return splusTemp;
}
int Baseline::returnPoints()
{
    if (forward)
    {
        BLChecks[0] = false;
    }
    else
    {
        BLChecks[BLChecks.size() - 1] = false;
    }

    int minIndex, maxIndex;
    double smallest, delta;
    double splus = 0, splusTemp, splusSigmaDiff = 999999;
    int low = 999999, high = -999999, lastChecked;
    int checkCount, counter, counterplus, smallIndex;
    std::vector<double> regressionLine, holderAVec;
    while (splus < scatter)
    {
        // DETERMINE THE LOWEST AND HIGHEST NON-REJECTED POINTS
        for (int k = 0; k < BLChecks.size(); k++)
        {
            if (BLChecks[k])
            {
                if (k < low)
                {
                    low = k;
                }
                if (k > high)
                {
                    high = k;
                }
            }
        }

        checkCount = sufficentPointCheck();
        regressionLine = autoWRegression();

        splus = 0.0;
        smallest = 999999;
        counter = 0;
        counterplus = 0;
        smallIndex = -1;

        minIndex = ((low - 1) > 0) ? (low - 1) : 0;
        maxIndex = ((high + 2) < BLChecks.size() ? (high + 2) : BLChecks.size());
        for (int k = minIndex; k < maxIndex; k++)
        {
            // CALCULATE VALUES FOR THE MODEL

            holderAVec = applyModel(BLAngDist[k], regressionLine); // I THINK THIS IS ZERO INDEX BECAUSE THERE WAS NO REGRESSION POINT TO GO THROUGH.
            delta = BLFlux[k] - holderAVec[0];
            // DETERMINE DEVIATION FROM MODEL AND SMALLEST *REJECTED* POINT
            if (BLChecks[k])
            {
                counter++;
                if (BLFlux[k] > holderAVec[0])
                {
                    splus = splus + delta * delta;
                    counterplus++;
                }
            }
            else
            {
                if (delta < smallest)
                {
                    smallest = delta;
                    smallIndex = k;
                }
            }
        }

        // IF MODEL HAS AT LEAST THREE POINTS, ADD VALUES BACK IN?

        splus = sqrt((splus / (counterplus)));
        if (((splus < scatter) || (counterplus == 0)) && (smallIndex > -1))
        {
            BLChecks[smallIndex] = true;
            splusSigmaDiff = std::abs(scatter - splus);
            lastChecked = smallIndex;
        }
        else
        {
            if (splusSigmaDiff < std::abs(scatter - splus))
            {
                BLChecks[lastChecked] = false;
            }
            splus = 999999;
        }
    }

    return counter;
}
void Baseline::setLocalModels()
{
    int low = 999999, high = -999999, counter = 0;
    double mean = 0, stDev = 0, kurtosis = 0, dumpSum = 0, bgWeight;
    std::vector<double> holderBVec, regressionLine;

    regressionLine = autoWRegression();

    for (int k = 0; k < BLChecks.size(); k++)
    {
        if (BLChecks[k])
        {
            if (k < low)
            {
                low = k;
            }
            if (k > high)
            {
                high = k;
            }
            counter++;
        }
    }

    for (int i = low; i < high + 1; i++)
    {
        if (BLChecks[i])
        {
            mean += BLDataDumps[i] * BLAngDist[i];
            dumpSum += BLDataDumps[i];
        }
    }
    mean = mean / dumpSum;

    for (int i = low; i < high + 1; i++)
    {
        if (BLChecks[i])
        {
            stDev += BLDataDumps[i] * pow(BLAngDist[i] - mean, 2);
            kurtosis += BLDataDumps[i] * pow(BLAngDist[i] - mean, 4);
        }
    }
    stDev = sqrt(stDev / dumpSum);
    kurtosis = pow(kurtosis / dumpSum, .25);

    for (int k = low; k < high + 1; k++)
    {
        holderBVec = applyModel(BLAngDist[k], regressionLine);

        if (!((holderBVec[0] > BLFlux[k]) && (BLChecks[k] == 0)))
        {
            if (holderBVec[0] == NAN || holderBVec[0] != holderBVec[0] || stDev == 0 || kurtosis == 0)
            {
                Debugger::print("Info", "Bad Local Model");
                bgWeight = NAN;
            }
            else if (counter == 1)
            {
                bgWeight = 1.0;
            }
            else
            {
                bgWeight = (dumpSum / (1 + pow((BLAngDist[k] - mean) / stDev, 2) + pow((BLAngDist[k] - mean) / kurtosis, 4)));
            }

            BLResults.push_back(k);
            BLResults.push_back(holderBVec[0]);
            BLResults.push_back(bgWeight);
        }
    }
}

BackgroundCUDA::BackgroundCUDA()
{
}
BackgroundCUDA::BackgroundCUDA(Scan &scan, bool lssElevation)
{
    size = scan.getSize();

    angDist = scan.getAngDist();
    if (lssElevation)
    {
        flux = scan.getElevation();
        scatter = scan.getElevationScatter();
    }
    else
    {
        flux = scan.getFlux();
        scatter = scan.getScatter();
    }

    dataDumps = scan.getDataDumps();
}
BackgroundCUDA::BackgroundCUDA(Spectra &spectra, int i)
{
    // Spectral Constructor
    baselineVec = spectra.getBaseline(i);
    angDist = spectra.getFreqDist(i);
    dataDumps = spectra.getWeight(i);
    scatter = spectra.getScatter(i);
    flux = spectra.getFlux(i);

    size = flux.size();
}
BackgroundCUDA::BackgroundCUDA(std::vector<double> baselines, std::vector<double> distances, std::vector<double> weights, std::vector<double> data, double scatter)
{
	// Spectral Constructor
	baselineVec = baselines;
	angDist = distances;
	dataDumps = weights;
	flux = data;
	this->scatter = scatter;

	size = flux.size();
}
BackgroundCUDA::~BackgroundCUDA()
{
}

std::vector<double> BackgroundCUDA::calculateBG(double baseline)
{
    Baseline baselineTemp;
    int start;
    std::vector<int> indicesForward, indicesBackward;
    std::vector<double> baselineLM, dubFiller, bg;
    std::vector<std::vector<double>> baselineResults;
    std::vector<std::vector<int>> indices;

    bg.resize(size, 999999);
    bgData.resize(size, dubFiller);
    bgWeights.resize(size, dubFiller);

    indices = findStartEndIndices(angDist, baseline);

    indicesForward = indices[0];
    indicesBackward = indices[1];

    baselineArray.reserve(indicesForward.size() / 2 + indicesBackward.size() / 2);
    for (int i = 0; i < bgData.size(); i++)
    {
        bgData[i].reserve(indicesForward.size() / 2 + indicesBackward.size() / 2);
        bgWeights[i].reserve(indicesForward.size() / 2 + indicesBackward.size() / 2);
    }

    buildBaselines(true, indicesForward);
    buildBaselines(false, indicesBackward);

    baselineResults.resize(baselineArray.size());

    std::vector<std::vector<double>> cleanBaselineArray;
    cleanBaselineArray.resize(baselineArray.size());

    for (int i = 0; i < baselineArray.size(); i++)
    {
        cleanBaselineArray[i] = clean(std::ref(baselineArray[i]));
    }
    for (int i = 0; i < indicesForward.size() / 2; i++)
    {
        start = indicesForward[2 * i];
        baselineLM = cleanBaselineArray[i];
        loadLocalModels(baselineLM, start);
    }
    for (int i = 0; i < indicesBackward.size() / 2; i++)
    {
        start = indicesBackward[2 * i];
        baselineLM = cleanBaselineArray[(indicesForward.size() / 2) + i];
        loadLocalModels(baselineLM, start);
    }

    bg = setBackground(bgData, bgWeights);

    return bg;
}
std::vector<double> BackgroundCUDA::calculateBGMulti(double baseline)
{
    Baseline baselineTemp;
    int start;
    std::vector<int> indicesForward, indicesBackward;
    std::vector<double> baselineLM, dubFiller, bg;
    std::vector<std::future<std::vector<double>>> futureVec;
    std::vector<std::vector<double>> baselineResults;
    std::vector<std::vector<int>> indices;

    bg.resize(size, 999999);
    bgData.resize(size, dubFiller);
    bgWeights.resize(size, dubFiller);
    baselineVec.resize(size, baseline);

    indices = findStartEndIndices(angDist, baseline);

    indicesForward = indices[0];
    indicesBackward = indices[1];

    baselineArray.reserve(indicesForward.size() / 2 + indicesBackward.size() / 2);
    for (int i = 0; i < bgData.size(); i++)
    {
        bgData[i].reserve(indicesForward.size() / 2 + indicesBackward.size() / 2);
        bgWeights[i].reserve(indicesForward.size() / 2 + indicesBackward.size() / 2);
    }

    buildBaselines(true, indicesForward);
    buildBaselines(false, indicesBackward);

    futureVec.resize(baselineArray.size());
    baselineResults.resize(baselineArray.size());

    for (int i = 0; i < baselineArray.size(); i++)
    {
        futureVec[i] = std::async(std::launch::async, clean, std::ref(baselineArray[i]));
    }
    for (int i = 0; i < indicesForward.size() / 2; i++)
    {
        start = indicesForward[2 * i];
        baselineLM = futureVec[i].get();
        loadLocalModels(baselineLM, start);
    }
    for (int i = 0; i < indicesBackward.size() / 2; i++)
    {
        start = indicesBackward[2 * i];
        baselineLM = futureVec[(indicesForward.size() / 2) + i].get();
        loadLocalModels(baselineLM, start);
    }

    bg = setBackground(bgData, bgWeights);

    return bg;
}
void BackgroundCUDA::buildBaselines(bool forward, std::vector<int> &indices)
{
    int start, end;
    double xAnchor, yAnchor;
    std::vector<bool> checksTemp;
    std::vector<double> angDistTemp, dumpsTemp;
    std::vector<double> zTemp;
    std::vector<double> results;
    Baseline baselineTemp;

    for (int i = 0; i < indices.size() / 2; i++)
    {
        start = indices[2 * i];
        end = indices[2 * i + 1];
        checksTemp.resize(end - start + 1, true);
        angDistTemp.resize(end - start + 1);
        zTemp.resize(end - start + 1);
        dumpsTemp.resize(end - start + 1);

        if (forward)
        {
            xAnchor = 0;
            yAnchor = flux[start];
        }
        else
        {
            xAnchor = angDist[end] - angDist[start];
            yAnchor = flux[end];
        }

        // STORE INFORMATION FOR DATA POINTS WITHIN DESIRED BASELINE
        for (int i = 0; i < angDistTemp.size(); i++)
        {
            angDistTemp[i] = angDist[start + i] - angDist[start];
            zTemp[i] = flux[start + i];
            dumpsTemp[i] = dataDumps[start + i];
        }

        baselineTemp = Baseline(forward, scatter, xAnchor, yAnchor, checksTemp, zTemp, dumpsTemp, angDistTemp);
        baselineArray.push_back(baselineTemp);
    }
}
void BackgroundCUDA::loadLocalModels(std::vector<double> &results, int start)
{
    int index;
    double LMValue, LMWeight;
    for (int i = 0; i < results.size() / 3; i++)
    {
        index = results[3 * i] + start;
        LMValue = results[3 * i + 1];
        LMWeight = results[3 * i + 2];
        if ((LMValue == LMValue) && (LMWeight == LMWeight))
        {
            bgData[index].push_back(LMValue);
            bgWeights[index].push_back(LMWeight);
        }
    }
}
std::vector<std::vector<int>> BackgroundCUDA::findStartEndIndices(std::vector<double> &angDist, double baseline)
{
    int size = angDist.size();
    int end = 0, start = 0;
    std::vector<int> intFiller;
    std::vector<std::vector<int>> indices;

    indices.resize(2, intFiller);
    indices[0].reserve(size);
    indices[1].reserve(size);

    for (; start < size - 2; start++)
    {
        while (end < size - 1 && angDist[end] - angDist[start] < baselineVec[start])
        {
            end++;
        }
        end--;
        if (end == size - 2 && angDist[end + 1] - angDist[start] < baselineVec[start])
        {
            end++;
        }
        // CALCULATE THE LOCAL MODEL VALUES
        if (end != start && baselineVec[start] != 0)
        {
            indices[0].push_back(start);
            indices[0].push_back(end);
        }
    }
    end = size - 1;
    start = size - 1;
    for (; end > 1; end--)
    {
        while (angDist[end] - angDist[start] < baselineVec[end] && start > 0)
        {
            start--;
        }
        start++;
        if (start == 1 && angDist[end] - angDist[start - 1] < baselineVec[end])
        {
            start--;
        }
        if (start != end && baselineVec[end] != 0)
        {
            indices[1].push_back(start);
            indices[1].push_back(end);
        }
    }

    return indices;
}
std::vector<double> BackgroundCUDA::setBackground(std::vector<std::vector<double>> &bgData, std::vector<std::vector<double>> &bgWeights)
{
    double low, high;
    int trueCount;
    int size = bgData.size();
    std::vector<double> data, weights, bg;
    bg.resize(size);

    RCR rcr = RCR(LS_MODE_DL);
    for (int i = 0; i < bgData.size(); i++)
    {
        trueCount = bgData[i].size();
        if (trueCount == 0)
        {
            bg[i] = 999999;
        }
        else if (trueCount == 1)
        {
            bg[i] = bgData[i][0];
        }
        else
        {
            rcr.performBulkRejection(bgWeights[i], bgData[i]);
            bg[i] = rcr.result.mu;
        }
        if (baselineVec[i] == 0)
        {
            bg[i] = 0.0;
        }
        /*
        else
        {
            for (int j = 0; j < bgData[i].size(); j++)
            {
                if (bgData[i][j] == NAN || bgData[i][j] != bgData[i][j] || bgWeights[i][j] == NAN || bgWeights[i][j] != bgWeights[i][j])
                {
                    trueCount--;
                }
                else
                {
                    data.push_back(bgData[i][j]);
                    weights.push_back(bgWeights[i][j]);
                }
            }
            if (trueCount == 0)
            {
                mean = 999999;
            }
            else if (trueCount == 1)
            {
                mean = bgData[i][0];
            }
            else
            {
                rcr.performBulkRejection(weights, data);
                mean = rcr.result.mu;
            }
            data.clear();
            weights.clear();
        }
        if (trueCount > 0)
        {
            bg[i] = mean;
        }
        else
        {
            bg[i] = 999999;
        }
        */
    }

    for (int k = 0; k < size; k++)
    {
        if (bg[k] == 999999)
        {
            low = Tools::max(0, k - 1);
            high = Tools::min(k + 1, size - 1);
            while (bg[low] == 999999 && low > 1)
            {
                low--;
            }
            while (bg[high] == 999999 && high < size - 1)
            {
                high++;
            }
            if (bg[low] != 999999 && bg[high] != 999999)
            {
                bg[k] = bg[low] + (bg[high] - bg[low]) / (high - low) * (k - low);
            }
            else if (bg[low] == 999999)
            {
                bg[low] = bg[high];
            }
            else if (bg[high] == 999999)
            {
                bg[high] = bg[low];
            }
        }
    }

    return bg;
}