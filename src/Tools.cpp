#define _USE_MATH_DEFINES
#include "Tools.h"
#include "RCR.h"
#include <fftw3.h>
#include <float.h>
#include <math.h>


//correlation
std::vector<std::vector<double> > ifft(std::vector<std::vector<double> > signal) {
    std::vector<std::vector<double> > toRet;
    std::vector<double> dubFiller;
    int size = signal.size();
    fftw_complex *data, *ifft_result;
    fftw_plan plan_backward;
    data = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * size);
    ifft_result = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * size);
    plan_backward = fftw_plan_dft_1d(size, data, ifft_result, FFTW_BACKWARD, FFTW_MEASURE);
    for (int i = 0; i < size; i++) {
        data[i][0] = signal[i][0];
        data[i][1] = signal[i][1];
    }
    fftw_execute(plan_backward);
    fftw_destroy_plan(plan_backward);
    fftw_free(data);
    for (int i = 0; i < size; i++) {
        toRet.push_back(dubFiller);
        toRet[i].push_back(ifft_result[i][0]);
        toRet[i].push_back(ifft_result[i][1]);
    }
    fftw_free(ifft_result);
    return toRet;
}

std::vector<std::vector<double> > fft(std::vector<double> signal) {
    int size = signal.size();
    std::vector<std::vector<double> > toRet;
    std::vector<double> dubFiller;
    double delta = 0;
    fftw_complex *data, *fft_result;
    fftw_plan plan_forward;
    data = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * size);
    fft_result = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * size);
    plan_forward = fftw_plan_dft_1d(size, data, fft_result, FFTW_FORWARD, FFTW_MEASURE);
    for (int i = 0; i < size; i++) {
        data[i][0] = signal[i];
        data[i][1] = 0.0;
    }
    fftw_execute(plan_forward);
    fftw_destroy_plan(plan_forward);
    fftw_free(data);
    toRet.resize(size, dubFiller);
    for (int i = 0; i < size; i++) {
        toRet[i].push_back(fft_result[i][0]);
        toRet[i].push_back(fft_result[i][1]);
    }
    fftw_free(fft_result);
    return toRet;
}

std::vector<std::vector<double> > correlate(std::vector<std::vector<double> > f, std::vector<std::vector<double> > g) {
    std::vector<double> fvec, gvec, dubFiller;
    std::vector<std::vector<double> > toRet;
    dubFiller.resize(2);
    toRet.resize(f.size(), dubFiller);
    for (int i = 0; i < f.size(); i++) {
        f[i][1] = -1.0 * f[i][1];
    }
    for (int i = 0; i < f.size(); i++) {
        fvec = f[i];
        gvec = g[i];
        toRet[i][0] = fvec[0] * gvec[0] - fvec[1] * gvec[1];
        toRet[i][1] = fvec[0] * gvec[1] + fvec[1] * gvec[0];
    }
    return toRet;
}

double shiftCalculator(std::vector<std::vector<double> > maxIfftInfo) {
    RCR rcr = RCR(LS_MODE_68);

    bool stop = false, first = true, alreadySorted = false;
    int maxIndex, counter, maxDevID;
    double medianAbove, medianBelow, stDevAbove, stDevBelow, belowAboveRatio, maxBelowAboveRatio = -999999;
    double maxBelowAboveIndex, max, median, stDev, maxDev, hold, shiftWeightSum = 0.0, weightSum = 0.0;
    std::vector<int> checks;
    std::vector<double> swap, shiftsAbove, shiftsBelow, weightTemps, shifts, weights, tempArray, tempWeights, tempArray2
            , tempWeights2;
    for (int j = 0; j < maxIfftInfo.size(); j++) {
        maxIfftInfo[j][0] = -1.0 * pow(-1.0, j) * maxIfftInfo[j][0];
    }
    for (int i = 0; i < maxIfftInfo.size() - 1; i++) {
        for (int j = 0; j < maxIfftInfo.size() - i - 1; j++) {
            if (maxIfftInfo[j][1] > maxIfftInfo[j + 1][1]) /* For decreasing order use < */
            {
                swap = maxIfftInfo[j];
                maxIfftInfo[j] = maxIfftInfo[j + 1];
                maxIfftInfo[j + 1] = swap;
            }
        }
    }
    for (int i = 2; i < maxIfftInfo.size() - 2; i++) {
        for (int j = 0; j < i; j++) {
            shiftsBelow.push_back(maxIfftInfo[j][0]);
        }
        for (int j = i + 1; j < maxIfftInfo.size(); j++) {
            shiftsAbove.push_back(maxIfftInfo[j][0]);
        }
        rcr.performRejection(shiftsBelow);
        medianBelow = rcr.result.mu;
        rcr.performRejection(shiftsAbove);
        medianAbove = rcr.result.mu;
        for (int j = 0; j < shiftsBelow.size(); j++) {
            shiftsBelow[j] -= medianBelow;
        }
        for (int j = 0; j < shiftsAbove.size(); j++) {
            shiftsAbove[j] -= medianAbove;
        }
        stDevBelow = Tools::get68th(shiftsBelow);
        stDevAbove = Tools::get68th(shiftsAbove);
        belowAboveRatio = stDevBelow / stDevAbove;
        if (belowAboveRatio > maxBelowAboveRatio) {
            maxBelowAboveRatio = belowAboveRatio;
            maxBelowAboveIndex = i;
        }
        shiftsAbove.clear();
        shiftsBelow.clear();
    }
    for (int j = maxBelowAboveIndex + 1; j < maxIfftInfo.size(); j++) {
        shiftsAbove.push_back(maxIfftInfo[j][0]);
        weightTemps.push_back(maxIfftInfo[j][1]);
    }

    rcr.performRejection(shiftsAbove);

    for (int i = 0; i < checks.size(); i++) {
        if (checks[i]) {
            shifts.push_back(shiftsAbove[i]);
            weights.push_back(weightTemps[i]);
        }
    }
    shifts = rcr.result.cleanY;
    for (int i = 0; i < shifts.size(); i++) {
        shiftWeightSum += shifts[i] * weights[i];
        weightSum += weights[i];
    }
    return shiftWeightSum / weightSum;
}


//sorting
void swapa(int a, int b, std::vector<double> &y) {
    double tmp;
    tmp = y[a];
    y[a] = y[b];
    y[b] = tmp;
}

void swapCircle(int a, int b, std::vector<double> &y) {
    double tmp;
    tmp = y[a];
    y[a] = y[b];
    y[b] = tmp;
}

void swapCircle(int a, int b, std::vector<int> &y) {
    double tmp;
    tmp = y[a];
    y[a] = y[b];
    y[b] = tmp;
}

void swapCircle(int a, int b, std::vector<bool> &y) {
    double tmp;
    tmp = y[a];
    y[a] = y[b];
    y[b] = tmp;
}

void QSa(int left, int right, std::vector<double> &y) {
    int i = left, j = right;
    double pivot = y[(left + right) / 2];

    while (i <= j) {
        while (y[i] < pivot) {
            i++;
        }
        while (y[j] > pivot) {
            j--;
        }
        if (i <= j) {
            swapa(i, j, y);
            i++;
            j--;
        }
    };

    if (left < j) {
        QSa(left, j, y);
    }
    if (i < right) {
        QSa(i, right, y);
    }
}


//linear algebra
std::vector<double> lowerTriangleSolver(int columnCount, std::vector<double> A, std::vector<double> b) {
    std::vector<double> x;
    double sub = 0.0;
    x.resize(columnCount, 0.0);
    for (int i = 0; i < columnCount; i++) {
        sub = 0.0;
        for (int j = 0; j < i + 1; j++) {
            sub += A[columnCount * i + j] * x[j];
        }
        x[i] = (b[i] - sub) / A[columnCount * i + i];
    }
    return x;
}

std::vector<double> performPivot(int columnCount, std::vector<double> A, std::vector<double> b) {
    int pivotIndex;
    double pivot, swap, multiplier;
    std::vector<double> coefficents;

    for (int i = columnCount - 1; i > 0; i--) {
        pivot = std::abs(A[columnCount * i + i]);
        pivotIndex = i;
        for (int j = i; j >= 0; j--) {
            if (std::abs(A[columnCount * j + i]) > pivot) {
                pivot = std::abs(A[columnCount * j + i]);
                pivotIndex = j;
            }
        }
        if (pivotIndex == i) {
            for (int j = 0; j < columnCount; j++) {
                A[columnCount * i + j] = A[columnCount * i + j] / pivot;
            }
            b[i] = b[i] / pivot;
        } else {
            for (int j = 0; j < columnCount; j++) {
                swap = A[columnCount * i + j];
                A[columnCount * i + j] = A[columnCount * pivotIndex + j] / pivot;
                A[columnCount * pivotIndex + j] = swap;
            }
            swap = b[i];
            b[i] = b[pivotIndex] / pivot;
            b[pivotIndex] = swap;
        }
        for (int j = i - 1; j >= 0; j--) {
            multiplier = A[j * columnCount + i] / A[i * columnCount + i];
            for (int k = 0; k < columnCount; k++) {
                A[j * columnCount + k] -= multiplier * A[i * columnCount + k];
            }
            b[j] -= multiplier * b[i];
        }
    }

    coefficents = lowerTriangleSolver(columnCount, A, b);

    return coefficents;
}


namespace Tools {
    //sorters
    void sort(std::vector<double> &y) {
        QSa(0, y.size() - 1, y);
    }

    void sortAll(int left, int right, std::vector<double> &a, std::vector<double> &b, std::vector<double> &c) {
        int i = left, j = right;
        double pivot = a[(left + right) / 2];

        while (i <= j) {
            while (a[i] < pivot) {
                i++;
            }
            while (a[j] > pivot) {
                j--;
            }
            if (i <= j) {
                swapCircle(i, j, a);
                swapCircle(i, j, b);
                swapCircle(i, j, c);
                i++;
                j--;
            }
        };

        if (left < j) {
            sortAll(left, j, a, b, c);
        }
        if (i < right) {
            sortAll(i, right, a, b, c);
        }
    }

    void sortAll(int left, int right, std::vector<double> &a, std::vector<int> &b, std::vector<int> &c) {
        int i = left, j = right;
        double pivot = a[(left + right) / 2];

        while (i <= j) {
            while (a[i] < pivot) {
                i++;
            }
            while (a[j] > pivot) {
                j--;
            }
            if (i <= j) {
                swapCircle(i, j, a);
                swapCircle(i, j, b);
                swapCircle(i, j, c);
                i++;
                j--;
            }
        };

        if (left < j) {
            sortAll(left, j, a, b, c);
        }
        if (i < right) {
            sortAll(i, right, a, b, c);
        }
    }

    void sortAll(int left, int right, std::vector<int> &a, std::vector<int> &b, std::vector<bool> &c) {
        int i = left, j = right;
        double pivot = a[(left + right) / 2];

        while (i <= j) {
            while (a[i] < pivot) {
                i++;
            }
            while (a[j] > pivot) {
                j--;
            }
            if (i <= j) {
                swapCircle(i, j, a);
                swapCircle(i, j, b);
                swapCircle(i, j, c);
                i++;
                j--;
            }
        };

        if (left < j) {
            sortAll(left, j, a, b, c);
        }
        if (i < right) {
            sortAll(i, right, a, b, c);
        }
    }

    // void bubbleSort2(std::vector<double> &sort, std::vector<int> &piggyBack1, std::vector<double> &piggyBack2, std::vector<int> &piggyBack3, std::vector<int> &piggyBack4, std::vector<double> &piggyBack5)
    // {
    // 	int intSwap;
    // 	double swap;
    // 	for (int i = 0; i < sort.size() - 1; i++)
    // 	{
    // 		for (int j = 0; j < sort.size() - i - 1; j++)
    // 		{
    // 			if (sort[j] > sort[j + 1]) // For decreasing order use <
    // 			{
    // 				intSwap = piggyBack1[j];
    // 				piggyBack1[j] = piggyBack1[j + 1];
    // 				piggyBack1[j + 1] = intSwap;
    //
    // 				intSwap = piggyBack3[j];
    // 				piggyBack3[j] = piggyBack3[j + 1];
    // 				piggyBack3[j + 1] = intSwap;
    //
    // 				intSwap = piggyBack4[j];
    // 				piggyBack4[j] = piggyBack4[j + 1];
    // 				piggyBack4[j + 1] = intSwap;
    //
    // 				swap = piggyBack2[j];
    // 				piggyBack2[j] = piggyBack2[j + 1];
    // 				piggyBack2[j + 1] = swap;
    //
    // 				swap = piggyBack5[j];
    // 				piggyBack5[j] = piggyBack5[j + 1];
    // 				piggyBack5[j + 1] = swap;
    //
    // 				swap = sort[j];
    // 				sort[j] = sort[j + 1];
    // 				sort[j + 1] = swap;
    // 			}
    // 		}
    // 	}
    // }
    // void bubbleSort3(std::vector<std::vector<double> > &sort)
    // {
    // 	int intSwap;
    // 	std::vector<double> swap;
    // 	for (int i = 0; i < sort.size() - 1; i++)
    // 	{
    // 		for (int j = 0; j < sort.size() - i - 1; j++)
    // 		{
    // 			if (sort[j][0] > sort[j + 1][0]) // For decreasing order use <
    // 			{
    // 				swap = sort[j];
    // 				sort[j] = sort[j + 1];
    // 				sort[j + 1] = swap;
    // 			}
    // 		}
    // 	}
    // }

    //coordinate conversions
    double convertToB(double ra, double dec) {
        double raInRad, decInRad, B;
        double toRad = M_PI / 180.0, toDeg = 180.0 / M_PI;
        double deltaG = 27.1284 * toRad;
        double alphaG = 192.8595 * toRad;

        // CONVERT RA,DEC TO RADIANS
        raInRad = ra * toRad; // IN RADIANS
        decInRad = dec * toRad; // IN RADIANS

        // CONVERT DEC TO B
        double x = sin(deltaG) * sin(decInRad);
        double y = cos(deltaG) * cos(decInRad) * cos(raInRad - alphaG);
        B = asin(x + y) * toDeg;

        return B; // IN DEGREES
    }

    double convertToL(double ra, double dec) {
        double raInRad, decInRad, L;
        double toRad = M_PI / 180.0, toDeg = 180.0 / M_PI;
        double deltaG = 27.1284 * toRad;
        double alphaG = 192.8595 * toRad;

        // CONVERT RA,DEC TO RADIANS
        raInRad = ra * toRad; // IN RADIANS
        decInRad = dec * toRad; // IN RADIANS

        // CONVERT RA TO L
        double numerator = cos(decInRad) * sin(raInRad - alphaG);
        double denominator = cos(deltaG) * sin(decInRad) - sin(deltaG) * cos(decInRad) * cos(raInRad - alphaG);
        // ATAN2 ADDS PI, BUT WE WANT TO SUBTRACT PI DUE TO HOW THE CONVERSION EQUATION WAS DEFINED.
        // YOU GET THE SAME ANGLE BUT MEASURED IN OPPOSITE DIRECTIONS. E.G., 99 DEG = -261 DEG BUT
        // WE NEED THE LATTER BY WAY OF DEFINITION, IN THIS CASE. NOT TRIVIAL!
        if (numerator > 0 && denominator < 0) {
            L = 122.9320 - (atan(numerator / denominator) - M_PI) * toDeg;
        } else {
            L = 122.9320 - atan2(numerator, denominator) * toDeg;
        }

        return L; // IN DEGREES
    }

    double convertToDec(double L, double B) {
        double lInRad, bInRad, DEC;
        double toRad = M_PI / 180.0, toDeg = 180.0 / M_PI;
        double deltaG = 27.1284 * toRad;

        // CONVERT L,B TO RADIANS
        lInRad = L * toRad; // IN RADIANS
        bInRad = B * toRad; // IN RADIANS

        // CONVERT B TO DEC
        double x = sin(deltaG) * sin(bInRad);
        double y = cos(deltaG) * cos(bInRad) * cos(122.9320 * toRad - lInRad);
        DEC = asin(x + y) * toDeg;

        return DEC; //IN DEGRES
    }

    double convertToRa(double L, double B) {
        double lInRad, bInRad, RA;
        double toRad = M_PI / 180.0, toDeg = 180.0 / M_PI;;
        double alphaG = 192.8595;
        double deltaG = 27.1284 * toRad;

        // CONVERT L,B TO RADIANS
        lInRad = L * toRad; // IN RADIANS
        bInRad = B * toRad; // IN RADIANS

        // CONVERT L TO RA
        double numerator = cos(bInRad) * sin(122.9320 * toRad - lInRad);
        double denominator = cos(deltaG) * sin(bInRad) - sin(deltaG) * cos(bInRad) * cos(122.9320 * toRad - lInRad);
        RA = atan2(numerator, denominator) * toDeg + alphaG;

        return RA; // IN DEGREES
    }

    //statistics tools
    double max(double a, double b) {
        return (a > b ? a : b);
    }

    double min(double a, double b) {
        return (a < b ? a : b);
    }

    double get68th(std::vector<double> &y) {
        double stDev = 0, totalSum = 0, runningSum, temp = 0, weightTemp = 0, sumCounter = 0;
        sort(y);
        for (int i = 0; i < y.size(); i++) {
            totalSum += 1.0;
        }
        if (y.size() > 1) {
            runningSum = 1.0 * .682689;
            while (runningSum < .682689 * totalSum) {
                sumCounter++;
                runningSum += 1.0 * .317311 + 1.0 * .682689;
            }
            if (sumCounter == 0) {
                stDev = y[0];
            } else {
                stDev = y[sumCounter - 1] + (.682689 * totalSum - (runningSum - (1.0 * .317311 + 1.0 * .682689))) / (
                            1.0 * .317311 + 1.0 * .682689) * (y[sumCounter] - y[sumCounter - 1]);
            }
        } else {
            stDev = y[0];
        }

        return stDev;
    }

    double getMean(std::vector<double> y) {
        double sum = 0;
        for (int i = 0; i < y.size(); i++) {
            sum += y[i];
        }
        return sum / y.size();
    }

    double getMean(std::vector<double> w, std::vector<double> y) {
        double wySum = 0, wSum = 0;
        for (int i = 0; i < y.size(); i++) {
            wySum += w[i] * y[i];
            wSum += w[i];
        }
        return wySum / wSum;
    }

    double getMedian(std::vector<double> &y) {
        sort(y);

        int high = floor(y.size() / 2);
        int low = high - 1;
        double runningSum = 0, median = 0;
        double totalSum = y.size();
        if (y.size() > 1) {
            if (y.size() % 2 == 0) {
                runningSum = y.size() / 2.0 + .5;
            } else {
                runningSum = y.size() / 2.0;
            }
            median = y[low] + (.5 * totalSum - runningSum + 1.0) * (y[high] - y[low]);
        } else {
            median = y[0];
        }
        return median;
    }

    double getStDev(double mu, double delta, std::vector<double> y) {
        double stDev = 0;
        for (int k = 0; k < y.size(); k++) {
            stDev = delta + pow(y[k] - mu, 2);
        }
        return sqrt(stDev / (y.size() - 1));
    }

    double getWStDev(double mu, double delta, std::vector<double> w, std::vector<double> y) {
        int size = w.size();
        double top = 0, wSum = 0, wSumSq = 0, weight;
        for (int i = 0; i < size; i++) {
            weight = w[i];
            top += weight * (y[i] - mu) * (y[i] - mu);
            wSum += weight;
            wSumSq += weight * weight;
        }
        return sqrt(top / (wSum - delta * wSumSq / wSum));
    }

    double erfc(double x) {
        x = x / sqrt(2.0);
        return 1.0 / pow(
                   (1 + x * (.0705230784 + x * (.0422820123 + x * (
                                                    .0092705272 + x * (
                                                        .0001520143 + x * (.0002765672 + .0000430638 * x)))))), 16);
    }


    //math tools
    double signedArcTan(double top, double bottom) {
        if (bottom > 0) {
            if (top > 0) {
                return atan(top / bottom);
            } else if (top < 0) {
                return atan(top / bottom) + 2.0 * M_PI;
            } else {
                return 0.0;
            }
        } else if (bottom < 0) {
            if (top > 0) {
                return atan(top / bottom) + M_PI;
            } else if (top < 0) {
                return atan(top / bottom) + M_PI;
            } else {
                return M_PI;
            }
        } else if (bottom == 0) {
            if (top > 0) {
                return M_PI / 2.0;
            } else if (top < 0) {
                return M_PI * 3.0 / 2.0;
            } else {
                //implies same point, must figure out clever way to handle this
            }
        }
        return 0.0;
    }

    double arcCos(double value) {
        if (value >= 1.0) {
            value = .99999999999999;
        } else if (value <= -1.0) {
            value = -.99999999999999;
        }
        return acos(value);
    }

    std::vector<std::vector<double> > pivotSystem(int columnCount, std::vector<double> &A, std::vector<double> &b) {
        int pivotIndex;
        double pivot, swap, multiplier;
        std::vector<std::vector<double> > tensorHold;

        for (int i = columnCount - 1; i > 0; i--) {
            pivot = std::abs(A[columnCount * i + i]);
            pivotIndex = i;
            for (int j = i; j >= 0; j--) {
                if (std::abs(A[columnCount * j + i]) > pivot) {
                    pivot = std::abs(A[columnCount * j + i]);
                    pivotIndex = j;
                }
            }
            if (pivotIndex == i) {
                for (int j = 0; j < columnCount; j++) {
                    A[columnCount * i + j] = A[columnCount * i + j] / pivot;
                }
                b[i] = b[i] / pivot;
            } else {
                for (int j = 0; j < columnCount; j++) {
                    swap = A[columnCount * i + j];
                    A[columnCount * i + j] = A[columnCount * pivotIndex + j] / pivot;
                    A[columnCount * pivotIndex + j] = swap;
                }
                swap = b[i];
                b[i] = b[pivotIndex] / pivot;
                b[pivotIndex] = swap;
            }
            for (int j = i - 1; j >= 0; j--) {
                multiplier = A[j * columnCount + i] / A[i * columnCount + i];
                for (int k = 0; k < columnCount; k++) {
                    A[j * columnCount + k] -= multiplier * A[i * columnCount + k];
                }
                b[j] -= multiplier * b[i];
            }
        }

        tensorHold.resize(2);
        tensorHold[0].resize(A.size());
        tensorHold[1].resize(b.size());

        tensorHold[0] = A;
        tensorHold[1] = b;

        return tensorHold;
    }

    std::vector<std::vector<double> > crossCorrelate(std::vector<double> f, std::vector<double> g) {
        return ifft(correlate(fft(f), fft(g)));
    }

    std::vector<double> matrixSolver(int columnCount, std::vector<double> A, std::vector<double> b) {
        std::vector<double> coefficents;
        coefficents = performPivot(columnCount, A, b);
        return coefficents;
    }

    /*
    //regression without pivoting
    std::vector<double> fixedRegression(int start, int end, std::vector<bool> &checks, std::vector<double> &weights, std::vector<double> &horizAxis, std::vector<double> &vertAxis, double horizAxisCoordinate, double vertAxisCoordinate)
    {
        double xx = 0, xy = 0;
        double m, b;
        double horizHold, vertHold;
        std::vector<double> ret;
        ret.resize(2);

        for (int i = start; i <= end; i++)
        {
            if (checks[i])
            {
                //horizHold = horizAxis[i] - horizAxisCoordinate;
                //vertHold = vertAxis[i] - vertAxisCoordinate;
                horizHold = horizAxis[i] - horizAxis[0];
                vertHold = vertAxis[i] - vertAxis[0];

                xx += weights[i] * horizHold*horizHold;
                xy += weights[i] * horizHold*vertHold;
            }
        }

        m = xy / xx;
        //b = vertAxisCoordinate - m*horizAxisCoordinate;
        b = vertAxis[0] - m*horizAxis[0];

        ret[0] = b;
        ret[1] = m;

        return ret;
    }
    std::vector<double> fixedQuadraticRegression(int start, int end, std::vector<bool> &checks, std::vector<double> &weights, std::vector<double> &horizAxis, std::vector<double> &vertAxis, double xAnchor, double yAnchor)
    {
        double w = 0;
        double xSum = 0, ySum = 0;
        double pivot, swap, multiplier;
        double finalAnswer;
        int pivotIndex;
        int columnCount = 3;
        std::vector<double> A;
        std::vector<double> b;
        A.resize(columnCount * columnCount, 0.0);
        b.resize(columnCount, 0.0);

        std::vector<double> ret;
        ret.resize(3);
        double x2Sum = 0, x3Sum = 0, x4Sum = 0, yxSum = 0, yx2Sum = 0, error = 0, horizHold, vertHold;
        double wHold, xHold, yHold;
        double wx = 0, wxx = 0, wxxx = 0, wxxxx = 0;
        double wy = 0, wyx = 0, wyxx = 0;
        for (int i = start; i <= end; i++)
        {
            if (checks[i])
            {
                wHold = weights[i];
                xHold = (horizAxis[i] - xAnchor);
                yHold = (vertAxis[i] - yAnchor);
                w += wHold;
                wx+= wHold * xHold;
                wxx += wHold * xHold*xHold;
                wxxx += wHold * xHold*xHold*xHold;
                wxxxx += wHold * xHold*xHold*xHold*xHold;

                wy += wHold * yHold;
                wyx += wHold * yHold*xHold;
                wyxx += wHold * yHold*xHold*xHold;
            }
        }

        double c = (wyx / wxx - wyxx / wxxx) / (wxxx / wxx - wxxxx / wxxx);
        double temp = (wyx - c * wxxx) / wxx;
        double f = temp - 2 * c* xAnchor;
        double a = yAnchor - temp * xAnchor + c * xAnchor*xAnchor;



        //double c = (yxSum / x2Sum - yx2Sum / x3Sum) / (x3Sum / x2Sum - x4Sum / x3Sum);
        //double temp = (yxSum - c * x3Sum) / x2Sum;
        //double f = temp - 2 * c*xAnchor;
        //double a = yAnchor - temp*xAnchor + c*xAnchor*xAnchor;

        ret[0] = a;
        ret[1] = f;
        ret[2] = c;
        return ret;
    }


    //regression with pivoting
    std::vector<double> regressionPivot(int start, int end, std::vector<bool> &checks, std::vector<double> &weights, std::vector<double> &horizAxis, std::vector<double> &vertAxis)
    {
        double xfirst = horizAxis[0], yfirst = vertAxis[0];
        double wxSum = 0, wySum = 0, wSum = 0, wxBar, wyBar, top = 0, bottom = 0, m, b;
        std::vector<double> toRet;
        toRet.resize(2, 0.0);
        for (int i = 0; i < horizAxis.size(); i++)
        {
            if (checks[i])
            {
                wySum += weights[i] * vertAxis[i];
                wxSum += weights[i] * horizAxis[i];
                wSum += weights[i];
            }
        }
        wyBar = wySum / wSum;
        wxBar = wxSum / wSum;
        for (int i = 0; i < horizAxis.size(); i++)
        {
            if (checks[i])
            {
                //top += weights[i] * (horizAxis[i] - wxBar) * (vertAxis[i] - wyBar);
                //bottom += weights[i] * (horizAxis[i] - wxBar) * (horizAxis[i] - wxBar);
                top += weights[i] * (horizAxis[i] - xfirst) * (vertAxis[i] - yfirst);
                bottom += weights[i] * (horizAxis[i] - xfirst) * (horizAxis[i] - xfirst);
            }
        }
        m = top / bottom;
        //b = wyBar;
        b = yfirst;

        toRet[0] = b;
        toRet[1] = m;
        //std::cout << "M: " << m << " B: " << b << "xBar: " << xBar << "\n";
        return toRet;
    }
    std::vector<double> quadraticRegressionPivot(int start, int end, std::vector<bool> &checks, std::vector<double> &weights, std::vector<double> &horizAxis, std::vector<double> &vertAxis)
    {
        std::vector<double> A;
        std::vector<double> b;
        std::vector<double> coef;
        int columnCount = 3;
        A.resize(columnCount * columnCount, 0.0);
        b.resize(columnCount, 0.0);
        double wHold, xHold, yHold, w = 0, wx = 0, wy = 0, wxx = 0, wxy = 0, wxxx = 0, wxxy = 0, wxxxx = 0;

        std::vector<double> ret;
        ret.resize(3);

        for (int i = start; i <= end; i++)
        {
            if (checks[i])
            {
                wHold = weights[i];
                xHold = horizAxis[i];
                yHold = vertAxis[i];
                w += wHold;
                wx += wHold*xHold;
                wy += wHold*yHold;
                wxx += wHold*xHold*xHold;
                wxy += wHold*xHold*yHold;
                wxxx += wHold*xHold*xHold*xHold;
                wxxy += wHold*xHold*xHold*yHold;
                wxxxx += wHold*xHold*xHold*xHold*xHold;



                //horizHold = horizAxis[i];
                //vertHold = vertAxis[i];
                //xSum += weights[i] * horizHold;
                //xSquareSum += weights[i] * horizHold * horizHold;
                //xCubeSum += weights[i] * horizHold * horizHold * horizHold;
                //xQuartSum += weights[i] * horizHold * horizHold * horizHold * horizHold;
                //prodSum += weights[i] * horizHold * vertHold;
                //ySum += weights[i] * vertHold;
                //xSquareYSum += weights[i] * horizHold * horizHold * vertHold;
                //length += weights[i];
            }
        }


        A[0] = w; A[1] = wx; A[2] = wxx;
        A[3] = wx; A[4] = wxx; A[5] = wxxx;
        A[6] = wxx; A[7] = wxxx; A[8] = wxxxx;
        b[0] = wy; b[1] = wxy; b[2] = wxxy;

        coef = performPivot(columnCount, A, b);

        ret[0] = coef[0];
        ret[1] = coef[1];
        ret[2] = coef[2];

        return ret;
    }
    std::vector<double> fixedRegressionPivot(int start, int end, std::vector<bool> &checks, std::vector<double> &weights, std::vector<double> &horizAxis, std::vector<double> &vertAxis, double horizAxisCoordinate, double vertAxisCoordinate)
    {
        double xx = 0, xy = 0;
        double m, b;
        double horizHold, vertHold;
        std::vector<double> ret;
        ret.resize(2);

        for (int i = start; i <= end; i++)
        {
            if (checks[i])
            {
                horizHold = horizAxis[i] - horizAxisCoordinate;
                vertHold = vertAxis[i] - vertAxisCoordinate;

                xx += weights[i] * horizHold*horizHold;
                xy += weights[i] * horizHold*vertHold;
            }
        }

        m = xy / xx;
        b = vertAxisCoordinate - m*horizAxisCoordinate;

        ret[0] = b;
        ret[1] = m;

        return ret;
    }
    std::vector<double> fixedQuadraticRegressionPivot(int start, int end, std::vector<bool> &checks, std::vector<double> &weights, std::vector<double> &horizAxis, std::vector<double> &vertAxis, double horizAxisCoordinate, double vertAxisCoordinate)
    {
        std::vector<double> coef;
        double w = 0;
        double xSum = 0, ySum = 0;
        double pivot, swap, multiplier;
        double finalAnswer;
        int pivotIndex;
        int columnCount = 3;
        std::vector<double> A;
        std::vector<double> b;
        A.resize(columnCount * columnCount, 0.0);
        b.resize(columnCount, 0.0);

        std::vector<double> ret;
        ret.resize(4);
        double x2Sum = 0, x3Sum = 0, x4Sum = 0, yxSum = 0, yx2Sum = 0, error = 0, horizHold, vertHold;
        for (int i = start; i <= end; i++)
        {
            if (checks[i])
            {
                horizHold = (horizAxis[i] - horizAxisCoordinate);
                vertHold = (vertAxis[i] - vertAxisCoordinate);
                w += weights[i];

                ySum += weights[i] * vertHold;
                xSum += weights[i] * horizHold;
                x2Sum += weights[i] * horizHold*horizHold;
                x3Sum += weights[i] * horizHold*horizHold*horizHold;
                x4Sum += weights[i] * horizHold*horizHold*horizHold*horizHold;
                yxSum += weights[i] * vertHold*horizHold;
                yx2Sum += weights[i] * vertHold*horizHold*horizHold;
            }
        }



        A[0] = w; A[1] = xSum; A[2] = x2Sum;
        A[3] = xSum; A[4] = x2Sum; A[5] = x3Sum;
        A[6] = x2Sum; A[7] = x3Sum; A[8] = x4Sum;
        b[0] = ySum; b[1] = yxSum; b[2] = yx2Sum;

        coef = performPivot(columnCount, A, b);

        ret[0] = coef[0];
        ret[1] = coef[1];
        ret[2] = coef[2];
        ret[3] = 0;

        return ret;
    }
    */

    //regression with pivoting
    std::vector<double> regressionPivot(std::vector<bool> &checks, std::vector<double> &weights,
                                        std::vector<double> &xAxis, std::vector<double> &yAxis) {
        /*
        This function uses standard linear regression and pivoting to satisfy a least squares fit of y = mx + b
        */
        double w = 0, wx = 0, wy = 0;
        double wxx = 0, wxy = 0;
        double wxBar, wyBar, top = 0, bottom = 0;
        double mPrime, bPrime, mCoef, bCoef;
        double wHold, xHold, yHold;
        double xShift = xAxis[0];
        double yShift = yAxis[0];
        std::vector<double> A(4), b(2);
        std::vector<double> coef;
        std::vector<double> toRet(2);

        for (int i = 0; i < xAxis.size(); i++) {
            if (checks[i]) {
                wHold = weights[i];
                xHold = xAxis[i] - xShift;
                //We subtract off the first value of the vector, NOT because it is an anchor. It is simply rescales the values for more well-behaved regression.
                yHold = yAxis[i] - yShift;

                w += wHold;
                wx += wHold * xHold;
                wy += wHold * yHold;
                wxx += wHold * xHold * xHold;
                wxy += wHold * xHold * yHold;
            }
        }

        A[0] = w;
        A[1] = wx;
        A[2] = wx;
        A[3] = wxx;
        b[0] = wy;
        b[1] = wxy;

        coef = performPivot(2, A, b);

        bPrime = coef[0];
        mPrime = coef[1];

        bCoef = bPrime + (yShift - mPrime * xShift);
        mCoef = mPrime;

        toRet[0] = bCoef;
        toRet[1] = mCoef;

        return toRet;
    }

    std::vector<double> quadraticRegressionPivot(std::vector<bool> &checks, std::vector<double> &weights,
                                                 std::vector<double> &xAxis, std::vector<double> &yAxis) {
        int columnCount = 3;
        std::vector<double> A(9);
        std::vector<double> b(3);
        std::vector<double> coef, coefPrime;
        coef.resize(3);
        double wHold, xHold, yHold;
        double xShift = xAxis[0];
        double yShift = yAxis[0];
        double w = 0, wx = 0, wy = 0, wxx = 0, wxy = 0, wxxx = 0, wxxy = 0, wxxxx = 0;

        for (int i = 0; i < checks.size(); i++) {
            if (checks[i]) {
                wHold = weights[i];
                xHold = xAxis[i] - xShift;
                yHold = yAxis[i] - yShift;
                w += wHold;
                wx += wHold * xHold;
                wxx += wHold * xHold * xHold;
                wxxx += wHold * xHold * xHold * xHold;
                wxxxx += wHold * xHold * xHold * xHold * xHold;

                wy += wHold * yHold;
                wxy += wHold * xHold * yHold;
                wxxy += wHold * xHold * xHold * yHold;
            }
        }

        A[0] = w;
        A[1] = wx;
        A[2] = wxx;
        A[3] = wx;
        A[4] = wxx;
        A[5] = wxxx;
        A[6] = wxx;
        A[7] = wxxx;
        A[8] = wxxxx;
        b[0] = wy;
        b[1] = wxy;
        b[2] = wxxy;

        coefPrime = performPivot(columnCount, A, b);

        double cPrime = coefPrime[2], bPrime = coefPrime[1], aPrime = coefPrime[0];
        double c_corr, b_corr, a_corr;

        c_corr = cPrime;
        b_corr = bPrime - 2 * cPrime * xShift;
        a_corr = yShift + aPrime - bPrime * xShift + cPrime * xShift * xShift;

        coef[0] = a_corr;
        coef[1] = b_corr;
        coef[2] = c_corr;

        return coef;
    }

    std::vector<double> fixedRegressionPivot(std::vector<bool> &checks, std::vector<double> &weights,
                                             std::vector<double> &xAxis, std::vector<double> &yAxis, double xAnchor,
                                             double yAnchor) {
        /*
        Linear Fixed Regression removes the b term from the standard y = mx + b form. As such, this algorithm does not require pivoting.
        */
        double wxx = 0, wxy = 0;
        double m, b;
        double wHold, xHold, yHold;
        std::vector<double> ret(2);


        for (int i = 0; i < checks.size(); i++) {
            if (checks[i]) {
                wHold = weights[i];
                xHold = xAxis[i] - xAnchor;
                yHold = yAxis[i] - yAnchor;

                wxx += wHold * xHold * xHold;
                wxy += wHold * xHold * yHold;
            }
        }

        //No pivot needed, matrix is a 1x1
        m = wxy / wxx;
        b = yAnchor - m * xAnchor;

        ret[0] = b;
        ret[1] = m;

        return ret;
    }

    std::vector<double> fixedQuadraticRegressionPivot(std::vector<bool> &checks, std::vector<double> &weights,
                                                      std::vector<double> &xAxis, std::vector<double> &yAxis,
                                                      double xAnchor, double yAnchor) {
        int columnCount = 2;
        std::vector<double> A(4);
        std::vector<double> b(2);
        std::vector<double> coef, coefPrime;
        coef.resize(3);

        double wHold, xHold, yHold;
        double w = 0, wx = 0, wxx = 0, wxxx = 0, wxxxx = 0, wy = 0, wyx = 0, wyxx = 0;
        for (int i = 0; i < checks.size(); i++) {
            if (checks[i]) {
                xHold = (xAxis[i] - xAnchor);
                yHold = (yAxis[i] - yAnchor);
                wHold = weights[i];;

                wxx += wHold * xHold * xHold;
                wxxx += wHold * xHold * xHold * xHold;
                wxxxx += wHold * xHold * xHold * xHold * xHold;

                wyx += wHold * yHold * xHold;
                wyxx += wHold * yHold * xHold * xHold;
            }
        }


        A[0] = wxx;
        A[1] = wxxx;
        A[2] = wxxx;
        A[3] = wxxxx;
        b[0] = wyx;
        b[1] = wyxx;

        coefPrime = performPivot(columnCount, A, b);

        //ADJUST FOR ANCHORS
        double cPrime = coefPrime[1], bPrime = coefPrime[0];
        double c_corr, b_corr, a_corr;

        c_corr = cPrime;
        b_corr = bPrime - 2 * cPrime * xAnchor;
        a_corr = yAnchor - bPrime * xAnchor + cPrime * xAnchor * xAnchor;

        coef[0] = a_corr;
        coef[1] = b_corr;
        coef[2] = c_corr;

        return coef;
    }


    //distance calculators
    double getGCDistance(double dec, double ra, double medianDec, double medianRa, double globalDecCenter) {
        double toRad = M_PI / 180.0;
        //double undoTransformRa = ra / cos(dec*toRad);
        //double undoTransformCenterRa = medianRa / cos(dec*toRad);

        double undoTransformRa = ra / cos(globalDecCenter * toRad);
        double undoTransformCenterRa = medianRa / cos(globalDecCenter * toRad);

        //double distance = acos(sin(dec*toRad) * sin(medianDec * toRad) + cos(dec*toRad) * cos(medianDec * toRad) * cos((ra - medianRa)*toRad));
        //double distance = acos(sin(dec*toRad) * sin(medianDec * toRad) + cos(dec*toRad) * cos(medianDec * toRad) * cos((undoTransformRa - undoTransformCenterRa)*toRad));
        double distance = acos(
            sin((dec + globalDecCenter) * toRad) * sin((medianDec + globalDecCenter) * toRad) +
            cos((dec + globalDecCenter) * toRad) * cos((medianDec + globalDecCenter) * toRad) * cos(
                (undoTransformRa - undoTransformCenterRa) * toRad));
        if (distance != distance) {
            distance = 0;
        }
        return distance;
    }

    double getModGCDistance(double dec, double ra, double medianDec, double medianRa) {
        double toRad = M_PI / 180.0;
        //double undoTransformRa = ra / cos(dec*toRad);
        //double undoTransformCenterRa = medianRa / cos(dec*toRad);

        double undoTransformRa = ra;
        double undoTransformCenterRa = medianRa;

        double distance = acos(
            sin(dec * toRad) * sin(medianDec * toRad) + cos(dec * toRad) * cos(medianDec * toRad) * cos(
                (ra - medianRa) * toRad));
        //double distance = acos(sin(dec*toRad) * sin(medianDec * toRad) + cos(dec*toRad) * cos(medianDec * toRad) * cos((undoTransformRa - undoTransformCenterRa)*toRad));
        //double distance = acos(sin((dec + globalDecCenter)*toRad) * sin((medianDec + globalDecCenter)* toRad) + cos((dec + globalDecCenter)*toRad) * cos((medianDec + globalDecCenter) * toRad) * cos((undoTransformRa - undoTransformCenterRa)*toRad));
        if (distance != distance) {
            distance = 0;
        }
        return distance;
    }

    double getPythDistance(double decRef, double raRef, double dec, double ra) {
        return std::abs(sqrt(pow((ra - raRef), 2) + pow((dec - decRef), 2)));
    }


    //misc
    double findPeriod(std::vector<double> &times, std::vector<double> &speeds) {
        //LOOK HERE: WE CHOSE H = 2 AND TIMES[1] BECAUSE WE WANTED TO PUSH SPEEDS BACK WITH A DUMMY ELEMENT IN THE FIRST BOX SO AS NOT TO HAVE TO CREATE A NEW TIMES ARRAY****************************
        int h, requiredPower = 1, maxBinIndex;
        double timeStep, timeHold, hold, maxBinVal = -999999, freq;
        std::vector<double> interpolated;
        std::vector<std::vector<double> > result;
        while (pow(2, requiredPower) < times.size()) {
            requiredPower++;
        }
        requiredPower += 3;
        interpolated.resize(pow(2, requiredPower), 0.0);
        timeStep = ((double) (times[times.size() - 1] - times[1])) / ((double) (pow(2, requiredPower)));
        timeHold = times[1];
        h = 2;
        for (int i = 0; i < interpolated.size(); i++) {
            while (h < times.size() && times[h] < timeHold) {
                h++;
            }
            if (times[h] == times[h - 1]) {
                interpolated[i] = speeds[h];
            } else {
                interpolated[i] = speeds[h - 1] + (timeHold - times[h - 1]) * (speeds[h] - speeds[h - 1]) / (
                                      times[h] - times[h - 1]);
            }
            timeHold += timeStep;
        }
        result = fft(interpolated);
        for (int i = 1; i < result.size() / 2.0; i++) {
            hold = sqrt(result[i][0] * result[i][0] + result[i][1] * result[i][1]);
            if (hold > maxBinVal) {
                maxBinVal = hold;
                maxBinIndex = i;
            }
        }
        freq = (double) maxBinIndex * (1.0 / timeStep) / (double) interpolated.size();
        return 1.0 / freq;
    }

    std::vector<double> daisyAngleBuilder(int i, int centerIndex, double centerDec, std::vector<double> ra,
                                          std::vector<double> dec) {
        double distance, negation, toRad = M_PI / 180.0, medianRa = ra[centerIndex], medianDec = dec[centerIndex];
        std::vector<double> angles;
        negation = -pow(-1.0, i);
        for (int j = 0; j < ra.size(); j++) {
            distance = getGCDistance(dec[j], ra[j], medianDec, medianRa, centerDec) * (180.0 / M_PI);
            if (distance != distance) {
                distance = 0;
            }
            if (j == centerIndex) {
                negation *= -1.0;
            }
            distance *= negation;
            angles.push_back(distance);
        }
        return angles;
    }

    void lasso(int pointsInRange, std::vector<int> &checks, std::vector<double> &alongScan,
               std::vector<double> &acrossScan, std::vector<double> &lows, std::vector<double> &highs) {
        int pointsNotRejected = 0;
        for (int i = 0; i < checks.size(); i++) {
            if (checks[i]) {
                pointsNotRejected++;
            }
        }
        if (pointsNotRejected > 1) {
            int pivot = 0, swapStart, swapEnd, whileCount = 0, jMem = -1;
            bool stop = false;
            double angle, angleHold = 999999, anglePrev = 0, top, bottom;
            std::vector<int> pivots;
            while (!stop) {
                for (int j = 0; j < pointsInRange; j++) {
                    if (j != pivot && checks[j]) {
                        top = acrossScan[j] - acrossScan[pivot];
                        bottom = alongScan[j] - alongScan[pivot];
                        angle = 180.0 * Tools::signedArcTan(top, bottom) / M_PI;
                        if (!(top == 0 && bottom == 0)) {
                            if (angle - anglePrev < 0) {
                                angle += 360;
                            }
                            if (angle - anglePrev < angleHold) {
                                angleHold = angle - anglePrev;
                                jMem = j;
                            }
                        }
                    }
                }
                anglePrev += angleHold;
                angleHold = 999999;
                pivot = jMem;
                if (anglePrev <= 360.0) {
                    pivots.push_back(pivot);
                } else {
                    stop = true;
                }
                whileCount++;
            }
            for (int i = 0; i < pivots.size() - 1; i++) {
                swapStart = pivots[i];
                swapEnd = pivots[i + 1];
                for (int k = 0; k < acrossScan.size(); k++) {
                    if ((alongScan[k] >= alongScan[swapStart] && alongScan[k] <= alongScan[swapEnd]) || (
                            alongScan[k] <= alongScan[swapStart] && alongScan[k] >= alongScan[swapEnd])) {
                        angleHold = acrossScan[swapStart] + (acrossScan[swapEnd] - acrossScan[swapStart]) * (
                                        alongScan[k] - alongScan[swapStart]) / (
                                        alongScan[swapEnd] - alongScan[swapStart]);
                        if (angleHold <= lows[k]) {
                            lows[k] = angleHold;
                        }
                        if (angleHold >= highs[k]) {
                            highs[k] = angleHold;
                        }
                    }
                }
            }
            swapStart = pivots[pivots.size() - 1];
            swapEnd = pivots[0];
            for (int k = 0; k < acrossScan.size(); k++) {
                if ((alongScan[k] >= alongScan[swapStart] && alongScan[k] <= alongScan[swapEnd]) || (
                        alongScan[k] <= alongScan[swapStart] && alongScan[k] >= alongScan[swapEnd])) {
                    angleHold = acrossScan[swapStart] + (acrossScan[swapEnd] - acrossScan[swapStart]) * (
                                    alongScan[k] - alongScan[swapStart]) / (alongScan[swapEnd] - alongScan[swapStart]);
                    if (angleHold <= lows[k]) {
                        lows[k] = angleHold;
                    }
                    if (angleHold >= highs[k]) {
                        highs[k] = angleHold;
                    }
                }
            }
        } else {
            //this code looks bugged, but pointsNotRejected cannot be 0, so it is OK.
            int indexNotRejected = 0;
            while (!checks[indexNotRejected]) {
                indexNotRejected++;
            }
            lows[indexNotRejected] = acrossScan[indexNotRejected];
            highs[indexNotRejected] = acrossScan[indexNotRejected];
        }
    }

    int determineNearestIndex(int index, std::vector<double> &vec1, std::vector<double> &vec2) {
        int j;
        if (vec1[vec1.size() - 1] > vec1[0]) {
            j = vec2.size() - 1;
            while (vec1[index] > vec2[j] && j != 0) {
                j--;
            }

            if (j == 0 && vec1[index] > vec2[j]) {
                return -1;
            }

            if (j < vec2.size() - 1) {
                if (std::abs(vec2[j] - vec1[index]) > std::abs(vec2[j + 1] - vec1[index])) {
                    j++;
                }
            }
            if (j == vec2.size() - 1 && vec2[j] > vec1[index]) {
                j = vec2.size();
            }
        } else {
            j = 0;
            while (vec1[index] > vec2[j] && j != vec2.size() - 1) {
                j++;
            }

            if (j == vec2.size() - 1 && vec1[index] > vec2[j]) {
                return vec2.size();
            }

            if (j < 0) {
                if (std::abs(vec2[j] - vec1[index]) > std::abs(vec2[j - 1] - vec1[index])) {
                    j--;
                }
            }

            if (j == 0 && vec2[j] > vec1[index]) {
                j = -1;
            }
        }
        return j;
    }

    //cartographer tools
    int determinePixel(double i, double iRef, double resolution) {
        return round((1.0 / resolution) * (i - iRef));
    }
}
