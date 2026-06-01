
#include <cstdio>
#include <vector>
#include <cmath>

std::vector<double> lowerTriangleSolver(int columnCount,
                                         std::vector<double> A,
                                         std::vector<double> b) {
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

std::vector<double> performPivot(int columnCount, std::vector<double> A,
                                  std::vector<double> b) {
    int pivotIndex;
    double pivot, swap, multiplier;
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
            for (int j = 0; j < columnCount; j++)
                A[columnCount * i + j] /= pivot;
            b[i] /= pivot;
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
            for (int k = 0; k < columnCount; k++)
                A[j * columnCount + k] -= multiplier * A[i * columnCount + k];
            b[j] -= multiplier * b[i];
        }
    }
    return lowerTriangleSolver(columnCount, A, b);
}

int main() {
    int cc;
    std::scanf("%d", &cc);
    std::vector<double> A(cc * cc), b(cc);
    for (int i = 0; i < cc * cc; i++) std::scanf("%lf", &A[i]);
    for (int i = 0; i < cc; i++) std::scanf("%lf", &b[i]);
    auto x = performPivot(cc, A, b);
    for (int i = 0; i < cc; i++) std::printf("%.17g\n", x[i]);
    return 0;
}
