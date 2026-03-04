#include "SpectralBackground.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

int main(int argc, char* argv[])
{
    if (argc != 3)
    {
        std::cerr << "Usage: spectral-cleaner <input_json> <output_json>\n";
        return 1;
    }

    std::ifstream ifile(argv[1]);
    if (!ifile.is_open())
    {
        std::cerr << "Failed to open input file.\n";
        return 1;
    }

    json j_in;
    ifile >> j_in;

    std::vector<double> freqs = j_in["frequencies"].get<std::vector<double>>();
    std::vector<double> flux = j_in["intensities"].get<std::vector<double>>();
    std::vector<double> weights;
    if (j_in.contains("weights")) {
        weights = j_in["weights"].get<std::vector<double>>();
    } else {
        weights.resize(freqs.size(), 1.0);
    }
    double scatter = j_in["scatter"].get<double>();
    double baseline = j_in["baseline"].get<double>();

    SpectralBackground bg(freqs, flux, weights, scatter);
    std::vector<double> modeled_spectrum = bg.calculateBGMulti(baseline);

    json j_out;
    j_out["modeled_spectrum"] = modeled_spectrum;

    std::ofstream ofile(argv[2]);
    if (!ofile.is_open())
    {
        std::cerr << "Failed to open output file.\n";
        return 1;
    }
    ofile << j_out.dump(4);
    
    return 0;
}
