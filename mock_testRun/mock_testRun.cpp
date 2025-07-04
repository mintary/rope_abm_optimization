#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <sstream>

// Helper to join values with commas
template<typename T>
std::string join(const std::vector<T>& values, const std::string& delim = ",") {
    std::ostringstream oss;
    for (size_t i = 0; i < values.size(); ++i) {
        if (i != 0) oss << delim;
        oss << values[i];
    }
    return oss.str();
}

int main(int argc, char* argv[]) {
    // Default values
    int numticks = 10;
    std::string inputfile = "config.txt";
    double wxw = 1.0, wyw = 1.0, wzw = 1.0;

    // Parse arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--numticks" && i + 1 < argc) {
            numticks = std::atoi(argv[++i]);
        } else if (arg == "--inputfile" && i + 1 < argc) {
            inputfile = argv[++i];
        } else if (arg == "--wxw" && i + 1 < argc) {
            wxw = std::atof(argv[++i]);
        } else if (arg == "--wyw" && i + 1 < argc) {
            wyw = std::atof(argv[++i]);
        } else if (arg == "--wzw" && i + 1 < argc) {
            wzw = std::atof(argv[++i]);
        }
    }

    // Read from Sample.txt
    std::ifstream samplefile("Sample.txt");
    if (samplefile.is_open()) {
        std::string line;
        std::cout << "Contents of Sample.txt:\n";
        while (std::getline(samplefile, line)) {
            std::cout << line << std::endl;
        }
        samplefile.close();
    } else {
        std::cerr << "Could not open Sample.txt\n";
    }

    // Output file
    std::ofstream outfile("output/Output_Biomarkers.csv");
    std::vector<std::string> headers = {
        "clock","TNF","TGF","FGF","IL6","IL8","IL10","Tropocollagen","Collagen","FragentedCollagen",
        "Tropoelastin","Elastin","FragmentedElastin","HA","FragmentedHA","Damage","ActivatedFibroblast",
        "Fibroblast","Elastic Mod (Pa)","Swelling Ratio","Mass Loss (%)"
    };
    outfile << join(headers) << "\n";

    // Write mock data
    for (int tick = 0; tick < numticks; ++tick) {
        std::vector<double> row = {
            static_cast<double>(tick),
            0.1 * tick,      // TNF
            0.2 * tick,      // TGF
            0.3 * tick,      // FGF
            0.4 * tick,      // IL6
            0.5 * tick,      // IL8
            0.6 * tick,      // IL10
            1.0 * tick,      // Tropocollagen
            2.0 * tick,      // Collagen
            0.5 * tick,      // FragentedCollagen
            1.1 * tick,      // Tropoelastin
            2.1 * tick,      // Elastin
            0.6 * tick,      // FragmentedElastin
            0.7 * tick,      // HA
            0.8 * tick,      // FragmentedHA
            0.9 * tick,      // Damage
            1.2 * tick,      // ActivatedFibroblast
            1.3 * tick,      // Fibroblast
            static_cast<double>(1000 + 10 * tick),// Elastic Mod (Pa)
            1.0 + 0.01 * tick,// Swelling Ratio
            0.5 * tick       // Mass Loss (%)
        };
        outfile << join(row) << "\n";
    }

    outfile.close();
    std::cout << "Mock test run complete. Output_Biomarkers.csv generated.\n";
    return 0;
}