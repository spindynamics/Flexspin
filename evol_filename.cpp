#include "objects.hpp"

#include <cmath>
#include <cstddef>
#include <iomanip>
#include <sstream>

namespace {

// Picks the largest SI prefix (n, u, m, k, M, ...) whose threshold the value
// clears, so filenames show e.g. "100.000uA" instead of "0.000" for 100e-6.
struct ScaledValue { double magnitude; const char* prefix; };
struct SiPrefix { double factor; const char* prefix; };

ScaledValue scaleMagnitude(double value) noexcept {
    static const SiPrefix steps[] = {
        {1e-9, "n"}, {1e-6, "u"}, {1e-3, "m"}, {1.0, ""},
        {1e3, "k"}, {1e6, "M"}, {1e9, "G"}
    };
    const double absValue = std::fabs(value);
    if (absValue == 0.0) return {0.0, ""};

    double chosenFactor = steps[0].factor;
    const char* chosenPrefix = steps[0].prefix;
    for (size_t i = 0; i < sizeof(steps) / sizeof(steps[0]); ++i) {
        if (absValue >= steps[i].factor) {
            chosenFactor = steps[i].factor;
            chosenPrefix = steps[i].prefix;
        }
    }
    return {value / chosenFactor, chosenPrefix};
}

// range1/range2 "type" is a path like "junction[x].junctionBias.dcCurrent.A";
// the last segment ("A", "H", "V", ...) is used as the unit label.
std::string unitSuffix(const std::string& typePath) {
    const size_t dot = typePath.find_last_of('.');
    return (dot == std::string::npos) ? "" : typePath.substr(dot + 1);
}

}  // namespace

std::string buildEvolFileName(const Simulation& simulation,
                              const std::vector<Junction>& junctions,
                              int simulationIndex,
                              int passIndex,
                              double fieldValue,
                              double voltageValue,
                              bool appendDirectionSuffix,
                              int directionIndex) noexcept {
    std::ostringstream fnameStream;
    fnameStream << std::fixed << std::setfill('0');
    fnameStream << std::setprecision(3);
    fnameStream << simulation.rootName << std::setw(2) << simulationIndex;
    if (simulation.range1.flag) {
        const ScaledValue scaled = scaleMagnitude(fieldValue);
        fnameStream << "_" << std::setw(6) << scaled.magnitude
                    << scaled.prefix << unitSuffix(simulation.range1.type);
    }
    if (simulation.range2.flag) {
        const ScaledValue scaled = scaleMagnitude(voltageValue);
        fnameStream << "_" << std::setw(6) << scaled.magnitude
                    << scaled.prefix << unitSuffix(simulation.range2.type);
    }
    fnameStream << "_" << std::setprecision(1) << std::setw(3)
                << (junctions.empty() ? 0.0 : junctions.front().temperature.T0);
    
                
    if (appendDirectionSuffix) {
        if (directionIndex == 0) {
            fnameStream << "K.evolup";
        }
        if (directionIndex == 1) {
            fnameStream << "K.evoldown";
        }
    }else{fnameStream << "K.evol";}

    return fnameStream.str();
}