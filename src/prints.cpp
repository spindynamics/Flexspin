#include "printing.hpp"

#include <iomanip>
#include <mpi.h>

using std::fixed;
using std::setprecision;
using std::setw;

namespace {

// MPI tags used to stream final-state rows from every non-root rank to rank 0.
// kFinalRowTag carries one packed row; kFinalDoneTag is an empty message a rank sends
// once it has no more rows coming.
constexpr int kFinalRowTag = 9001;
constexpr int kFinalDoneTag = 9002;

std::vector<double> packFinalStateRow(const std::vector<Junction>& junctions,
                                       double fieldValue, double voltageValue, int pass) {
    std::vector<double> buf;
    buf.reserve(finalStateRowSize(junctions));
    buf.push_back(fieldValue);
    buf.push_back(voltageValue);
    for (const auto& junction : junctions) {
        for (const auto& layer : junction.layers) {
            buf.push_back(layer.m.x);
            buf.push_back(layer.m.y);
            buf.push_back(layer.m.z);
        }
    }
    buf.push_back(static_cast<double>(pass));
    for (const auto& junction : junctions) {
        buf.push_back(junction.T);
        buf.push_back(junction.R);
    }
    return buf;
}

// Mirrors writeFinalStateRow()'s column order/formatting, but reads values out of a flat
// buffer (received over MPI from a non-root rank) instead of live Junction objects.
void writeFinalStateRowFromBuffer(std::ostream& out, const std::vector<Junction>& junctions,
                                   const std::vector<double>& buf) {
    size_t idx = 0;
    out << fixed << setprecision(10) << setw(13) << buf[idx++];
    out << "\t" << fixed << setprecision(10) << setw(13) << buf[idx++];
    for (const auto& junction : junctions) {
        for (size_t li = 0; li < junction.layers.size(); ++li) {
            out << "\t" << fixed << setprecision(10) << setw(13) << buf[idx++];
            out << "\t" << fixed << setprecision(10) << setw(13) << buf[idx++];
            out << "\t" << fixed << setprecision(10) << setw(13) << buf[idx++];
        }
    }
    const int pass = static_cast<int>(buf[idx++]);
    out << "\t" << setw(13) << (pass == 0 ? "up" : "down");
    for (size_t ji = 0; ji < junctions.size(); ++ji) {
        out << "\t" << fixed << setprecision(10) << setw(13) << buf[idx++];
        out << "\t" << fixed << setprecision(10) << setw(13) << buf[idx++];
    }
    out << "\n";
    out.flush();
}

} // namespace

///////////////////////////////////////////////////////////////////////////////////////
// evol file                                                                          //
///////////////////////////////////////////////////////////////////////////////////////

void writeEvolHeader(std::ostream& evol, const std::vector<std::string>& labels) {
    evol << setw(13) << "t_ns";
    for (const auto& label : labels) evol << "\t" << setw(13) << label;
    evol << "\n";
}

void writeEvolRow(std::ostream& evol, const DoubleRefs& values, double t_ns) {
    evol << fixed << setprecision(10) << setw(13) << t_ns;
    for (double v : values) evol << "\t" << fixed << setprecision(10) << setw(13) << v;
    evol << "\n";
    evol.flush(); // Push each written row to disk immediately, instead of waiting for
                  // evolFile to close at the end of the rangepoint (or being lost
                  // entirely if the run is killed/crashes mid-rangepoint).
}

///////////////////////////////////////////////////////////////////////////////////////
// 1-State-Diagram.txt                                                                //
///////////////////////////////////////////////////////////////////////////////////////

void writeFinalStateHeader(std::ostream& out, const std::vector<Junction>& junctions) {
    out << setw(13) << "H" << "\t" << setw(13) << "V";
    for (size_t ji = 0; ji < junctions.size(); ++ji) {
        const std::string j = "_J" + std::to_string(ji);
        for (size_t li = 0; li < junctions[ji].layers.size(); ++li) {
            const std::string jl = j + "L" + std::to_string(li);
            out << "\t" << setw(13) << ("mx" + jl);
            out << "\t" << setw(13) << ("my" + jl);
            out << "\t" << setw(13) << ("mz" + jl);
        }
    }
    out << "\t" << setw(13) << "direction";
    for (size_t ji = 0; ji < junctions.size(); ++ji) {
        const std::string j = "_J" + std::to_string(ji);
        out << "\t" << setw(13) << ("T" + j);
        out << "\t" << setw(13) << ("R" + j);
    }
    out << "\n";
}

void writeFinalStateRow(std::ostream& out, const std::vector<Junction>& junctions,
                         double fieldValue, double voltageValue, int pass) {
    out << fixed << setprecision(10) << setw(13) << fieldValue;
    out << "\t" << fixed << setprecision(10) << setw(13) << voltageValue;
    for (const auto& junction : junctions) {
        for (const auto& layer : junction.layers) {
            out << "\t" << fixed << setprecision(10) << setw(13) << layer.m.x;
            out << "\t" << fixed << setprecision(10) << setw(13) << layer.m.y;
            out << "\t" << fixed << setprecision(10) << setw(13) << layer.m.z;
        }
    }
    out << "\t" << setw(13) << (pass == 0 ? "up" : "down");
    for (const auto& junction : junctions) {
        out << "\t" << fixed << setprecision(10) << setw(13) << junction.T;
        out << "\t" << fixed << setprecision(10) << setw(13) << junction.R;
    }
    out << "\n";
    out.flush();
}

size_t finalStateRowSize(const std::vector<Junction>& junctions) {
    size_t n = 2; // H, V
    for (const auto& junction : junctions) n += 3 * junction.layers.size();
    n += 1; // pass
    n += 2 * junctions.size(); // T, R per junction
    return n;
}

void drainFinalStateRows(std::ostream& out, const std::vector<Junction>& junctions, size_t rowSize) {
    int flag = 0;
    MPI_Status status;
    while (true) {
        MPI_Iprobe(MPI_ANY_SOURCE, kFinalRowTag, MPI_COMM_WORLD, &flag, &status);
        if (!flag) break;
        std::vector<double> buf(rowSize);
        MPI_Recv(buf.data(), static_cast<int>(rowSize), MPI_DOUBLE, status.MPI_SOURCE,
                 kFinalRowTag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        writeFinalStateRowFromBuffer(out, junctions, buf);
    }
}

// MPI's non-overtaking guarantee (messages from a given sender are matched in the order
// they were sent) means a rank's still-queued rows always get matched here before that
// same rank's later "done" message does.
void finalizeFinalStateRows(std::ostream& out, const std::vector<Junction>& junctions,
                             size_t rowSize, int total_task) {
    int doneCount = 0;
    while (doneCount < total_task - 1) {
        MPI_Status status;
        MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (status.MPI_TAG == kFinalRowTag) {
            std::vector<double> buf(rowSize);
            MPI_Recv(buf.data(), static_cast<int>(rowSize), MPI_DOUBLE, status.MPI_SOURCE,
                     kFinalRowTag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            writeFinalStateRowFromBuffer(out, junctions, buf);
        } else { // kFinalDoneTag
            char dummy = 0;
            MPI_Recv(&dummy, 0, MPI_CHAR, status.MPI_SOURCE, kFinalDoneTag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            ++doneCount;
        }
    }
}

void sendFinalStateRow(const std::vector<Junction>& junctions,
                        double fieldValue, double voltageValue, int pass, size_t rowSize) {
    std::vector<double> buf = packFinalStateRow(junctions, fieldValue, voltageValue, pass);
    MPI_Send(buf.data(), static_cast<int>(rowSize), MPI_DOUBLE, 0, kFinalRowTag, MPI_COMM_WORLD);
}

void sendFinalStateDone() {
    char dummy = 0;
    MPI_Send(&dummy, 0, MPI_CHAR, 0, kFinalDoneTag, MPI_COMM_WORLD);
}
