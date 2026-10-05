#include "printing.hpp"

#include <cmath>
#include <fstream>
#include <stdexcept>

// Splits range1 across MPI tasks, then runs: for each iteration, for each hysteresis
// pass (forward, and backward if hysteresis is set), sweep range1 x range2.
void Simulation::range(std::vector<Junction> &junctions, std::ostream &log, int task_id, int total_task)
{
    // Split range1 into total_task contiguous slices of N1 steps; this task sweeps [R1_start, R1_stop]
    // If a range's flag is off, it's held fixed at its start value (N=0, single point) instead of sweeping.
    int N1 = range1.flag ? (range1.step / total_task) : 0;
    double step1 = range1.flag ? (range1.stop - range1.start) / range1.step : 0.0;
    int start_step = task_id * N1;
    int stop_step = (task_id + 1) * N1;
    double R1_start = range1.start + start_step * step1;
    double R1_stop = range1.start + stop_step * step1;

    // Range1 is the one cut though MPI, so it can't be hysteretic or continuous, these parameter concern only Range2

    int N2 = range2.flag ? range2.step : 0;

    // Resolved once: every live field range1/range2's swept value gets written into
    // (more than one if range1.type/range2.type uses a "[x]" wildcard)
    DoubleRefs targets1 = resolvePath(*this, junctions, range1.type);
    DoubleRefs targets2 = resolvePath(*this, junctions, range2.type);

    // .evol columns: resolved once here (not per row) from "saving columns" in the input.
    // Each entry must resolve to exactly one field -- no [x] wildcards, since each column
    // needs a single unambiguous header/value pair.
    DoubleRefs savingRefs;
    savingRefs.reserve(savingColumns.size());
    for (const std::string &columnPath : savingColumns)
    {
        DoubleRefs resolved = resolvePath(*this, junctions, columnPath);
        if (resolved.size() != 1)
        {
            throw std::runtime_error("simulation.saving columns entry '" + columnPath +
                                     "' must resolve to exactly one field (no [x] wildcards)");
        }
        savingRefs.push_back(resolved.front());
    }

    // One file for the whole sweep: each rangepoint appends a single row once its
    // temporal loop finishes, instead of the per-rangepoint evol files above. Only rank 0
    // ever opens/writes this file directly; every other rank streams its rows to rank 0
    // over MPI instead (see printing.hpp for why).
    const size_t finalRowSize = finalStateSaving ? finalStateRowSize(junctions) : 0;
    std::ofstream finalStateFile;
    if (finalStateSaving && task_id == 0)
    {
        const std::string finalStateFileName = "1-State-Diagram.txt";
        log << "final state file name = " << finalStateFileName << "\n";
        finalStateFile.open(finalStateFileName);
        writeFinalStateHeader(finalStateFile, junctions);
    }
    else if (finalStateSaving)
    {
        log << "final state rows will be streamed to rank 0 over MPI\n";
    }

    //////////////////////////////////////////////////////
    //                  ITERATIONS
    //////////////////////////////////////////////////////
    for (int it = 0; it < iterations; ++it)
    {
        int passes = hysteresis ? 2 : 1;

        //////////////////////////////////////////////////////
        //                    HYSTERESIS
        //////////////////////////////////////////////////////
        for (int pass = 0; pass < passes; ++pass)
        {
            // pass 0: range2.start -> range2.stop ; pass 1 (hysteresis only): range2.stop -> range2.start
            double R2_start = (pass == 0) ? range2.start : range2.stop;
            double R2_stop = (pass == 0) ? range2.stop : range2.start;

            if (pass == 0)
            {
                reset_mag(junctions, 1);
            } // Reset magnetization to initial value at the start of the first pass
            else if (pass == 1)
            {
                reset_mag(junctions, -1);
            } // Reset magnetization to initial value at the start of the second pass if hysteresis is not enabled

            //////////////////////////////////////////////////////
            //                    Range 1 (MPI)
            //////////////////////////////////////////////////////
            for (int i1 = 0; i1 <= N1; ++i1)
            {
                double v1 = getValue(i1, R1_start, R1_stop, N1 + 1);
                // Only overwrite the target when range1 is actually swept; if the flag is off,
                // leave whatever value was set at initialization (e.g. junction Transport.Rp) untouched.
                if (range1.flag)
                    for (double &ref : targets1)
                        ref = v1;

                //////////////////////////////////////////////////////
                //                    Range 2 (Continuous)
                //////////////////////////////////////////////////////
                for (int i2 = 0; i2 <= N2; ++i2)
                {
                    double v2 = getValue(i2, R2_start, R2_stop, N2 + 1);
                    if (range2.flag)
                        for (double &ref : targets2)
                            ref = v2;

                    // Transport/Demag initialisation for every junction and layer (at each rangepoint).
                    // Diagnostic dump (Gp/P/aDL + demag/dipolar tensors) only on the very first
                    // rangepoint of the whole run: those quantities don't change sweep to sweep.
                    const bool logDetails = (it == 0 && pass == 0 && i1 == 0 && i2 == 0);
                    for (size_t ji = 0; ji < junctions.size(); ++ji)
                    {
                        junctions[ji].initializeRangePoint(log, ji, !continuous, logDetails);
                        // Start every rangepoint assuming the ceiling dt for CalcHeff's thermal-field
                        // amplitude; adaptiveTimestep narrows it down from the first step onward.
                        for (auto &layer : junctions[ji].layers)
                            layer.dt = dt;
                    }

                    std::ofstream evolFile;
                    if (saving)
                    {
                        const std::string evolFileName = buildEvolFileName(
                            *this, junctions, it, pass, v1, v2, hysteresis, pass);

                        log << "evol file name = " << evolFileName << "\n";

                        evolFile.open(evolFileName);
                        writeEvolHeader(evolFile, savingColumns);
                    }

                    //////////////////////////////////////////////////////
                    //////////////////////////////////////////////////////
                    //                   Time loop
                    //////////////////////////////////////////////////////
                    //////////////////////////////////////////////////////

                    const int totalSteps = static_cast<int>(std::llround(tmax / dt));
                    double t = 0.0;
                    int step = 0;
                    do
                    {
                        // dtUsed == dt unless adaptiveTimestep narrowed it down this step.
                        const double dtUsed = timeStep(junctions, evolFile, t, step, totalSteps, targets1, targets2, v1, v2, savingRefs);
                        t += dtUsed;
                        ++step;
                    } while (t <= tmax + dt); // Temporal loop

                    if (finalStateSaving)
                    {
                        if (task_id == 0)
                        {
                            writeFinalStateRow(finalStateFile, junctions, v2, v1, pass);
                            drainFinalStateRows(finalStateFile, junctions, finalRowSize); // pick up any rows other ranks sent meanwhile
                        }
                        else
                        {
                            sendFinalStateRow(junctions, v2, v1, pass, finalRowSize);
                        }
                    }

                    // Range hceck
                    // log << "iteration " << it << ", pass " << pass << ", range1 step " << i1 << "/" << N1 << " = " << v1 << ", range2 step " << i2 << "/" << N2 << " = " << v2 << "\n";

                } // Range 2
            } // Range1
        } // Up and down
    } // Iterations

    // Every rank's rangepoints are done: rank 0 keeps receiving until it has heard from
    // all the others; every other rank tells rank 0 it's finished sending rows.
    if (finalStateSaving)
    {
        if (task_id == 0)
        {
            finalizeFinalStateRows(finalStateFile, junctions, finalRowSize, total_task);
        }
        else
        {
            sendFinalStateDone();
        }
    }
}
