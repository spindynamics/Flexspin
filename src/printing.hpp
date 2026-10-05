#pragma once

#include "objects.hpp"

#include <ostream>
#include <vector>

///////////////////////////////////////////////////////////////////////////////////////
// evol file: one row per sampled timestep, per rangepoint. Columns are whatever       //
// "saving columns" in the input resolved to (see resolvePath in objects.hpp) --       //
// labels and refs are resolved once, outside the time loop, and reused every row.      //
///////////////////////////////////////////////////////////////////////////////////////

// t_ns column, then one column per label (the resolvePath path string itself, used as-is).
void writeEvolHeader(std::ostream& evol, const std::vector<std::string>& labels);

// t_ns column, then the current value of each pre-resolved ref, in the same order as
// the labels writeEvolHeader was given.
void writeEvolRow(std::ostream& evol, const DoubleRefs& values, double t_ns);

///////////////////////////////////////////////////////////////////////////////////////
// 1-State-Diagram.txt: one row per rangepoint, appended once its temporal loop        //
// settles. Only rank 0 ever opens/writes this file directly: MPI ranks can run on     //
// different nodes over a shared/network filesystem, where append-at-EOF isn't atomic  //
// across separate client handles the way local O_APPEND is, so every rank writing it  //
// directly could interleave mid-write and corrupt it. Every other rank instead streams //
// its row to rank 0 over MPI via sendFinalStateRow()/sendFinalStateDone(); rank 0      //
// receives with drainFinalStateRows() (non-blocking, between its own rangepoints) and  //
// finalizeFinalStateRows() (blocking, once its own sweep is done).                    //
///////////////////////////////////////////////////////////////////////////////////////

// Column names must match, in order, whatever writeFinalStateRow()/the MPI-received rows
// write: H, V, per-layer m={mx,my,mz}, direction ("up"/"down"), then per-junction T, R.
void writeFinalStateHeader(std::ostream& out, const std::vector<Junction>& junctions);

// Rank 0 only: writes its own row directly.
void writeFinalStateRow(std::ostream& out, const std::vector<Junction>& junctions,
                         double fieldValue, double voltageValue, int pass);

// Number of doubles in one row as sent over MPI. Fixed for a given run, since every rank
// parses the same junctions from the same input file.
size_t finalStateRowSize(const std::vector<Junction>& junctions);

// Rank 0 only. Non-blocking: writes any rows that have already arrived from other ranks
// without waiting for more, so incoming rows get appended promptly (between rank 0's own
// rangepoints) instead of only being picked up at the very end.
void drainFinalStateRows(std::ostream& out, const std::vector<Junction>& junctions, size_t rowSize);

// Rank 0 only, called once its own rangepoints are all done. Keeps receiving rows/"done"
// signals from every other rank until all total_task-1 of them have reported done, so no
// row sent late by a slower rank is lost.
void finalizeFinalStateRows(std::ostream& out, const std::vector<Junction>& junctions,
                             size_t rowSize, int total_task);

// Non-root ranks only: packs and sends one row to rank 0.
void sendFinalStateRow(const std::vector<Junction>& junctions,
                        double fieldValue, double voltageValue, int pass, size_t rowSize);

// Non-root ranks only: tells rank 0 this rank has no more rows coming.
void sendFinalStateDone();
