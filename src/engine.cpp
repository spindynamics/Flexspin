#include "printing.hpp"

#include <algorithm>
#include <cmath>
#include "Constants.hpp"

// Advances every junction/layer by one step: pulse update, effective-field recompute, and
// one LLG solver step per free layer. If saving is enabled, appends one row of every
// saving-columns value at time t to evol, but only every samplingTime worth of physical
// time (plus always the first and last step of the run). Returns the step size actually
// used (== dt unless adaptiveTimestep is set).
double Simulation::timeStep(std::vector<Junction>& junctions, std::ostream& evol, double t, int step, int totalSteps,
                           DoubleRefs& targets1, DoubleRefs& targets2, double v1, double v2, DoubleRefs& savingRefs) {

    //Initialise pulse for range values (no shard change in A/V/H etc)
    initializePulse(t, range1, range2, targets1, targets2, v1, v2);

    //Total conductance for current distribution (reduced to systemBias.Resistance = 1/G at the end)
    double G_all_junction_next = 0.0;





    //For every junction: effective field + interface resistance (doesn't depend on the step size)
    for (size_t ji = 0; ji < junctions.size(); ++ji) {
        Junction& junction = junctions[ji];
        const int tot_layers = static_cast<int>(junction.layers.size());

        const double sharedCurrent = (systemBias.Resistance > 0.0 && junction.R > 0.0)
            ? systemBias.dcCurrent.A * systemBias.Resistance / junction.R
            : 0.0;
        junction.Itot = junction.Iapp + sharedCurrent;
        junction.Vtot = junction.Vapp + systemBias.dcVoltage.V + junction.Itot * junction.R;

        junction.R = 0.0;
        //For every layer
        for (size_t li = 0; li < junction.layers.size(); ++li) {
            const int nlayer = static_cast<int>(li);

            CalcHeff(junction, nlayer, tot_layers);

            if (li + 1 < junction.layers.size()) {
                junction.layers[li].updateResistance(junction.layers[li + 1]);
                junction.R += junction.layers[li].R;
            }
        }//Layer Heff
    }//Junction Heff

    // Torque-based adaptive step: cap the precession angle of the fastest free layer to
    // thetaMax this step, bounded to [dtMin, dt] (dt is the ceiling here). One dt shared
    // by every layer this step, so the most-restrictive (largest-torque) layer decides it.
    double dt_step = dt;
    if (adaptiveTimestep) {
        double omega_max = 0.0;
        for (const auto& junction : junctions) {
            for (const auto& layer : junction.layers) {
                if (!layer.flags.free) continue;
                const Vec3 torque = layer.m * layer.Heff; 
                omega_max = std::max(omega_max, gamma0 * std::sqrt(dot(torque, torque)));
            }
        }
        if (omega_max > 0.0) dt_step = std::min(std::max(thetaMax / omega_max, dtMin), dt);

        // Feeds CalcHeff's thermal-field amplitude on the *next* step (one step behind,
        // same approximation already made for systemBias.Resistance below).
        for (auto& junction : junctions)
            for (auto& layer : junction.layers)
                layer.dt = dt_step;
    }

    //For every junction: LLG step + self-heating (both use this step's dt_step)
    for (size_t ji = 0; ji < junctions.size(); ++ji) {
        Junction& junction = junctions[ji];
        const int tot_layers = static_cast<int>(junction.layers.size());

        for (size_t li = 0; li < junction.layers.size(); ++li) {
            if (!junction.layers[li].flags.free) continue;
            const int nlayer = static_cast<int>(li);

            if (solver == "Heun1") {LLG_Heun1(junction, nlayer, tot_layers, dt_step);}
            else if (solver == "Heun2") {LLG_Heun2(junction, nlayer, tot_layers, dt_step);}
            else if (solver == "RK4") {LLG_RK4(junction, nlayer, tot_layers, dt_step);}
            //else if (solver == "Symplectic") {LLG_Symplectic(junction, nlayer, tot_layers, dt_step);}

            //Callen-Callen Not included for now
            //CalcCallenCallen(junction, nlayer, tot_layers);
        }//Layer Temporal

        junction.T += ((dt_step * junction.Vtot * junction.Vtot) / (junction.temperature.C * junction.R)) - (dt_step * junction.temperature.Q * (junction.T - junction.temperature.T0) / junction.temperature.C);
        if (junction.R > 0.0) G_all_junction_next += 1.0 / junction.R;

    }//Junction loop
    systemBias.Resistance = (G_all_junction_next > 0.0) ? 1.0 / G_all_junction_next : 0.0;

    // Only every samplingTime worth of physical time (always the first and last step).
    bool shouldSave;
    if (adaptiveTimestep) {
        const long long idxBefore = static_cast<long long>(t / samplingTime);
        const long long idxAfter  = static_cast<long long>((t + dt_step) / samplingTime);
        const bool isLastStep = (t + dt_step >= tmax);
        shouldSave = saving && (step == 0 || isLastStep || idxAfter > idxBefore);
    } else {
        const int sampleEvery = std::max(1, static_cast<int>(std::llround(samplingTime / dt)));
        shouldSave = saving && (step == 0 || step == totalSteps || step % sampleEvery == 0);
    }
    if (shouldSave) writeEvolRow(evol, savingRefs, t * 1e9); //Time in ns

    return dt_step;
}

void Simulation::reset_mag(std::vector<Junction>& junctions, int sign) {
    for (auto& junction : junctions) {
        for (auto& layer : junction.layers) {
            layer.m = sign * layer.initial_m; // Reset magnetization to initial value
        }
    }
}
