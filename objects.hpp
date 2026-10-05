#pragma once

#include <string>
#include <vector>
#include <stdexcept>
#include <ostream>
#include <functional>

#include <json-c/json.h>

using json = json_object;

///////////////////////////////////////////////////////////////////////////////////////
// Small json-c convenience layer (replaces nlohmann's .at()/.get<T>()/.contains())   //
///////////////////////////////////////////////////////////////////////////////////////

inline json* jsonAt(json* obj, const std::string& key) {
    json_object* value = nullptr;
    if (!obj || !json_object_object_get_ex(obj, key.c_str(), &value)) {
        throw std::out_of_range("Missing JSON key: " + key);
    }
    return value;
}

inline bool jsonContains(json* obj, const std::string& key) {
    return obj && json_object_object_get_ex(obj, key.c_str(), nullptr);
}

inline double jsonDouble(json* obj) { return json_object_get_double(obj); }
inline int jsonInt(json* obj) { return json_object_get_int(obj); }
inline bool jsonBool(json* obj) { return json_object_get_boolean(obj) != 0; }
inline std::string jsonString(json* obj) { return json_object_get_string(obj); }
inline size_t jsonSize(json* obj) { return static_cast<size_t>(json_object_array_length(obj)); }
inline json* jsonIndex(json* obj, size_t idx) { return json_object_array_get_idx(obj, static_cast<int>(idx)); }

///////////////////////////////////////////////////////////////////////////////////////
//                  Mathematical quantities (vector, surface, etc.)                  //
///////////////////////////////////////////////////////////////////////////////////////

struct Vec3{
    double x;
    double y;
    double z;
};

inline Vec3 makeVec(double x, double y, double z) {
    return {x, y, z};
}

inline Vec3 operator+(const Vec3& a, const Vec3& b) {
    return makeVec(a.x + b.x, a.y + b.y, a.z + b.z);
}

inline Vec3 operator-(const Vec3& a, const Vec3& b) {
    return makeVec(a.x - b.x, a.y - b.y, a.z - b.z);
}

inline Vec3& operator+=(Vec3& a, const Vec3& b) {
    a = a + b;
    return a;
}

inline Vec3& operator-=(Vec3& a, const Vec3& b) {
    a = a - b;
    return a;
}

inline Vec3 operator*(double s, const Vec3& v) {
    return makeVec(s * v.x, s * v.y, s * v.z);
}

inline Vec3 operator*(const Vec3& v, double s) {
    return s * v;
}

inline Vec3 operator/(const Vec3& v, double s) {
    return makeVec(v.x / s, v.y / s, v.z / s);
}

inline double dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return makeVec(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x);
}

inline Vec3 operator*(const Vec3& a, const Vec3& b) {
    return cross(a, b);
}

struct tensor {
    double xx;
    double xy;
    double xz;
    double yx;
    double yy;
    double yz;
    double zx;
    double zy;
    double zz;
};

inline Vec3 tensorTimesVector(const tensor& N, const Vec3& v) {
    return makeVec(
        N.xx * v.x + N.xy * v.y + N.xz * v.z,
        N.yx * v.x + N.yy * v.y + N.yz * v.z,
        N.zx * v.x + N.zy * v.y + N.zz * v.z);
}

inline Vec3 diagonalTensorTimesVector(const tensor& N, const Vec3& v) {
    return makeVec(N.xx * v.x, N.yy * v.y, N.zz * v.z);
}


///////////////////////////////////////////////////////////////////////////////////////
//                                 Program entities                                  //
///////////////////////////////////////////////////////////////////////////////////////

struct Flags {
    bool free = false;
    bool dipolar = false;
    bool callenCallen = false;
    bool tmrTempDep = false;
    bool sttTempDep = false;
};

struct STTTerm {
    double beta = 0.0;
    double gamma = 0.0;
};

struct RKKY {
    std::vector<int> layerIDFrom;
    std::vector<double> constant;
};

struct Layer{
    double Ms;
    double alpha;
    double Ku1;
    double Ku2;
    double Ki1;
    double Ki2;
    Vec3 mca_axis;
    Vec3 size;
    Vec3 center;
    Vec3 initial_m;
    Vec3 m;
    Vec3 exchange_bias;
    double vcmaKuCoeff;
    double vcmaKu2Coeff;
    double Rp;
    double TMR;
    double Gp;
    double Gap;
    double R;
    double P;
    double aDL_0;
    double aDL;
    STTTerm sttDL;
    STTTerm sttFL;
    RKKY rkky;
    Flags flags;
    std::vector<tensor> ND;

    Vec3 Hanis;
    Vec3 Hdem;
    Vec3 Hdip;
    Vec3 Hrkky;
    Vec3 Hst;
    Vec3 Hfl;
    Vec3 Hdl;
    Vec3 Hth; // Thermal fluctuation field (Langevin), recomputed every CalcHeff call
    Vec3 Heff;

    // Persistent ran2/gasdev seed for this layer's thermal fluctuation field.
    // Must be seeded once (to a negative value) before the first gasdev() call,
    // then left alone: ran2 evolves it in place across calls.
    long rngState = 0;

    // Cloned once from Simulation (constant for the whole run), so CalcHeff can see
    // them without needing them threaded through every LLG solver's signature.
    double dt = 0.0;
    bool thermalField = false;

    // Refreshes R (this layer's interface resistance to the next layer) from Gp/Gap
    // and the angle between this layer's and next layer's magnetization.
    void updateResistance(const Layer& next);

    static Layer fromJson(json* j, std::ostream& log, size_t junctionNumber, size_t layerIndex);
};

struct Temperature{
    double T0;
    double C;
    double Q;
    double alpha;
    double Tb;
    double Vh;

    static Temperature fromJson(json* j, std::ostream& log, size_t junctionNumber);
};

struct Pulse{
    int type;
    bool flag;
    double Tdelay;
    double Trise;
    double Tpulse;
    double Tfall;

    static Pulse fromJson(json* j, std::ostream& log, const std::string& label);
};

struct BiasChannel{
    bool flag;
    double As;
    double A;
    double V;
    double phi;
    double theta;
    double H;
    double F;

    static BiasChannel fromJson(json* j, std::ostream& log, const std::string& label);
};

struct SystemBias{
    BiasChannel dcCurrent;
    BiasChannel dcVoltage;
    BiasChannel dcField;
    BiasChannel acCurrent;
    BiasChannel acVoltage;
    BiasChannel acField;

    // Equivalent resistance of every junction wired in parallel, i.e. 1 / sum(1/R_i).
    // Only meaningful on Simulation::systemBias (updated at the end of each timeStep()
    // from that step's junction.R values, then read back at the start of the *next* call
    // to split systemBias.dcCurrent.A across junctions -- one dt behind, negligible given
    // how small dt is here). Unused on Junction::junctionBias, which shares this type.
    double Resistance = 0.0;

    static SystemBias fromJson(json* j, std::ostream& log, const std::string& label);
};

struct SweepRange{
    bool flag;
    std::string type;
    double start;
    double tempo_start;
    double stop;
    double tempo_stop;
    int step;
    Pulse pulse;

    static SweepRange fromJson(json* j, std::ostream& log, const std::string& label);
};

struct Junction{
    std::vector<Layer> layers;
    Temperature temperature;
    double T;
    SystemBias junctionBias;
    bool hasJunctionBias = false;

    // Cached values evaluated at each range point.
    double Iapp = 0.0;
    double Iac = 0.0;
    double f_I = 0.0;
    Vec3 Happ{0.0, 0.0, 0.0};
    Vec3 Htot{0.0, 0.0, 0.0};
    double Vapp = 0.0;
    double Vac = 0.0;
    double f_V = 0.0;
    double f_H = 0.0;
    Vec3 H_ac{0.0, 0.0, 0.0};
    double R;

    // Recomputed every timeStep() (not just once per range point): Vapp/Iapp above are
    // this junction's own DC bias, Vtot/Itot fold in the shared systemBias applied across
    // every junction wired in parallel (same node voltage added to all; systemBias's
    // current split across junctions in proportion to each one's conductance 1/R).
    double Vtot = 0.0;
    double Itot = 0.0;


    // Per-range-point initialization for junction and layer derived parameters.
    // logDetails gates the one-time diagnostic dump (Gp/P/aDL + demag/dipolar tensors),
    // so it's printed once for the whole run instead of on every swept point.
    void initializeRangePoint(std::ostream& log, size_t junctionIndex, bool resetMagnetization, bool logDetails);

    // j = the full "junction" array; blockIndex is advanced past the 1-3 blocks (layer [+ Temperature] [+ junction-bias]) this junction consumes
    static Junction fromJson(json* j, std::ostream& log, size_t& blockIndex, size_t junctionNumber);
};

using DoubleRefs = std::vector<std::reference_wrapper<double>>;

struct Simulation{
    double dt; // fixed step when !adaptiveTimestep; the step's ceiling when adaptiveTimestep is set
    bool adaptiveTimestep;
    double thetaMax; // max precession angle allowed per step (rad), only used if adaptiveTimestep
    double dtMin;    // floor for the adaptive step (s), only used if adaptiveTimestep
    double samplingTime;
    double tmax;
    bool saving;
    std::vector<std::string> savingColumns; // object paths (resolvePath grammar), e.g. "junction[0].layers[0].m.z"
    bool finalStateSaving;
    std::string rootName;
    int iterations;
    bool continuous;
    bool hysteresis;
    bool thermalField;
    std::string solver;
    SystemBias systemBias;
    SweepRange range1;
    SweepRange range2;

    static Simulation fromJson(json* j, std::ostream& log);
    static void reset_mag(std::vector<Junction>& junctions, int sign);

    // Splits range1 across MPI tasks, then runs: for each iteration, for each
    // hysteresis pass (forward, and backward if hysteresis is set), sweep range1 x range2.
    // If finalStateSaving is set, only rank 0 ever touches "1-State-Diagram.txt": every
    // other rank streams its rows to rank 0 over MPI as they're computed, instead of all
    // ranks writing the file directly (which isn't safe on a shared/network filesystem).
    void range(std::vector<Junction>& junctions, std::ostream& log, int task_id, int total_task);

    // Advances every junction/layer by one step: pulse update, effective-field recompute,
    // and one LLG solver step per free layer. If adaptiveTimestep is set, the step size is
    // chosen so no free layer precesses by more than thetaMax this step (bounded to
    // [dtMin, dt]); otherwise it's always exactly dt. Returns the step size actually used,
    // so the caller can advance its own time accumulator by the right amount.
    // If saving is enabled, appends one row of savingRefs' current values at time t to evol,
    // but only every samplingTime worth of physical time (plus always the first and last step).
    double timeStep(std::vector<Junction>& junctions, std::ostream& evol, double t, int step, int totalSteps,
                  DoubleRefs& targets1, DoubleRefs& targets2, double v1, double v2, DoubleRefs& savingRefs);
};

// Generic object-path resolver: e.g. "simulation.systemBias.dcCurrent.A" or "junction[1].layers[2].size.z"
// "[x]" instead of a numeric index is a wildcard meaning "every element" (fans out and
// concatenates results), e.g. "junction[x].layers[x].size.z" returns every layer's size.z,
// in every junction.
DoubleRefs resolvePath(Simulation& simulation, std::vector<Junction>& junctions, const std::string& path);

///////////////////////////////////////////////////////////////////////////////////////
//                       List parsers (loop over a json array)                       //
///////////////////////////////////////////////////////////////////////////////////////

std::vector<Layer> parseLayerList(json* j, std::ostream& log, size_t junctionNumber);
std::vector<Junction> parseJunctionList(json* j, std::ostream& log);

double getValue(int i, double Val_Start, double Val_Stop, int NVal) noexcept;
double pulse(int typepulse, double t, double Tdelai, double Trise, double Tpulse, double Tfall) noexcept;
tensor coeffdemag(const Layer& sourceLayer, const Layer& targetLayer);

// L'Ecuyer ran2 uniform deviate + Box-Muller gasdev Gaussian deviate (Numerical Recipes style).
// idum must be seeded to a negative value before the first call, then left alone across calls.
double ran2(long* idum);
double gasdev(long* idum);

// Effective field and LLG solvers (one layer, indexed by nlayer, within a junction)
void CalcHeff(Junction& junction, int nlayer, int tot_layers);
void LLG_Heun1(Junction& junction, int nlayer, int tot_layers, double dt);
void LLG_Heun2(Junction& junction, int nlayer, int tot_layers, double dt);
void LLG_RK4(Junction& junction, int nlayer, int tot_layers, double dt);
void LLG_Symplectic(Junction& junction, int nlayer1, int nlayer2, int tot_layers, double dt);
void initializePulse(double t,
                     const SweepRange& range1,
                     const SweepRange& range2,
                     DoubleRefs& targets1,
                     DoubleRefs& targets2,
                     double value1,
                     double value2);
std::string buildEvolFileName(const Simulation& simulation,const std::vector<Junction>& junctions,
                              int simulationIndex,
                              int passIndex,
                              double fieldValue,
                              double voltageValue,
                              bool appendDirectionSuffix,
                              int directionIndex) noexcept;
