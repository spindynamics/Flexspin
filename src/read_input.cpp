#include <fstream>
#include <cstdlib>
#include <stdexcept>

#include "objects.hpp"

///////////////////////////////////////////////////////////////////////////////////////
//                  Mathematical quantities (vector, surface, etc.)                  //
///////////////////////////////////////////////////////////////////////////////////////

Vec3 readVec3(
    json* j,
    const std::string& kx,
    const std::string& ky,
    const std::string& kz,
    std::ostream& log,
    const std::string& label)
{
    Vec3 value{jsonDouble(jsonAt(j, kx)),
               jsonDouble(jsonAt(j, ky)),
               jsonDouble(jsonAt(j, kz))};
    log << label << " = {" << value.x << ", " << value.y << ", " << value.z << "}\n";
    return value;
}

///////////////////////////////////////////////////////////////////////////////////////
//                                     Read functions                                //
///////////////////////////////////////////////////////////////////////////////////////

Layer Layer::fromJson(json* j, std::ostream& log, size_t junctionNumber, size_t layerIndex) {
    Layer l{};
    l.Ms = jsonDouble(jsonAt(j, "Ms (A/m)"));
    log << "junction[" << junctionNumber << "].layer[" << layerIndex << "].Ms (A/m) = " << l.Ms << "\n";
    l.alpha = jsonDouble(jsonAt(j, "alpha"));
    log << "junction[" << junctionNumber << "].layer[" << layerIndex << "].alpha = " << l.alpha << "\n";
    l.Ku1 = jsonDouble(jsonAt(j, "Ku (J/m3)"));
    log << "junction[" << junctionNumber << "].layer[" << layerIndex << "].Ku1 (J/m3) = " << l.Ku1 << "\n";
    l.Ku2 = jsonDouble(jsonAt(j, "Ku2 (J/m3)"));
    log << "junction[" << junctionNumber << "].layer[" << layerIndex << "].Ku2 (J/m3) = " << l.Ku2 << "\n";
    l.Ki1 = jsonDouble(jsonAt(j, "Ki (J/m2)"));
    log << "junction[" << junctionNumber << "].layer[" << layerIndex << "].Ki1 (J/m2) = " << l.Ki1 << "\n";
    l.Ki2 = jsonDouble(jsonAt(j, "Ki2 (J/m2)"));
    log << "junction[" << junctionNumber << "].layer[" << layerIndex << "].Ki2 (J/m2) = " << l.Ki2 << "\n";
    l.mca_axis = readVec3(
        jsonAt(j, "MCA axis"),
        "ukx",
        "uky",
        "ukz",
        log,
        "junction[" + std::to_string(junctionNumber) + "].layer[" + std::to_string(layerIndex) + "].MCA axis");
    l.size = 1e-9 * readVec3(
        jsonAt(j, "Layer Size (nm)"),"Lx","Ly","Lz",
        log,
        "junction[" + std::to_string(junctionNumber) + "].layer[" + std::to_string(layerIndex) + "].Layer Size (nm)");
    l.center = 1e-9 * readVec3(
        jsonAt(j, "Layer center (nm)"),"Cx","Cy","Cz",
        log,
        "junction[" + std::to_string(junctionNumber) + "].layer[" + std::to_string(layerIndex) + "].Layer center (nm)");
    l.initial_m = readVec3(
        jsonAt(j, "Initial magnetization"),"mx","my","mz",
        log,
        "junction[" + std::to_string(junctionNumber) + "].layer[" + std::to_string(layerIndex) + "].Initial magnetization");
    l.exchange_bias = readVec3(
        jsonAt(j, "Exchange bias field (A/m)"),"Hx","Hy","Hz",
        log,
        "junction[" + std::to_string(junctionNumber) + "].layer[" + std::to_string(layerIndex) + "].Exchange bias field (A/m)");
    l.vcmaKuCoeff =
        jsonDouble(jsonAt(jsonAt(j, "VCMA (J/(m2*V))"), "KuCoeff"));
    log << "junction[" << junctionNumber << "].layer[" << layerIndex << "].VCMA (J/(m2*V)).KuCoeff = " << l.vcmaKuCoeff << "\n";
    l.vcmaKu2Coeff =
        jsonDouble(jsonAt(jsonAt(j, "VCMA (J/(m2*V))"), "Ku2Coeff"));
    log << "junction[" << junctionNumber << "].layer[" << layerIndex << "].VCMA (J/(m2*V)).Ku2Coeff = " << l.vcmaKu2Coeff << "\n";

    json* sttDamping = jsonAt(j, "STT Damping Like");
    l.sttDL.beta = jsonDouble(jsonAt(sttDamping, "beta"));
    l.sttDL.gamma = jsonDouble(jsonAt(sttDamping, "gamma"));
    log << "junction[" << junctionNumber << "].layer[" << layerIndex << "].STT Damping Like.beta = " << l.sttDL.beta << "\n";
    log << "junction[" << junctionNumber << "].layer[" << layerIndex << "].STT Damping Like.gamma = " << l.sttDL.gamma << "\n";

    json* sttField = jsonAt(j, "STT Field Like");
    l.sttFL.beta = jsonDouble(jsonAt(sttField, "beta"));
    l.sttFL.gamma = jsonDouble(jsonAt(sttField, "gamma"));
    log << "junction[" << junctionNumber << "].layer[" << layerIndex << "].STT Field Like.beta = " << l.sttFL.beta << "\n";
    log << "junction[" << junctionNumber << "].layer[" << layerIndex << "].STT Field Like.gamma = " << l.sttFL.gamma << "\n";

    json* rkky = jsonAt(j, "RKKY (J/m2)");
    json* rkkyIds = jsonAt(rkky, "layer_ID_from");
    json* rkkyConstants = jsonAt(rkky, "Constant");
    size_t rkkyIdCount = jsonSize(rkkyIds);
    size_t rkkyConstCount = jsonSize(rkkyConstants);
    size_t rkkyCount = (rkkyIdCount < rkkyConstCount) ? rkkyIdCount : rkkyConstCount;
    l.rkky.layerIDFrom.reserve(rkkyCount);
    l.rkky.constant.reserve(rkkyCount);
    for (size_t idx = 0; idx < rkkyCount; ++idx) {
        l.rkky.layerIDFrom.push_back(jsonInt(jsonIndex(rkkyIds, idx)));
        l.rkky.constant.push_back(jsonDouble(jsonIndex(rkkyConstants, idx)));
        log << "junction[" << junctionNumber << "].layer[" << layerIndex << "].RKKY (J/m2)[" << idx
            << "] = {layer_ID_from=" << l.rkky.layerIDFrom[idx]
            << ", Constant=" << l.rkky.constant[idx] << "}\n";
    }

    json* flags = jsonAt(j, "Flags");
    l.flags.free = jsonBool(jsonAt(flags, "Free"));
    l.flags.dipolar = jsonBool(jsonAt(flags, "Dipolar"));
    l.flags.callenCallen = jsonBool(jsonAt(flags, "Callen-Callen"));
    l.flags.tmrTempDep = jsonBool(jsonAt(flags, "TMR_Temp_dep"));
    l.flags.sttTempDep = jsonBool(jsonAt(flags, "STT_Temp_dep"));
    log << "junction[" << junctionNumber << "].layer[" << layerIndex << "].Flags = {Free=" << l.flags.free
        << ", Dipolar=" << l.flags.dipolar
        << ", Callen-Callen=" << l.flags.callenCallen
        << ", TMR_Temp_dep=" << l.flags.tmrTempDep
        << ", STT_Temp_dep=" << l.flags.sttTempDep << "}\n";

    l.Rp =
        jsonDouble(jsonAt(jsonAt(j, "Transport"), "Rp"));
    log << "junction[" << junctionNumber << "].layer[" << layerIndex << "].Transport.Rp = " << l.Rp << "\n";
    l.TMR =
        jsonDouble(jsonAt(jsonAt(j, "Transport"), "TMR"));
    log << "junction[" << junctionNumber << "].layer[" << layerIndex << "].Transport.TMR = " << l.TMR << "\n";
    return l;
}

Temperature Temperature::fromJson(json* j, std::ostream& log, size_t junctionNumber) {
    Temperature t{};
    t.T0 = jsonDouble(jsonAt(j, "T0 (K)"));
    log << "junction[" << junctionNumber << "].Temperature.T0 (K) = " << t.T0 << "\n";
    t.C = jsonDouble(jsonAt(j, "C"));
    log << "junction[" << junctionNumber << "].Temperature.C = " << t.C << "\n";
    t.Q = jsonDouble(jsonAt(j, "Q"));
    log << "junction[" << junctionNumber << "].Temperature.Q = " << t.Q << "\n";
    t.alpha = jsonDouble(jsonAt(j, "alpha"));
    log << "junction[" << junctionNumber << "].Temperature.alpha = " << t.alpha << "\n";
    t.Tb = jsonDouble(jsonAt(j, "Tb (K)"));
    log << "junction[" << junctionNumber << "].Temperature.Tb (K) = " << t.Tb << "\n";
    t.Vh = jsonDouble(jsonAt(j, "Vh (V)"));
    log << "junction[" << junctionNumber << "].Temperature.Vh (V) = " << t.Vh << "\n";
    return t;
}

Pulse Pulse::fromJson(json* j, std::ostream& log, const std::string& label) {
    Pulse p{};
    p.flag = jsonBool(jsonAt(j, "flag"));
    p.type = jsonInt(jsonAt(j, "type"));
    p.Tdelay = jsonDouble(jsonAt(j, "Tdelay"));
    p.Trise = jsonDouble(jsonAt(j, "Trise"));
    p.Tpulse = jsonDouble(jsonAt(j, "Tpulse"));
    p.Tfall = jsonDouble(jsonAt(j, "Tfall"));
    log << label << " = {type=" << p.type
        << ", Tdelay=" << p.Tdelay
        << ", Trise=" << p.Trise
        << ", Tpulse=" << p.Tpulse
        << ", Tfall=" << p.Tfall
        << "}\n";
    return p;
}

BiasChannel BiasChannel::fromJson(json* j, std::ostream& log, const std::string& label) {
    BiasChannel b{};
    b.flag = jsonBool(jsonAt(j, "flag"));
    log << label << ".flag = " << b.flag << "\n";

    if (jsonContains(j, "A")) {
        b.A = jsonDouble(jsonAt(j, "A"));
        log << label << ".A = " << b.A << "\n";
    }
    if (jsonContains(j, "V")) {
        b.V = jsonDouble(jsonAt(j, "V"));
        log << label << ".V = " << b.V << "\n";
    }
    if (jsonContains(j, "phi")) {
        b.phi = jsonDouble(jsonAt(j, "phi"));
        log << label << ".phi = " << b.phi << "\n";
    }
    if (jsonContains(j, "theta")) {
        b.theta = jsonDouble(jsonAt(j, "theta"));
        log << label << ".theta = " << b.theta << "\n";
    }
    if (jsonContains(j, "H")) {
        b.H = jsonDouble(jsonAt(j, "H"));
        log << label << ".H = " << b.H << "\n";
    }
    if (jsonContains(j, "F")) {
        b.F = jsonDouble(jsonAt(j, "F"));
        log << label << ".F = " << b.F << "\n";
    }

    return b;
}

SystemBias SystemBias::fromJson(json* j, std::ostream& log, const std::string& label) {
    SystemBias b{};
    b.dcCurrent = BiasChannel::fromJson(jsonAt(j, "DC current"), log, label + ".DC current");
    b.dcVoltage = BiasChannel::fromJson(jsonAt(j, "DC voltage"), log, label + ".DC voltage");
    b.dcField = BiasChannel::fromJson(jsonAt(j, "DC field"), log, label + ".DC field");
    b.acCurrent = BiasChannel::fromJson(jsonAt(j, "AC current"), log, label + ".AC current");
    b.acVoltage = BiasChannel::fromJson(jsonAt(j, "AC voltage"), log, label + ".AC voltage");
    b.acField = BiasChannel::fromJson(jsonAt(j, "AC field"), log, label + ".AC field");
    return b;
}

SweepRange SweepRange::fromJson(json* j, std::ostream& log, const std::string& label) {
    SweepRange r{};
    r.flag = jsonBool(jsonAt(j, "flag"));
    r.type = jsonString(jsonAt(j, "type"));
    r.start = jsonDouble(jsonAt(j, "start"));
    r.stop = jsonDouble(jsonAt(j, "stop"));
    r.step = jsonInt(jsonAt(j, "step"));
    log << label << " = {flag=" << r.flag
        << ", type=" << r.type
        << ", start=" << r.start
        << ", stop=" << r.stop
        << ", step=" << r.step
        << "}\n";
    r.pulse = Pulse::fromJson(jsonAt(j, "pulse (ns)"), log, label + ".pulse (ns)");
    return r;
}

Simulation Simulation::fromJson(json* j, std::ostream& log) {
    Simulation s{};
    s.dt = jsonDouble(jsonAt(j, "dt (s)"));
    log << "simulation.dt (s) = " << s.dt << "\n";
    s.adaptiveTimestep = jsonBool(jsonAt(j, "adapatative"));
    log << "simulation.adapatative = " << s.adaptiveTimestep << "\n";
    s.thetaMax = jsonDouble(jsonAt(j, "theta max (rad)"));
    log << "simulation.theta max (rad) = " << s.thetaMax << "\n";
    s.dtMin = jsonDouble(jsonAt(j, "dt min (s)"));
    log << "simulation.dt min (s) = " << s.dtMin << "\n";
    s.samplingTime = jsonDouble(jsonAt(j, "sampling time (s)"));
    log << "simulation.sampling time (s) = " << s.samplingTime << "\n";
    s.tmax = jsonDouble(jsonAt(j, "tmax (s)"));
    log << "simulation.tmax (s) = " << s.tmax << "\n";
    s.saving = jsonBool(jsonAt(j, "saving"));
    log << "simulation.saving = " << s.saving << "\n";
    json* savingColumnsArr = jsonAt(j, "saving columns");
    size_t savingColumnCount = jsonSize(savingColumnsArr);
    s.savingColumns.reserve(savingColumnCount);
    for (size_t idx = 0; idx < savingColumnCount; ++idx) {
        s.savingColumns.push_back(jsonString(jsonIndex(savingColumnsArr, idx)));
        log << "simulation.saving columns[" << idx << "] = " << s.savingColumns.back() << "\n";
    }
    s.finalStateSaving = jsonBool(jsonAt(j, "final-state saving"));
    log << "simulation.final-state saving = " << s.finalStateSaving << "\n";
    s.rootName = jsonString(jsonAt(j, "root name"));
    log << "simulation.root name = " << s.rootName << "\n";
    s.iterations = jsonInt(jsonAt(j, "iterations"));
    log << "simulation.iterations = " << s.iterations << "\n";
    s.continuous = jsonBool(jsonAt(j, "Continuous"));
    log << "simulation.Continuous = " << s.continuous << "\n";
    s.hysteresis = jsonBool(jsonAt(j, "hysteresis"));
    log << "simulation.hysteresis = " << s.hysteresis << "\n";
    s.thermalField = jsonBool(jsonAt(j, "thermal field"));
    log << "simulation.thermal field = " << s.thermalField << "\n";
    s.solver = jsonString(jsonAt(j, "solver (Heun1 ,Heun2, RK4, Symplectic(!))"));
    log << "simulation.solver (Heun1 ,Heun2, RK4, Symplectic(!)) = " << s.solver << "\n";
    s.systemBias = SystemBias::fromJson(jsonAt(j, "system-bias"), log, "simulation.system-bias");
    s.range1 = SweepRange::fromJson(jsonAt(j, "range1"), log, "simulation.range1");
    s.range2 = SweepRange::fromJson(jsonAt(j, "range2"), log, "simulation.range2");
    return s;
}

///////////////////////////////////////////////////////////////////////////////////////
//                       List parsers (loop over a json array)                       //
///////////////////////////////////////////////////////////////////////////////////////

// Parses every element of a "layer" array until none are left.
std::vector<Layer> parseLayerList(json* j, std::ostream& log, size_t junctionNumber) {
    std::vector<Layer> layers;
    for (size_t layerIndex = 0; layerIndex < jsonSize(j); ++layerIndex) {
        layers.push_back(Layer::fromJson(jsonIndex(j, layerIndex), log, junctionNumber, layerIndex));
    }
    return layers;
}

// One junction = one mandatory "layer" block, optionally followed by a "Temperature"
// block and a "junction-bias" block. blockIndex is advanced past whatever was consumed.
Junction Junction::fromJson(json* j, std::ostream& log, size_t& blockIndex, size_t junctionNumber) {
    Junction junction{};
    
    junction.layers = parseLayerList(jsonAt(jsonIndex(j, blockIndex), "layer"), log, junctionNumber);
    ++blockIndex;
    junction.temperature = Temperature::fromJson(jsonAt(jsonIndex(j, blockIndex), "Temperature"), log, junctionNumber);
    ++blockIndex;    
    junction.junctionBias = SystemBias::fromJson(jsonAt(jsonIndex(j, blockIndex), "junction-bias"), log,
     "junction[" + std::to_string(junctionNumber) + "].junction-bias");
    junction.hasJunctionBias = true;
    ++blockIndex;    

    return junction;
}

// Repeatedly calls Junction::fromJson, each call consuming the next 1-3 blocks of the
// "junction" array, until the whole array has been consumed.
std::vector<Junction> parseJunctionList(json* j, std::ostream& log) {
    std::vector<Junction> junctions;
    size_t blockIndex = 0;
    while (blockIndex < jsonSize(j)) {
        junctions.push_back(Junction::fromJson(j, log, blockIndex, junctions.size()));
    }
    return junctions;
}
