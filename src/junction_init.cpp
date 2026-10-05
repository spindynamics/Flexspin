#include "objects.hpp"

#include <cmath>

#include "Constants.hpp"

namespace {

Vec3 sphericalToCartesian(double amplitude, double theta, double phi) {
	Vec3 v{};
	v.x = amplitude * sin(theta) * cos(phi);
	v.y = amplitude * sin(theta) * sin(phi);
	v.z = amplitude * cos(theta);
	return v;
}

} // namespace

void Layer::updateResistance(const Layer& next) {
	R = 2.0 / (Gap + Gp + (Gp - Gap) * dot(m, next.m));
}

void Junction::initializeRangePoint(std::ostream& log, size_t junctionIndex, bool resetMagnetization, bool logDetails) {

    // Reset the self-heating temperature to the ambient bath value at the start of
    // every rangepoint sweep (T is otherwise never initialized before first use).
    T = temperature.T0;

    //////////////////////////////////////////////////////
    //          Voltage, Current, Field
    //////////////////////////////////////////////////////
    // DC current / voltage / field
	Iapp = junctionBias.dcCurrent.A;
	Vapp = junctionBias.dcVoltage.V;
	// dcField.H/acField.H are entered in Tesla; convert B->H (H = B/mu0) here so Happ/H_ac
	// end up in A/m, matching every other field CalcHeff sums into Heff.
	Happ = sphericalToCartesian(
		junctionBias.dcField.H / mu0,
		junctionBias.dcField.theta,
		junctionBias.dcField.phi);

	// AC current / voltage / field
	Iac = junctionBias.acCurrent.A;
	f_I = junctionBias.acCurrent.F;
	Vac = junctionBias.acVoltage.V;
	f_V = junctionBias.acVoltage.F;
	f_H = junctionBias.acField.F;
	H_ac = sphericalToCartesian(
		junctionBias.acField.H / mu0,
		junctionBias.acField.theta,
		junctionBias.acField.phi);


    //////////////////////////////////////////////////////
    //          Transport layer per layer
    //////////////////////////////////////////////////////
	for (size_t layerIndex = 0; layerIndex < layers.size(); ++layerIndex) {
		Layer& layer = layers[layerIndex];

		layer.Gp = (layer.Rp != 0.0) ? (1.0 / layer.Rp) : 0.0;
		layer.Gap = layer.Gp / (1.0 + layer.TMR); // Rap = Rp * (1 + TMR)
		layer.P = sqrt(layer.TMR * (layer.TMR + 2.0) / (2.0 * (layer.TMR + 1.0)));

		// Damping-like STT torque asymmetry coefficient
		const double volume = layer.size.x * layer.size.y * layer.size.z;
		layer.aDL_0 = (h_bar * layer.P * layer.Gp) / ((-2.0 * e_charge) * layer.Ms * mu0 * volume);
		layer.aDL = layer.aDL_0;


        ////////////////////////////////////////////////////// 
        //          Demag tensor layer to layer
        //////////////////////////////////////////////////////
		layer.ND.clear();
		layer.ND.reserve(layers.size());
		for (size_t j2 = 0; j2 < layers.size(); ++j2) {
			layer.ND.push_back(coeffdemag(layer, layers[j2]));
		}

        // One-time diagnostic dump: transport/STT coefficients plus the self-demag
        // and mutual-dipolar tensors, none of which change from one swept point to
        // the next, so printing them on every rangepoint would just be noise.
        if (logDetails) {
			log << "junction[" << junctionIndex << "].layer[" << layerIndex << "]"
				<< " Gp=" << layer.Gp
				<< " P=" << layer.P
				<< " aDL=" << layer.aDL
				<< "\n";

			for (size_t j2 = 0; j2 < layer.ND.size(); ++j2) {
				const tensor& N = layer.ND[j2];
				log << "junction[" << junctionIndex << "].layer[" << layerIndex << "]."
					<< (j2 == layerIndex ? "Ndemag" : "Ndipolar") << "[" << j2 << "]"
					<< " xx=" << N.xx << " yy=" << N.yy << " zz=" << N.zz
					<< " xy=" << N.xy << " xz=" << N.xz << " yz=" << N.yz
					<< "\n";
			}
		}
	}
}
void initializePulse(double t,
	                 const SweepRange& range1,
	                 const SweepRange& range2,
	                 DoubleRefs& targets1,
	                 DoubleRefs& targets2,
	                 double value1,
	                 double value2) {

	const double gate1 = pulse(range1.pulse.type, t, range1.pulse.Tdelay, range1.pulse.Trise, range1.pulse.Tpulse, range1.pulse.Tfall);
	const double gate2 = pulse(range2.pulse.type, t, range2.pulse.Tdelay, range2.pulse.Trise, range2.pulse.Tpulse, range2.pulse.Tfall);

	const double appliedGate1 = range1.pulse.flag ? gate1 : 1.0;
	const double appliedGate2 = range2.pulse.flag ? gate2 : 1.0;

	for (double& ref : targets1) {
		ref = appliedGate1 * value1;
	}
	for (double& ref : targets2) {
		ref = appliedGate2 * value2;
	}
}
