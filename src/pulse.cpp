#include "objects.hpp"
#include <cmath>

#include "Constants.hpp"

double getValue(int i, double Val_Start, double Val_Stop, int NVal) noexcept {
	if (NVal <= 1) return Val_Start;
	double step = (Val_Stop - Val_Start) / (NVal - 1);
	return Val_Start + i * step;
}

static inline double safe_div(double num, double den) noexcept {
	return (den == 0.0) ? ( (num>0.0) ? INFINITY : (num<0.0 ? -INFINITY : 0.0) ) : (num/den);
}

double pulse(int typepulse, double t, double Tdelai, double Trise, double Tpulse, double Tfall) noexcept
{
	const double period = Tdelai + Trise + Tpulse + Tfall;
	if (period <= 0.0) return 0.0;

	double cycles = std::floor(t / period);
	double tp = t - cycles * period; // position inside current period (>=0)

	switch (typepulse) {
	case 1: // linear variation
		if (tp <= Tdelai) return 0.0;
		if (tp <= Tdelai + Trise) return safe_div(tp - Tdelai, Trise);
		if (tp > Tdelai + Trise + Tpulse) return safe_div(period - tp, Tfall);
		return 1.0;

	case 2: // sinusoidal variation type 1
		if (tp <= Tdelai) return 0.0;
		if (tp <= Tdelai + Trise) return std::sin((M_PI/2.0) * safe_div(tp - Tdelai, Trise));
		if (tp > Tdelai + Trise + Tpulse) return std::sin((M_PI/2.0) * safe_div(period - tp, Tfall));
		return 1.0;

	case 3: // sinusoidal variation type 2
		if (tp <= Tdelai) return 0.0;
		if (tp <= Tdelai + Trise) return 0.5 * (1.0 - std::sin(M_PI * safe_div(tp - Tdelai + Trise/2.0, Trise)));
		if (tp > Tdelai + Trise + Tpulse) return 0.5 * (1.0 + std::sin(M_PI * safe_div(tp - (Tdelai + Trise + Tpulse) + Tfall/2.0, Tfall)));
		return 1.0;

	case 4: // variant with explicit zero outside pulse (robust handling)
		if (tp <= Tdelai) return 0.0;
		if (tp <= Tdelai + Trise) return 0.5 * (1.0 - std::sin(M_PI * safe_div(tp - Tdelai + Trise/2.0, Trise)));
		if (tp <= Tdelai + Trise + Tpulse) return 1.0;
		if (tp <= period) return 0.5 * (1.0 + std::sin(M_PI * safe_div(tp - (Tdelai + Trise + Tpulse) + Tfall/2.0, Tfall)));
		return 0.0;

	default:
		return 0.0;
	}
}
