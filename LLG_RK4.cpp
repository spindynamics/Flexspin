#include "objects.hpp"

#include <cmath>

const double gamma0 = 2.210173e5; /*  m/(As)  */

/*******************************************************/
void LLG_RK4(Junction& junction, int nlayer, int tot_layers, double dt)
{
    auto& l = junction.layers[nlayer];
    const double alpha = l.alpha;
    const double denom = 1.0 + alpha * alpha;

    const Vec3 mp = l.m;

    /******** step 1 **************/
    const Vec3 f1 = (-gamma0 / denom) * (mp*l.Heff + alpha * mp*(mp*l.Heff));
    const Vec3 m1 = mp + 0.5 * dt * f1;

    /******** step 2 **************/
    l.m = m1;
    CalcHeff(junction, nlayer, tot_layers); // Update Heff with intermediate m
    const Vec3 f2 = (-gamma0 / denom) * (m1*l.Heff + alpha * m1*(m1*l.Heff));
    const Vec3 m2 = mp + 0.5 * dt * f2;

    /******** step 3 **************/
    l.m = m2;
    CalcHeff(junction, nlayer, tot_layers); // Update Heff with intermediate m
    const Vec3 f3 = (-gamma0 / denom) * (m2*l.Heff + alpha * m2*(m2*l.Heff));
    const Vec3 m3 = mp + dt * f3;

    /******** step 4 **************/
    l.m = m3;
    CalcHeff(junction, nlayer, tot_layers); // Update Heff with intermediate m
    const Vec3 f4 = (-gamma0 / denom) * (m3*l.Heff + alpha * m3*(m3*l.Heff));

    // Update m final
    const Vec3 mf = mp + (dt / 6.0) * (f1 + 2.0 * f2 + 2.0 * f3 + f4);

    // Renormalization |m| = 1
    l.m = mf / sqrt(dot(mf, mf));
}
