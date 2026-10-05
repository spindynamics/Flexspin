#include "objects.hpp"

#include <cmath>

#include "Constants.hpp"

/*******************************************************/
void LLG_Heun1(Junction &junction, int nlayer, int tot_layers, double dt)
{
    auto &l = junction.layers[nlayer];
    const double alpha = l.alpha;
    const double denom = 1.0 + alpha * alpha;

    const Vec3 mp = l.m;

    /******** predictor step ***************/
    const Vec3 fp = (-gamma0 / denom) * (mp * l.Heff + alpha * mp * (mp * l.Heff));

    /******** corrector step ***************/
    const Vec3 mc = (mp + dt * fp);
    const Vec3 fc = (-gamma0 / denom) * (mc * l.Heff + alpha * mc * (mc * l.Heff));

    // Update m final
    const Vec3 mf = mp + 0.5 * dt * (fp + fc);

    // Renormalization |m| = 1
    l.m = mf / sqrt(dot(mf, mf));
}
