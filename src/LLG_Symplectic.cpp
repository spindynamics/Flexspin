#include "objects.hpp"
#include <cmath>
#include "Constants.hpp"

/********************************************************/
void AdvanceIndividualSpin(Junction &junction, int nlayer, double timestep)
{
    auto &l = junction.layers[nlayer];
    const double alpha = l.alpha;
    const double fact = gamma0 / (1.0 + alpha * alpha);

    const Vec3 m = l.m;
    const Vec3 Htot = l.Heff;

    // Calculate omega
    const Vec3 mXH = m * Htot;
    const Vec3 omega = fact * (Htot + alpha * mXH);

    // Calculate intermediate variables
    const double norm_omega2 = dot(omega, omega);
    const double omegam = dot(omega, m);
    const Vec3 omegaXm = omega * m;
    const double denom = 1.0 + timestep * timestep * 0.25 * norm_omega2;

    // Update m
    l.m = (m + timestep * omegaXm +
           (timestep * timestep * 0.25) * (2.0 * omega * omegam - norm_omega2 * m)) /
          denom;
}

/********************************************************/
void LLG_Symplectic(Junction &junction, int nlayer1, int nlayer2, int tot_layers, double dt)
{
    // Update the effective fields for both layers
    CalcHeff(junction, nlayer1, tot_layers);
    CalcHeff(junction, nlayer2, tot_layers);

    /*************************************/
    /* Advance M2 by dt/2 */
    /*************************************/
    AdvanceIndividualSpin(junction, nlayer2, dt / 2.0);
    CalcHeff(junction, nlayer2, tot_layers); // Update the effective field for layer nlayer2

    /*************************************/
    /* Advance M1 by dt */
    /*************************************/
    AdvanceIndividualSpin(junction, nlayer1, dt);
    CalcHeff(junction, nlayer1, tot_layers); // Update the effective field for layer nlayer1

    /*************************************/
    /* Advance M2 by dt/2 */
    /*************************************/
    AdvanceIndividualSpin(junction, nlayer2, dt / 2.0);
    CalcHeff(junction, nlayer2, tot_layers); // Update the effective field for layer nlayer2
}
