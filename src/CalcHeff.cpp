#include "objects.hpp"
#include <cmath>
#include "Constants.hpp"

void CalcHeff(Junction &junction, int nlayer, int tot_layers)
{
    auto &l = junction.layers[nlayer];
    auto &j = junction;

    //////////////////////////////////////////////////////////////////////
    // Thermal fluctuation field (Langevin), drawn fresh every call
    //////////////////////////////////////////////////////////////////////
    l.Hfl = makeVec(0.0, 0.0, 0.0); // Field-like STT torque: not implemented yet

    if (l.thermalField)
    {
        const double volume = l.size.x * l.size.y * l.size.z;
        const double Hth_ampl = sqrt(2.0 * l.alpha * kB * j.T / (gamma0 * mu0 * volume * l.Ms * l.dt));
        l.Hth = makeVec(Hth_ampl * gasdev(&l.rngState),
                        Hth_ampl * gasdev(&l.rngState),
                        Hth_ampl * gasdev(&l.rngState));
    }
    else
    {
        l.Hth = makeVec(0.0, 0.0, 0.0);
    }

    const double uk_m = dot(l.mca_axis, l.m);
    const double uk_m3 = uk_m * uk_m * uk_m;
    const double volDen = mu0 * l.Ms;

    //////////////////////////////////////////////////////////////////////
    // Internal fields
    //////////////////////////////////////////////////////////////////////
    const double bulkAnis = (2.0 * l.Ku1 * uk_m + 4.0 * l.Ku2 * uk_m3) / volDen;
    const double intAnis = (2.0 * (l.Ki1 / l.size.z) * uk_m + 4.0 * (l.Ki2 / l.size.z) * uk_m3) / volDen;
    l.Hanis = (bulkAnis + intAnis) * l.mca_axis;

    l.Hdem = makeVec(0.0, 0.0, 0.0);
    l.Hdip = makeVec(0.0, 0.0, 0.0);
    l.Hrkky = makeVec(0.0, 0.0, 0.0);

    //////////////////////////////////////////////////////////////////////
    // External fields (except demag)
    //////////////////////////////////////////////////////////////////////
    for (int index = 0; index < tot_layers; ++index)
    {
        const auto &layer_j = junction.layers[index];

        // If the layer is itself it is demag
        if (index == nlayer)
        {
            l.Hdem += (-layer_j.Ms) * diagonalTensorTimesVector(l.ND[index], layer_j.m);
            continue;
        }

        // If the layer is not itself it is dipolar
        if (l.flags.dipolar == 1)
        {
            l.Hdip += (-layer_j.Ms) * tensorTimesVector(l.ND[index], layer_j.m);
        }

        // RKKY is checking is the layer tested is in iteraction as put in put file
        // and using the corresponding constant
        double Jrkky = 0.0;
        for (size_t rkkyIndex = 0; rkkyIndex < l.rkky.layerIDFrom.size(); ++rkkyIndex)
        {
            if (l.rkky.layerIDFrom[rkkyIndex] == index)
            {
                Jrkky = l.rkky.constant[rkkyIndex];
                break;
            }
        }
        if (Jrkky != 0.0)
        {
            l.Hrkky += (Jrkky / (mu0 * l.Ms * l.size.z)) * layer_j.m;
        }
    }

    //////////////////////////////////////////////////////////////////////
    // Spin Transfer Torque (STT) #and Field Like Torque (FLT) (maybe)
    //////////////////////////////////////////////////////////////////////

    l.Hst = makeVec(0.0, 0.0, 0.0);

    if (nlayer > 0)
    {
        // Direct STT
        l.Hst += l.aDL * (j.Vtot * l.sttDL.beta - l.sttDL.gamma * j.Vtot * j.Vtot) * l.m * j.layers[nlayer - 1].m;
    }

    if (nlayer < tot_layers - 1)
    {
        // Back-scatter STT
        l.Hst -= l.aDL * (j.Vtot * l.sttDL.beta - l.sttDL.gamma * j.Vtot * j.Vtot) * l.m * j.layers[nlayer + 1].m;
    }

    //////////////////////////////////////////////////////////////////////
    // Adding every fields
    //////////////////////////////////////////////////////////////////////
    l.Heff = l.Hdip + l.Hst + l.Hfl + l.Hdem + j.Happ + j.H_ac + l.Hanis + l.Hrkky + l.Hth;
}
