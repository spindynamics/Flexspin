#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "objects.hpp"

// Generic "object path" resolver: walks a dot-separated path of real C++ member names
// (with [index] or [x] for vectors) down the Simulation/Junction object graph, and
// returns every live double& the path reaches.
// "[x]" is a wildcard meaning "every element": it fans out and the results of each
// branch get concatenated, so e.g. "junction[x].layers[x].size.z" returns every
// layer's size.z, in every junction.
// Examples: "simulation.systemBias.dcCurrent.A", "junction[1].layers[2].size.z"

namespace
{

    struct PathToken
    {
        std::string name;
        bool hasIndex = false;
        bool isWildcard = false;
        size_t index = 0;
    };

    std::vector<PathToken> parsePath(const std::string &path)
    {
        std::vector<PathToken> tokens;
        std::stringstream ss(path);
        std::string segment;
        while (std::getline(ss, segment, '.'))
        {
            PathToken tok;
            size_t bracket = segment.find('[');
            if (bracket != std::string::npos)
            {
                size_t close = segment.find(']', bracket);
                tok.name = segment.substr(0, bracket);
                std::string indexText = segment.substr(bracket + 1, close - bracket - 1);
                if (indexText == "x")
                {
                    tok.isWildcard = true;
                }
                else
                {
                    tok.index = static_cast<size_t>(std::stoul(indexText));
                    tok.hasIndex = true;
                }
            }
            else
            {
                tok.name = segment;
            }
            tokens.push_back(tok);
        }
        return tokens;
    }

    [[noreturn]] void fail(const std::string &path, const std::string &reason)
    {
        throw std::runtime_error("Invalid path '" + path + "': " + reason);
    }

    DoubleRefs resolve(Vec3 &v, const std::string &path, const std::vector<PathToken> &tokens, size_t pos)
    {
        if (pos != tokens.size() - 1)
            fail(path, "'" + tokens[pos].name + "' is not a leaf field");
        const std::string &name = tokens[pos].name;
        if (name == "x")
            return {v.x};
        if (name == "y")
            return {v.y};
        if (name == "z")
            return {v.z};
        fail(path, "Vec3 has no field '" + name + "'");
    }

    DoubleRefs resolve(BiasChannel &b, const std::string &path, const std::vector<PathToken> &tokens, size_t pos)
    {
        if (pos != tokens.size() - 1)
            fail(path, "'" + tokens[pos].name + "' is not a leaf field");
        const std::string &name = tokens[pos].name;
        if (name == "As")
            return {b.As};
        if (name == "A")
            return {b.A};
        if (name == "V")
            return {b.V};
        if (name == "phi")
            return {b.phi};
        if (name == "theta")
            return {b.theta};
        if (name == "H")
            return {b.H};
        if (name == "F")
            return {b.F};
        fail(path, "BiasChannel has no double field '" + name + "'");
    }

    DoubleRefs resolve(SystemBias &sb, const std::string &path, const std::vector<PathToken> &tokens, size_t pos)
    {
        const std::string &name = tokens[pos].name;
        bool isLeaf = (pos == tokens.size() - 1);
        if (name == "Resistance")
        {
            if (!isLeaf)
                fail(path, "'Resistance' is a leaf field");
            return {sb.Resistance};
        }
        if (isLeaf)
            fail(path, "'" + name + "' is not a leaf field");
        size_t next = pos + 1;
        if (name == "dcCurrent")
            return resolve(sb.dcCurrent, path, tokens, next);
        if (name == "dcVoltage")
            return resolve(sb.dcVoltage, path, tokens, next);
        if (name == "dcField")
            return resolve(sb.dcField, path, tokens, next);
        if (name == "acCurrent")
            return resolve(sb.acCurrent, path, tokens, next);
        if (name == "acVoltage")
            return resolve(sb.acVoltage, path, tokens, next);
        if (name == "acField")
            return resolve(sb.acField, path, tokens, next);
        fail(path, "SystemBias has no field '" + name + "'");
    }

    DoubleRefs resolve(Temperature &t, const std::string &path, const std::vector<PathToken> &tokens, size_t pos)
    {
        if (pos != tokens.size() - 1)
            fail(path, "'" + tokens[pos].name + "' is not a leaf field");
        const std::string &name = tokens[pos].name;
        if (name == "T0")
            return {t.T0};
        if (name == "C")
            return {t.C};
        if (name == "Q")
            return {t.Q};
        if (name == "alpha")
            return {t.alpha};
        if (name == "Tb")
            return {t.Tb};
        if (name == "Vh")
            return {t.Vh};
        fail(path, "Temperature has no field '" + name + "'");
    }

    DoubleRefs resolve(Pulse &p, const std::string &path, const std::vector<PathToken> &tokens, size_t pos)
    {
        if (pos != tokens.size() - 1)
            fail(path, "'" + tokens[pos].name + "' is not a leaf field");
        const std::string &name = tokens[pos].name;
        if (name == "Tdelay")
            return {p.Tdelay};
        if (name == "Trise")
            return {p.Trise};
        if (name == "Tpulse")
            return {p.Tpulse};
        if (name == "Tfall")
            return {p.Tfall};
        fail(path, "Pulse has no double field '" + name + "'");
    }

    DoubleRefs resolve(SweepRange &r, const std::string &path, const std::vector<PathToken> &tokens, size_t pos)
    {
        const std::string &name = tokens[pos].name;
        bool isLeaf = (pos == tokens.size() - 1);
        if (name == "start")
        {
            if (!isLeaf)
                fail(path, "'start' is a leaf field");
            return {r.start};
        }
        if (name == "stop")
        {
            if (!isLeaf)
                fail(path, "'stop' is a leaf field");
            return {r.stop};
        }
        if (name == "pulse")
        {
            if (isLeaf)
                fail(path, "'pulse' is not a leaf field");
            return resolve(r.pulse, path, tokens, pos + 1);
        }
        fail(path, "SweepRange has no field '" + name + "'");
    }

    DoubleRefs resolve(Layer &l, const std::string &path, const std::vector<PathToken> &tokens, size_t pos)
    {
        const std::string &name = tokens[pos].name;
        bool isLeaf = (pos == tokens.size() - 1);
        if (name == "Ms")
        {
            if (!isLeaf)
                fail(path, "'Ms' is a leaf field");
            return {l.Ms};
        }
        if (name == "alpha")
        {
            if (!isLeaf)
                fail(path, "'alpha' is a leaf field");
            return {l.alpha};
        }
        if (name == "Ku1")
        {
            if (!isLeaf)
                fail(path, "'Ku1' is a leaf field");
            return {l.Ku1};
        }
        if (name == "Ku2")
        {
            if (!isLeaf)
                fail(path, "'Ku2' is a leaf field");
            return {l.Ku2};
        }
        if (name == "Ki1")
        {
            if (!isLeaf)
                fail(path, "'Ki1' is a leaf field");
            return {l.Ki1};
        }
        if (name == "Ki2")
        {
            if (!isLeaf)
                fail(path, "'Ki2' is a leaf field");
            return {l.Ki2};
        }
        if (name == "mca_axis")
        {
            if (isLeaf)
                fail(path, "'mca_axis' is not a leaf field");
            return resolve(l.mca_axis, path, tokens, pos + 1);
        }
        if (name == "size")
        {
            if (isLeaf)
                fail(path, "'size' is not a leaf field");
            return resolve(l.size, path, tokens, pos + 1);
        }
        if (name == "center")
        {
            if (isLeaf)
                fail(path, "'center' is not a leaf field");
            return resolve(l.center, path, tokens, pos + 1);
        }
        if (name == "initial_m")
        {
            if (isLeaf)
                fail(path, "'initial_m' is not a leaf field");
            return resolve(l.initial_m, path, tokens, pos + 1);
        }
        if (name == "m")
        {
            if (isLeaf)
                fail(path, "'m' is not a leaf field");
            return resolve(l.m, path, tokens, pos + 1);
        }
        if (name == "exchange_bias")
        {
            if (isLeaf)
                fail(path, "'exchange_bias' is not a leaf field");
            return resolve(l.exchange_bias, path, tokens, pos + 1);
        }
        if (name == "vcmaKuCoeff")
        {
            if (!isLeaf)
                fail(path, "'vcmaKuCoeff' is a leaf field");
            return {l.vcmaKuCoeff};
        }
        if (name == "vcmaKu2Coeff")
        {
            if (!isLeaf)
                fail(path, "'vcmaKu2Coeff' is a leaf field");
            return {l.vcmaKu2Coeff};
        }
        if (name == "Rp")
        {
            if (!isLeaf)
                fail(path, "'Rp' is a leaf field");
            return {l.Rp};
        }
        if (name == "TMR")
        {
            if (!isLeaf)
                fail(path, "'TMR' is a leaf field");
            return {l.TMR};
        }
        if (name == "Gp")
        {
            if (!isLeaf)
                fail(path, "'Gp' is a leaf field");
            return {l.Gp};
        }
        if (name == "P")
        {
            if (!isLeaf)
                fail(path, "'P' is a leaf field");
            return {l.P};
        }
        fail(path, "Layer has no field '" + name + "'");
    }

    DoubleRefs resolve(Junction &j, const std::string &path, const std::vector<PathToken> &tokens, size_t pos)
    {
        const std::string &name = tokens[pos].name;
        bool isLeaf = (pos == tokens.size() - 1);
        if (name == "layers")
        {
            if (!tokens[pos].hasIndex && !tokens[pos].isWildcard)
                fail(path, "'layers' needs an index or [x], e.g. layers[0] or layers[x]");
            if (isLeaf)
                fail(path, "'layers[...]' is not a leaf field");
            DoubleRefs results;
            if (tokens[pos].isWildcard)
            {
                for (size_t idx = 0; idx < j.layers.size(); ++idx)
                {
                    DoubleRefs sub = resolve(j.layers[idx], path, tokens, pos + 1);
                    results.insert(results.end(), sub.begin(), sub.end());
                }
            }
            else
            {
                if (tokens[pos].index >= j.layers.size())
                    fail(path, "layer index out of range");
                results = resolve(j.layers[tokens[pos].index], path, tokens, pos + 1);
            }
            return results;
        }
        if (name == "temperature")
        {
            if (isLeaf)
                fail(path, "'temperature' is not a leaf field");
            return resolve(j.temperature, path, tokens, pos + 1);
        }
        if (name == "junctionBias")
        {
            if (isLeaf)
                fail(path, "'junctionBias' is not a leaf field");
            return resolve(j.junctionBias, path, tokens, pos + 1);
        }
        if (name == "R")
        {
            if (!isLeaf)
                fail(path, "'R' is a leaf field");
            return {j.R};
        }
        if (name == "T")
        {
            if (!isLeaf)
                fail(path, "'T' is a leaf field");
            return {j.T};
        }
        if (name == "Itot")
        {
            if (!isLeaf)
                fail(path, "'Itot' is a leaf field");
            return {j.Itot};
        }
        if (name == "Vtot")
        {
            if (!isLeaf)
                fail(path, "'Vtot' is a leaf field");
            return {j.Vtot};
        }
        if (name == "Iapp")
        {
            if (!isLeaf)
                fail(path, "'Iapp' is a leaf field");
            return {j.Iapp};
        }
        if (name == "Vapp")
        {
            if (!isLeaf)
                fail(path, "'Vapp' is a leaf field");
            return {j.Vapp};
        }
        if (name == "Happ")
        {
            if (isLeaf)
                fail(path, "'Happ' is not a leaf field");
            return resolve(j.Happ, path, tokens, pos + 1);
        }
        fail(path, "Junction has no field '" + name + "'");
    }

    DoubleRefs resolve(Simulation &s, const std::string &path, const std::vector<PathToken> &tokens, size_t pos)
    {
        const std::string &name = tokens[pos].name;
        bool isLeaf = (pos == tokens.size() - 1);
        if (name == "dt")
        {
            if (!isLeaf)
                fail(path, "'dt' is a leaf field");
            return {s.dt};
        }
        if (name == "samplingTime")
        {
            if (!isLeaf)
                fail(path, "'samplingTime' is a leaf field");
            return {s.samplingTime};
        }
        if (name == "tmax")
        {
            if (!isLeaf)
                fail(path, "'tmax' is a leaf field");
            return {s.tmax};
        }
        if (name == "systemBias")
        {
            if (isLeaf)
                fail(path, "'systemBias' is not a leaf field");
            return resolve(s.systemBias, path, tokens, pos + 1);
        }
        if (name == "range1")
        {
            if (isLeaf)
                fail(path, "'range1' is not a leaf field");
            return resolve(s.range1, path, tokens, pos + 1);
        }
        if (name == "range2")
        {
            if (isLeaf)
                fail(path, "'range2' is not a leaf field");
            return resolve(s.range2, path, tokens, pos + 1);
        }
        fail(path, "Simulation has no double field '" + name + "'");
    }

} // namespace

DoubleRefs resolvePath(Simulation &simulation, std::vector<Junction> &junctions, const std::string &path)
{
    std::vector<PathToken> tokens = parsePath(path);
    if (tokens.empty())
        fail(path, "empty path");

    const std::string &root = tokens[0].name;
    if (root == "simulation")
    {
        if (tokens.size() == 1)
            fail(path, "'simulation' is not a leaf field");
        return resolve(simulation, path, tokens, 1);
    }
    if (root == "junction")
    {
        if (!tokens[0].hasIndex && !tokens[0].isWildcard)
            fail(path, "'junction' needs an index or [x], e.g. junction[0] or junction[x]");
        if (tokens.size() == 1)
            fail(path, "'junction[...]' is not a leaf field");
        DoubleRefs results;
        if (tokens[0].isWildcard)
        {
            for (Junction &j : junctions)
            {
                DoubleRefs sub = resolve(j, path, tokens, 1);
                results.insert(results.end(), sub.begin(), sub.end());
            }
        }
        else
        {
            if (tokens[0].index >= junctions.size())
                fail(path, "junction index out of range");
            results = resolve(junctions[tokens[0].index], path, tokens, 1);
        }
        return results;
    }
    fail(path, "unknown root '" + root + "' (expected 'simulation' or 'junction[...]')");
}
