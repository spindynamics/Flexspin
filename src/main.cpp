#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <mpi.h>

#include "objects.hpp"
#include "Constants.hpp"

using namespace std;
using namespace std::chrono;

/*****************************************************************************/
void version(void)
{
    printf("===================================================\n");
    printf("Flexspin : Stochastic LLG program for multiple junction and layers \n");
    printf("May 2026\n");
    printf("Arthur Courberand\n");
    printf("===================================================\n");
}

int main(int argc, char **argv)
{

    // MPI initialization
    MPI_Init(&argc, &argv);
    int task_id, total_task;
    MPI_Comm_rank(MPI_COMM_WORLD, &task_id);    // Get task id (which task THIS FILE will be)
    MPI_Comm_size(MPI_COMM_WORLD, &total_task); // Get total number of task
    if (task_id == 0)
    { // Print the world size to be sure
        cout << " Node size: " << total_task << endl;
    }
    // So on a 10 nodes standard Q, total_task will be 10, and task_id will be 0 to 9
    MPI_Barrier(MPI_COMM_WORLD);

    ////////////////////////////////////////////////////////////////////////////////////
    //               Read input file and parse it into structured data                //
    ////////////////////////////////////////////////////////////////////////////////////
    std::string filePath = (argc > 1) ? std::string(argv[1]) : std::string("input.json"); // Either a specified file or "input.json"
    std::ifstream file(filePath);                                                         // Open the file
    std::ostringstream fileContent;                                                       // Take its content
    fileContent << file.rdbuf();                                                          // This write the content somewhere
    json *root = json_tokener_parse(fileContent.str().c_str());                           // This reconstruct the json structure (this is json-c library)

    // Open the log file to check values read (only on task 0, so every MPI rank
    // doesn't fight over the same file and stomp on each other's writes)
    std::ofstream log;
    if (task_id == 0)
    {
        log.open("out.evol");
        log << std::unitbuf; // Auto-flush on every write, so out.evol is live/complete
                             // even while the run is still going or gets killed mid-run.
    }

    // Make the import
    Simulation simulation = Simulation::fromJson(jsonAt(root, "simulation"), log);      // Import from the json the simulation parameters (dt, tmax, sampling time, etc...)
    std::vector<Junction> junctions = parseJunctionList(jsonAt(root, "junction"), log); // Import from the json every junction (layers, temperature, bias, etc...)

    json_object_put(root); // Done with the json tree, free it

    size_t N_junction = junctions.size();

    // Number of layers in each junction: N_layer[i] == junctions[i].layers.size()
    std::vector<size_t> N_layer(N_junction);
    for (size_t i = 0; i < N_junction; ++i)
    {
        N_layer[i] = junctions[i].layers.size();
    }

    // Seed each layer's thermal-noise generator once (decorrelated across MPI ranks
    // via task_id*3), same scheme as before. Also clone dt/thermalField from Simulation
    // onto each layer (constant for the whole run), so CalcHeff can see them without
    // needing them threaded through every LLG solver's signature.
    for (auto &junction : junctions)
    {
        for (auto &layer : junction.layers)
        {
            srand(time(NULL) + task_id * 3);
            layer.rngState = -rand();
            layer.dt = simulation.dt;
            layer.thermalField = simulation.thermalField;
        }
    }

    /////////////////////////////////////////////////////////////////////////////
    //                             Program                                     //
    /////////////////////////////////////////////////////////////////////////////

    simulation.range(junctions, log, task_id, total_task);

    MPI_Finalize();
    return 0;
}