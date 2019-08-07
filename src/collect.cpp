#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <fstream>
#include <string>

#define FMT_HEADER_ONLY
#include "fmt/format.h"

#include "generic_lattice.h"

#ifndef L // typically will be defined from the shell when compiling
#define L 16
#endif

constexpr unsigned int nt = 32;
constexpr unsigned int nx = L, ny = L;

const double beta = 4.;
const double field_min = 0.5,
      field_max = 5.,
      field_step = .1;

ising::generic<nx, ny, nt,
    ising::lattice_structure::square> lattice(beta, field_min);


void out(std::ostream &file_out, const std::string &str) {
    std::cout << str;
    file_out << str;
}


void thermalize(std::string file_name) {

    std::string data_path = std::move(file_name);
    std::ofstream data_file(data_path);

    constexpr auto thermalize_steps = L * 50;
    double dur_sum = 0.;

    out(data_file, fmt::format("{:<10} {:>20} {:>20} {:>20}\n",
                "step", "energy density", "mag density", "duration"));
    for (int i = 0; i < thermalize_steps; ++i) {
        double dur;

        timing::push_time_point();
        lattice.cluster_update();
        dur = timing::pop_time_point();
        dur_sum += dur;

        if ( (i%(thermalize_steps/100)) == 0 ) {
            out(data_file,
                    fmt::format("{:<10} {:> 20.7E} {:> 20.7E} {:> 20.7E}\n",
                        i, lattice.energy_density(),
                        lattice.magnetization_density(), dur));
        }
    }
    out(data_file, fmt::format("\ntotal thermalization time is {:< 20.7E}\n\n",
                dur_sum));
}

auto measure(int autocorr_time, int n_measure) {
    std::vector<double> energies, mags;
    double ene, mag;
    for (int i = 0; i < n_measure; ++i) {
        ene = mag = 0.;
        for (int t = 0; t < autocorr_time; ++t) {
            lattice.cluster_update();
            ene += lattice.energy_density();
            mag += lattice.magnetization_density();
        }
        ene /= autocorr_time;
        mag /= autocorr_time;
        energies.push_back(ene);
        mags.push_back(mag);
    }

    return std::make_tuple(stats::mean(energies),
            stats::standard_error(energies),
            stats::mean(mags),
            stats::standard_error(mags),
            stats::binder_cumulant(energies),
            stats::binder_cumulant(mags));
}


void collect(std::string file_name) {

    std::string data_path = std::move(file_name);
    std::ofstream data_file(data_path);

    out(data_file, 
            fmt::format("{:<16} {:>16} {:>16} {:>16} {:>16} {:>16} {:>16}\n",
                "field", "energy", "energy error", "magnetization",
                "mag error", "ene binder", "mag binder"));

    for (double field = field_min; field < field_max; field += field_step) {
        lattice.set_params(beta, field);

        // re-thermalize
        constexpr auto rethermalize_steps = 10 * L;
        for (int i = 0; i < rethermalize_steps; ++i) {
            lattice.cluster_update();
        }

        constexpr auto autocorr_time = std::max(20, L/2);
        constexpr auto n_measure = 200;
        double energy, energy_err, mag, mag_err, binder_energy, binder_mag;

        std::tie(energy, energy_err, mag, mag_err, binder_energy, binder_mag)
            = measure(autocorr_time, n_measure);

        out(data_file,
                fmt::format("{:< 16.7E} {:> 16.7E} {:> 16.7E} {:> 16.7E} "
                    "{:> 16.7E} {:> 16.7E} {:> 16.7E}\n", field, energy,
                    energy_err, mag, mag_err, binder_energy, binder_mag));
    }
}

int main() {
    thermalize(fmt::format("data/thermalization-L{}.sq.txt", L));
    collect(fmt::format("data/collecting-L{}.sq.txt", L));

    return 0;
}

