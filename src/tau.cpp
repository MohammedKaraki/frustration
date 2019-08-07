#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <ratio>
// #include <filesystem>
#include <algorithm>

#include "generic_lattice.h"

#define FMT_HEADER_ONLY
#include "fmt/format.h"



constexpr int nx = 32, ny = 32, nt = 32;
constexpr double beta = 4.;
constexpr double field_min = 1.;
constexpr double field_max = 5.;
constexpr double field_step = .1;

ising::generic<nx, ny, nt,
    ising::lattice_structure::hexagonal> lattice(beta, field_min);


std::chrono::high_resolution_clock::time_point time_point;
void print_start(std::string process)
{
    std::cerr << process;
    time_point = std::chrono::high_resolution_clock::now();
}
void print_done()
{
    auto t2 = std::chrono::high_resolution_clock::now();

    double dur = std::chrono::duration<double, std::milli>
        (t2 - time_point).count();

    std::cerr << fmt::format("done: {:>10.4G} ms\n", dur);
}


// std::string make_unused_path()
// {
//     constexpr auto N_files_allowed = 100'000;
//     for (int i = 1; i < N_files_allowed; ++i) {
//         std::string path = fmt::format("./data/{:d}.txt", i);
//
//         if (!std::filesystem::exists(path)) {
//             return path;
//         }
//     }
//     std::cerr << "Cannot find unused path\n";
//     exit(1);
// }


double auto_corr_func(const std::vector<double> &data, int t)
{
    const int N = data.size();
    double mu = std::accumulate(data.begin(), data.end(), 0.) / double(N);

    double ret = 0.;
    for (int i = 0; i < N - t; ++i) {
        ret += (data[i] - mu) * (data[i+t] - mu);
    }

    ret /= N - t;
    return ret;
}

// ++(beta, field, N)
// # i    u(i)    acf(i)    acf(i)/acf0     tau(i) 
// --(tau)
int main()
{
    std::string acf_path = "acf_data_2.txt";
    std::ofstream acf_file(acf_path);

    constexpr auto n_thermalize_steps = nx * ny * nt;
    constexpr auto n_rethermalize_steps = nx * ny * nt;

    print_start("thermalizing ... ");
    for (int i = 0; i < n_thermalize_steps; ++i) {
        lattice.cluster_update();
    }
    print_done();

    for (double field = field_min; field <= field_max; field += field_step) {

        print_start(fmt::format("rethermalizing for field = {} ... ", field));
        for (int i = 0; i < n_rethermalize_steps; ++i) {
            lattice.cluster_update();
        }
        print_done();

        constexpr auto N_min = 1'000;
        constexpr auto N_max = 2'000'000;

        bool found_tau = false;
        int n_collected_points = 0;
        std::vector<double> data(N_max);
        int N = N_min;
        while (N <= N_max && !found_tau) {

            data.resize(N);

            print_start(fmt::format("collecting for N = {} ... ", N));
            while (n_collected_points < N) {
                lattice.cluster_update();
                double energy = lattice.energy_density();
                data[n_collected_points] = energy;
                n_collected_points++;
            }
            print_done();

            double tau = 1.;
            double acf0 = auto_corr_func(data, 0);
            const int M_max = N / 100; // make sure M << N
            bool N_was_updated = false;
            for (int M = 0; M <= M_max; ++M) {
                double acf = auto_corr_func(data, M);
                double acf_norm = acf / acf0;

                tau += 2. * acf_norm;

                if (5. * tau < M) {
                    if (N > 900. * tau) {
                        std::cout << fmt::format("found tau={} at (N, M) = "
                                "({}, {})", tau, M, N) << std::endl;
                        found_tau = true;
                    } else {
                        if (N >= N_max) {
                            std::cerr << "failed to find tau\n";
                            break;
                        }
                        N = std::min(N_max, int(1000. * tau));
                        N_was_updated = true;
                        break;
                    }
                }
            }

            tau = 1.;
            acf_file << fmt::format("++({}, {}, {})\n", beta, field, N);
            for (int M = 0; M <= M_max; ++M) {
                double acf = auto_corr_func(data, M);
                double acf_norm = acf / acf0;

                tau += 2. * acf_norm;
                acf_file << fmt::format("{:<10} {:< 20.10E} {:< 20.10E} "
                        "{:< 20.10E} {:< 20.10E}\n",
                        M, data[M], acf, acf_norm, tau);
            }
            if (!N_was_updated) {
                N *= 10;
            }
            acf_file << "\n\n" << std::flush;
        }
    }


//     std::string data_path = make_unused_path();
//     std::string log_path = data_path + ".log";
//     std::ofstream data_file(data_path);
//     std::ofstream log_file(log_path);
//
//     log_file << fmt::format("{} {} {} {} {} {}"
//
//     constexpr int N_thermalize = 20 * 2 * nx * ny * nt;
//     constexpr int N_collect =100'000;
//     constexpr int N_tofile = 1000;
//
//     state("thermalizing");
//     for (int i = 0; i < N_thermalize; ++i) {
//         lattice.cluster_update(nullptr);
//     }
//     done();
//
//
//     std::vector<double> vals;
//     bool found_good_tau = false;
//     for (int N = 1000; N < 10'000'000; N *= 10) {
//         std::cout << "N = " << N << ':' << std::endl;
//
//         vals.resize(N);
//         for (int i = 0; i < N; ++i) {
//             lattice.cluster_update(nullptr);
//             vals[i] = lattice.energy();
//         }
//
//         double acf0 = lattice.acf(vals, 0);
//         double tau = 1.;
//         for (int M = 1; M < N/10; ++M) {
//             std::cout << "\tM = " << M << ':' << std::endl;
//
//             double acf_norm = lattice.acf(vals, M) / acf0;
//             tau += 2. * acf_norm;
//
//             if (M > 5. * tau) {
//                 std::cout << "possibly found" << std::endl;
//                 if (N >= int(900. * tau)) {
//                     found_good_tau = true;
//                     std::cout << "good tau = " << tau << std::endl;
//                     goto EXIT_TAU_SEARCH;
//                 } else {
//                     N = int(1000. * tau);
//                     break;
//                 }
//             }
//         }
//     }
// EXIT_TAU_SEARCH:
//     if (found_good_tau) {
//
//         double acf0 = lattice.acf(vals, 0);
//         double tau = 1.;
//         for (int t = 0; t < N_tofile; ++t) {
//             double acf = lattice.acf(vals, t);
//             double acf_norm = acf / acf0;
//
//             if (t != 0) {
//                 tau += 2. * acf_norm;
//             }
//
//             // data_file << t << " \t" << acf << " \t" << acf_norm <<
//             // " \t" << acf_norm_sum << std::endl;
//             data_file << 
//                 fmt::format("{:<10d} {:< 20.10E} {:< 20.10E} {:< 20.10E}\n",
//                     t, vals[t], acf_norm, tau);
//         }
//     } else {
//     }
// return 0;
//
//
//
//     state("collecting data");
//     std::vector<double> data(N_collect);
//     for (int i = 0; i < N_collect; ++i) {
//         double energy_val;
//
//         lattice.cluster_update(nullptr);
//
//         energy_val = lattice.energy();
//         data[i] = energy_val;
//     }
//     done();
//
//     state("writing data to file");
//     double acf0 = lattice.acf(data, 0);
//     double tau = 1.;
//     for (int t = 0; t < N_tofile; ++t) {
//         double acf = lattice.acf(data, t);
//         double acf_norm = acf / acf0;
//
//         if (t != 0) {
//             tau += 2. * acf_norm;
//         }
//
//         // data_file << t << " \t" << acf << " \t" << acf_norm <<
//             // " \t" << acf_norm_sum << std::endl;
//         data_file << fmt::format("{:<10d} {:< 20.10E} {:< 20.10E} {:< 20.10E}\n",
//                 t, data[t], acf_norm, tau);
//     }
//     done();
//
//
//     return 0;
}
