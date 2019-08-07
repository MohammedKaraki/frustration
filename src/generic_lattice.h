#ifndef GENERIC_LATTICE
#define GENERIC_LATTICE

// modify this line if you want to use mersenne-twister
#define rand_uint random_gen::xorshift

#include <sstream>
#include <vector>
#include <array>
#include <cmath>
#include <random>
#include <chrono>


namespace random_gen {

    // fast random number generator (xorshift algorithm)
    unsigned int xorshift() {
        static unsigned int x=123456789, y=362436069, z=521288629;
        unsigned int t;
        x ^= x << 16;
        x ^= x >> 5;
        x ^= x << 1;
        t = x;
        x = y;
        y = z;
        z = t ^ x ^ y;
        return z;
    }

    // mersenne twister
    static constexpr unsigned int mt_seed = 1239u;
    std::mt19937 mt{mt_seed};
}


// easy-to-use timing functionality
namespace timing {

    std::vector<std::chrono::high_resolution_clock::time_point> time_points;

    void push_time_point() {
        time_points.push_back(std::chrono::high_resolution_clock::now());
    }

    double last_elapsed() {
        auto t2 = std::chrono::high_resolution_clock::now();
        auto t1 = time_points.back();
        return std::chrono::duration<double>(t2 - t1).count();
    }

    double pop_time_point() {
        double ret = last_elapsed();
        time_points.pop_back();
        return ret;
    }
}


namespace stats {

    template<typename T>
    double mean(const T &v) {
        return std::accumulate(v.begin(), v.end(), 0.) / v.size();
    }

    template<typename T>
    double standard_error(const T &v) {
        double sample_mean = mean(v);
        double sq_sum = 0.;
        for (auto it = v.begin(); it != v.end(); ++it) {
            double val = *it;
            sq_sum += val*val;
        }
        double N = v.size();
        double sq_mean = sq_sum / N;
        double tmp = (sq_mean - sample_mean * sample_mean) / (N-1);

        // tmp could be slightly negative, when all elements are equal,
        // due to the finite precision
        if (tmp < 0.) 
            tmp = 0.;

        return std::sqrt(tmp);
    } 

    template<typename T>
    double binder_cumulant(const T &v) {
        double s2_sum = 0.;
        double s4_sum = 0.;

        for (auto s : v) {
            auto s2 = s * s;
            s2_sum += s2;
            s4_sum += s2 * s2;
        }

        double N = v.size();
        double s2_bar = s2_sum / N;
        double s4_bar = s4_sum / N;

        return 1. - (s4_bar / (s2_bar * s2_bar)) / 3.;

    }
}

namespace ising {

    enum class lattice_structure {hexagonal, square};

    template<int nx, int ny, int nt, lattice_structure structure>
    class generic {
    private:
        static constexpr unsigned int n_sites = nx * ny * nt;
        std::array<char, n_sites> state;

    private:
        double beta, field;
        double eta; // function of beta and field (see set_params()).

        // used in cluster_update() and determine the probability of adding
        // a link to the cluster
        double p_xy;
        double p_t;

    private:

        class coords {
        public:
            coords(unsigned int index) {
                t = index / (nx*ny);
                index %= (nx*ny);
                y = index / nx;
                x = index % nx;
            }

            coords(unsigned char x_p, unsigned char y_p, unsigned char t_p)
                : x{x_p}, y{y_p}, t{t_p} { }

            unsigned int get_index() {
                return x + y*nx + t*(nx*ny);
            }

            int get_nbr_index(char delta_x, char delta_y, char delta_t) {
                unsigned char nbr_x = (nx + delta_x + x) % nx;
                unsigned char nbr_y = (ny + delta_y + y) % ny;
                unsigned char nbr_t = (nt + delta_t + t) % nt;

                return nbr_x + nbr_y*nx + nbr_t*(nx*ny);
            }

            unsigned char x, y, t;
        };

    private:

        // return true with probability p
        bool prob_true(double p) {
            // constexpr auto max = decltype(rand_gen)::max();
            constexpr unsigned int max = 0xffffffff;
            return (rand_uint() < p * max);
        }

    public:
        generic(double beta_p, double field_p) {
            set_params(beta_p, field_p);
            generate_random_state();
        }


        void generate_random_state() {
            for (auto &spin : state) {
                spin = 1-2*(rand_uint()%2);
            }
        }

        void set_params(double beta_p, double field_p) {
            beta = beta_p;
            field = field_p;

            double delta = beta / nt;

            eta = -.5 * std::log(std::tanh(delta * field));
            double Gamma = eta / delta;

            p_xy = 1. - std::exp(-2. * delta);
            p_t = 1. - std::exp(-2. * delta * Gamma);
        }


        int cluster_update() {
            std::vector<unsigned int> cluster, old_set, new_set;

            // index of the first site on the cluster
            unsigned int seed_index = rand_uint() % n_sites;

            cluster.push_back(seed_index);
            old_set.push_back(seed_index);

            // using this array, it takes O(1) to check if a given site 
            // belongs to the cluster; no need to do a search every time.
            std::array<bool, n_sites> on_cluster;
            for (auto &it : on_cluster) {
                it = false;
            }
            on_cluster[seed_index] = true;

            while (!old_set.empty()) {
                new_set.clear();
                for (auto my_index : old_set) {
                    coords my_coords{my_index};
                    char my_spin_state = state[my_index];

                    auto try_add = [&](int nbr_index, double p, char sign) {
                        if (prob_true(p) == true) {
                        if (on_cluster[nbr_index] == false) {
                        if (my_spin_state * state[nbr_index] == sign) {
                            new_set.push_back(nbr_index);
                            cluster.push_back(nbr_index);
                            on_cluster[nbr_index] = true;
                        }
                        }
                        }
                    };

                    // temporal-bond: common to all lattice structures
                    try_add(my_coords.get_nbr_index(0,0,+1), p_t, +1);
                    try_add(my_coords.get_nbr_index(0,0,-1), p_t, +1);

                    // lattice structure-dependent bonds
                    if constexpr (structure == lattice_structure::hexagonal) {
                        if (my_coords.y%2 == 0) {
                            try_add(my_coords.get_nbr_index(0,+1,0), p_xy, -1);
                            try_add(my_coords.get_nbr_index(0,-1,0), p_xy, +1);
                        } else {
                            try_add(my_coords.get_nbr_index(0,+1,0), p_xy, +1);
                            try_add(my_coords.get_nbr_index(0,-1,0), p_xy, -1);
                        }

                        char y_rem = my_coords.y % 4;
                        switch (y_rem) {
                        case 0:
                            try_add(my_coords.get_nbr_index(+1,-1,0), p_xy, +1);
                            break;
                        case 1:
                            try_add(my_coords.get_nbr_index(+1,+1,0), p_xy, +1);
                            break;
                        case 2:
                            try_add(my_coords.get_nbr_index(-1,-1,0), p_xy, +1);
                            break;
                        case 3:
                            try_add(my_coords.get_nbr_index(-1,+1,0), p_xy, +1);
                            break;
                        }
                    } else if constexpr
                        (structure == lattice_structure::square) {
                        try_add(my_coords.get_nbr_index(+1,0,0), p_xy, +1);
                        try_add(my_coords.get_nbr_index(-1,0,0), p_xy, +1);

                        if (my_coords.y%2 == 0) {
                            try_add(my_coords.get_nbr_index(0,+1,0), p_xy, -1);
                            try_add(my_coords.get_nbr_index(0,-1,0), p_xy, -1);
                        } else {
                            try_add(my_coords.get_nbr_index(0,+1,0), p_xy, +1);
                            try_add(my_coords.get_nbr_index(0,-1,0), p_xy, +1);
                        }
                    }

                }
                old_set = new_set;
            }

            for (auto index : cluster) {
                state[index] *= -1;
            }

            return cluster.size();
        }

        double energy_density() {
            double V = 0., U = 0.;

            for (unsigned int t = 0; t < nt; ++t) {
            for (unsigned int x = 0; x < nx; ++x) {
            for (unsigned int y = 0; y < ny; ++y) {

                coords my_coords(x, y, t);
                char my_spin = state[my_coords.get_index()];

                if constexpr (structure == lattice_structure::hexagonal) {
                    char nbr1 = state[my_coords.get_nbr_index(0,+1,0)];
                    if (y%2 == 0) {
                        nbr1 *= -1;
                    }
                    char nbr2;
                    switch (y%4) {
                        case 0:
                            nbr2 = state[my_coords.get_nbr_index(+1,-1,0)];
                            break;
                        case 1:
                            nbr2 = state[my_coords.get_nbr_index(+1,+1,0)];
                            break;
                        default:
                            nbr2 = 0;
                    }

                    U += (-1) * my_spin * (nbr1 + nbr2);
                } else if constexpr (structure == lattice_structure::square) {
                    char nbr1 = state[my_coords.get_nbr_index(0,+1,0)];
                    char nbr2 = state[my_coords.get_nbr_index(+1,0,0)];

                    if (x%2 == 0) {
                        nbr1 *= -1;
                    }

                    U += (-1) * my_spin * (nbr1 + nbr2);
                }

                int t_nbr_spin = state[my_coords.get_nbr_index(0,0,+1)];
                V += std::exp((-2.) * eta * my_spin * t_nbr_spin);
            }
            }
            }

            V = V * (-field) / n_sites;
            U = U / n_sites;

            return (U + V);
        }

        double magnetization_density() {
            int sum = 0;
            for (auto spin : state) {
                sum += spin;
            }
            return double(sum) / n_sites;
        }

    };
}

#endif
