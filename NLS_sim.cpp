#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>
#include <cmath>

using namespace std;

const complex<double> I(0.0, 1.0);
// Class to store simulation and initial data parameters
class Config {
    // Simulation Parameters
public:
    double x_range;
    double dx;
    double dt;
    double smoothing;
    int domain_size;
    double t_max;
    double g;
    double exponent;

    // Initial Data Parameters
    double L;
    double amplitude;
    double exterior_density;
    double c;

    Config(){
    }

    Config(std::string filename) {
        std::ifstream in(filename);
        if (!in) {
            throw std::runtime_error("Could not open config file: " + filename);
        }

        std::string line;
        while (std::getline(in, line)) {
            // Remove leading and trailing whitespaces
            std::stringstream ss(line);
            std::string key;
            double value;

            // Skip comments or empty lines
            if (line.empty() || line[0] == '#') {
                continue;
            }

            // Read the key-value pair
            if (ss >> key >> value) {
                if (key == "x_range") x_range = value;
                else if (key == "dx") dx = value;
                else if (key == "dt") dt = value;
                else if (key == "domain_size") domain_size = static_cast<int>(value);
                else if (key == "t_max") t_max = value;
                else if (key == "g") g = value;
                else if (key == "exponent") exponent = value;
                else if (key == "smoothing") smoothing = value;
                else if (key == "L") L = value;
                else if (key == "amplitude") amplitude = value;
                else if (key == "exterior_density") exterior_density = value;
                else if (key == "c") c = value;
                else {
                    std::cerr << "Warning: Unknown config key '" << key << "'\n";
                }
            }
        }

        // Derived values
        domain_size = static_cast<int>(2 * x_range / dx);
    }
};

// Helper functions
auto IdxToCoord(const Config& config, int idx) {
    return -config.x_range + idx * config.dx;
}

auto SolitonId(const Config& config, double x) {
    return pow(1. / cosh(x), 2);
}

auto StepDensity(const Config& config, double x) {
    double absx = std::abs(x);
    double eps = config.dx * config.smoothing;
    if (absx >= config.L + eps) {
        return config.exterior_density;
    } else if (absx <= config.L - eps) {
        return config.amplitude + config.exterior_density;
    } else {
        // Corrected z to map linearly from L+eps to L-eps
        double z = (config.L + eps - absx) / (2.0 * eps);  // z goes from 0 at L+eps to 1 at L-eps
        // Smooth tanh transition between L-eps and L+eps
        return config.exterior_density + (config.amplitude) * 0.5 * (1 + tanh(10 * (z - 0.5)));  // Smooth transition
    }
}

auto StepPhase(const Config& config, double x) {
    double absx = std::abs(x);
    double eps = config.smoothing * config.dx;

    if (absx <= config.L - eps) {
        return 0.0;
    } else if (absx >= config.L + eps) {
        return config.c * (absx - (config.L + eps));
    } else {
        // Arctangent smoothing in transition region
        double z = (absx - (config.L - eps)) / (2 * eps);  // z ∈ [0, 1]
        double smooth = (1.0 / M_PI) * atan(10 * (z - 0.5)) + 0.5;  // Smoothstep from 0 to 1
        return config.c * smooth * (absx - (config.L + eps));
    }
}

auto toComplex(double density, double phase) {
    return  pow(density, 0.5) * exp(I * phase);
}

auto InitialData(const Config& config) {
    auto id = vector<complex<double>>{};
    for (int i = 0; i < config.domain_size; i++) {
        auto x = IdxToCoord(config, i);
        auto density = StepDensity(config, x);
        auto phase = StepPhase(config, x);
        id.emplace_back(toComplex(density, phase));
    }
    return id;
}

auto DoubleDeriv(const Config& config, const vector<complex<double>>& psi) {
    auto gradient = vector<complex<double>>(config.domain_size);

    for (int i = 0; i < config.domain_size; i++) {
        int left  = (i - 1 + config.domain_size) % config.domain_size;
        int right = (i + 1) % config.domain_size;

        gradient[i] = (psi[left] - complex<double>(2.0) * psi[i] + psi[right]) 
                      / (config.dx * config.dx);
    }

    return gradient;
}

auto linear_comb(const vector<complex<double>>& v_1, const vector<complex<double>>& v_2, complex<double> coefficient) {
    auto sum = vector<complex<double>>{};
    for (int i = 0; i < v_1.size(); i++) {
        auto added = v_1[i] + coefficient * v_2[i];
        sum.emplace_back(added);
    }
    return sum;
}

auto F(const Config& config, const vector<complex<double>>& psi) {
    auto result = vector<complex<double>>{};
    auto gradient = DoubleDeriv(config, psi);
    for (int i = 0; i < config.domain_size; i++) {
        auto absval = std::norm(psi[i]);
        auto term_1 = 0.5 * gradient[i] - config.g * psi[i] * std::pow(absval, config.exponent - 1.);
        auto final_t = I * term_1;
        result.emplace_back(final_t);
    }
    return result;
}

// Psi class for evolving the state using Runge-Kutta
class Psi {
    double time;
    vector<complex<double>> data;
    Config config;

public:
    Psi(Config config_) {
        time = 0.;
        config = config_;
        data = InitialData(config_);
    }

    void EvolveRK() {
        auto k_1 = F(config, data);
        auto k_2 = F(config, linear_comb(data, k_1, config.dt / 2.));
        auto k_3 = F(config, linear_comb(data, k_2, config.dt / 2.));
        auto k_4 = F(config, linear_comb(data, k_3, config.dt));
        auto evolved = vector<complex<double>>{};

        for (int i = 0; i < config.domain_size; i++) {
            auto evolved_i = data[i] + config.dt / 6. * (k_1[i] + 2. * k_2[i] + 2. * k_3[i] + k_4[i]);
            evolved.emplace_back(evolved_i);
        }
        data = evolved;
        time = time + config.dt;
    }

    void Write(const std::string& filename) const {
        std::ofstream out(filename, std::ios::app);  // append mode
        if (!out) {
            throw std::runtime_error("Could not open file: " + filename);
        }

        out << time;
        for (const auto& c : data) {
            out << "," << c.real() << "," << c.imag();
        }
        out << "\n";
    }
};

// Main function to run the simulation
int main() {
    string filename = "evolution.csv";

    // Load configuration from the file
    Config config("config.txt");
    auto step = static_cast<int>(0.1/config.dt);
    // Calculate timesteps
    auto timesteps = static_cast<int>(config.t_max / config.dt);
    auto state = Psi(config);
    // Evolve system and save data to file
    for (int t = 0; t < timesteps; t++) {
        if(t%step == 0){
	cout << "t = " << config.dt*t << "\n";
	state.Write(filename);
	}
	state.EvolveRK();
    }

    return 0;
}

