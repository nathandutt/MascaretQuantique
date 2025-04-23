
#include <fstream>
#include <complex>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>
#include <cmath>

#include <iostream>

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
    double save_every;

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
                else if (key == "save_every") save_every = value;

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

auto SlabDensity(const Config& config, double x){
    auto alpha = 1./(config.smoothing*config.dx); 
    return config.exterior_density + (config.amplitude - config.exterior_density) * 0.5 * ( std::tanh(alpha * (x + config.L))
                 - std::tanh(alpha * (x - config.L)) );
}

auto StepDensity(const Config& config, double x){
    return 0.5 + 0.5*( std::tanh((1./(config.smoothing*config.dx)) * (x + config.L)));}
                 
auto toComplex(double density, double phase) {
   // return  pow(density, 0.5) * exp(I * phase);i
   return  complex(pow(density, 0.5));

}

auto InitialData(const Config& config) {
    auto id = vector<complex<double>>{};
    for (int i = 0; i < config.domain_size; i++) {
        auto x = IdxToCoord(config, i);
        auto density = SlabDensity(config, x);
        auto phase = 0.;
        id.emplace_back(toComplex(density, phase));
    }
    return id;
}

auto DoubleDeriv(const Config& config, const vector<complex<double>>& psi) {
    auto gradient = vector<complex<double>>(config.domain_size);

    for (int i = 0; i < config.domain_size; i++) {
        if(i==0 || i == config.domain_size -1){
            gradient.emplace_back(0.);
            continue;
        }
        gradient[i] = (psi[i-1] - complex<double>(2.0) * psi[i] + psi[i+1]) 
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

auto GetCurrent(const Config& config,
                const std::vector<std::complex<double>>& psi)
{
    const int N   = config.domain_size;
    const double dx = config.dx;

    std::vector<double> j(N);

    // Bulk points: central difference for ∂ψ/∂x
    for(int i = 1; i < N-1; ++i) {
        std::complex<double> grad = (psi[i+1] - psi[i-1]) / (2.0 * dx);
        j[i] = std::imag(std::conj(psi[i]) * grad);
    }

    // Boundaries: one‐sided difference
    {
        auto grad0 = (psi[1]   - psi[0])   / dx;
        j[0]      = std::imag(std::conj(psi[0]) * grad0);

        auto gradN = (psi[N-1] - psi[N-2]) / dx;
        j[N-1]    = std::imag(std::conj(psi[N-1]) * gradN);
    }

    return j;
}

// Psi class for evolving the state using Runge-Kutta
class Psi {
    double time;
    vector<complex<double>> data;
    Config config;
    vector<double> phase_gradient;

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

    void GetPhaseGradient(){
       auto current = GetCurrent(config, data);
       auto u = vector<double>{};
       for(int i = 0; i<config.domain_size; i++){
           if(abs(data[i]) < pow(10.,-3)){
                   u.emplace_back(0.);
                   continue;

            }
            
            auto velocity = current[i]/norm(data[i]);
            u.emplace_back(velocity);
       }
       phase_gradient = u;
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
        out << time;
        for (const auto& c : phase_gradient){
            out << "," << c;
        }
        out << "\n";
    }
};

// Main function to run the simulation
int main() {
    string filename = "evolution.csv";

    // Load configuration from the file
    Config config("config.txt");
    auto step = static_cast<int>(config.save_every/config.dt);
    // Calculate timesteps
    auto timesteps = static_cast<int>(config.t_max / config.dt);
    auto state = Psi(config);
    // Evolve system and save data to file
    for (int t = 0; t < timesteps; t++) {
        if(t%step == 0){
	cout << "t = " << config.dt*t << "\n";
    state.GetPhaseGradient();
	state.Write(filename);
	}
	state.EvolveRK();
    }

    return 0;
}

