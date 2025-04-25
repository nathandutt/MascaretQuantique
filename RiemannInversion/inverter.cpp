#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

double nan2 = std::numeric_limits<double>::quiet_NaN(); 

struct Config{
    double dx = 0.1;
    double dt = 0.1;
    double x_min = -10.;
    double x_max = 10.;
    double t_max = 20;
    double lambda_bound = 5;
    double dl = 0.1;
}

//Just to get closest idx to given cooridnates
auto findClosestX(const Config& config, double x){
    return static_cast<int>((x-config.x_min)/config.dx);
}


auto findClosestT(const Config& config, double t){
    return static_cast<int>(t/config.dt);
}


//Here you put the analytically obtained expressions for x, t
auto x(double l_plus, double l_minus){
    return 0.0;
}

auto t(double l_plus, double l_minus){
    return 0.0;
}


//Lambda class to do the inverting
class Lambdas{
    vector<vector<double>> lambda_plus;
    vector<vector<double>> lambda_minus;
    Config config;
    
    Lambdas(Config config_){
        config = config_
        auto time_size = static_cast<int>(config_.t_max/config.dt);
        auto x_size = static_cast<int>((config.x_max-config.x_min)/config.dx);

        lambda_plus = vector<vector<double>>(time_size, vector<double>(x_size, nan2));
        lambda_minus = vector<vector<double>>(time_size, vector<double>(x_size, nan2));

    }

    void Invert(){  
        auto max_idx = static_cast<int>(2.*config.lambda_bound/config.dl);
        for(int i_plus = 0; i_plus < max_idx; i_plus++){
            for(int i_minus = 0; i_minus < max_idx; i_minus++){
                auto l_plus = -l_bound + static_cast<double>(i_plus)*dl;
                auto l_minus = -l_bound + static_cast<double>(i_minus)*dl;
                auto x = x(l_plus, l_minus);
                auto t = t(l_plus, l_minus);
                if(x>config.x_max || x< config.x_min || t < 0. || t> config.t_max){
                    continue;
                }
                auto x_idx = findClosestX(config, x);
                auto t_idx = findClosestT(config, t);
                lambda_plus[t_idx][x_idx] = l_plus;
                lambda_minus[t_idx][x_idx] = l_minus;
            }
        }
    }


   void Write(const std::string& filename) const {
        WriteMatrixToCSV(lambda_plus, filename + "_plus.csv");
        WriteMatrixToCSV(lambda_minus, filename + "_minus.csv");
    }

private:
    void WriteMatrixToCSV(const std::vector<std::vector<double>>& matrix, const std::string& out_filename) const {
        std::ofstream out(out_filename);
        if (!out) {
            throw std::runtime_error("Failed to open file: " + out_filename);
        }

        out << std::fixed << std::setprecision(8); // consistent decimal format

        for (size_t t_idx = 0; t_idx < matrix.size(); ++t_idx) {
            double t = t_idx * config.dt;
            out << t; // first column = time

            const auto& row = matrix[t_idx];
            for (double val : row) {
                out << ",";
                if (std::isnan(val)) {
                    out << "NaN";
                } else {
                    out << val;
                }
            }
            out << "\n";
        }
    }
};

