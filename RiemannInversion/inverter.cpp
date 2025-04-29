#include <iomanip>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>

using namespace std;

double nan2 = std::numeric_limits<double>::quiet_NaN(); 

struct Config{
    double dx = 0.1;
    double dt = 0.1;
    double x_min = -250.;
    double x_max = 250.;
    double t_max = 60;
    double lambda_bound = -500;
    double dl = 0.01;
    
    Config() = default;
    Config(const string& filename){
    }

};

//Just to get closest idx to given cooridnates
auto findClosestX(const Config& config, double x){
    return static_cast<int>((x-config.x_min)/config.dx);
}


auto findClosestT(const Config& config, double t){
    auto t_t = static_cast<int>(t/config.dt);
    return t_t;
}


//Here you put the analytically obtained expressions for x, t
auto x(double l_plus, double l_minus){
    return l_plus - l_minus;
}

auto t(double l_plus, double l_minus){
    return l_plus + l_minus;
}


//Lambda class to do the inverting
class Lambdas{
    vector<vector<double>> lambda_plus;
    vector<vector<double>> lambda_minus;
    Config config;
public:    
    Lambdas(Config config_){
        config = config_;
        auto time_size = static_cast<int>(config_.t_max/config.dt);
        auto x_size = static_cast<int>((config.x_max-config.x_min)/config.dx);

        lambda_plus = vector<vector<double>>(time_size, vector<double>(x_size, nan2));
        lambda_minus = vector<vector<double>>(time_size, vector<double>(x_size, nan2));

    }

    void Invert(){  
        auto max_idx = static_cast<int>(2.*config.lambda_bound/config.dl);
        auto x_size = static_cast<int>((config.x_max-config.x_min)/config.dx);
        auto time_size = static_cast<int>(config.t_max/config.dt);


        for(int i_plus = 0; i_plus < max_idx; i_plus++){
            for(int i_minus = 0; i_minus < max_idx; i_minus++){
                auto l_plus = -config.lambda_bound + static_cast<double>(i_plus)*config.dl;
                auto l_minus = -config.lambda_bound + static_cast<double>(i_minus)*config.dl;

                auto x_val = x(l_plus, l_minus);
                auto t_val = t(l_plus, l_minus);
                if(x_val>=config.x_max || x_val<= config.x_min || t_val <= 0. || t_val>= config.t_max){
                    continue;
                }
                auto x_idx = findClosestX(config, x_val);
                auto t_idx = findClosestT(config, t_val);
                if(x_idx >= x_size || t_idx >= time_size){std::cout << "We gone overfloooow \n";}
                lambda_plus[t_idx][x_idx] = l_plus;
                lambda_minus[t_idx][x_idx] = l_minus;
            }
        }
    }


   void Write(const std::string& filename) const {
        WriteMatrixToCSV(lambda_plus, filename + "_plus.csv");
        WriteMatrixToCSV(lambda_minus, filename + "_minus.csv");
    }
double write_every = 2.;
int write_step = static_cast<int>(write_every/config.dt);
private:
    void WriteMatrixToCSV(const std::vector<std::vector<double>>& matrix, const std::string& out_filename) const {
        std::ofstream out(out_filename);
        if (!out) {
            throw std::runtime_error("Failed to open file: " + out_filename);
        }

        out << std::fixed << std::setprecision(8); // consistent decimal format

        for (size_t t_idx = 0; t_idx < matrix.size(); ++t_idx) {
            if(t_idx % write_step != 0){
                continue;
            }
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

int main(){
    auto config = Config();
    auto lambdas = Lambdas(config);
    std::cout << "iNverting" << "\n" ;
    lambdas.Invert();
    std::cout << "Done inverting" << "\n";
    string output = "inverted";
    lambdas.Write(output);
}
