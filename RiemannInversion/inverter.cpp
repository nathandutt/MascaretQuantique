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
    double dt = 0.05;
    double x_min = -70.;
    double x_max = 70.;
    double t_max = 100.;
    double lambda_bound = 100;
    double dl = 0.01;
    
    Config() = default;
    Config(const string& filename){
    }

};
int inserted = 0;
//Just to get closest idx to given cooridnates
auto findClosestX(const Config& config, double x){
    return static_cast<int>((x-config.x_min)/config.dx);
}


auto findClosestT(const Config& config, double t){
    auto t_t = static_cast<int>(t/config.dt);
    return t_t;
}

//Physical Functions
//
auto V_plus(double l_plus, double l_minus){
    //return 1./3. * l_plus + 2./3. * l_minus;
    return l_plus;
}

auto V_minus(double l_plus, double l_minus){
    //return 1./3. * l_minus + 2./3. * l_plus;
    return l_minus;
}

const auto l_0 = sqrt(6.0);
const auto a_0 = sqrt(6.0)/10.;
auto F(double l_plus){
    return 1./(2.*a_0)*(l_plus+l_0)*pow(l_plus-l_0, 2); 
 }
auto G(double l_minus){
    return  1./(2.*a_0)* (l_0 - l_minus)*pow(l_minus + l_0, 2);
}

auto F_prime(double l_plus){
    return  1./(2.*a_0) *(l_plus+l_0)*(l_plus-l_0)*2. + pow(l_plus - l_0, 2);
}

auto G_prime(double l_minus){
    return 2.* 1./(2.*a_0) *(l_minus + l_0)*(l_0 - l_minus) - pow(l_minus + l_0, 2);
}

auto chi(double l_plus, double l_minus){
    return 1./(l_plus - l_minus) * (F(l_plus)+G(l_minus));
}

auto chi_1(double l_plus, double l_minus){
    return (l_plus - sqrt(2.))/(2./10. * sqrt(2.));
    //return 1./(l_plus - l_minus) * (F_prime(l_plus) - chi(l_plus, l_minus)); 
}
auto chi_2(double l_plus, double l_minus){
    return (l_minus+sqrt(2.))/(2/10. * sqrt(2.));
    //return 1./(l_plus - l_minus) * (G_prime(l_plus) + chi(l_plus, l_minus)); 
}

//Final expressions!
auto x(double l_plus, double l_minus){
    return (V_minus(l_plus, l_minus)*chi_1(l_plus, l_minus) - V_plus(l_plus, l_minus)*chi_2(l_plus, l_minus))/(V_minus(l_plus, l_minus) - V_plus(l_plus, l_minus));
}

auto t(double l_plus, double l_minus){
    return (chi_1(l_plus, l_minus) - chi_2(l_plus, l_minus))/(V_minus(l_plus, l_minus) - V_plus(l_plus, l_minus));
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
        //cout << max_idx << "\n";
        auto x_size = static_cast<int>((config.x_max-config.x_min)/config.dx);
        auto time_size = static_cast<int>(config.t_max/config.dt);


        for(int i_plus = 0; i_plus < max_idx; i_plus++){
            for(int i_minus = 0; i_minus < i_plus; i_minus++){
                auto l_plus = -config.lambda_bound + static_cast<double>(i_plus)*config.dl;
                auto l_minus = -config.lambda_bound + static_cast<double>(i_minus)*config.dl;
                                auto x_val = x(l_plus, l_minus);
                auto t_val = t(l_plus, l_minus);
                if(x_val>=config.x_max || x_val<= config.x_min || t_val <= 0. || t_val>= config.t_max){
                    continue;
                }
                cout << l_plus << " " << l_minus << "\n"; 
                auto x_idx = findClosestX(config, x_val);
                auto t_idx = findClosestT(config, t_val);
                cout << "x " << x_idx << " t " << t_idx << " t max " << time_size << " x max " << x_size << "\n";
                if(x_idx >= x_size || t_idx >= time_size){std::cout << "We gone overfloooow \n";}
                lambda_plus[t_idx][x_idx] = l_plus;
                lambda_minus[t_idx][x_idx] = l_minus;
                inserted++;
            }
        }
    }


   void Write(const std::string& filename) const {
        WriteMatrixToCSV(lambda_plus, filename + "_plus.csv");
        WriteMatrixToCSV(lambda_minus, filename + "_minus.csv");
    }
double write_every = 0.25;
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
    cout << inserted << " insertions \n";
}
