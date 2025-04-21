// Want to do Runge Kuta
#include <fstream>
#include <complex>
#include <vector>
#include <string>
using namespace std;

double x_range =  20;
double dx = 0.1;
double dt = 0.08 * dx*dx; //RK condition, maybe can be more precise
int domain_size = 2*x_range/dx;
double t_max = 5;
double g = 0.01; //g>1 is repulsive here
double exponent = 2.;
double epsilon = 10*dx;
auto L = 3.; //half length of initial profile step
double amplitude = 1.0;
double exterior = 0.3;
double c = 1.;
double a =1;

const complex<double> I(0.0, 1.0);

auto IdxToCoord(int idx){
	return -x_range + idx*dx;
}
auto SolitonId(double x, double L, double eps){
	return 1./cosh(x);
}
auto StepDensityC2(double x, double L, double eps){
    double absx = std::abs(x);

    if (absx >= L + eps) {
        return exterior;
    } else if (absx <= L - eps) {
        return amplitude;
    } else {
        double z = (L + eps - absx) / (2.0 * eps);  //map to [0,1]
        return exterior+amplitude*(6*z*z*z*z*z - 15*z*z*z*z + 10*z*z*z); // This function gives a C2 transition
    }
}

auto StepDensity(double x, double L, double eps){
    double absx = std::abs(x);

    if (absx >= L + eps) {
        return exterior;
    } else if (absx <= L - eps) {
        return amplitude + exterior;
    } else {
        // Corrected z to map linearly from L+eps to L-eps
        double z = (L + eps - absx) / (2.0 * eps);  // z goes from 0 at L+eps to 1 at L-eps
        // Smooth tanh transition between L-eps and L+eps
        return exterior + (amplitude) * 0.5 * (1 + tanh(10 * (z - 0.5)));  // Smooth transition
    }
}


auto StepPhase(double x, double L, double eps){
    double absx = std::abs(x);
    
    if (absx <= L - eps) {
        return 0.0;
    } else if (absx >= L + eps) {
        return c * (absx - (L + eps));
    } else {
        // In the smoothing zone: [L - eps, L + eps]
        // z goes from 0 to 1
        double z = (absx - (L - eps)) / (2 * eps);
        double smooth = z * z * (3 - 2 * z); // C¹ smoothstep (Hermite polynomial)
        return c * smooth * (absx - (L + eps));
    }
}
auto toComplex(double density, double phase){
	return pow(density, 0.5) * exp(I*phase);
}
auto InitialData(){
	auto id = vector<complex<double>>{};
	for(int i =0; i< domain_size; i++){
		auto x = IdxToCoord(i);
		auto density = StepDensity(x, L, epsilon);
		auto phase = StepPhase(x, L, epsilon);
		id.emplace_back(toComplex(density, phase));
	}
	return id;
}
auto DoubleDeriv(vector<complex<double>>& psi){
	auto gradient = vector<complex<double>>{};
	for(int i=0; i<domain_size; i++){
		if(i ==0||(i+1 == domain_size)){
			gradient.emplace_back(0);
		}
		else{
			auto term = (psi[i-1]- complex<double>(2) * psi[i]+ psi[i+1]) / pow(complex<double>(dx), 2);
			gradient.emplace_back(term);
		}
	}
	return gradient;
}

auto linear_comb(vector<complex<double>>& v_1, vector<complex<double>>& v_2, complex<double> coefficient){
	auto sum = vector<complex<double>>{};
	for(int i = 0; i < domain_size; i++){
		auto added = v_1[i] + coefficient * v_2[i];
		sum.emplace_back(added);
	}
	return sum;
}

auto F(vector<complex<double>> psi){
	auto result = vector<complex<double>>{};
	auto gradient = DoubleDeriv(psi);
	for(int i = 0; i < domain_size; i++){
		auto power = 2. * (exponent - double(1.));
		auto term_1 =  -a*gradient[i] + g*psi[i]*pow(abs(psi[i]),power);
		auto final_t = (1./I) * term_1;
		result.emplace_back(final_t);
	}
	return result;
}

class Psi{
	double time;
	vector<complex<double>> data;
	public: 
	
	Psi(){
		time = 0.;
		data = InitialData();
	}
	void EvolveRK(){
		auto k_1 = F(data);
		auto k_2 = F(linear_comb(data, k_1, dt/2.));
		auto k_3 = F(linear_comb(data, k_2, dt/2.));
		auto k_4 = F(linear_comb(data, k_3, dt));
		auto evolved = vector<complex<double>>{};

		for(int i = 0; i < domain_size; i++){
			auto evolved_i = data[i] + dt/6. * (k_1[i] + 2.*k_2[i]+ 2.*k_3[i]+k_4[i]);
			evolved.emplace_back(evolved_i);
		}
		data = evolved;
		time = time+dt;
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

int main(){
	string filename = "evolution.csv";
	auto timesteps = t_max/dt;
	auto state = Psi();

	for(int t = 0; t < timesteps; t++){
		state.Write(filename);
		state.EvolveRK();
	}
	return 0;
}