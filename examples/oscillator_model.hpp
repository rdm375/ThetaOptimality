#pragma once
#include <theta_opt/theta_opt.hpp>
#include <complex>
#include <numbers>

namespace oscillator {
using namespace theta_opt;
inline constexpr double t0=2.,tf=10.,horizon=tf-t0;
inline constexpr double mass=10.,decay=.2,frequency=2.;
inline constexpr double damping=2*mass*decay,stiffness=mass*(decay*decay+frequency*frequency);
inline constexpr double sigma=.0005;
inline constexpr Vec<2> truth{.02,.12}; // metres, metres/second
inline constexpr Mat<2,2> A{{{0.,1.},{-stiffness/mass,-damping/mass}}};
inline constexpr std::array<double,6> elapsed{.5,1.,2.,3.,4.,6.};
inline constexpr std::array<double,6> perturbation{.4,-.6,.2,.8,-.3,.5};
inline Mat<2,2> propagator(double tau){
    double c=std::cos(frequency*tau),s=std::sin(frequency*tau),e=std::exp(-decay*tau);
    return {{{e*(c+decay*s/frequency),e*s/frequency},{-e*(decay*decay+frequency*frequency)*s/frequency,e*(c-decay*s/frequency)}}};
}
inline Vec<2> exact(double tau,const Vec<2>& xi){return mv(propagator(tau),xi);}
inline double observation(std::size_t j,bool noisy=true){return exact(elapsed[j],truth)[0]+(noisy?sigma*perturbation[j]:0.);}
inline Mat<2,2> gramian(double duration){
    // Exact integral of Phi(t)^T C^T C Phi(t), C=(1,0).
    double I0=-std::expm1(-2*decay*duration)/(2*decay);
    std::complex<double> z(-2*decay,2*frequency),integral=(std::exp(z*duration)-1.)/z;
    double cc=.5*(I0+integral.real()),ss=.5*(I0-integral.real()),cs=.5*integral.imag(),r=decay/frequency;
    return {{{cc+2*r*cs+r*r*ss,(cs+r*ss)/frequency},{(cs+r*ss)/frequency,ss/(frequency*frequency)}}};
}
inline double energy(const Vec<2>& x){return .5*(stiffness*x[0]*x[0]+mass*x[1]*x[1]);}
inline Mat<2,2> outer(const Vec<2>& v){return {{{v[0]*v[0],v[0]*v[1]},{v[1]*v[0],v[1]*v[1]}}};}
inline Mat<2,2> inverse(const Mat<2,2>& H){Mat<2,2> inv{};LU<2> lu(H);for(int j=0;j<2;++j){Vec<2> e{};e[j]=1.;auto col=lu.solve(e);for(int i=0;i<2;++i)inv[i][j]=col[i];}return inv;}
inline Mat<2,2> information(){Mat<2,2> H{};for(double tau:elapsed)H=add(H,scaled(1/(sigma*sigma),outer(propagator(tau)[0])));return H;}
inline Vec<2> analytical_fit(bool noisy=true){Vec<2> rhs{};for(std::size_t j=0;j<elapsed.size();++j)rhs=rhs+(observation(j,noisy)/(sigma*sigma))*propagator(elapsed[j])[0];return LU<2>(information()).solve(rhs);}
inline std::size_t node_at(const Grid& g,double tau){for(std::size_t j=0;j<g.time.size();++j)if(std::abs(g.time[j]-t0-tau)<1e-12)return j;throw std::invalid_argument("Oscillator observation must coincide with a node");}
inline Problem<2,2> problem(const Grid& g,bool running=false,bool noisy=true){
    Problem<2,2> p;
    p.initial.value=[](const Vec<2>& xi){return xi;};
    p.initial.jacobian=[](const Vec<2>&){return identity<2>();};
    p.initial.directional_jacobian=[](const Vec<2>&,const Vec<2>&){return Mat<2,2>{};};
    p.dynamics.value=[](double,const Vec<2>& x,const Vec<2>&){return mv(A,x);};
    p.dynamics.derivatives=[](double,const Vec<2>&,const Vec<2>&){Derivatives<2,2> d;d.x=A;return d;};
    p.dynamics.directional_derivatives=[](double,const Vec<2>&,const Vec<2>&,const Vec<2>&,const Vec<2>&){return Derivatives<2,2>{};};
    if(running){
        // A separate verification functional, not synthetic extra sensor data:
        // J_run=1/2 integral q(t)^2 dt. Used to test quadrature and costates.
        p.running.evaluate=[](double,const Vec<2>& x,const Vec<2>&){ScalarValue<2,2> v;v.value=.5*x[0]*x[0];v.x[0]=x[0];return v;};
        p.running.directional_gradient=[](double,const Vec<2>&,const Vec<2>&,const Vec<2>& s,const Vec<2>&){ScalarDirection<2,2> v;v.x[0]=s[0];return v;};
    }else for(std::size_t j=0;j<elapsed.size();++j){auto cost=ScalarFunction<2,2>::zero();double y=observation(j,noisy);
        cost.evaluate=[y](double,const Vec<2>& x,const Vec<2>&){ScalarValue<2,2> v;double r=x[0]-y;v.value=.5*r*r/(sigma*sigma);v.x[0]=r/(sigma*sigma);return v;};
        cost.directional_gradient=[](double,const Vec<2>&,const Vec<2>&,const Vec<2>& s,const Vec<2>&){ScalarDirection<2,2> v;v.x[0]=s[0]/(sigma*sigma);return v;};
        p.nodal.emplace_back(node_at(g,elapsed[j]),cost);
    }return p;
}
inline System<2,2> system(std::size_t steps,double theta,bool running=false,bool noisy=true){auto g=Grid::uniform(t0,tf,steps,theta);return System<2,2>(problem(g,running,noisy),g);}
inline Optimum<2> fit(const System<2,2>& s){return minimize(s,Vec<2>{0.,0.},OptimizeOptions{1e-8,5,40});}
inline double exact_cost(const Vec<2>& xi){double J=0.;for(std::size_t j=0;j<elapsed.size();++j){double r=exact(elapsed[j],xi)[0]-observation(j);J+=.5*r*r/(sigma*sigma);}return J;}
inline Vec<2> exact_gradient(const Vec<2>& xi){Vec<2> g{};for(std::size_t j=0;j<elapsed.size();++j)g=g+((exact(elapsed[j],xi)[0]-observation(j))/(sigma*sigma))*propagator(elapsed[j])[0];return g;}
inline Vec<2> nodal_costate(double tau,const Vec<2>& xi){Vec<2> lambda{};for(std::size_t j=0;j<elapsed.size();++j)if(elapsed[j]+1e-12>=tau)lambda=lambda+((exact(elapsed[j],xi)[0]-observation(j))/(sigma*sigma))*propagator(elapsed[j]-tau)[0];return lambda;}
inline Vec<2> running_costate(double tau,const Vec<2>& xi){return mv(gramian(horizon-tau),exact(tau,xi));}
inline double eigen_min(const Mat<2,2>& H){double big=.5*(H[0][0]+H[1][1]+std::hypot(H[0][0]-H[1][1],2*H[0][1]));return (H[0][0]*H[1][1]-H[0][1]*H[1][0])/big;}
}
