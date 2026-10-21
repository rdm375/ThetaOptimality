// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Richard Myers
#include "peng_robinson_model.hpp"
#include <array>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

// Well-mixed vessel with a varying ethane feed. Two independent species
// fractions are the state; one feed-composition schedule is the control.
namespace tvg {
using namespace theta_opt;
inline constexpr double rate=.7, n_feed=.05, sigma_rel=.001;
inline constexpr std::array<double,3> knots{0.,5.,10.};
inline constexpr Vec<5> truth{.08,.04,.06,.18,.08}; // e0,n0,u(0),u(5),u(10)
using Schedule=PiecewiseLinearControl<1,3>;
inline double feed(double t,const Vec<5>& z){return Schedule(knots).value(t,Vec<3>{z[2],z[3],z[4]})[0];}
inline Vec<2> exact(double t,const Vec<5>& z){
    Vec<2> x{z[0],z[1]}; double a=0.;
    while(a<t-1e-13){
        double b=t;for(double k:knots)if(k>a+1e-13)b=std::min(b,k);
        double h=b-a,u0=feed(a,z),u1=feed(b,z),s=(u1-u0)/h;
        x[0]=u0+s*(h-1./rate)+(x[0]-u0+s/rate)*std::exp(-rate*h);
        x[1]=n_feed+(x[1]-n_feed)*std::exp(-rate*h);a=b;
    } return x;
}
struct Obs{double t;thermogpu::State thermo;double rho,sigma;};
inline thermogpu::State conditions(std::size_t j){
    constexpr std::array<double,4> p{1e6,5e6,1e7,2e7},t{300.,300.,350.,400.};
    return {p[j%4],t[j%4]};
}
inline std::vector<Obs> data(){
    std::vector<Obs> out;
    for(std::size_t j=0;j<=20;++j){double t=.5*j;auto c=conditions(j);
        double y=gas_inverse::density(c,exact(t,truth)).value;out.push_back({t,c,y,sigma_rel*y});}
    return out;
}
inline Problem<2,5> problem(const std::vector<Obs>& obs,std::size_t N){
    Problem<2,5> p;
    p.initial.value=[](const Vec<5>& z){return Vec<2>{z[0],z[1]};};
    p.initial.jacobian=[](const Vec<5>&){Mat<2,5> J{};J[0][0]=J[1][1]=1.;return J;};
    p.initial.directional_jacobian=[](const Vec<5>&,const Vec<5>&){return Mat<2,5>{};};
    Dynamics<2,1> f;
    f.value=[](double,const Vec<2>& x,const Vec<1>& u){return Vec<2>{rate*(u[0]-x[0]),rate*(n_feed-x[1])};};
    f.derivatives=[](double,const Vec<2>&,const Vec<1>&){Derivatives<2,1> d;d.x[0][0]=d.x[1][1]=-rate;d.control[0][0]=rate;return d;};
    f.directional_derivatives=[](double,const Vec<2>&,const Vec<1>&,const Vec<2>&,const Vec<1>&){return Derivatives<2,1>{};};
    auto base=Schedule(knots).parameterization();ControlParameterization<1,5> u;
    u.value=[base](double t,const Vec<5>& z){return base.value(t,Vec<3>{z[2],z[3],z[4]});};
    u.jacobian=[base](double t,const Vec<5>& z){Mat<1,5> J{};auto B=base.jacobian(t,Vec<3>{z[2],z[3],z[4]});for(int j=0;j<3;++j)J[0][j+2]=B[0][j];return J;};
    u.directional_jacobian=[](double,const Vec<5>&,const Vec<5>&){return Mat<1,5>{};};
    p.dynamics=compose_dynamics(std::move(f),std::move(u));
    for(const auto& o:obs){auto n=static_cast<std::size_t>(std::llround(o.t*N/10.));auto g=ScalarFunction<2,5>::zero();
        g.evaluate=[o](double,const Vec<2>& x,const Vec<5>&){auto q=gas_inverse::density(o.thermo,x);double r=q.value-o.rho;ScalarValue<2,5> v;v.value=.5*r*r/(o.sigma*o.sigma);v.x=(r/(o.sigma*o.sigma))*q.gradient;return v;};
        g.directional_gradient=[o](double,const Vec<2>& x,const Vec<5>&,const Vec<2>& dx,const Vec<5>&){auto q=gas_inverse::density(o.thermo,x);double r=q.value-o.rho;ScalarDirection<2,5> v;v.x=(1/(o.sigma*o.sigma))*(dot(q.gradient,dx)*q.gradient+r*mv(q.hessian,dx));return v;};
        p.nodal.emplace_back(n,std::move(g));}
    return p;
}
inline System<2,5> system(std::size_t N,double theta,const std::vector<Obs>& o){
    return System<2,5>(problem(o,N),Grid::uniform(0.,10.,N,theta),NewtonOptions{2e-13,0.,10,40});
}
inline double trajectory_error(const System<2,5>& s,const Vec<5>& z){
    auto y=s.forward(z);double e=0.;for(std::size_t n=0;n<y.state.size();++n)e=std::max(e,norm(y.state[n]-exact(s.grid().time[n],z)));return e;
}
inline PsdSpectrum<5> information(const System<2,5>& s,const Vec<5>& z,
                                  const Trajectory<2>& y,const std::vector<Obs>& obs){
    std::array<std::vector<Vec<2>>,5> tangent;
    for(std::size_t j=0;j<5;++j){Vec<5>d{};d[j]=1.;tangent[j]=s.tangent(z,y,d);}
    Mat<5,5> I{};
    for(const auto& o:obs){auto n=static_cast<std::size_t>(std::llround(o.t*s.grid().steps()/10.));
        auto q=gas_inverse::density(o.thermo,y.state[n]);
        for(int i=0;i<5;++i)for(int j=0;j<5;++j)
            I[i][j]+=dot(q.gradient,tangent[i][n])*dot(q.gradient,tangent[j][n])/(o.sigma*o.sigma);}
    return psd_spectrum(I);
}
inline double grad_error(const System<2,5>& s,const Vec<5>& z){
    Vec<5> probe=z;probe[0]+=.003;probe[1]-=.002;probe[2]+=.004;probe[3]-=.003;probe[4]+=.002;
    auto g=s.evaluate(probe).gradient;double e=0.;
    for(int j=0;j<5;++j){Vec<5>d{};d[j]=1.;double h=1e-4;
        auto c=[&](Vec<5> q){auto y=s.forward(q);return s.cost(q,y);};
        double fd=(c(probe+h*d)-c(probe-h*d))/(2*h);e=std::max(e,std::abs(fd-g[j])/std::max(1.,std::abs(fd)));}
    return e;
}
}
int main(int argc,char** argv){try{
    using namespace theta_opt;std::filesystem::path out=argc>1?argv[1]:"results/time_varying_gas";std::filesystem::create_directories(out);
    auto obs=tvg::data();std::ofstream conv(out/"convergence.csv");conv<<"theta,N,state_error,parameter_error,adjoint_gradient_error\n";
    std::ofstream pars(out/"parameter_refinement.csv");pars<<"theta,N,name,truth,estimate,error\n";
    std::array<std::ofstream,2> state;
    state[0].open(out/"state_theta_050.csv");state[1].open(out/"state_theta_065.csv");
    for(auto& f:state)f<<"time,ethane_true,nitrogen_true,methane_true,ethane_fit,nitrogen_fit,methane_fit,feed_ethane_true,feed_ethane_fit\n";
    std::ofstream residual(out/"density_residuals.csv");residual<<"theta,time,pressure_Pa,temperature_K,observed,predicted,residual,sigma\n";
    std::ofstream spectra(out/"information_spectrum.csv");spectra<<"theta,rank,condition,eigenvalue_0,eigenvalue_1,eigenvalue_2,eigenvalue_3,eigenvalue_4\n";
    constexpr std::array<std::size_t,4> Ns{80,160,320,640};
    constexpr std::array<const char*,5> names{"initial_ethane","initial_nitrogen","feed_ethane_0","feed_ethane_5","feed_ethane_10"};
    std::cout<<"DYNAMIC GAS MIXING-VESSEL STATE FINDING\nOne time-varying feed-composition schedule; Peng-Robinson density observations.\n";
    for(double theta:{.5,.65}){std::array<double,4> err{};Vec<5> finest{};
        for(std::size_t k=0;k<Ns.size();++k){auto s=tvg::system(Ns[k],theta,obs);
            auto fit=minimize_box_constrained(s,{.09,.035,.08,.15,.10},{.01,.01,.01,.01,.01},{.25,.15,.25,.25,.25},OptimizeOptions{1e-8,40,60});
            err[k]=tvg::trajectory_error(s,tvg::truth);if(k==3)finest=fit.control;
            conv<<theta<<','<<Ns[k]<<','<<std::setprecision(12)<<err[k]<<','<<norm(fit.control-tvg::truth)<<','<<tvg::grad_error(s,fit.control)<<'\n';
            for(int j=0;j<5;++j)pars<<theta<<','<<Ns[k]<<','<<names[j]<<','<<tvg::truth[j]<<','<<fit.control[j]<<','<<fit.control[j]-tvg::truth[j]<<'\n';
            if(k==3){auto path=s.forward(fit.control);
                auto info=tvg::information(s,fit.control,path,obs);
                spectra<<theta<<','<<info.rank<<','<<(info.rank==5?info.eigenvalues[0]/info.eigenvalues[4]:std::numeric_limits<double>::infinity());
                for(double e:info.eigenvalues)spectra<<','<<e;
                spectra<<'\n';
                auto& states=state[theta==.5?0:1];
                for(std::size_t n=0;n<path.state.size();++n){double t=s.grid().time[n];auto x=tvg::exact(t,tvg::truth),y=path.state[n];
                    states<<t<<','<<x[0]<<','<<x[1]<<','<<1-x[0]-x[1]<<','<<y[0]<<','<<y[1]<<','<<1-y[0]-y[1]<<','<<tvg::feed(t,tvg::truth)<<','<<tvg::feed(t,fit.control)<<'\n';}
                for(const auto& o:obs){auto n=static_cast<std::size_t>(std::llround(o.t*Ns[k]/10.));double pred=gas_inverse::density(o.thermo,path.state[n]).value;
                    residual<<theta<<','<<o.t<<','<<o.thermo.pressure_Pa<<','<<o.thermo.temperature_K<<','<<o.rho<<','<<pred<<','<<pred-o.rho<<','<<o.sigma<<'\n';}
                std::cout<<"  theta="<<theta<<" fit (e0,n0,u0,u5,u10)=(";
                for(int j=0;j<5;++j)std::cout<<(j?", ":"")<<fit.control[j];
                std::cout<<"), information rank="<<info.rank<<"/5, condition="
                         <<(info.rank==5?info.eigenvalues[0]/info.eigenvalues[4]:std::numeric_limits<double>::infinity())
                         <<", gradient check="<<tvg::grad_error(s,fit.control)<<'\n';
                if(info.rank!=5)throw std::runtime_error("The transient density samples do not locally identify the five fitted coefficients");
                if(tvg::grad_error(s,fit.control)>3e-5)throw std::runtime_error("Discrete adjoint gradient check failed");}
        }
        double order=std::log(err[1]/err[2])/std::log(2.);std::cout<<"  theta="<<theta<<" exact-state convergence order="<<order<<'\n';
        if(theta==.5&&std::abs(order-2.)>.2)throw std::runtime_error("Expected second-order midpoint convergence");
        if(theta==.65&&std::abs(order-1.)>.25)throw std::runtime_error("Expected first-order off-centered convergence");
    }
    std::cout<<"Results: "<<out<<'\n';return 0;
}catch(const std::exception& e){std::cerr<<"time_varying_gas: "<<e.what()<<'\n';return 1;}}
