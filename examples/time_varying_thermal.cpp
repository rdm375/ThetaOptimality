// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Richard Myers

#include <theta_opt/theta_opt.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

// A two-node thermal system with one fluid thermometer and an unknown,
// piecewise-linear heater schedule. The five fitted coefficients are the two
// initial temperatures and heater power at three schedule knots.
namespace tvt {
using namespace theta_opt;
inline constexpr Vec<2> equilibrium{30.,20.};
inline constexpr Mat<2,2> A{{{-2.,1.},{.25,-.5}}};
inline constexpr Vec<2> heater_direction{.15,0.};
inline constexpr std::array<double,3> knot_time{0.,5.,10.};
inline constexpr Vec<5> truth{35.,15.,.2,1.2,.1};
inline constexpr double sigma=.08;
using Schedule=PiecewiseLinearControl<1,3>;

inline double control(double t,const Vec<5>& z){
    Vec<3> c{z[2],z[3],z[4]};return Schedule(knot_time).value(t, c)[0];
}
inline double i0(double rate,double h){return std::expm1(rate*h)/rate;}
inline double i1(double rate,double h){return (std::expm1(rate*h)-rate*h)/(rate*rate);}
inline Vec<2> advance_exact(const Vec<2>& x,double t,double h,const Vec<5>& z){
    const double mid=t+.5*h;
    const auto S=scaled(1./( (-5.+std::sqrt(13.))/4.-(-5.-std::sqrt(13.))/4.),
                        add(A,scaled(-((-5.-std::sqrt(13.))/4.),identity<2>())));
    const auto F=add(identity<2>(),scaled(-1.,S));
    const double u0=control(t,z),slope=(control(t+h,z)-u0)/h;
    const auto y=x-equilibrium;
    const double ls=(-5.+std::sqrt(13.))/4.,lf=(-5.-std::sqrt(13.))/4.;
    auto result=mv(S,y);result=std::exp(ls*h)*result;
    result=result+std::exp(lf*h)*mv(F,y);
    const double qs=u0*i0(ls,h)+slope*i1(ls,h),qf=u0*i0(lf,h)+slope*i1(lf,h);
    result=result+qs*mv(S,heater_direction)+qf*mv(F,heater_direction);
    (void)mid;return equilibrium+result;
}
inline Vec<2> exact(double t,const Vec<5>& z){
    Vec<2> x{z[0],z[1]};double now=0.;
    while(now<t-1e-14){double end=t;for(double knot:knot_time)if(knot>now+1e-14)end=std::min(end,knot);
        x=advance_exact(x,now,end-now,z);now=end;}
    return x;
}
inline Problem<2,5> problem(const std::vector<std::pair<double,double>>& observations,std::size_t steps){
    Problem<2,5> p;
    p.initial.value=[](const Vec<5>& z){return Vec<2>{z[0],z[1]};};
    p.initial.jacobian=[](const Vec<5>&){Mat<2,5> J{};J[0][0]=1.;J[1][1]=1.;return J;};
    p.initial.directional_jacobian=[](const Vec<5>&,const Vec<5>&){return Mat<2,5>{};};
    Dynamics<2,1> model;
    model.value=[](double,const Vec<2>& x,const Vec<1>& u){return mv(A,x-equilibrium)+u[0]*heater_direction;};
    model.derivatives=[](double,const Vec<2>&,const Vec<1>&){Derivatives<2,1> d;d.x=A;d.control[0][0]=heater_direction[0];d.control[1][0]=heater_direction[1];return d;};
    model.directional_derivatives=[](double,const Vec<2>&,const Vec<1>&,const Vec<2>&,const Vec<1>&){return Derivatives<2,1>{};};
    auto base=Schedule(knot_time).parameterization();
    ControlParameterization<1,5> scheduled;
    scheduled.value=[base](double t,const Vec<5>& z){return base.value(t,Vec<3>{z[2],z[3],z[4]});};
    scheduled.jacobian=[base](double t,const Vec<5>& z){Mat<1,5> J{};auto q=base.jacobian(t,Vec<3>{z[2],z[3],z[4]});for(std::size_t j=0;j<3;++j)J[0][j+2]=q[0][j];return J;};
    scheduled.directional_jacobian=[](double,const Vec<5>&,const Vec<5>&){return Mat<1,5>{};};
    p.dynamics=compose_dynamics(std::move(model),std::move(scheduled));
    for(const auto& [t,y]:observations){
        const auto node=static_cast<std::size_t>(std::llround(t*static_cast<double>(steps)/10.));
        auto f=ScalarFunction<2,5>::zero();
        f.evaluate=[y](double,const Vec<2>& x,const Vec<5>&){double r=x[0]-y;ScalarValue<2,5> v;v.value=.5*r*r/(sigma*sigma);v.x[0]=r/(sigma*sigma);return v;};
        f.directional_gradient=[](double,const Vec<2>&,const Vec<5>&,const Vec<2>& dx,const Vec<5>&){ScalarDirection<2,5> d;d.x[0]=dx[0]/(sigma*sigma);return d;};
        p.nodal.emplace_back(node,std::move(f));
    }
    return p;
}
inline System<2,5> make_system(std::size_t steps,double theta,const std::vector<std::pair<double,double>>& observations){
    return System<2,5>(problem(observations,steps),Grid::uniform(0.,10.,steps,theta),NewtonOptions{2e-13,0.,10,40});
}
inline std::vector<std::pair<double,double>> data(){
    std::vector<std::pair<double,double>> out;
    for(int k=0;k<=20;++k){double t=.5*k;out.emplace_back(t,exact(t,truth)[0]);}return out;
}
inline double order(double coarse,double fine){return std::log(coarse/fine)/std::log(2.);}
inline double trajectory_error(const System<2,5>& sys,const Vec<5>& z){
    auto tr=sys.forward(z);double e=0.;for(std::size_t n=0;n<tr.state.size();++n){auto x=exact(sys.grid().time[n],z);e=std::max(e,norm(tr.state[n]-x));}return e;
}
inline PsdSpectrum<5> information(const System<2,5>& sys,const Vec<5>& z,
                                  const Trajectory<2>& path,
                                  const std::vector<std::pair<double,double>>& obs){
    std::array<std::vector<Vec<2>>,5> sensitivity;
    for(std::size_t j=0;j<5;++j){Vec<5> direction{};direction[j]=1.;sensitivity[j]=sys.tangent(z,path,direction);}
    Mat<5,5> matrix{};
    for(const auto& [t,value]:obs){(void)value;const auto n=static_cast<std::size_t>(std::llround(t*static_cast<double>(sys.grid().steps())/10.));
        for(std::size_t i=0;i<5;++i)for(std::size_t j=0;j<5;++j)
            matrix[i][j]+=sensitivity[i][n][0]*sensitivity[j][n][0]/(sigma*sigma);}
    return psd_spectrum(matrix);
}
inline double gradient_fd(const System<2,5>& sys,const Vec<5>& z){
    auto g=sys.evaluate(z).gradient;double err=0.;for(std::size_t j=0;j<5;++j){auto e=Vec<5>{};e[j]=1.;double h=2e-5;
        auto cost=[&](const Vec<5>& q){auto path=sys.forward(q);return sys.cost(q,path);};
        double fd=(cost(z+h*e)-cost(z-h*e))/(2*h);err=std::max(err,std::abs(fd-g[j])/std::max(1.,std::abs(fd)));}return err;
}
}

int main(int argc,char** argv){
    try{
        using namespace theta_opt;
        const std::filesystem::path dir=argc>1?argv[1]:"results/time_varying_thermal";std::filesystem::create_directories(dir);
        const auto observations=tvt::data();
        std::cout<<"TIME-VARYING TWO-TEMPERATURE RECONSTRUCTION\n"
                 <<"One fluid thermometer; estimate initial fluid/wall temperatures and three heater knots.\n";
        std::ofstream convergence(dir/"convergence.csv");convergence<<"theta,steps,trajectory_error,fit_error,gradient_fd_error\n";
        std::ofstream parameters(dir/"parameter_refinement.csv");parameters<<"theta,steps,parameter,truth,estimate,error\n";
        std::ofstream residuals(dir/"observation_residuals.csv");residuals<<"theta,time,observation,prediction,residual,measurement_scale\n";
        std::ofstream spectra(dir/"information_spectrum.csv");spectra<<"theta,steps,rank,condition,eigenvalue_0,eigenvalue_1,eigenvalue_2,eigenvalue_3,eigenvalue_4\n";
        std::ofstream uncertainty(dir/"local_parameter_uncertainty.csv");uncertainty<<"theta,steps,parameter,local_standard_deviation\n";
        for(double theta:{.5,.65}){
            std::array<std::size_t,4> meshes{80,160,320,640};std::array<Vec<5>,4> fit{};std::array<double,4> err{};
            for(std::size_t k=0;k<meshes.size();++k){auto sys=tvt::make_system(meshes[k],theta,observations);
                auto opt=minimize_box_constrained(sys,{30.,25.,.5,.5,.5},{0.,0.,0.,0.,0.},{60.,60.,2.,2.,2.},OptimizeOptions{2e-9,30,50});
                fit[k]=opt.control;err[k]=tvt::trajectory_error(sys,tvt::truth);
                convergence<<theta<<','<<meshes[k]<<','<<std::setprecision(12)<<err[k]<<','<<norm(fit[k]-tvt::truth)<<','<<tvt::gradient_fd(sys,fit[k])<<'\n';
                static constexpr std::array<const char*,5> names{"fluid_initial","wall_initial","heater_knot_0","heater_knot_5","heater_knot_10"};
                for(std::size_t j=0;j<5;++j)parameters<<theta<<','<<meshes[k]<<','<<names[j]<<','<<tvt::truth[j]<<','<<opt.control[j]<<','<<opt.control[j]-tvt::truth[j]<<'\n';
                if(k+1==meshes.size()){
                    auto path=sys.forward(opt.control);std::ofstream states(dir/(theta==.5?"state_theta_050.csv":"state_theta_065.csv"));
                    states<<"time,fluid_reference,wall_reference,fluid_fit,wall_fit,heater_truth,heater_fit\n";
                    for(std::size_t n=0;n<path.state.size();++n){double t=sys.grid().time[n];auto ref=tvt::exact(t,tvt::truth);states<<t<<','<<ref[0]<<','<<ref[1]<<','<<path.state[n][0]<<','<<path.state[n][1]<<','<<tvt::control(t,tvt::truth)<<','<<tvt::control(t,opt.control)<<'\n';}
                    const auto spectrum=tvt::information(sys,opt.control,path,observations);
                    spectra<<theta<<','<<meshes[k]<<','<<spectrum.rank<<','<<(spectrum.rank==5?spectrum.eigenvalues[0]/spectrum.eigenvalues[4]:std::numeric_limits<double>::infinity());
                    for(double eigenvalue:spectrum.eigenvalues)spectra<<','<<eigenvalue;
                    spectra<<'\n';
                    static constexpr std::array<const char*,5> names{"fluid_initial","wall_initial","heater_knot_0","heater_knot_5","heater_knot_10"};
                    for(std::size_t j=0;j<5;++j){double variance=0.;for(std::size_t mode=0;mode<5;++mode)
                            variance+=spectrum.directions[mode][j]*spectrum.directions[mode][j]/spectrum.eigenvalues[mode];
                        uncertainty<<theta<<','<<meshes[k]<<','<<names[j]<<','<<std::sqrt(variance)<<'\n';}
                    for(const auto& [t,y]:observations){const auto node=static_cast<std::size_t>(std::llround(t*static_cast<double>(meshes[k])/10.));
                        const double predicted=path.state[node][0];residuals<<theta<<','<<t<<','<<y<<','<<predicted<<','<<predicted-y<<','<<tvt::sigma<<'\n';}
                    std::cout<<"  theta="<<theta<<" fitted (Tfluid0,Twall0,u0,u5,u10)=(";
                    for(std::size_t j=0;j<5;++j)std::cout<<(j?", ":"")<<opt.control[j];
                    std::cout<<") gradient="<<norm(sys.evaluate(opt.control).gradient)
                             <<", information rank="<<spectrum.rank<<"/5, condition="
                             <<(spectrum.rank==5?spectrum.eigenvalues[0]/spectrum.eigenvalues[4]:std::numeric_limits<double>::infinity())<<"\n";
                    if(spectrum.rank!=5)throw std::runtime_error("The sampled experiment lost local parameter identifiability");
                    if(tvt::gradient_fd(sys,opt.control)>2e-6)throw std::runtime_error("Discrete adjoint gradient check failed");
                }
            }
            if(theta==.5){double p=tvt::order(err[1],err[2]);std::cout<<"  theta=0.5 analytical trajectory convergence order="<<p<<"\n";if(std::abs(p-2.)>.15)throw std::runtime_error("Expected second-order theta trajectory convergence");}
            else {double p=tvt::order(err[1],err[2]);std::cout<<"  theta=0.65 analytical trajectory convergence order="<<p<<"\n";if(std::abs(p-1.)>.2)throw std::runtime_error("Expected first-order off-centered theta convergence");}
        }
        std::cout<<"Results: \""<<dir.string()<<"\"\n";return 0;
    }catch(const std::exception& e){std::cerr<<"time_varying_thermal: "<<e.what()<<'\n';return 1;}
}
