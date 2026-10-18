// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include <theta_opt/theta_opt.hpp>
#include <thermogpu/peng_robinson.hpp>
#include <algorithm>

namespace gas_inverse {
using namespace theta_opt;
using theta_opt::operator+;
using theta_opt::operator-;
using theta_opt::operator*;
inline constexpr Vec<2> truth{.10,.05}; // ethane, nitrogen; methane is the remainder
inline const std::vector<thermogpu::Component> components{
    {"methane",190.564,4599200.,.01142,.01604246},
    {"ethane",305.322,4872200.,.0995,.03006904},
    {"nitrogen",126.192,3395800.,.0372,.0280134}};
inline bool admissible(const Vec<2>& x){return x[0]>0&&x[1]>0&&x[0]+x[1]<1;}
inline thermogpu::Mixture mixture(const Vec<2>& x){return {components,{1-x[0]-x[1],x[0],x[1]},std::vector<double>(9,0.)};}
inline std::array<double,3> fractions(const Vec<2>& x){return {1-x[0]-x[1],x[0],x[1]};}
// Small second-order forward AD for the derivative adapter only. All state
// values and root selection come from ThermoGPU, not this adapter.
struct Jet {double v{};Vec<2> g{};Mat<2,2> h{};Jet()=default;Jet(double value):v(value){};};
inline Jet operator+(Jet a,const Jet& b){a.v+=b.v;a.g=a.g+b.g;a.h=add(a.h,b.h);return a;}
inline Jet operator*(const Jet& a,const Jet& b){Jet c(a.v*b.v);c.g=b.v*a.g+a.v*b.g;for(int i=0;i<2;++i)for(int j=0;j<2;++j)c.h[i][j]=a.h[i][j]*b.v+a.v*b.h[i][j]+a.g[i]*b.g[j]+b.g[i]*a.g[j];return c;}
inline Jet operator-(Jet a,const Jet& b){return a+Jet(-1.)*b;}
inline Jet reciprocal(const Jet& a){Jet b(1/a.v);b.g=(-1/(a.v*a.v))*a.g;for(int i=0;i<2;++i)for(int j=0;j<2;++j)b.h[i][j]=2*a.g[i]*a.g[j]/(a.v*a.v*a.v)-a.h[i][j]/(a.v*a.v);return b;}
inline Jet operator/(const Jet& a,const Jet& b){return a*reciprocal(b);}
inline Jet variable(double v,int j){Jet a(v);a.g[j]=1;return a;}
struct Property {double value{},Z{},residual{},root_slope{};Vec<2> gradient{};Mat<2,2> hessian{};};
inline Property density(const thermogpu::State& state,const Vec<2>& control){
    auto d=thermogpu::evaluate_mixture_diagnostics(mixture(control),state);auto r=d.result;
    if(r.root_selection_ambiguous||r.root_classification!=thermogpu::RootClassification::one_real)throw std::domain_error("Experiment requires a unique PR root");
    std::array<Jet,3> x{Jet(1.)-variable(control[0],0)-variable(control[1],1),variable(control[0],0),variable(control[1],1)};
    Jet a,b,m;
    for(int i=0;i<3;++i){b=b+x[i]*Jet(d.b[i]);m=m+x[i]*Jet(components[i].molar_mass_kg_per_mol);for(int j=0;j<3;++j)a=a+x[i]*x[j]*Jet(std::sqrt(d.a_alpha[i]*d.a_alpha[j]));}
    const double RT=thermogpu::gas_constant*state.temperature_K;
    Jet A=a*Jet(state.pressure_Pa/(RT*RT)),B=b*Jet(state.pressure_Pa/RT),z(r.Z);
    const double slope=3*r.Z*r.Z+2*(r.B-1)*r.Z+r.A-3*r.B*r.B-2*r.B;
    if(std::abs(slope)<1e-5)throw std::domain_error("Near-multiple PR root");
    // At fixed Z, F carries total composition partials. F_Z carries its
    // composition partials too. Implicit differentiation supplies Z_x and Z_xx.
    Jet F=z*z*z+(B-Jet(1))*z*z+(A-Jet(3)*B*B-Jet(2)*B)*z-A*B+B*B+B*B*B;
    Jet Fz=Jet(3)*z*z+Jet(2)*(B-Jet(1))*z+A-Jet(3)*B*B-Jet(2)*B;
    z.g=(-1/slope)*F.g;
    double Fzz=6*r.Z+2*(r.B-1);
    for(int i=0;i<2;++i)for(int j=0;j<2;++j)z.h[i][j]=-(F.h[i][j]+Fz.g[i]*z.g[j]+Fz.g[j]*z.g[i]+Fzz*z.g[i]*z.g[j])/slope;
    Jet rho=Jet(state.pressure_Pa/RT)*m/z;
    return {r.density_kg_per_m3,r.Z,r.cubic_residual,slope,rho.g,rho.h};
}
struct Observation {thermogpu::State state;double value{},sigma{},bound{};};
inline std::vector<Observation> observations(bool noisy=true,double relative_precision=.001){
    std::array<thermogpu::State,4> states{{{1e6,300.},{5e6,300.},{1e7,350.},{2e7,400.}}};
    std::array<double,4> perturbation{.4,-.6,.2,.8};std::vector<Observation> result;
    for(int j=0;j<4;++j){double rho=density(states[j],truth).value,scale=relative_precision*rho;result.push_back({states[j],rho+(noisy?perturbation[j]*scale:0.),scale,scale});}return result;
}
inline ScalarFunction<2,2> measurement(const Observation& o){
    auto f=ScalarFunction<2,2>::zero();
    f.evaluate=[o](double,const Vec<2>& x,const Vec<2>&){auto p=density(o.state,x);double r=p.value-o.value;ScalarValue<2,2> v;v.value=.5*r*r/(o.sigma*o.sigma);v.x=(r/(o.sigma*o.sigma))*p.gradient;return v;};
    f.directional_gradient=[o](double,const Vec<2>& x,const Vec<2>&,const Vec<2>& s,const Vec<2>&){auto p=density(o.state,x);ScalarDirection<2,2> v;v.x=(1/(o.sigma*o.sigma))*(dot(p.gradient,s)*p.gradient+(p.value-o.value)*mv(p.hessian,s));return v;};return f;
}
inline System<2,2> system(const std::vector<Observation>& obs,double theta=.5){
    // Static specialization: composition is unchanged between experiments.
    // Node coordinates are observation indices, not elapsed physical time.
    Problem<2,2> p;p.initial.value=[](const Vec<2>& x){return x;};p.initial.jacobian=[](const Vec<2>&){return identity<2>();};p.initial.directional_jacobian=[](const Vec<2>&,const Vec<2>&){return Mat<2,2>{};};
    p.dynamics.value=[](double,const Vec<2>&,const Vec<2>&){return Vec<2>{};};p.dynamics.derivatives=[](double,const Vec<2>&,const Vec<2>&){return Derivatives<2,2>{};};p.dynamics.directional_derivatives=[](double,const Vec<2>&,const Vec<2>&,const Vec<2>&,const Vec<2>&){return Derivatives<2,2>{};};
    for(std::size_t j=0;j<obs.size();++j)p.nodal.emplace_back(j,measurement(obs[j]));
    auto grid=Grid::uniform(0.,static_cast<double>(obs.size()-1),obs.size()-1,theta);return System<2,2>(p,grid);
}
inline Optimum<2> fit(const System<2,2>& s,Vec<2> guess={.07,.08}){return minimize(s,guess,OptimizeOptions{1e-7,60,50,1e-11},admissible);}
inline Mat<2,2> information(const System<2,2>& sys,const std::vector<Observation>& obs,const Vec<2>& xi){
    auto traj=sys.forward(xi);Mat<2,2> I{};
    for(std::size_t j=0;j<obs.size();++j){auto p=density(obs[j].state,traj.state[j]);auto g=sys.goal_gradient(xi,traj,j,p.gradient);for(int a=0;a<2;++a)for(int b=0;b<2;++b)I[a][b]+=g[a]*g[b]/(obs[j].sigma*obs[j].sigma);}return I;
}
}
