// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Richard Myers

#include "oscillator_model.hpp"
#include "oscillator_qualification_model.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

// Continuously deform a well-spread oscillator sampling design into the
// half-period design that is blind to initial velocity. Each fit uses the
// actual midpoint theta-discrete objective; exact formulas independently
// measure physical information and bounded-data uncertainty.
namespace {
using namespace oscillator;
namespace oq=oscillator::qualification;
using Schedule=std::array<double,4>;
inline constexpr Schedule proper{.3,.9,1.8,3.2};
inline constexpr Schedule blind{std::numbers::pi/2.,std::numbers::pi,
                                3*std::numbers::pi/2.,2*std::numbers::pi};
inline constexpr std::array<double,4> noise{.4,-.6,.2,.8};
inline constexpr Vec<2> scale{.02,.12}, prior{4.,2.5};

Schedule schedule(double alpha){Schedule t{};for(std::size_t j=0;j<4;++j)t[j]=(1-alpha)*proper[j]+alpha*blind[j];return t;}
double datum(std::size_t j,double tau){return exact(tau,truth)[0]+sigma*noise[j];}
std::size_t node(const Grid& g,double tau){for(std::size_t n=0;n<g.time.size();++n)if(std::abs(g.time[n]-t0-tau)<2e-12)return n;throw std::logic_error("sample time missing from mesh");}
Grid mesh(std::size_t n,double theta,const Schedule& t){std::vector<double> x;for(std::size_t j=0;j<=n;++j)x.push_back(t0+horizon*static_cast<double>(j)/n);for(double tau:t)x.push_back(t0+tau);std::sort(x.begin(),x.end());x.erase(std::unique(x.begin(),x.end(),[](double a,double b){return std::abs(a-b)<2e-14;}),x.end());const auto steps=x.size()-1;return Grid(std::move(x),std::vector<double>(steps,theta));}
System<2,2> make_system(std::size_t n,double theta,const Schedule& t){auto g=mesh(n,theta,t);Problem<2,2> p;
    p.initial.value=[](const Vec<2>& u){return Vec<2>{scale[0]*u[0],scale[1]*u[1]};};
    p.initial.jacobian=[](const Vec<2>&){return Mat<2,2>{{{scale[0],0.},{0.,scale[1]}}};};
    p.initial.directional_jacobian=[](const Vec<2>&,const Vec<2>&){return Mat<2,2>{};};
    p.dynamics.value=[](double,const Vec<2>& x,const Vec<2>&){return mv(A,x);};
    p.dynamics.derivatives=[](double,const Vec<2>&,const Vec<2>&){Derivatives<2,2> d;d.x=A;return d;};
    p.dynamics.directional_derivatives=[](double,const Vec<2>&,const Vec<2>&,const Vec<2>&,const Vec<2>&){return Derivatives<2,2>{};};
    for(std::size_t j=0;j<t.size();++j){const double y=datum(j,t[j]);auto c=ScalarFunction<2,2>::zero();
        c.evaluate=[y](double,const Vec<2>& x,const Vec<2>&){ScalarValue<2,2> v;double r=x[0]-y;v.value=.5*r*r/(sigma*sigma);v.x[0]=r/(sigma*sigma);return v;};
        c.directional_gradient=[](double,const Vec<2>&,const Vec<2>&,const Vec<2>& s,const Vec<2>&){ScalarDirection<2,2> v;v.x[0]=s[0]/(sigma*sigma);return v;};p.nodal.emplace_back(node(g,t[j]),c);}
    return System<2,2>(std::move(p),std::move(g));}
Mat<2,2> exact_information(const Schedule& t){Mat<2,2> H{};for(double tau:t){auto r=propagator(tau)[0];Vec<2> g{r[0]*scale[0],r[1]*scale[1]};H=add(H,scaled(1/(sigma*sigma),outer(g)));}return H;}
ConvexPolygon2 feasible_set(const Schedule& t){auto p=rectangle_polygon(-1.*prior,prior);for(std::size_t j=0;j<t.size();++j){auto r=propagator(t[j])[0];Vec<2> g{r[0]*scale[0],r[1]*scale[1]};double y=datum(j,t[j]);p=clip_halfspace(std::move(p),g,y+sigma);p=clip_halfspace(std::move(p),-1.*g,-y+sigma);}return p;}
Vec<2> centroid(const ConvexPolygon2& p){Vec<2> c{};for(auto x:p)c=c+(1./static_cast<double>(p.size()))*x;return c;}
double halfwidth(const ConvexPolygon2& p,int j){auto c=centroid(p);double w=0.;for(auto x:p)w=std::max(w,std::abs((x[j]-c[j])*scale[j]));return w;}
double information_sd(const Mat<2,2>& H,int j){auto e=psd_spectrum_2d(H,1e-13);if(e.rank<2)return std::numeric_limits<double>::infinity();auto C=inverse(H);return std::sqrt(C[j][j])*scale[j];}
}

int main(int argc,char** argv){try{
    const std::filesystem::path out=argc>1?argv[1]:"results/oscillator_information_sweep";std::filesystem::create_directories(out);
    std::ofstream csv(out/"sweep.csv");csv<<std::setprecision(17)<<"alpha,first_sample_s,last_sample_s,analytic_rank,analytic_lambda_small,analytic_lambda_large,condition_number,discrete_lambda_small,discrete_lambda_large,q0_fit_m,v0_fit_m_s,gradient_norm,q0_bounded_halfwidth_m,v0_bounded_halfwidth_m_s,q0_local_sd_m,v0_local_sd_m_s,q0_exact_qualified_t0,v0_exact_qualified_t0,support_force_exact_qualified_t0,joint_exact_qualification_time_s\n";
    std::ofstream samples(out/"sampling_schedules.csv");samples<<std::setprecision(17)<<"alpha,sample,elapsed_s\n";
    std::size_t checks=0;auto check=[&](const std::string& name,double e,double tol){++checks;if(!std::isfinite(e)||e>tol)throw std::runtime_error(name+" failed: "+std::to_string(e));};
    double last_width=-1.;
    for(int k=0;k<=20;++k){double a=k/20.;auto t=schedule(a);for(std::size_t j=0;j<4;++j)samples<<a<<','<<j<<','<<t[j]<<'\n';
        auto H=exact_information(t);auto spectrum=psd_spectrum_2d(H,1e-11);auto set=feasible_set(t);if(set.size()<3)throw std::runtime_error("empty bounded-data feasible set");
        auto sys=make_system(512,.5,t);auto opt=minimize(sys,Vec<2>{0.,0.},OptimizeOptions{1e-8,8,50});double grad=norm(sys.evaluate(opt.control).gradient);auto physical=Vec<2>{scale[0]*opt.control[0],scale[1]*opt.control[1]};Mat<2,2> Hd{};for(int j=0;j<2;++j){Vec<2> d{};d[j]=1.;auto col=sys.hessian_vector({0.,0.},d);for(int i=0;i<2;++i)Hd[i][j]=col[i];}auto ds=psd_spectrum_2d(Hd,1e-14);
        double qw=halfwidth(set,0),vw=halfwidth(set,1),qsd=information_sd(H,0),vsd=information_sd(H,1);
        auto c=centroid(set);std::array<double,3> width{};for(std::size_t j=0;j<3;++j){Vec<2> gu{scale[0]*oq::goals[j][0],scale[1]*oq::goals[j][1]};for(auto v:set)width[j]=std::max(width[j],std::abs(dot(gu,v-c)));}
        std::array<bool,3> qualified{};for(std::size_t j=0;j<3;++j)qualified[j]=width[j]<=oq::tolerances[j];
        ConvexPolygon2 physical_set;for(auto v:set)physical_set.push_back({scale[0]*v[0],scale[1]*v[1]});std::array<std::optional<double>,3> qt{};auto center=Vec<2>{scale[0]*c[0],scale[1]*c[1]};for(std::size_t j=0;j<3;++j)qt[j]=oq::analytical_persistent_time(physical_set,center,oq::goals[j],oq::tolerances[j]);
        std::optional<double> joint;for(auto q:qt)if(!q){joint.reset();break;}else joint=std::max(joint.value_or(0.),*q);
        csv<<a<<','<<t.front()<<','<<t.back()<<','<<spectrum.rank<<','<<spectrum.eigenvalues[1]<<','<<spectrum.eigenvalues[0]<<','<<(spectrum.eigenvalues[1]>0?spectrum.eigenvalues[0]/spectrum.eigenvalues[1]:std::numeric_limits<double>::infinity())<<','<<ds.eigenvalues[1]<<','<<ds.eigenvalues[0]<<','<<physical[0]<<','<<physical[1]<<','<<grad<<','<<qw<<','<<vw<<','<<qsd<<','<<vsd<<','<<qualified[0]<<','<<qualified[1]<<','<<qualified[2]<<',';if(joint)csv<<*joint;else csv<<"nan";csv<<'\n';
        if(k==0){check("proper schedule full analytical rank",spectrum.rank==2?0.:1.,0.);check("proper schedule identifies bounded velocity",vw<.03?0.:1.,0.);}if(k==20){check("half-period schedule rank one",spectrum.rank==1?0.:1.,0.);check("blind endpoint retains velocity prior",std::abs(vw-scale[1]*prior[1]),2e-11);check("blind endpoint cannot qualify velocity",qualified[1]?1.:0.,0.);check("midpoint phase leakage creates false curvature",ds.rank==2&&ds.eigenvalues[1]>1e-5?0.:1.,0.);check("stationary blind fit does not recover velocity",std::abs(physical[1]-truth[1])>1.?0.:1.,0.);check("blind fit is nevertheless stationary",grad,2e-7);last_width=vw;}
        check("discrete optimizer stationarity",grad,2e-7);
        if(k==10){Vec<2> u{.1,-.2},d{.03,-.04};auto ev=sys.evaluate(u);auto tangent=sys.tangent(u,ev.trajectory,d);double pair=0.;for(std::size_t j=0;j<t.size();++j){auto n=node(sys.grid(),t[j]);pair+=(ev.trajectory.state[n][0]-datum(j,t[j]))*tangent[n][0]/(sigma*sigma);}check("sweep tangent-adjoint pairing",std::abs(pair-dot(ev.gradient,d)),2e-8);double eps=1e-5;auto ep=sys.evaluate(u+eps*d),em=sys.evaluate(u-eps*d);check("sweep discrete objective gradient",std::abs((ep.cost-em.cost)/(2*eps)-dot(ev.gradient,d)),2e-5);check("sweep incremental-adjoint Hessian",norm(sys.hessian_vector(u,d)-(1/(2*eps))*(ep.gradient-em.gradient)),2e-5);}
    }
    if(!(last_width>0.))throw std::runtime_error("blind-design uncertainty did not remain positive");
    std::cout<<"WEAK-INFORMATION OSCILLATOR DESIGN SWEEP\n"
             <<"Interpolates four displacement sample times from a full-rank schedule to exact half-periods.\n"
             <<"Every point uses a midpoint theta-discrete fit; exact information and bounded-data sets assess identifiability.\n"
             <<"Qualification columns use the exact physical observation model with bounded errors and no discretization allowance.\n"
             <<"At the blind endpoint: exact rank 1/2 and velocity remains unqualified despite optimizer stationarity.\n"
             <<"Checks passed: "<<checks<<". Results: \""<<out.string()<<"\"\n";
    return 0;
}catch(const std::exception& e){std::cerr<<"oscillator_information_sweep: "<<e.what()<<'\n';return 1;}}
