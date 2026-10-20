#pragma once
#include "optimize.hpp"

namespace theta_opt {
// Spectral diagnostics in the coordinates/scales supplied by the caller.
// Small two-control PSD utility; no diagonal regularization or null-space prior.
struct PsdSpectrum2 {
    Vec<2> eigenvalues{};
    std::array<Vec<2>,2> directions{{{1.,0.},{0.,1.}}};
    std::size_t rank{};
    double cutoff{};
    Vec<2> solve_compatible(const Vec<2>& rhs,double relative_tolerance=1e-10) const {
        (void)norm(rhs);
        if(!std::isfinite(relative_tolerance)||relative_tolerance<=0)throw std::invalid_argument("Positive compatibility tolerance required");
        Vec<2> result{};double size=std::hypot(rhs[0],rhs[1]);
        for(std::size_t j=0;j<2;++j){double component=dot(directions[j],rhs);if(j<rank)result=result+(component/eigenvalues[j])*directions[j];else if(std::abs(component)>relative_tolerance*size)throw std::domain_error("Right hand side has unresolved null-space component");}
        return result;
    }
};
inline PsdSpectrum2 psd_spectrum_2d(const Mat<2,2>& H,double relative_cutoff=1e-12){
    if(!std::isfinite(relative_cutoff)||relative_cutoff<=0||relative_cutoff>=1)throw std::invalid_argument("Rank cutoff must lie in (0,1)");
    double scale=0.;for(auto row:H){(void)norm(row);for(double value:row)scale=std::max(scale,std::abs(value));}
    PsdSpectrum2 result;if(scale==0)return result;
    if(std::abs(H[0][1]-H[1][0])/scale>relative_cutoff)throw std::domain_error("Symmetric information matrix required");
    double a=H[0][0]/scale,d=H[1][1]/scale,b=.5*(H[0][1]/scale+H[1][0]/scale);
    double large=.5*(a+d+std::hypot(a-d,2*b));
    if(large<=0)throw std::domain_error("Positive semidefinite matrix required");
    double small=(a*d-b*b)/large;
    if(small < -relative_cutoff*large)throw std::domain_error("Indefinite information matrix");
    Vec<2> v{b,large-a},other{large-d,b};if(std::hypot(other[0],other[1])>std::hypot(v[0],v[1]))v=other;
    double length=std::hypot(v[0],v[1]);if(length>0)v=(1/length)*v;else v={1.,0.};
    if(v[0]<0||(v[0]==0&&v[1]<0))v=-1.*v;
    result.directions={v,Vec<2>{-v[1],v[0]}};
    result.eigenvalues={scale*large,scale*std::max(0.,small)};result.cutoff=relative_cutoff*result.eigenvalues[0];
    if(!std::isfinite(result.eigenvalues[0]))throw std::overflow_error("Information spectrum overflow");
    result.rank=1+(result.eigenvalues[1]>result.cutoff?1:0);return result;
}
struct QuadraticOptimum2 {Optimum<2> optimum;PsdSpectrum2 spectrum;};
// For unconstrained convex QUADRATIC reduced costs only. Preserve the seed's
// Euclidean null-space component; this choice is a gauge, not measured data.
// Throws if the supplied problem is incompatible, indefinite, nonquadratic at
// the tested endpoints, or not stationary after the projected correction.
template<std::size_t N>
QuadraticOptimum2 minimize_quadratic_psd_2d(const System<N,2>& system,const Vec<2>& seed,double gradient_tolerance=1e-8,double relative_cutoff=1e-12){
    if(!std::isfinite(gradient_tolerance)||gradient_tolerance<=0)throw std::invalid_argument("Positive gradient tolerance required");
    auto initial=system.evaluate(seed);Mat<2,2> H{};for(int j=0;j<2;++j){Vec<2> e{};e[j]=1.;auto h=system.hessian_vector(seed,e);for(int i=0;i<2;++i)H[i][j]=h[i];}
    QuadraticOptimum2 result;result.spectrum=psd_spectrum_2d(H,relative_cutoff);
    auto step=result.spectrum.solve_compatible(-1.*initial.gradient);auto control=seed+step;auto final=system.evaluate(control);
    double change=0.,hscale=std::max(norm(H[0]),norm(H[1]));for(int j=0;j<2;++j){Vec<2> e{};e[j]=1.;auto h=system.hessian_vector(control,e);for(int i=0;i<2;++i)change=std::max(change,std::abs(h[i]-H[i][j]));}
    if(change>1e-10*hscale)throw std::domain_error("Quadratic optimizer requires constant Hessian");
    double slack=64*std::numeric_limits<double>::epsilon()*std::max(1.,std::abs(initial.cost));
    if(norm(final.gradient)>gradient_tolerance||final.cost>initial.cost+slack)throw std::runtime_error("Projected quadratic correction did not reach a minimum");
    result.optimum.control=control;result.optimum.trace={{seed,initial.cost,norm(initial.gradient),1.},{control,final.cost,norm(final.gradient),0.}};return result;
}
}
