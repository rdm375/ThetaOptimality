#include "../examples/decay_model.hpp"
#include <iomanip>
#include <iostream>
#include <string>
#include <theta_opt/theta_opt.hpp>

using namespace theta_opt;
int checks=0;
void check(bool ok,const std::string& name){++checks;if(!ok)throw std::runtime_error("FAIL: "+name);}
void close(double a,double b,double tol,const std::string& name){check(std::isfinite(a)&&std::isfinite(b)&&std::abs(a-b)<tol,name);}
template<class F> void rejects(F f,const std::string& name){bool caught=false;try{f();}catch(const std::exception&){caught=true;}check(caught,name);}

Problem<2,2> coupled(){
    Problem<2,2> p;
    p.initial.value=[](const Vec<2>& v){return Vec<2>{v[0]*v[0]+.5*v[1],std::sin(v[1])};};
    p.initial.jacobian=[](const Vec<2>& v){return Mat<2,2>{{{2*v[0],.5},{0.,std::cos(v[1])}}};};
    p.initial.directional_jacobian=[](const Vec<2>& v,const Vec<2>& d){return Mat<2,2>{{{2*d[0],0.},{0.,-std::sin(v[1])*d[1]}}};};
    p.dynamics.value=[](double t,const Vec<2>& x,const Vec<2>& v){return Vec<2>{-v[1]*x[0]+std::sin(t)*x[1]+v[0]*v[0],x[0]*x[1]-.3*x[1]+v[0]*std::sin(t)};};
    p.dynamics.derivatives=[](double t,const Vec<2>& x,const Vec<2>& v){return Derivatives<2,2>{{{{-v[1],std::sin(t)},{x[1],x[0]-.3}}},{{{2*v[0],-x[0]},{std::sin(t),0.}}}};};
    p.dynamics.directional_derivatives=[](double,const Vec<2>&,const Vec<2>&,const Vec<2>& s,const Vec<2>& d){return Derivatives<2,2>{{{{-d[1],0.},{s[1],s[0]}}},{{{2*d[0],-s[0]},{0.,0.}}}};};
    p.running.evaluate=[](double t,const Vec<2>& x,const Vec<2>& v){double r=x[0]-v[0],c=std::cos(t);return ScalarValue<2,2>{.5*r*r+.7*x[1]*x[1]+.1*v[1]*v[1]+.2*v[1]*x[1]*c,{r,1.4*x[1]+.2*v[1]*c},{-r,.2*v[1]+.2*x[1]*c}};};
    p.running.directional_gradient=[](double t,const Vec<2>&,const Vec<2>&,const Vec<2>& s,const Vec<2>& d){double c=std::cos(t);return ScalarDirection<2,2>{{s[0]-d[0],1.4*s[1]+.2*d[1]*c},{d[0]-s[0],.2*d[1]+.2*s[1]*c}};};
    p.terminal.evaluate=[](double,const Vec<2>& x,const Vec<2>& v){double r=x[1]+v[0]*v[1];return ScalarValue<2,2>{.5*r*r+.1*x[0]*v[1],{.1*v[1],r},{r*v[1],r*v[0]+.1*x[0]}};};
    p.terminal.directional_gradient=[](double,const Vec<2>& x,const Vec<2>& v,const Vec<2>& s,const Vec<2>& d){double r=x[1]+v[0]*v[1],dr=s[1]+d[0]*v[1]+v[0]*d[1];return ScalarDirection<2,2>{{.1*d[1],dr},{dr*v[1]+r*d[1],dr*v[0]+r*d[0]+.1*s[0]}};};
    p.penalty.value=[](const Vec<2>& v){return .05*(v[0]*v[0]+2*v[1]*v[1]);};
    p.penalty.gradient=[](const Vec<2>& v){return Vec<2>{.1*v[0],.2*v[1]};};
    p.penalty.hessian_vector=[](const Vec<2>&,const Vec<2>& d){return Vec<2>{.1*d[0],.2*d[1]};};return p;
}

int main(){try{
    {
        Box2 box{{-2.,-1.},{2.,1.}};Quadratic2 square{0.,{},{{{1.,0.},{0.,1.}}}};
        auto range=quadratic_range(square,box);check(range.lower<=0&&range.upper>=5.,"quadratic box range encloses corners");
        Quadratic2 linear{1.,{2.,-3.},{}};auto lr=quadratic_range(linear,box);check(lr.lower<=-6.&&lr.upper>=8.,"linear box range exact up to slack");
        check(contains(box,{2.,1.})&&!contains(box,{3.,1.}),"box membership boundaries");
        auto cover=subdivide_feasible_boxes(box,{.03,.03},[](const Box2& b){Quadratic2 q{0.,{},{{{1.,0.},{0.,0.}}}};auto r=quadratic_range(q,b);return r.lower<=1.01&&r.upper>=.99;});
        check(cover.resolution_complete&&!cover.boxes.empty(),"disconnected feasible cover resolves");
        check(cover_contains(cover,{1.,0.})&&cover_contains(cover,{-1.,0.}),"both solution branches retained");
        auto union_range=cover_quadratic_range(cover,Quadratic2{0.,{1.,0.},{}});check(union_range.lower<=-1.&&union_range.upper>=1.,"union range retains opposite signs");
        auto whole_union=cover_polynomial_quadratic_range(cover,{Quadratic2{0.,{1.,0.},{}},Quadratic2{1.,{},{}}});check(whole_union.lower<=-1.&&whole_union.upper>=2.,"whole interval union retains all branches");
        auto limited=subdivide_feasible_boxes(box,{.001,.001},[](const Box2&){return true;},1);check(!limited.resolution_complete&&cover_contains(limited,{1.7,.8}),"budget preserves unresolved feasible cells");
        auto empty=subdivide_feasible_boxes(box,{1.,1.},[](const Box2&){return false;});check(empty.boxes.empty(),"empty cover");
        rejects([&]{cover_quadratic_range(empty,square);},"empty union has no goal enclosure");
        auto pr=polynomial_quadratic_range({Quadratic2{1.,{},{}},Quadratic2{-2.,{},{}},Quadratic2{3.,{},{}}},box);check(pr.lower<=1&&pr.upper>=2.,"whole-time polynomial range");
        rejects([&]{subdivide_feasible_boxes(box,{0.,1.},[](const Box2&){return true;});},"reject invalid subdivision resolution");
        rejects([]{quadratic_range(Quadratic2{},Box2{{1.,0.},{0.,1.}});},"reject inverted box");
        rejects([&]{polynomial_quadratic_range({},box);},"reject empty temporal polynomial");
    }

    {
        auto rank1=psd_spectrum_2d(Mat<2,2>{{{2.,-2.},{-2.,2.}}});
        check(rank1.rank==1,"rank-one PSD information");
        close(rank1.eigenvalues[0],4.,1e-14,"nonzero information eigenvalue");
        close(norm(rank1.solve_compatible({2.,-2.})-Vec<2>{.5,-.5}),0.,1e-14,"compatible pseudoinverse solve");
        close(std::abs(dot(rank1.directions[1],Vec<2>{1.,1.})),std::sqrt(2.),1e-14,"missing common direction");
        rejects([&]{rank1.solve_compatible({1.,1.});},"reject incompatible null-space right hand side");
        auto full=psd_spectrum_2d(Mat<2,2>{{{3.,1.},{1.,2.}}});
        check(full.rank==2,"full-rank information");
        close(norm(mv(Mat<2,2>{{{3.,1.},{1.,2.}}},full.solve_compatible({2.,3.}))-Vec<2>{2.,3.}),0.,1e-14,"full rank spectral solve");
        auto zero=psd_spectrum_2d(Mat<2,2>{});check(zero.rank==0,"zero information rank");
        close(norm(zero.solve_compatible({0.,0.})),0.,1e-15,"zero compatible solve");
        rejects([&]{zero.solve_compatible({1.,0.});},"zero incompatible solve");
        rejects([]{psd_spectrum_2d(Mat<2,2>{{{1.,0.},{0.,-1.}}});},"reject indefinite spectrum");
        rejects([]{psd_spectrum_2d(Mat<2,2>{{{1.,1.},{0.,1.}}});},"reject nonsymmetric spectrum");
        rejects([]{psd_spectrum_2d(identity<2>(),0.);},"reject invalid rank cutoff");
        auto small=psd_spectrum_2d(Mat<2,2>{{{1e-100,0.},{0.,1e-105}}});check(small.rank==2,"rank is relative to information scale");
        auto weak=psd_spectrum_2d(Mat<2,2>{{{1.,0.},{0.,1e-14}}});check(weak.rank==1,"explicit numerical rank threshold");
        rejects([]{psd_spectrum_2d(Mat<2,2>{{{NAN,0.},{0.,1.}}});},"reject nonfinite information");
        {
            const Mat<3,3> H{{{2.,-2.,0.},{-2.,2.,0.},{0.,0.,3.}}};
            const auto spectrum=psd_spectrum(H);
            check(spectrum.rank==2,"generic spectrum detects rank-two three-control information");
            close(spectrum.eigenvalues[0],4.,1e-13,"generic spectrum largest eigenvalue");
            close(spectrum.eigenvalues[1],3.,1e-13,"generic spectrum second eigenvalue");
            close(norm(mv(H,spectrum.solve_compatible({2.,-2.,3.}))-Vec<3>{2.,-2.,3.}),0.,1e-12,
                  "generic compatible spectral solve");
            rejects([&]{spectrum.solve_compatible({1.,1.,0.});},"generic spectrum rejects null-space rhs");
            const auto diagonal=psd_spectrum(Mat<4,4>{{{1.,0.,0.,0.},{0.,2.,0.,0.},{0.,0.,3.,0.},{0.,0.,0.,4.}}});
            check(diagonal.rank==4,"generic spectrum supports four-control full rank");
            close(diagonal.eigenvalues[0],4.,1e-14,"generic spectrum sorts eigenvalues");
            rejects([]{psd_spectrum(Mat<3,3>{{{1.,0.,0.},{0.,1.,0.},{0.,0.,-1.}}});},"generic spectrum rejects indefinite matrix");
            rejects([]{psd_spectrum(Mat<3,3>{{{1.,1.,0.},{0.,1.,0.},{0.,0.,1.}}});},"generic spectrum rejects nonsymmetric matrix");
        }
        auto quadratic_problem=[](){Problem<1,2> p;p.initial.value=[](auto){return Vec<1>{0.};};p.initial.jacobian=[](auto){return Mat<1,2>{};};p.initial.directional_jacobian=[](auto,auto){return Mat<1,2>{};};p.dynamics.value=[](double,auto,auto){return Vec<1>{};};p.dynamics.derivatives=[](double,auto,auto){return Derivatives<1,2>{};};p.dynamics.directional_derivatives=[](double,auto,auto,auto,auto){return Derivatives<1,2>{};};return p;};
        auto empty_system=System<1,2>(quadratic_problem(),Grid::uniform(2.,3.,2,.5));
        close(norm(minimize_quadratic_psd_2d(empty_system,{2.,3.}).optimum.control-Vec<2>{2.,3.}),0.,1e-15,"zero-information optimizer retains representative");
        rejects([&]{minimize_quadratic_psd_2d(empty_system,{2.,3.},0.);},"reject invalid projected optimizer tolerance");
        auto linear=quadratic_problem();linear.penalty.value=[](auto x){return x[0];};linear.penalty.gradient=[](auto){return Vec<2>{1.,0.};};
        rejects([&]{minimize_quadratic_psd_2d(System<1,2>(linear,Grid::uniform(2.,3.,2,.5)),{2.,3.});},"reject unbounded null-direction linear cost");
        auto nonlinear=quadratic_problem();nonlinear.penalty.value=[](auto x){return .25*std::pow(x[0],4)+.5*x[1]*x[1];};nonlinear.penalty.gradient=[](auto x){return Vec<2>{std::pow(x[0],3),x[1]};};nonlinear.penalty.hessian_vector=[](auto x,auto d){return Vec<2>{3*x[0]*x[0]*d[0],d[1]};};
        rejects([&]{minimize_quadratic_psd_2d(System<1,2>(nonlinear,Grid::uniform(2.,3.,2,.5)),{2.,3.});},"reject changing Hessian in quadratic optimizer");

    }

    {
        auto bounded_quadratic=[](double target){Problem<1,1> p;
            p.initial.value=[](auto){return Vec<1>{0.};};p.initial.jacobian=[](auto){return Mat<1,1>{};};p.initial.directional_jacobian=[](auto,auto){return Mat<1,1>{};};
            p.dynamics.value=[](double,auto,auto){return Vec<1>{0.};};p.dynamics.derivatives=[](double,auto,auto){return Derivatives<1,1>{};};p.dynamics.directional_derivatives=[](double,auto,auto,auto,auto){return Derivatives<1,1>{};};
            p.penalty.value=[target](const Vec<1>& x){return .5*(x[0]-target)*(x[0]-target);};p.penalty.gradient=[target](const Vec<1>& x){return Vec<1>{x[0]-target};};p.penalty.hessian_vector=[](const Vec<1>&,const Vec<1>& d){return d;};return p;};
        const auto grid=Grid::uniform(0.,1.,2,.5);
        auto upper=minimize_box_constrained(System<1,1>(bounded_quadratic(2.),grid),{.5},{0.},{1.});
        close(upper.control[0],1.,1e-14,"box optimizer upper boundary");close(upper.gradient[0],-1.,1e-14,"upper-bound raw gradient sign");close(upper.upper_multiplier[0],1.,1e-14,"upper-bound KKT multiplier");check(upper.active_upper[0]&&upper.kkt_residual<1e-14,"upper-bound KKT stationarity");
        auto lower=minimize_box_constrained(System<1,1>(bounded_quadratic(-2.),grid),{.5},{0.},{1.});
        close(lower.control[0],0.,1e-14,"box optimizer lower boundary");close(lower.gradient[0],2.,1e-14,"lower-bound raw gradient sign");close(lower.lower_multiplier[0],2.,1e-14,"lower-bound KKT multiplier");check(lower.active_lower[0]&&lower.kkt_residual<1e-14,"lower-bound KKT stationarity");
        auto interior=minimize_box_constrained(System<1,1>(bounded_quadratic(.25),grid),{.5},{0.},{1.});close(interior.control[0],.25,1e-14,"box optimizer interior solution");check(!interior.active_lower[0]&&!interior.active_upper[0]&&interior.kkt_residual<1e-14,"interior KKT stationarity");
        rejects([&]{minimize_box_constrained(System<1,1>(bounded_quadratic(.25),grid),{1.5},{0.},{1.});},"reject box optimizer start outside bounds");
        rejects([&]{minimize_box_constrained(System<1,1>(bounded_quadratic(.25),grid),{.5},{1.},{0.});},"reject inverted optimizer bounds");

        Problem<1,1> path_problem;
        path_problem.initial.value=[](const Vec<1>&){return Vec<1>{0.};};
        path_problem.initial.jacobian=[](const Vec<1>&){return Mat<1,1>{};};
        path_problem.initial.directional_jacobian=[](const Vec<1>&,const Vec<1>&){return Mat<1,1>{};};
        path_problem.dynamics.value=[](double,const Vec<1>&,const Vec<1>& z){return Vec<1>{z[0]};};
        path_problem.dynamics.derivatives=[](double,const Vec<1>&,const Vec<1>&){
            Derivatives<1,1> d;d.control[0][0]=1.;return d;};
        path_problem.dynamics.directional_derivatives=[](double,const Vec<1>&,const Vec<1>&,
                                                          const Vec<1>&,const Vec<1>&){return Derivatives<1,1>{};};
        path_problem.terminal.evaluate=[](double,const Vec<1>& x,const Vec<1>&){
            const double r=x[0]-1.;return ScalarValue<1,1>{.5*r*r,Vec<1>{r},{}};};
        path_problem.terminal.directional_gradient=[](double,const Vec<1>&,const Vec<1>&,
                                                       const Vec<1>& dx,const Vec<1>&){
            return ScalarDirection<1,1>{Vec<1>{dx[0]}, {}};};
        auto path=ScalarFunction<1,1>::zero();
        path.evaluate=[](double,const Vec<1>& x,const Vec<1>&){return ScalarValue<1,1>{x[0]-.4,Vec<1>{1.},{}};};
        path.directional_gradient=[](double,const Vec<1>&,const Vec<1>&,const Vec<1>&,const Vec<1>&){return ScalarDirection<1,1>{};};
        path_problem.path_inequalities.push_back({1,path});
        auto path_system=System<1,1>(path_problem,grid);
        const auto path_result=minimize_path_constrained(path_system,{.2},{-1.},{2.});
        close(path_result.control[0],.8,2e-8,"path KKT constrained optimum");
        close(path_result.path_multiplier[0],.4,2e-7,"path KKT multiplier");
        check(path_result.path_value[0]<=1e-9&&path_result.kkt_residual<1e-8,
              "path KKT primal feasibility, stationarity, and complementarity");
        auto path_augmented=path_system.with_path_lagrangian({.4},12.);
        const auto aug_e=path_augmented.evaluate({.7});
        const auto aug_p=path_augmented.evaluate({.7+1e-5});
        const auto aug_m=path_augmented.evaluate({.7-1e-5});
        close(aug_e.gradient[0],(aug_p.cost-aug_m.cost)/2e-5,2e-8,
              "augmented path cost uses discrete adjoint gradient");
        close(path_augmented.hessian_vector({.7},{1.})[0],
              (aug_p.gradient[0]-aug_m.gradient[0])/2e-5,2e-8,
              "augmented path cost uses incremental-adjoint Hessian");

        SoftmaxCompositionControl<3,3> composition({0.,.4,1.});
        Vec<6> logits{-.4,.1,.6,.3,-.2,.2},direction{.2,-.1,.3,.1,.2,-.2};
        const double time=.63,step=1e-5;const auto fraction=composition.value(time,logits);
        close(fraction[0]+fraction[1]+fraction[2],1.,2e-15,"scheduled composition sums to one");
        check(fraction[0]>0.&&fraction[1]>0.&&fraction[2]>0.,"scheduled composition remains positive");
        const auto J=composition.jacobian(time,logits),dJ=composition.directional_jacobian(time,logits,direction);
        const auto Jp=composition.jacobian(time,logits+step*direction);
        const auto Jm=composition.jacobian(time,logits-step*direction);
        const auto yp=composition.value(time,logits+step*direction),ym=composition.value(time,logits-step*direction);
        close(norm(mv(J,direction)-(.5/step)*(yp-ym)),0.,2e-9,"softmax schedule Jacobian");
        for(std::size_t i=0;i<3;++i)for(std::size_t j=0;j<6;++j)
            close(dJ[i][j],(Jp[i][j]-Jm[i][j])/(2.*step),2e-9,
                  "softmax schedule directional Jacobian");
        auto raw_observation=ScalarFunction<1,3>::zero();
        raw_observation.evaluate=[](double,const Vec<1>&,const Vec<3>& u){
            const double r=u[1]-.35;return ScalarValue<1,3>{.5*r*r,{},Vec<3>{0.,r,0.}};};
        raw_observation.directional_gradient=[](double,const Vec<1>&,const Vec<3>&,
                                                  const Vec<1>&,const Vec<3>& du){
            return ScalarDirection<1,3>{{},{0.,du[1]}};};
        const auto scheduled_observation=compose_scalar_function(raw_observation,composition.parameterization());
        const auto ov=scheduled_observation.evaluate(time,Vec<1>{},logits);
        const auto od=scheduled_observation.directional_gradient(time,Vec<1>{},logits,Vec<1>{},direction);
        const auto ovp=scheduled_observation.evaluate(time,Vec<1>{},logits+step*direction);
        const auto ovm=scheduled_observation.evaluate(time,Vec<1>{},logits-step*direction);
        close(dot(ov.control,direction),(ovp.value-ovm.value)/(2.*step),2e-9,
              "scheduled scalar gradient chain rule");
        for(std::size_t j=0;j<6;++j)
            close(od.control[j],(ovp.control[j]-ovm.control[j])/(2.*step),2e-8,
                  "scheduled scalar Hessian chain rule");
        PiecewiseLinearControl<2,3> linear_schedule({0.,.5,1.});
        const Vec<6> control_knots{1.,2.,3.,4.,5.,6.};
        const auto scheduled_value=linear_schedule.value(.25,control_knots);
        close(scheduled_value[0],1.5,1e-15,"piecewise-linear control interpolation channel one");
        close(scheduled_value[1],4.5,1e-15,"piecewise-linear control interpolation channel two");
        close(norm(mv(linear_schedule.jacobian(.25,control_knots),direction)-
                   (.5/step)*(linear_schedule.value(.25,control_knots+step*direction)-
                              linear_schedule.value(.25,control_knots-step*direction))),0.,2e-10,
              "piecewise-linear schedule Jacobian");
        Dynamics<1,3> mixture_model;
        mixture_model.value=[](double,const Vec<1>&,const Vec<3>& u){return Vec<1>{u[0]+2.*u[1]-u[2]};};
        mixture_model.derivatives=[](double,const Vec<1>&,const Vec<3>&){
            Derivatives<1,3> d;d.control[0]={1.,2.,-1.};return d;};
        mixture_model.directional_derivatives=[](double,const Vec<1>&,const Vec<3>&,
                                                  const Vec<1>&,const Vec<3>&){return Derivatives<1,3>{};};
        const auto scheduled_dynamics=compose_dynamics(mixture_model,composition.parameterization());
        Problem<1,6> scheduled_problem;scheduled_problem.dynamics=scheduled_dynamics;
        scheduled_problem.initial.value=[](const Vec<6>&){return Vec<1>{0.};};
        scheduled_problem.initial.jacobian=[](const Vec<6>&){return Mat<1,6>{};};
        scheduled_problem.initial.directional_jacobian=[](const Vec<6>&,const Vec<6>&){return Mat<1,6>{};};
        scheduled_problem.terminal.evaluate=[](double,const Vec<1>& x,const Vec<6>&){
            const double r=x[0]-.25;return ScalarValue<1,6>{.5*r*r,Vec<1>{r},{}};};
        scheduled_problem.terminal.directional_gradient=[](double,const Vec<1>&,const Vec<6>&,
                                                             const Vec<1>& dx,const Vec<6>&){
            return ScalarDirection<1,6>{Vec<1>{dx[0]},{}};};
        System<1,6> scheduled_system(scheduled_problem,Grid::uniform(0.,1.,12,.5));
        const auto se=scheduled_system.evaluate(logits);
        const auto st=scheduled_system.tangent(logits,se.trajectory,direction);
        const auto sep=scheduled_system.evaluate(logits+step*direction);
        const auto sem=scheduled_system.evaluate(logits-step*direction);
        for(std::size_t n=0;n<st.size();++n)
            close(st[n][0],(sep.trajectory.state[n][0]-sem.trajectory.state[n][0])/(2.*step),2e-8,
                  "scheduled dynamics tangent finite difference");
        close(dot(se.gradient,direction),(sep.cost-sem.cost)/(2.*step),2e-8,
              "scheduled dynamics adjoint gradient finite difference");
        close(scheduled_system.hessian_vector(logits,direction)[0],
              (sep.gradient[0]-sem.gradient[0])/(2.*step),2e-7,
              "scheduled dynamics incremental-adjoint Hessian");
    }

    close(worst_amplification_2d(identity<2>()).gain,1.,1e-14,"identity amplification");
    close(worst_amplification_2d(Mat<2,2>{}).gain,0.,1e-14,"zero amplification");
    close(worst_amplification_2d(Mat<2,2>{{{1.,1.},{0.,1.}}}).gain,(1+std::sqrt(5.))/2,1e-14,"shear singular gain");
    auto weighted=worst_amplification_2d(Mat<2,2>{{{2.,0.},{0.,3.}}},{3.,2.},{2.,3.});
    close(weighted.gain,3.,1e-14,"weighted amplification");close(norm(weighted.input_direction-Vec<2>{3.,0.}),0.,1e-14,"physical worst direction");
    close(std::hypot(weighted.state_response[0]/2,weighted.state_response[1]/3),3.,1e-14,"weighted gain attained");
    rejects([]{worst_amplification_2d(identity<2>(),{0.,1.});},"zero amplification scale");
    rejects([]{worst_amplification_2d(Mat<2,2>{{{std::nan(""),0.},{0.,1.}}});},"nonfinite amplification");
    auto polygon=rectangle_polygon({-1.,-1.},{1.,1.});polygon=clip_halfspace(polygon,{1.,0.},.5);
    check(polygon_contains(polygon,{0.,0.}),"polygon contains interior");check(!polygon_contains(polygon,{.75,0.}),"polygon excludes clipped point");
    check(!polygon_contains({{0.,0.},{1.,0.},{2.,0.}},{3.,0.}),"degenerate polygon is not a full-dimensional set");
    close(polygon_goal_bound(polygon,{1.,0.},{-.25,0.}),.75,1e-14,"polygon support about center");
    check(clip_halfspace(polygon,{1.,0.},-2.).empty(),"empty polygon intersection");
    rejects([]{rectangle_polygon({1.,0.},{0.,1.});},"invalid rectangle");rejects([]{polygon_goal_bound({}, {1.,0.},{0.,0.});},"empty polygon goal bound");
    close(propagate_defect_bound(1.,2.,.5,0.),2.,1e-14,"neutral defect propagation");
    close(propagate_defect_bound(1.,0.,1.,.2),std::exp(.2),1e-14,"amplifying defect propagation");
    close(propagate_defect_bound(1.,0.,1.,-.2),std::exp(-.2),1e-14,"contracting defect propagation");
    rejects([]{propagate_defect_bound(-1.,0.,1.,0.);},"negative error allowance");
    close(quadratic_absolute_maximum(0.,4.,-4.),1.,1e-14,"interior quadratic peak");
    close(quadratic_absolute_maximum(0.,-4.,4.),1.,1e-14,"interior absolute quadratic peak");
    auto intervals=qualify_intervals<2>({2.,3.,4.,5.},{{{.5,.5}},{{1.5,.2}},{{.6,.2}}},{1.,1.});
    check(intervals.first_below[0]==2.&&intervals.persistent[0]==4.&&intervals.persistent[1]==2.&&intervals.joint==4.,"whole-interval rebound and joint time");
    auto failed=qualify_intervals<1>({0.,1.,2.},{{{.5}},{{2.}}},{1.});check(!failed.joint&&!failed.persistent[0],"whole-interval unqualified suffix");
    rejects([]{qualify_intervals<1>({0.,0.},{{{1.}}},{1.});},"invalid interval edges");
    rejects([]{qualify_intervals<1>({0.,1.},{{{1.}}},{0.});},"invalid goal tolerance");
    rejects([]{qualify_intervals<1>({0.,1.},{{{-1.}}},{1.});},"negative interval bound");
    close(norm(linear_quadratic_reconstruction(identity<2>(),Vec<2>{1.,2.},Vec<2>{2.,3.},.1,1.)-Vec<2>{2.,3.}),0.,1e-14,"quadratic reconstruction endpoint");
    rejects([]{linear_quadratic_reconstruction(identity<2>(),Vec<2>{},Vec<2>{},.1,2.);},"invalid reconstruction fraction");
    rejects([]{linear_theta_quadratic_defect_bound(identity<2>(),Vec<2>{},Vec<2>{},.1,2.,[](auto v){return norm(v);});},"invalid defect theta");

    rejects([]{Grid({0.,0.},{.5});},"duplicate times");
    rejects([]{Grid({0.,1.},{1.1});},"invalid theta");
    rejects([]{Grid({0.,1.},{std::nan("")});},"NaN theta");
    rejects([]{Grid::uniform(0.,1.,0,.5);},"zero steps");
    rejects([]{LU<2> lu(Mat<2,2>{{{1.,2.},{2.,4.}}});},"singular LU");
    // A pivot exchange is necessary in both direct and transpose systems.
    Mat<2,2> aa{{{0.,2.},{3.,4.}}};Vec<2> expected{1.,-2.};
    close(norm(LU<2>(aa).solve(mv(aa,expected))-expected),0.,1e-14,"pivoted LU");
    close(norm(LU<2>(transpose(aa)).solve(mv(transpose(aa),expected))-expected),0.,1e-14,"transpose LU");

    Grid variable({.2,.23,.31,.48,.7},{0.,.5,.8,1.});auto ww=variable.weights();
    close(ww[0],.03,1e-14,"first quadrature endpoint");
    close(ww[2],.08*.5+.17*.2,1e-14,"nonuniform interior quadrature");
    close(ww.back(),.22,1e-14,"last quadrature endpoint");
    double sum=0;for(double w:ww)sum+=w;close(sum,.5,1e-14,"quadrature sum");
    auto cp=coupled();System<2,2> sys(cp,variable);Vec<2> xi{.6,.8},d{.3,-.4};double eps=1e-5;
    auto e=sys.evaluate(xi),plus=sys.evaluate(xi+eps*d),minus=sys.evaluate(xi-eps*d);auto s=sys.tangent(xi,e.trajectory,d);auto jacobians=control_to_state_jacobians(sys,xi,e.trajectory);for(std::size_t n=0;n<s.size();++n)close(norm(mv(jacobians[n],d)-s[n]),0.,1e-13,"assembled control-state Jacobian");auto hv=sys.hessian_vector(xi,d);
    for(std::size_t n=0;n<s.size();++n)close(norm(s[n]-(.5/eps)*(plus.trajectory.state[n]-minus.trajectory.state[n])),0.,2e-8,"coupled tangent FD");
    close(dot(e.gradient,d),(plus.cost-minus.cost)/(2*eps),2e-8,"coupled gradient directional FD");
    close(norm(hv-(.5/eps)*(plus.gradient-minus.gradient)),0.,2e-7,"coupled incremental adjoint FD");
    auto h0=sys.hessian_vector(xi,{1.,0.}),h1=sys.hessian_vector(xi,{0.,1.});close(h0[1],h1[0],2e-12,"Hessian symmetry");
    for(std::size_t j=0;j<2;++j){Vec<2> unit{};unit[j]=1.;auto ep=sys.evaluate(xi+eps*unit),em=sys.evaluate(xi-eps*unit);close(e.gradient[j],(ep.cost-em.cost)/(2*eps),3e-8,"direct control cost terms");}
    double assembled=cp.penalty.value(xi);for(std::size_t n=0;n<ww.size();++n)assembled+=ww[n]*cp.running.evaluate(variable.time[n],e.trajectory.state[n],xi).value;
    assembled+=cp.terminal.evaluate(variable.time.back(),e.trajectory.state.back(),xi).value;close(assembled,e.cost,1e-14,"explicit cost assembly");
    close(sys.cost(xi,e.trajectory),e.cost,1e-14,"standalone cost API");
    close(norm(sys.adjoint(xi,e.trajectory).gradient-e.gradient),0.,1e-14,"standalone adjoint API");

    // Independent all-at-once global Jacobian, including nonlinear initial map.
    Mat<10,10> F{};Mat<10,2> rhs{};Vec<10> jz{};
    auto initial_jac=cp.initial.jacobian(xi);
    for(std::size_t i=0;i<2;++i){F[i][i]=1.;rhs[i]=initial_jac[i];}
    for(std::size_t n=0;n<4;++n){double h=variable.time[n+1]-variable.time[n],th=variable.theta[n];
        auto old=cp.dynamics.derivatives(variable.time[n],e.trajectory.state[n],xi),next=cp.dynamics.derivatives(variable.time[n+1],e.trajectory.state[n+1],xi);
        auto B=add(identity<2>(),scaled(-h*th,next.x)),C=add(identity<2>(),scaled(h*(1-th),old.x));auto K=add(scaled(h*(1-th),old.control),scaled(h*th,next.control));
        for(std::size_t i=0;i<2;++i){rhs[2*(n+1)+i]=K[i];for(std::size_t j=0;j<2;++j){F[2*(n+1)+i][2*n+j]=-C[i][j];F[2*(n+1)+i][2*(n+1)+j]=B[i][j];}}
    }
    auto direct=cp.penalty.gradient(xi);
    for(std::size_t n=0;n<5;++n){auto c=cp.running.evaluate(variable.time[n],e.trajectory.state[n],xi);direct=direct+ww[n]*c.control;for(std::size_t i=0;i<2;++i)jz[2*n+i]=ww[n]*c.x[i];}
    auto term=cp.terminal.evaluate(variable.time.back(),e.trajectory.state.back(),xi);direct=direct+term.control;for(std::size_t i=0;i<2;++i)jz[8+i]+=term.x[i];
    auto adj=LU<10>(transpose(F)).solve(jz);auto global_grad=direct+mv(transpose(rhs),adj);auto global_s=LU<10>(F).solve(mv(rhs,d));
    close(norm(global_grad-e.gradient),0.,3e-12,"global KKT versus swept gradient");
    for(std::size_t n=0;n<5;++n)for(std::size_t i=0;i<2;++i)close(global_s[2*n+i],s[n][i],3e-12,"global tangent versus sweep");

    auto first_order=cp;first_order.dynamics.directional_derivatives={};System<2,2> first(first_order,variable);
    close(first.evaluate(xi).cost,e.cost,1e-14,"second derivatives optional for first order");
    rejects([&]{first.hessian_vector(xi,d);},"missing Hessian callback rejected");
    rejects([&]{sys.tangent(xi,Trajectory<2>{},d);},"wrong trajectory size");
    rejects([]{System<2,1> s0(decay::problem(),Grid::uniform(0.,1.,2,.5),NewtonOptions{-1.,0.,20,20});},"invalid Newton options");
    rejects([]{System<2,1> s0(decay::problem(),Grid::uniform(0.,1.,2,.5));s0.forward({-1.});},"inadmissible initial state");

    // Mixed running and unweighted nodal costs, including duplicate endpoints
    // and direct state/control cross derivatives.
    auto nodal_problem=cp;auto nodal_cost=ScalarFunction<2,2>::zero();
    nodal_cost.evaluate=[](double,const Vec<2>& x,const Vec<2>& p){double r=x[0]-p[1];return ScalarValue<2,2>{.5*r*r,{r,0.},{0.,-r}};};
    nodal_cost.directional_gradient=[](double,const Vec<2>&,const Vec<2>&,const Vec<2>& s,const Vec<2>& d){double dr=s[0]-d[1];return ScalarDirection<2,2>{{dr,0.},{0.,-dr}};};
    for(auto node:{0U,2U,4U,4U})nodal_problem.nodal.emplace_back(node,nodal_cost);
    System<2,2> nodal(nodal_problem,variable);auto ne=nodal.evaluate(xi),np=nodal.evaluate(xi+eps*d),nm=nodal.evaluate(xi-eps*d);
    double addition=0.;for(auto node:{0U,2U,4U,4U})addition+=nodal_cost.evaluate(variable.time[node],e.trajectory.state[node],xi).value;
    close(ne.cost,e.cost+addition,2e-14,"nodal observations have no quadrature weights");
    close(nodal.cost(xi,ne.trajectory),ne.cost,2e-14,"mixed cost API");
    close(dot(ne.gradient,d),(np.cost-nm.cost)/(2*eps),3e-8,"nodal direct-control gradient FD");
    close(norm(nodal.hessian_vector(xi,d)-(.5/eps)*(np.gradient-nm.gradient)),0.,2e-7,"nodal mixed incremental adjoint FD");
    // Cross-time correlated residuals contribute r^T W r/2 and their shared
    // adjoint/Hessian actions must match finite differences of that same cost.
    auto correlated_problem=cp;CorrelatedObservationTerm<2,2> correlated;
    correlated.node={0,4};correlated.precision={{2.,.5},{.5,3.}};
    for(double observed:{.15,.35}){
        auto residual=ScalarFunction<2,2>::zero();
        residual.evaluate=[observed](double,const Vec<2>& x,const Vec<2>& p){
            return ScalarValue<2,2>{x[0]-.2*p[1]-observed,{1.,0.},{0.,-.2}};
        };
        residual.directional_gradient=[](double,const Vec<2>&,const Vec<2>&,const Vec<2>&,const Vec<2>&){
            return ScalarDirection<2,2>{};
        };
        correlated.residual.push_back(residual);
    }
    correlated_problem.correlated_observations.push_back(correlated);
    System<2,2> correlated_system(correlated_problem,variable);
    auto ce=correlated_system.evaluate(xi),cep=correlated_system.evaluate(xi+eps*d),cem=correlated_system.evaluate(xi-eps*d);
    double correlated_addition=0.;
    for(std::size_t i=0;i<2;++i)for(std::size_t j=0;j<2;++j){
        double ri=ce.trajectory.state[correlated.node[i]][0]-.2*xi[1]-(i==0?.15:.35);
        double rj=ce.trajectory.state[correlated.node[j]][0]-.2*xi[1]-(j==0?.15:.35);
        correlated_addition+=.5*ri*correlated.precision[i][j]*rj;
    }
    close(ce.cost,e.cost+correlated_addition,2e-14,"correlated observation quadratic cost");
    close(correlated_system.cost(xi,ce.trajectory),ce.cost,2e-14,"correlated standalone cost API");
    close(dot(ce.gradient,d),(cep.cost-cem.cost)/(2*eps),3e-8,"correlated observation adjoint gradient FD");
    close(norm(correlated_system.hessian_vector(xi,d)-(.5/eps)*(cep.gradient-cem.gradient)),0.,2e-7,
          "correlated observation incremental adjoint FD");
    auto nonsymmetric=correlated_problem;nonsymmetric.correlated_observations[0].precision[0][1]=.6;
    rejects([&]{System<2,2> invalid(nonsymmetric,variable);},"nonsymmetric observation precision rejected");
    auto indefinite=correlated_problem;indefinite.correlated_observations[0].precision={{1.,2.},{2.,1.}};
    rejects([&]{System<2,2> invalid(indefinite,variable);},"indefinite observation precision rejected");
    auto goal=sys.goal_gradient(xi,e.trajectory,4,{1.,0.});
    close(dot(goal,d),s.back()[0],2e-12,"goal adjoint versus tangent with direct controls");
    close(dot(goal,d),(plus.trajectory.state.back()[0]-minus.trajectory.state.back()[0])/(2*eps),2e-8,"goal adjoint versus FD");
    rejects([&]{sys.goal_gradient(xi,e.trajectory,5,{1.,0.});},"goal node bounds");
    close(local_goal_variance(Vec<2>{1.,2.},Mat<2,2>{{{2.,.5},{.5,3.}}}),16.,1e-14,"correlated local goal variance");
    close(box_support(Vec<2>{-2.,3.},Vec<2>{.1,.2}),.8,1e-14,"box support");
    auto tq=persistent_node_qualification({0.,1.,2.,3.},{.5,1.2,.9,.8});
    check(tq&&*tq==2.,"persistent qualification ignores early temporary crossing");
    check(!persistent_node_qualification({0.,1.},{.5,1.1}),"unqualified at horizon");
    rejects([]{persistent_node_qualification({0.,1.},{.5,std::nan("")});},"invalid qualification ratios");
    double gaussian_radius=independent_gaussian_box_radius(.03,9,.95);
    close(std::pow(std::erf(gaussian_radius/(.03*std::sqrt(2.))),9),.95,1e-13,"joint Gaussian noise-box probability");
    rejects([]{independent_gaussian_box_radius(.03,0,.95);},"Gaussian measurement count");
    auto joint_time=all_goals_qualification_time<2>({.5,2.});
    check(joint_time&&*joint_time==2.,"all goals qualify from latest persistent time");
    check(!all_goals_qualification_time<2>({.5,std::nullopt}),"one unqualified goal blocks joint qualification");
    auto initially_qualified=all_goals_qualification_time<2>({0.,0.});
    check(initially_qualified&&*initially_qualified==0.,"zero qualification time is a qualified result");

    // End-to-end decay: compare all quantities at a fixed nonstationary input.
    double analytical_min=decay::root([](double v){return decay::reference(v)[1];});
    for(double th:{.5,.65}){
        Vec<7> previous{},latest{};
        for(std::size_t steps:{20,40,80,160,320,640}){
            System<2,1> ds(decay::problem(),Grid::uniform(decay::t0,decay::tf,steps,th));Vec<1> v{2.4};auto de=ds.evaluate(v);auto st=ds.tangent(v,de.trajectory,{1.});auto lam=ds.cost_to_go_adjoint(v,de);auto H=ds.hessian_vector(v,{1.});
            auto fit=minimize(ds,Vec<1>{1.},OptimizeOptions{},[](const Vec<1>& z){return z[0]>0;});auto re=decay::reference(v[0]);double es=0,et=0,ea=0;
            for(std::size_t n=0;n<=steps;++n){auto ex=decay::exact(ds.grid().time[n]-decay::t0,v[0]);es=std::max(es,norm(de.trajectory.state[n]-ex.state));et=std::max(et,norm(st[n]-ex.tangent));ea=std::max(ea,norm(lam[n]-decay::costate(ds.grid().time[n],v[0])));}
            previous=latest;latest={es,std::abs(de.cost-re[0]),et,ea,std::abs(de.gradient[0]-re[1]),std::abs(H[0]-re[2]),std::abs(fit.control[0]-analytical_min)};
            check(de.trajectory.max_residual<5e-13,"decay forward residual");close(ds.evaluate(fit.control).gradient[0],0.,1e-11,"decay optimized stationarity");
            auto weights=ds.grid().weights();double tangent_gradient=decay::beta*(v[0]-decay::prior);
            for(std::size_t n=0;n<=steps;++n)tangent_gradient+=weights[n]*decay::kappa*de.trajectory.state[n][0]*st[n][0];
            tangent_gradient+=decay::gamma*(de.trajectory.state.back()[1]-decay::target())*st.back()[1];close(tangent_gradient,de.gradient[0],2e-12,"decay tangent-adjoint pairing");
            // Independent scalar quadratic-root reference for every theta step.
            double x=v[0],q=0.;
            for(std::size_t n=0;n<steps;++n){double h=ds.grid().time[n+1]-ds.grid().time[n],r=x-h*(1-th)*(decay::a*x+decay::b*x*x),c=1+h*th*decay::a,coef=h*th*decay::b;double next=2*r/(c+std::sqrt(c*c+4*coef*r));q+=h*((1-th)*x+th*next);x=next;}
            close(norm(de.trajectory.state.back()-Vec<2>{x,q}),0.,2e-11,"Newton versus independent quadratic roots");
            if(th==.5&&steps==160){double discrete_min=decay::root([&](double z){return ds.evaluate({z}).gradient[0];});std::vector<double> errors;for(auto it:fit.trace)errors.push_back(std::abs(it.control[0]-discrete_min));
                check(errors.size()>=5,"Newton iterations available");double p=std::log(errors[4]/errors[3])/std::log(errors[3]/errors[2]);close(p,2.,.01,"local quadratic Newton convergence");
                close(fit.trace[3].step_length,1.,1e-15,"full local Newton step");
                auto dp=ds.evaluate({v[0]+eps}),dm=ds.evaluate({v[0]-eps});close(de.gradient[0],(dp.cost-dm.cost)/(2*eps),3e-8,"decay cost gradient FD");close(H[0],(dp.gradient[0]-dm.gradient[0])/(2*eps),3e-8,"decay Hessian FD");
            }
        }
        for(std::size_t j=0;j<7;++j){double order=std::log2(previous[j]/latest[j]);close(order,th==.5?2.:1.,.09,"temporal convergence component "+std::to_string(j));}
    }
    std::cout<<checks<<" validation checks passed\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
