#include "../examples/decay_model.hpp"
#include <iomanip>
#include <iostream>
#include <string>
#include <theta_opt/qualification.hpp>

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
    auto e=sys.evaluate(xi),plus=sys.evaluate(xi+eps*d),minus=sys.evaluate(xi-eps*d);auto s=sys.tangent(xi,e.trajectory,d);auto hv=sys.hessian_vector(xi,d);
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
