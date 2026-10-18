#include "decay_model.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

int main(int argc,char** argv){try{
    using namespace theta_opt;namespace fs=std::filesystem;
    fs::path output=argc>1?argv[1]:"results";fs::create_directories(output);
    std::ofstream csv(output/"convergence.csv");if(!csv)throw std::runtime_error("Cannot open CSV");csv<<std::setprecision(17);
    csv<<"theta,steps,h,state,cost,tangent,adjoint,gradient,hessian,initial_state\n";
    double true_min=decay::root([](double xi){return decay::reference(xi)[1];});
    for(double th:{.5,.65})for(std::size_t n:{20,40,80,160,320,640}){
        System<2,1> sys(decay::problem(),Grid::uniform(decay::t0,decay::tf,n,th));Vec<1> xi{2.4};
        auto e=sys.evaluate(xi);auto s=sys.tangent(xi,e.trajectory,{1.});auto lam=sys.cost_to_go_adjoint(xi,e);auto hess=sys.hessian_vector(xi,{1.});
        auto opt=minimize(sys,Vec<1>{1.},OptimizeOptions{},[](const Vec<1>& p){return p[0]>0;});
        double state=0,tangent=0,adjoint=0;
        for(std::size_t k=0;k<=n;++k){auto ref=decay::exact(sys.grid().time[k]-decay::t0,xi[0]);state=std::max(state,norm(e.trajectory.state[k]-ref.state));tangent=std::max(tangent,norm(s[k]-ref.tangent));adjoint=std::max(adjoint,norm(lam[k]-decay::costate(sys.grid().time[k],xi[0])));}
        auto ref=decay::reference(xi[0]);csv<<th<<','<<n<<','<<(decay::tf-decay::t0)/static_cast<double>(n)<<','<<state<<','<<std::abs(e.cost-ref[0])<<','<<tangent<<','<<adjoint<<','<<std::abs(e.gradient[0]-ref[1])<<','<<std::abs(hess[0]-ref[2])<<','<<std::abs(opt.control[0]-true_min)<<'\n';
        if(th==.5&&n==160){
            double discrete_root=decay::root([&](double p){return sys.evaluate({p}).gradient[0];});
            std::ofstream trace(output/"newton.csv");trace<<std::setprecision(17)<<"iteration,xi,cost,gradient_norm,alpha,error\n";
            for(std::size_t k=0;k<opt.trace.size();++k){const auto& it=opt.trace[k];trace<<k<<','<<it.control[0]<<','<<it.cost<<','<<it.gradient_norm<<','<<it.step_length<<','<<std::abs(it.control[0]-discrete_root)<<'\n';}
            std::ofstream states(output/"trajectory.csv");states<<std::setprecision(17)<<"time,x,q,exact_x,exact_q,tlm_x,exact_tlm_x,lambda_x,exact_lambda_x,raw_lambda_x\n";
            for(std::size_t k=0;k<=n;++k){double t=sys.grid().time[k];auto ex=decay::exact(t-decay::t0,xi[0]);states<<t<<','<<e.trajectory.state[k][0]<<','<<e.trajectory.state[k][1]<<','<<ex.state[0]<<','<<ex.state[1]<<','<<s[k][0]<<','<<ex.tangent[0]<<','<<lam[k][0]<<','<<decay::costate(t,xi[0])[0]<<','<<e.nodal_adjoint[k][0]<<'\n';}
        }
    }
    std::cout<<std::setprecision(15)<<"Analytical-objective minimizer: "<<true_min<<"\nResults: "<<output<<'\n';
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
