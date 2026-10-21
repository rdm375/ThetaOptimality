#include "oscillator_model.hpp"
#include "oscillator_qualification_model.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace oscillator;
namespace oq=oscillator::qualification;

namespace {
struct Checks {
    std::ofstream csv;
    std::size_t count=0;
    explicit Checks(const std::filesystem::path& path):csv(path){
        csv<<"check,error,tolerance,passed\n"<<std::setprecision(17);
    }
    void operator()(const std::string& name,double error,double tolerance){
        ++count;
        const bool passed=std::isfinite(error)&&error<=tolerance;
        csv<<name<<','<<error<<','<<tolerance<<','<<passed<<'\n';
        if(!passed)throw std::runtime_error(name+" failed: "+std::to_string(error));
    }
};

struct Case {
    const char* name;
    std::size_t steps;
    double theta;
    bool simultaneous_gaussian;
};
}

int main(int argc,char** argv){try{
    const std::filesystem::path output=argc>1?argv[1]:"results/goal_dependence";
    std::filesystem::create_directories(output);
    Checks checks(output/"checks.csv");
    std::ofstream decisions(output/"decisions.csv"),envelopes(output/"envelopes.csv"),
        states(output/"state_reconstruction.csv");
    decisions<<std::setprecision(17)
        <<"case,theta,steps,uncertainty_model,goal,units,tolerance,qualified,persistent_time_s,joint_qualified,joint_time_s,blocking_goals\n";
    envelopes<<std::setprecision(17)
        <<"case,theta,elapsed_s,interval_end_s,goal,units,tolerance,node_bound,whole_interval_bound,normalized_whole_interval_bound\n";
    states<<std::setprecision(17)
        <<"case,theta,steps,elapsed_s,true_displacement_m,analytical_fit_displacement_m,"
        <<"discrete_fit_displacement_m,displacement_bound_m,true_velocity_m_s,"
        <<"analytical_fit_velocity_m_s,discrete_fit_velocity_m_s,velocity_bound_m_s,"
        <<"true_support_force_N,analytical_fit_support_force_N,"
        <<"discrete_fit_support_force_N,support_force_bound_N\n";
    std::cout<<"GOAL-DEPENDENCE CROSS-CUTTING TEST\n"
             <<"Each case uses one fitted state and one shared observation-consistent initial-state set.\n"
             <<"Qualification is goal-specific; joint qualification requires every goal to pass.\n";

    constexpr std::array<Case,2> cases{{
        {"different_goal_times",1024,.5,true},
        {"displacement_blocker",4096,.65,false}
    }};
    for(const auto& item:cases){
        const auto model=system(item.steps,item.theta);
        const double radius=item.simultaneous_gaussian
            ?independent_gaussian_box_radius(sigma,elapsed.size(),.95):sigma;
        const auto result=oq::experiment(model,radius,
            [&](const std::string& name,double error,double tolerance){
                checks(std::string(item.name)+": "+name,error,tolerance);
            });

        std::array<bool,3> qualified{};
        std::array<double,3> times{};
        std::string blockers;
        for(std::size_t k=0;k<oq::goals.size();++k){
            qualified[k]=static_cast<bool>(result.continuous_time[k]);
            if(qualified[k])times[k]=*result.continuous_time[k];
            else {if(!blockers.empty())blockers+='+';blockers+=oq::names[k];}
        }
        if(blockers.empty())blockers="none";

        if(std::string(item.name)=="different_goal_times"){
            checks("all three goals qualify in different_goal_times",
                   std::count(qualified.begin(),qualified.end(),true)==3?0.:1.,0.);
            checks("all-goals decision exists when every goal qualifies",
                   result.joint?0.:1.,0.);
            const double distinctness=std::min({std::abs(times[0]-times[1]),
                std::abs(times[0]-times[2]),std::abs(times[1]-times[2])});
            checks("goal-specific persistent times are distinct",
                   distinctness>1e-8?0.:1.,0.);
            checks("joint time is the latest individual goal time",
                   std::abs(*result.joint-*std::max_element(times.begin(),times.end())),2e-14);
        }else{
            checks("displacement is the sole unqualified goal",
                   !qualified[0]&&qualified[1]&&qualified[2]?0.:1.,0.);
            checks("one unqualified goal blocks joint qualification",
                   result.joint?1.:0.,0.);
            checks("reported blocking goal is displacement",
                   blockers=="displacement"?0.:1.,0.);
        }

        std::cout<<"\n"<<item.name<<": theta="<<item.theta<<", N="<<item.steps
                 <<", "<<(item.simultaneous_gaussian
                    ?"95% simultaneous Gaussian noise box":"bounded measurement errors")<<'\n';
        for(std::size_t k=0;k<oq::goals.size();++k){
            decisions<<item.name<<','<<item.theta<<','<<item.steps<<','
                <<(item.simultaneous_gaussian?"simultaneous_gaussian_95":"bounded")<<','
                <<oq::names[k]<<','<<oq::units[k]<<','<<oq::tolerances[k]<<','<<qualified[k]<<',';
            if(qualified[k])decisions<<times[k];
            decisions<<','<<static_cast<bool>(result.joint)<<',';
            if(result.joint)decisions<<*result.joint;
            decisions<<','<<blockers<<'\n';
            std::cout<<"  "<<oq::names[k]<<": ";
            if(qualified[k])std::cout<<"QUALIFIED persistently from "<<times[k]<<" s\n";
            else std::cout<<"NOT QUALIFIED over the horizon\n";
        }
        std::cout<<"  joint: ";
        if(result.joint)std::cout<<"QUALIFIED from "<<*result.joint<<" s\n";
        else std::cout<<"NOT QUALIFIED; blocking goal(s): "<<blockers<<'\n';

        const auto& grid=model.grid();
        const auto reconstruction=model.forward(result.center);
        const auto analytical_center=analytical_fit();
        for(std::size_t n=0;n<grid.time.size();++n){
            const double tau=grid.time[n]-t0;
            const auto truth_state=exact(tau,truth);
            const auto analytical_state=exact(tau,analytical_center);
            const auto discrete_state=reconstruction.state[n];
            const auto support_force=[](const Vec<2>& x){return stiffness*x[0]+damping*x[1];};
            states<<item.name<<','<<item.theta<<','<<item.steps<<','<<tau<<','
                <<truth_state[0]<<','<<analytical_state[0]<<','<<discrete_state[0]<<','
                <<result.nodes[n][0]<<','<<truth_state[1]<<','<<analytical_state[1]<<','
                <<discrete_state[1]<<','<<result.nodes[n][1]<<','
                <<support_force(truth_state)<<','<<support_force(analytical_state)<<','
                <<support_force(discrete_state)<<','<<result.nodes[n][2]<<'\n';
            for(std::size_t k=0;k<oq::goals.size();++k){
                envelopes<<item.name<<','<<item.theta<<','<<tau<<',';
                if(n<grid.steps())envelopes<<grid.time[n+1]-t0;
                else envelopes<<tau;
                envelopes<<','<<oq::names[k]<<','<<oq::units[k]<<','<<oq::tolerances[k]<<','
                    <<result.nodes[n][k]<<',';
                if(n<grid.steps())envelopes<<result.intervals[n][k];
                envelopes<<',';
                if(n<grid.steps())envelopes<<result.intervals[n][k]/oq::tolerances[k];
                envelopes<<'\n';
            }
        }
    }
    std::cout<<"\nChecks passed: "<<checks.count<<". Goal dependence confirmed.\n"
             <<"Results: \""<<output.string()<<"\"\n";
    return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
