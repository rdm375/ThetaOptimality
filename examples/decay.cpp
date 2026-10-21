// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Richard Myers

#include "decay_model.hpp"

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

namespace {
struct MeshResult {
    double theta{};
    std::size_t steps{};
    double h{};
    std::array<double,7> error{};
    double fitted_initial{};
};

void write_control_convergence(const std::filesystem::path& output,
                               const std::vector<MeshResult>& rows,
                               double analytical_minimizer) {
    std::ofstream csv(output/"control_convergence.csv");
    if(!csv) throw std::runtime_error("Cannot open control convergence CSV");

    csv << std::setprecision(17)
        << "theta,steps,h,discrete_minimizer,analytical_minimizer,absolute_error,observed_order\n";

    std::cout << "Observed orders from the finest mesh pair\n"
              << "theta       state       cost     tangent     adjoint    gradient     Hessian    initial-state\n";
    for(double theta : {.5,.65}) {
        std::array<const MeshResult*,6> group{};
        std::size_t count=0;
        for(const auto& row : rows)
            if(row.theta==theta) group.at(count++)=&row;
        if(count!=group.size()) throw std::runtime_error("Unexpected mesh count");

        std::array<double,7> order{};
        for(std::size_t j=0;j<order.size();++j)
            order[j]=std::log(group[4]->error[j]/group[5]->error[j])/std::log(2.);

        std::cout << std::fixed << std::setprecision(4) << std::setw(5) << theta;
        for(double p : order) std::cout << std::setw(12) << p;
        std::cout << '\n';

        for(std::size_t i=0;i<group.size();++i) {
            const auto& row=*group[i];
            const double control_error=row.error[6];
            const double p=i ? std::log(group[i-1]->error[6]/control_error)/std::log(2.)
                             : std::numeric_limits<double>::quiet_NaN();
            csv << row.theta << ',' << row.steps << ',' << row.h << ','
                << row.fitted_initial << ',' << analytical_minimizer << ','
                << control_error << ',' << p << '\n';
        }
    }
}
}

int main(int argc,char** argv) {
    try {
        using namespace theta_opt;
        namespace fs=std::filesystem;
        const fs::path output=argc>1 ? argv[1] : "results";
        fs::create_directories(output);

        std::ofstream convergence(output/"convergence.csv");
        if(!convergence) throw std::runtime_error("Cannot open convergence CSV");
        convergence << std::setprecision(17)
                    << "theta,steps,h,state,cost,tangent,adjoint,gradient,hessian,initial_state\n";

        // Analytical minimizer: reference only. The discrete optimization below
        // uses the library's adjoint and incremental-adjoint Hessian.
        const double analytical_minimizer=
            decay::root([](double xi){ return decay::reference(xi)[1]; });
        std::vector<MeshResult> mesh_results;

        for(double theta : {.5,.65}) {
            for(std::size_t steps : {20U,40U,80U,160U,320U,640U}) {
                const Grid grid=Grid::uniform(decay::t0,decay::tf,steps,theta);
                const System<2,1> system(decay::problem(),grid);
                const Vec<1> reference_control{2.4};

                // Each object comes from the reusable C++20 library.
                const auto trajectory=system.forward(reference_control);
                const double discrete_cost=system.cost(reference_control,trajectory);
                const auto tangent=system.tangent(reference_control,trajectory,{1.});
                const auto evaluation=system.adjoint(reference_control,trajectory);
                const auto costate=system.cost_to_go_adjoint(reference_control,evaluation);
                const auto hessian=system.hessian_vector(reference_control,{1.});

                // Optimize this mesh's discrete objective, starting away from
                // the minimizer. No analytical derivatives enter this solve.
                const auto optimum=minimize(system,Vec<1>{1.},OptimizeOptions{},
                    [](const Vec<1>& trial){ return trial[0]>0.; });

                double state_error=0.,tangent_error=0.,adjoint_error=0.;
                for(std::size_t n=0;n<=steps;++n) {
                    const double time=grid.time[n];
                    const auto exact=decay::exact(time-decay::t0,reference_control[0]);
                    state_error=std::max(state_error,norm(trajectory.state[n]-exact.state));
                    tangent_error=std::max(tangent_error,norm(tangent[n]-exact.tangent));
                    adjoint_error=std::max(adjoint_error,
                        norm(costate[n]-decay::costate(time,reference_control[0])));
                }
                const auto exact_objective=decay::reference(reference_control[0]);
                const std::array<double,7> error{
                    state_error,
                    std::abs(discrete_cost-exact_objective[0]),
                    tangent_error,
                    adjoint_error,
                    std::abs(evaluation.gradient[0]-exact_objective[1]),
                    std::abs(hessian[0]-exact_objective[2]),
                    std::abs(optimum.control[0]-analytical_minimizer)};

                const double h=(decay::tf-decay::t0)/static_cast<double>(steps);
                convergence << theta << ',' << steps << ',' << h;
                for(double e : error) convergence << ',' << e;
                convergence << '\n';
                mesh_results.push_back({theta,steps,h,error,optimum.control[0]});

                // Keep the complete iteration history and pointwise state,
                // tangent, and both adjoint conventions for a representative grid.
                if(theta==.5 && steps==160) {
                    const double discrete_minimizer=decay::root(
                        [&](double xi){ return system.evaluate({xi}).gradient[0]; });
                    std::ofstream trace(output/"newton.csv");
                    if(!trace) throw std::runtime_error("Cannot open Newton CSV");
                    trace << std::setprecision(17)
                          << "iteration,xi,cost,gradient_norm,alpha,error\n";
                    for(std::size_t k=0;k<optimum.trace.size();++k) {
                        const auto& row=optimum.trace[k];
                        trace << k << ',' << row.control[0] << ',' << row.cost << ','
                              << row.gradient_norm << ',' << row.step_length << ','
                              << std::abs(row.control[0]-discrete_minimizer) << '\n';
                    }
                    std::ofstream states(output/"trajectory.csv");
                    if(!states) throw std::runtime_error("Cannot open trajectory CSV");
                    states << std::setprecision(17)
                           << "time,x,q,exact_x,exact_q,tlm_x,exact_tlm_x,"
                              "lambda_x,exact_lambda_x,raw_lambda_x\n";
                    for(std::size_t n=0;n<=steps;++n) {
                        const double time=grid.time[n];
                        const auto exact=decay::exact(time-decay::t0,reference_control[0]);
                        const double exact_lambda=decay::costate(time,reference_control[0])[0];
                        states << time << ',' << trajectory.state[n][0] << ','
                               << trajectory.state[n][1] << ',' << exact.state[0] << ','
                               << exact.state[1] << ',' << tangent[n][0] << ','
                               << exact.tangent[0] << ',' << costate[n][0] << ','
                               << exact_lambda << ',' << evaluation.nodal_adjoint[n][0] << '\n';
                    }
                }
            }
        }

        convergence.close();
        write_control_convergence(output,mesh_results,analytical_minimizer);
        std::cout << "Analytical minimizer: " << std::setprecision(15)
                  << analytical_minimizer << "\nResults: " << output << '\n';
        return 0;
    } catch(const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
