// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Richard Myers
#pragma once
#include "peng_robinson_model.hpp"
#include <limits>

namespace gas_inverse {
// Outward rounding of arithmetic operations. ThermoGPU's floating-point
// component coefficients are held fixed: uncertainty in EOS/data is excluded.
struct Interval {double lo{},hi{};Interval()=default;Interval(double x):lo(x),hi(x){};Interval(double a,double b):lo(a),hi(b){};};
inline double down(double v){return std::nextafter(v,-std::numeric_limits<double>::infinity());}
inline double up(double v){return std::nextafter(v,std::numeric_limits<double>::infinity());}
inline Interval operator+(Interval a,Interval b){return {down(a.lo+b.lo),up(a.hi+b.hi)};}
inline Interval operator-(Interval a,Interval b){return {down(a.lo-b.hi),up(a.hi-b.lo)};}
inline Interval operator*(Interval a,Interval b){std::array<double,4> v{a.lo*b.lo,a.lo*b.hi,a.hi*b.lo,a.hi*b.hi};return {down(*std::min_element(v.begin(),v.end())),up(*std::max_element(v.begin(),v.end()))};}
inline Interval operator/(Interval a,Interval b){if(b.lo<=0&&b.hi>=0)throw std::domain_error("Interval denominator crosses zero");return a*Interval(down(1/b.hi),up(1/b.lo));}
struct Cell {Vec<2> lo{},hi{};};
inline bool possible(const Cell& cell,const std::vector<Observation>& obs){
    Interval e(cell.lo[0],cell.hi[0]),n(cell.lo[1],cell.hi[1]);
    if(e.lo+n.lo>1)return false;
    Interval methane=Interval(1)-e-n;methane.lo=std::max(0.,methane.lo);
    std::array<Interval,3> x{methane,e,n};
    for(const auto& o:obs){Interval a,b,m;std::array<double,3> ai{};
        for(int i=0;i<3;++i)ai[i]=thermogpu::a_parameter(components[i])*thermogpu::alpha(components[i],o.state.temperature_K);
        for(int i=0;i<3;++i){b=b+x[i]*Interval(thermogpu::b_parameter(components[i]));m=m+x[i]*Interval(components[i].molar_mass_kg_per_mol);for(int j=0;j<3;++j)a=a+x[i]*x[j]*Interval(std::sqrt(ai[i]*ai[j]));}
        Interval rho(down(o.value-o.bound),up(o.value+o.bound));
        Interval v=m/rho,den=v-b;
        // Keep boxes where a denominator is unresolved, rather than excluding
        // a possible physical solution. Small accepted boxes have positive ones.
        if(den.lo<=0)continue;
        Interval attraction_den=v*(v+b)+b*(v-b);if(attraction_den.lo<=0)continue;
        Interval pressure=Interval(thermogpu::gas_constant*o.state.temperature_K)/den-a/attraction_den;
        if(pressure.hi<o.state.pressure_Pa||pressure.lo>o.state.pressure_Pa)return false;
    }return true;
}
struct OuterSet {std::vector<Cell> cells;std::size_t visited{};};
inline OuterSet outer_set(const std::vector<Observation>& obs,double resolution=.0005){
    // An explicit deterministic prior domain, not a regularization penalty.
    std::vector<Cell> queue{{{0.,0.},{.30,.20}}};OuterSet result;
    while(!queue.empty()){
        Cell c=queue.back();queue.pop_back();++result.visited;
        if(result.visited>4000000)throw std::runtime_error("Qualification subdivision budget exceeded");
        if(!possible(c,obs))continue;
        int j=(c.hi[0]-c.lo[0]>=c.hi[1]-c.lo[1])?0:1;
        if(c.hi[j]-c.lo[j]<=resolution){result.cells.push_back(c);continue;}
        double mid=.5*(c.lo[j]+c.hi[j]);Cell right=c;c.hi[j]=mid;right.lo[j]=mid;queue.push_back(right);queue.push_back(c);
    }return result;
}
inline std::array<Vec<2>,3> ranges(const OuterSet& set){
    std::array<Vec<2>,3> r{{{1.,0.},{1.,0.},{1.,0.}}};
    for(auto c:set.cells){auto methane=Interval(1)-Interval(c.lo[0],c.hi[0])-Interval(c.lo[1],c.hi[1]);r[0][0]=std::min(r[0][0],methane.lo);r[0][1]=std::max(r[0][1],methane.hi);for(int j=0;j<2;++j){r[j+1][0]=std::min(r[j+1][0],c.lo[j]);r[j+1][1]=std::max(r[j+1][1],c.hi[j]);}}return r;
}
}
