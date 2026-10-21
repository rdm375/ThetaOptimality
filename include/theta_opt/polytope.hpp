// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Richard Myers

#pragma once
#include "linalg.hpp"
#include <vector>
namespace theta_opt {
// Counterclockwise convex polygons, ordinary floating-point geometry.
using ConvexPolygon2=std::vector<Vec<2>>;
inline ConvexPolygon2 rectangle_polygon(const Vec<2>& lower,const Vec<2>& upper){
    (void)norm(lower);(void)norm(upper);for(int j=0;j<2;++j)if(!(lower[j]<upper[j]))throw std::invalid_argument("Nonempty rectangle required");return {{lower[0],lower[1]},{upper[0],lower[1]},{upper[0],upper[1]},{lower[0],upper[1]}};
}
inline ConvexPolygon2 clip_halfspace(ConvexPolygon2 polygon,const Vec<2>& normal,double bound){
    (void)norm(normal);if(!std::isfinite(bound))throw std::invalid_argument("Finite halfspace bound required");ConvexPolygon2 out;if(polygon.empty())return out;auto previous=polygon.back();(void)norm(previous);double fp=dot(normal,previous)-bound;
    for(auto current:polygon){(void)norm(current);double fc=dot(normal,current)-bound;if((fc<=0)!=(fp<=0)){double alpha=fp/(fp-fc);out.push_back(previous+alpha*(current-previous));}if(fc<=0)out.push_back(current);previous=current;fp=fc;}return out;
}
inline bool polygon_contains(const ConvexPolygon2& polygon,const Vec<2>& x,double tolerance=1e-13){
    (void)norm(x);if(!std::isfinite(tolerance)||tolerance<0)throw std::invalid_argument("Nonnegative finite geometry tolerance required");if(polygon.size()<3)return false;
    double area=0.;for(std::size_t j=0;j<polygon.size();++j){auto a=polygon[j],b=polygon[(j+1)%polygon.size()];area+=a[0]*b[1]-a[1]*b[0];}if(!(area>0))return false;
    for(std::size_t j=0;j<polygon.size();++j){auto a=polygon[j],b=polygon[(j+1)%polygon.size()];(void)norm(a);(void)norm(b);if((b[0]-a[0])*(x[1]-a[1])-(b[1]-a[1])*(x[0]-a[0])< -tolerance)return false;}return true;
}
inline double polygon_goal_bound(const ConvexPolygon2& polygon,const Vec<2>& gradient,const Vec<2>& center){
    if(polygon.empty())throw std::invalid_argument("Goal bound needs a nonempty polygon");
    (void)norm(gradient);(void)norm(center);double bound=0.;for(auto vertex:polygon){(void)norm(vertex);bound=std::max(bound,std::abs(dot(gradient,vertex-center)));}return bound;
}
}
