#pragma once
#include "linalg.hpp"
#include <vector>
#include <functional>
namespace theta_opt {
// Ordinary double-precision conservative ranges, not outward-rounded intervals.
struct Range {double lower{},upper{};};
struct Box2 {Vec<2> lower{},upper{};};
inline void validate_box(const Box2& b){(void)norm(b.lower);(void)norm(b.upper);for(int j=0;j<2;++j)if(b.lower[j]>b.upper[j])throw std::invalid_argument("Ordered box endpoints required");}
inline bool contains(const Box2& b,const Vec<2>& x,double tolerance=1e-12){validate_box(b);(void)norm(x);if(!std::isfinite(tolerance)||tolerance<0)throw std::invalid_argument("Invalid box tolerance");for(int j=0;j<2;++j)if(x[j]<b.lower[j]-tolerance||x[j]>b.upper[j]+tolerance)return false;return true;}
struct Quadratic2 {double constant{};Vec<2> linear{};Mat<2,2> quadratic{};};
inline Quadratic2 operator+(Quadratic2 a,const Quadratic2& b){a.constant+=b.constant;a.linear=a.linear+b.linear;a.quadratic=add(a.quadratic,b.quadratic);return a;}
inline Quadratic2 operator*(double s,Quadratic2 a){a.constant*=s;a.linear=s*a.linear;a.quadratic=scaled(s,a.quadratic);return a;}
inline double evaluate(const Quadratic2& q,const Vec<2>& x){return q.constant+dot(q.linear,x)+dot(x,mv(q.quadratic,x));}
inline Range quadratic_range(const Quadratic2& q,const Box2& b){validate_box(b);(void)norm(q.linear);for(auto row:q.quadratic)(void)norm(row);if(!std::isfinite(q.constant))throw std::invalid_argument("Finite quadratic required");Vec<2> m{},r{};for(int j=0;j<2;++j){m[j]=.5*b.lower[j]+.5*b.upper[j];r[j]=.5*b.upper[j]-.5*b.lower[j];}auto g=q.linear+mv(add(q.quadratic,transpose(q.quadratic)),m);double radius=0.;for(int i=0;i<2;++i){radius+=std::abs(g[i])*r[i];for(int j=0;j<2;++j)radius+=std::abs(q.quadratic[i][j])*r[i]*r[j];}double center=evaluate(q,m);double slack=128*std::numeric_limits<double>::epsilon()*(1+std::abs(center)+radius);if(!std::isfinite(radius)||!std::isfinite(center))throw std::overflow_error("Quadratic range overflow");return {center-radius-slack,center+radius+slack};}
inline double absolute_bound(Range r){return std::max(std::abs(r.lower),std::abs(r.upper));}
// q(s,x)=sum s^j q_j(x), s in [0,1]. No temporal node sampling.
inline Range polynomial_quadratic_range(const std::vector<Quadratic2>& coefficients,const Box2& b){if(coefficients.empty())throw std::invalid_argument("Polynomial coefficients required");auto r=quadratic_range(coefficients[0],b);for(std::size_t j=1;j<coefficients.size();++j){auto term=quadratic_range(coefficients[j],b);r.lower+=std::min(0.,term.lower);r.upper+=std::max(0.,term.upper);}return r;}
struct BoxCover2 {std::vector<Box2> boxes;std::size_t visited{},rejected{};bool resolution_complete=true;};
// possibly_feasible must return true whenever its box may contain a feasible
// point. A budget limit retains unresolved cells; it never silently drops them.
inline BoxCover2 subdivide_feasible_boxes(const Box2& prior,const Vec<2>& target_width,const std::function<bool(const Box2&)>& possibly_feasible,std::size_t budget=500000){validate_box(prior);for(double h:target_width)if(!std::isfinite(h)||h<=0)throw std::invalid_argument("Positive finite target widths required");if(!budget||!possibly_feasible)throw std::invalid_argument("Box cover arguments");BoxCover2 result;std::vector<Box2> pending{prior};while(!pending.empty()){auto b=pending.back();pending.pop_back();++result.visited;if(!possibly_feasible(b)){++result.rejected;continue;}int axis=(b.upper[1]-b.lower[1])/target_width[1]>(b.upper[0]-b.lower[0])/target_width[0]?1:0;bool small=(b.upper[axis]-b.lower[axis])<=target_width[axis];if(small||result.visited>=budget){result.boxes.push_back(b);if(!small)result.resolution_complete=false;continue;}double mid=.5*b.lower[axis]+.5*b.upper[axis];if(mid==b.lower[axis]||mid==b.upper[axis]){result.boxes.push_back(b);result.resolution_complete=false;continue;}auto left=b,right=b;left.upper[axis]=mid;right.lower[axis]=mid;pending.push_back(right);pending.push_back(left);}return result;}
inline bool cover_contains(const BoxCover2& cover,const Vec<2>& x){for(auto b:cover.boxes)if(contains(b,x))return true;return false;}
// Union bounds must keep all retained components, including unselected optima.
inline Range cover_quadratic_range(const BoxCover2& cover,const Quadratic2& q){if(cover.boxes.empty())throw std::domain_error("Empty feasible cover");Range r{std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity()};for(auto b:cover.boxes){auto s=quadratic_range(q,b);r.lower=std::min(r.lower,s.lower);r.upper=std::max(r.upper,s.upper);}return r;}
inline Range cover_polynomial_quadratic_range(const BoxCover2& cover,const std::vector<Quadratic2>& coefficients){if(cover.boxes.empty())throw std::domain_error("Empty feasible cover");Range r{std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity()};for(auto b:cover.boxes){auto s=polynomial_quadratic_range(coefficients,b);r.lower=std::min(r.lower,s.lower);r.upper=std::max(r.upper,s.upper);}return r;}

}
