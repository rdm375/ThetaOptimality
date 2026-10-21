// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Richard Myers

#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace theta_opt {
template<std::size_t N> using Vec = std::array<double,N>;
template<std::size_t R,std::size_t C> using Mat = std::array<Vec<C>,R>;
template<std::size_t N> Vec<N> operator+(Vec<N> a,const Vec<N>& b){for(std::size_t i=0;i<N;++i)a[i]+=b[i];return a;}
template<std::size_t N> Vec<N> operator-(Vec<N> a,const Vec<N>& b){for(std::size_t i=0;i<N;++i)a[i]-=b[i];return a;}
template<std::size_t N> Vec<N> operator*(double s,Vec<N> a){for(auto& v:a)v*=s;return a;}
template<std::size_t N> double dot(const Vec<N>& a,const Vec<N>& b){double v=0;for(std::size_t i=0;i<N;++i)v+=a[i]*b[i];return v;}
template<std::size_t N> double norm(const Vec<N>& a){double v=0;for(auto x:a){if(!std::isfinite(x))throw std::runtime_error("Nonfinite vector");v=std::max(v,std::abs(x));}return v;}
template<std::size_t R,std::size_t C> Vec<R> mv(const Mat<R,C>& a,const Vec<C>& x){Vec<R> y{};for(std::size_t i=0;i<R;++i)y[i]=dot(a[i],x);return y;}
template<std::size_t R,std::size_t C> Mat<C,R> transpose(const Mat<R,C>& a){Mat<C,R> b{};for(std::size_t i=0;i<R;++i)for(std::size_t j=0;j<C;++j)b[j][i]=a[i][j];return b;}
template<std::size_t N> Mat<N,N> identity(){Mat<N,N> a{};for(std::size_t i=0;i<N;++i)a[i][i]=1;return a;}
template<std::size_t R,std::size_t C> Mat<R,C> scaled(double s,Mat<R,C> a){for(auto& row:a)for(auto& x:row)x*=s;return a;}
template<std::size_t R,std::size_t C> Mat<R,C> add(Mat<R,C> a,const Mat<R,C>& b){for(std::size_t i=0;i<R;++i)for(std::size_t j=0;j<C;++j)a[i][j]+=b[i][j];return a;}

// Dense partial-pivot LU, deliberately small and dependency-free.
template<std::size_t N> class LU {
    Mat<N,N> a_;std::array<std::size_t,N> piv_{};
public:
    explicit LU(Mat<N,N> a):a_(a){
        double scale=0;
        for(const auto& row:a_)for(double x:row){if(!std::isfinite(x))throw std::runtime_error("Nonfinite matrix");scale=std::max(scale,std::abs(x));}
        if(scale==0)throw std::runtime_error("Singular matrix");
        for(std::size_t k=0;k<N;++k){
            auto p=k;for(std::size_t i=k+1;i<N;++i)if(std::abs(a_[i][k])>std::abs(a_[p][k]))p=i;
            if(std::abs(a_[p][k])<=32*std::numeric_limits<double>::epsilon()*scale)throw std::runtime_error("Singular or poorly scaled matrix");
            piv_[k]=p;std::swap(a_[k],a_[p]);
            for(std::size_t i=k+1;i<N;++i){a_[i][k]/=a_[k][k];for(std::size_t j=k+1;j<N;++j)a_[i][j]-=a_[i][k]*a_[k][j];}
        }
    }
    Vec<N> solve(Vec<N> b) const {
        for(std::size_t k=0;k<N;++k)std::swap(b[k],b[piv_[k]]);
        for(std::size_t i=0;i<N;++i)for(std::size_t j=0;j<i;++j)b[i]-=a_[i][j]*b[j];
        for(std::size_t k=N;k-->0;){for(std::size_t j=k+1;j<N;++j)b[k]-=a_[k][j]*b[j];b[k]/=a_[k][k];}
        (void)norm(b);return b;
    }
};
}
