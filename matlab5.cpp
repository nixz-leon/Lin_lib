#include"Lin_alg.hpp"

template<typename T>
inline vec<T> fx(T x){
    vec<T> temp(3);
    temp(0) =  3.0;
    temp(1) = 2.0*x;
    temp(2) = x*x;
    return temp;
}
template <typename T>
inline vec<T> fx_1(T x, vec<T> inputs){
    vec<T> temp(2);
    temp(1) =  (inputs(1)*(T)3) + ((T)2 * x *inputs(0)) + (x*x);
    temp(0) = inputs(1);
    return temp;
}
template <typename T>
inline vec<T> fx_2(T x, vec<T> inputs){
    vec<T> temp(2);
    temp(1) =  (inputs(1)*(T)3) + ((T)2 * x *inputs(0));
    temp(0) = inputs(1);
    return temp;
}






int main(){
    vec<double>a(linear_shooting(&fx_1, &fx_2, 0.25, 0.0, 1.0, 1.0, 4.0));
    vec<double>b(linear_shooting(&fx_1, &fx_2, 0.1, 0.0, 1.0, 1.0, 4.0));
    vec<double>c(linear_shooting(&fx_1, &fx_2, 0.05, 0.0, 1.0, 1.0, 4.0));
    vec<double>a2(finite_differences(&fx, 0.25, 0.0, 1.0, 1.0, 4.0));
    vec<double>b2(finite_differences(&fx, 0.1, 0.0, 1.0, 1.0, 4.0));
    vec<double>c2(finite_differences(&fx, 0.05, 0.0, 1.0, 1.0, 4.0));
    (a2-a).printout();
    (b2-b).printout();
    (c2-c).printout();

}
