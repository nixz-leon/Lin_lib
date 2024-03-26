#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include<cstdlib>
#include <chrono>

//next step is to add a second optional arg to the template, which i will expand on in a different file for 2by2 3by3 4by4 matricies, I won't add the implimentation just yet, however 

template<class T>
struct _init_list_with_square_brackets {
    const std::initializer_list<T>& list;
    _init_list_with_square_brackets(const std::initializer_list<T>& _list): list(_list) {}
    T operator[](unsigned int index) {
        return *(list.begin() + index);
    }
    int size(){
        return list.size();
    }

};


template <typename T> class matrix{
    public:
        matrix(int n, int m);
        matrix(matrix<T> &a);
        matrix(_init_list_with_square_brackets<T> l, int rows, int col);
        inline void printout();
        inline void set_diag(T c){for(int i = 0; i < rows; i++){data[(i*cols)+i]=c;}};
        inline void row_swap(int row1, int row2);
        inline void row_add(int target, int op, T fac);
        inline matrix<T> operator=(matrix<T> a) noexcept
        {
            matrix<T> temp(a.cols,a.rows);
            for(int i = 0; i < a.rows*a.cols; i ++){
                temp.data[i] = a.data[i];
            }
            return temp;
        };
        
        /*inline matrix<T> operator=(_init_list_with_square_brackets<T> l) noexcept{
            if(l.size() != rows*cols){
                std::cerr << "list wrong size" << '\n';
                std::abort;
            }
            matrix<T> temp;
            temp.data = new T[l.size()];
            for(int i = 0; i < l.size(); i ++){
                temp.data[i] = l[i];
            }
            return temp;
        }
        */
        inline friend matrix<T> operator*(matrix<T> a, matrix<T> b){
            matrix<T> temp(a.rows, b.cols);
            if((a.cols!=b.rows)){
                std::cerr << "dim mis match" << '\n';
                return temp;
            }
            for(int i = 0; i < a.rows; i++){
                for(int j = 0; j < b.cols; j++){
                    for(int k=0; k < b.rows; k++ ){
                        temp(i,j) += a(i,k)*b(k,j);
                    }
                }
            }
            return temp;
        };
        inline matrix operator *=(T c){
            for(int i = 0; i < this->rows;i++){
                for(int j = 0; j < this->cols;j++){
                    data[(i*cols)+j] *= c;
                }
            }
            return *this;
        };
        inline matrix operator/=(T c){
            for(int i = 0; i < this->rows;i++){
                for(int j = 0; j < this->cols;j++){
                    data[(i*cols)+j] /= c;
                }
            }
            return *this;
        }
        inline friend matrix<T> operator +(matrix<T> a, matrix<T> b){
            matrix<T> temp(a.rows, b.cols);
            if((a.cols!=b.cols)||(a.rows!=b.rows)){
                std::cerr << "dim mis match" << '\n';
                return temp;
            }
            for(int i= 0; i < a.rows; i++){
                for(int j=0; j < a.cols; j++){
                    temp(i,j) += a(i,j) + b(i,j);
                }
            }
            return temp;
        };
        inline T& operator()(int a, int b) {
            if(a > rows-1 || b > cols-1){
                std::cerr << "dim do not match" << '\n';
                std::cout << "a\n";
                return data[0];
            }
            return data[(a*cols)+b];
        };
        inline const T&operator()(int a, int b) const {
            if(a > rows-1 || b > cols-1){
                std::cerr << "dim do not match" << '\n';
                std::cout << "b\n";
                return data[0];
            }
            return data[(a*cols)+b];
        };
        ~matrix(){delete[] data;rows =0; cols =0;}
        int rows;
        int cols;
    private:
        T *data;

};
template <typename T>
inline matrix<T>::matrix(_init_list_with_square_brackets<T> l, int r, int c){
    data = new T[r*c]();
    rows = r;
    cols = c;
    for(int i =0; i < r*c; i++){
        data[i] = l[i];
    }
}
template <typename T>
inline matrix<T>::matrix(int n, int m)
{
    //data = (T*)calloc(n*m,sizeof(T));
    data = new T[n*m]();
    rows = n;
    cols = m;
};
template <typename T>
inline matrix<T>::matrix(matrix<T> &a){
    data = new T[a.cols + a.rows]();
    rows = a.rows;
    cols = a.cols;
    for(int i =0; i < rows*cols; i++){
        data[i] = a.data[i];
    }
}

template<class T> 
inline void matrix<T>::printout(){
    std::cout << '\n';
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols; j++){
           std::cout << data[(cols*i)+j] << ' ';
        }
       std::cout << '\n';
    }

};
template <typename T>
inline void matrix<T>::row_swap(int row1, int row2)
{
    if((row1>rows-1)||(row2>rows-1)){
        std::cerr << "out of bounds row" << " row1: " << row1 << " row2: "<< row2 << '\n';
        return;
    }else if(row1==row2){
        return;
    }
    for(int i =0; i < cols; i++){
        T temp = data[(row1*cols)+i];
        data[(row1*cols)+i] = data [(row2*cols)+i];
        data[(row2*cols)+i] = temp;
    } 
};
template <typename T>
inline void matrix<T>::row_add(int target, int op, T fac)
{
    if((target>rows-1)||(op>rows-1)){
        std::cerr << "out of bounds row" << " row1:" << target << " row2:"<<  op<<'\n';
        return;
    }else if((target==op)){
        return;
    }
    for(int i =0; i < cols; i++){
        data[(target*cols)+i] += (fac)*data[(op*cols)+i];
    } 

};


template <typename T> class vec{
    public:
        vec(int n);
        vec(vec<T> &a);
        vec(_init_list_with_square_brackets<T> l);
        inline void printout();
        inline void row_add(int target, int op, T fac){data[target] += (data[op]*fac);};
        inline void row_swap(int target, int op){T temp = data[target]; data[target]= data[op];data[op]=temp;};
        T* return_data(){return data;};
        void set(int a, T d){data[a]=d;};
        inline vec<T> operator=(vec<T> a) noexcept
        {
            for(int i = 0; i < a.size; i ++){
                this->data[i] = a.data[i];
            }
            return *this;
        };
        inline vec<T> operator=(_init_list_with_square_brackets<T> l) noexcept{
            vec<T> temp(l.size());
            for(int i = 0; i < l.size(); i ++){
                temp(i) = l[i];
            }
            return temp;
        };
        inline friend vec<T> operator+(vec<T> a, vec<T> b){
            vec<T> temp(a.size);
            if(a.size!=b.size){
                std::cerr << "dim do not match" << '\n';
                std::cout << "c\n";
                return temp;
            }
            for(int i = 0; i < a.size;i++){
                temp(i) = (a(i)+b(i));
            }
            return temp;
        };
        inline friend vec<T> operator-(vec<T> a, vec<T> b){
            vec<T> temp(a.size);
            if(a.size!=b.size){
                std::cerr << "dim do not match" << '\n';
                std::cout << "d\n";
                return temp;
            }
            for(int i = 0; i < a.size;i++){
                temp(i) = (a(i)-b(i));
            }
            return temp;
        };
        inline vec<T> operator+=(T c){for(int i = 0; i < this->size;i++){data[i]+=c;}return *this;};
        inline vec<T> operator+=(vec<T> a){for(int i =0; i < this->size;i++){this->data[i] += a(i);}return *this;};
        inline vec<T> operator-=(T c){for(int i = 0; i < this->size;i++){data[i]-=c;}return *this;};
        inline vec<T> operator-=(vec<T> a){for(int i =0; i < this->size;i++){this->data[i] -= a(i);}return *this;};
        inline friend vec<T> operator*(vec<T> a,T c){
            vec<T> temp(a.size);
            for(int i =0; i < a.size; i++)
                temp(i) = (a(i) * c);
            return temp;
        };
        inline friend vec<T> operator*(matrix<T> a, vec<T> b){
            vec<T> temp(b.size);
            if(b.size != a.cols){
                std::cerr << "Dim do not match" << '\n';
                return temp;
            }
            T sum  = 0;
            for(int i = 0; i<a.rows; i++){
                sum = 0;
                for(int j = 0; j < a.cols; j++){
                    sum += a(i,j) * b(j);
                }
                temp(i) = sum;
            }
            return temp;
        };
        inline vec<T> operator*=(T scalar){
            for(int i =0; i < this->size; i++){
                this->data[i] *= scalar;
            }
            return *this;
        };
        // x = A/b; system's solver;
        inline T& operator()(int a) {
            if(a > size-1){
                std::cerr << "Tried to access elm " << a << ", But size is " << size-1 << '\n';
                exit(0);
                return data[0];
            }else{
                return data[a];    
            }
            
        };
        inline const T&operator()(int a) const {
            if(a > size-1){
                std::cerr << "Tried to access elm " << a << ", But size is " << size-1 << '\n';
                exit(0);
                return data[0];
            }else{
                return data[a];
            }
        };
        inline const T back(){
            return data[size-1];
        }
        int size;
        ~vec(){delete[] data;size=0;}   
    private:
        T *data;
};

template <typename T>
inline vec<T>::vec(int n)
{
    data = new T[n]();
    size = n;
}
template <typename T>
inline vec<T>::vec(vec<T> &a){
    data = new T[a.size]();
    //data = a.data;
    size = a.size;
    for(int i =0; i < size;i++){
        data[i] = a.data[i];
    }
    
}
template <typename T>
inline vec<T>::vec(_init_list_with_square_brackets<T> l) {
    size = l.size();
    data = new T[size]();
    for(int i =0; i < size; i++){
        data[i] = l[i];
    }
};
template <typename T>
inline void vec<T>::printout(){
    for(int i=0; i < size; i++){
       std::cout << data[i] << '\n';
    }
   std::cout << '\n';
};


template <typename T>
constexpr T norm1(vec<T> a){
    T sum = 0;
    for(int i =0; i < a.size; i++){
        sum += a(i);
    }
    return sum;
};

template <typename T> 
inline void gausian_reduction(matrix<T> &a){
    int biggest_row;
    for(int i =0; i < a.cols;i++){
        biggest_row = i;
        for(int j = i+1; j < a.rows;j++){
            if(a(j,i) > a(biggest_row,i)){
                biggest_row = j;
            }
        }
        a.row_swap(biggest_row, i);
    }
    //gaussian elimination to ref
    for(int i = 0; i < a.rows;i++ ){
        for(int j = i+1; j <a.cols;j++){
            T m = (a(j,i)/a(i,i))*(-1.f);
            a.row_add(j,i,m);
        }
    }
}

template <typename T> 
inline vec<T> lin_solver(matrix<T> a, vec<T> b){
    matrix<T> tempa(a.rows, a.cols);
    vec<T> tempb(b.size), x(b.size); 
    //Two if else statements for size matching
    if(a.rows!=a.cols){
        std::cerr << "matrix not square" << '\n';
        return x;
    }
    if(a.rows!=b.size){
        std::cerr << "dim of A and b do not match" << '\n';
        return x;
    }
    //partial pivoting
    //iterate through each coloum
    tempa = a;
    tempb = b;
    int biggest_row;
    for(int i =0; i < tempa.cols;i++){
        biggest_row = i;
        for(int j = i+1; j < tempa.rows;j++){
            if(tempa(j,i) > tempa(biggest_row,i)){
                biggest_row = j;
            }
        }
        tempa.row_swap(biggest_row, i);
        tempb.row_swap(biggest_row, i);
    }
    //gaussian elimination to ref
    for(int i = 0; i < tempa.rows;i++ ){
        for(int j = i+1; j <tempa.cols;j++){
            T m = (tempa(j,i)/tempa(i,i))*(-1.f);
            tempa.row_add(j,i,m);
            tempb.row_add(j,i,m);
        }
    }
    //back sub
    T sum = 0;
    for(int i = tempa.rows-1; i>-1; i--){
        for(int j = tempa.rows-1; j>i-1; j--){
            if(i == tempa.rows-1){
                sum += 0;
            }else{
                sum += tempa(i,j)*x(j);
            }
        }
        x(i) = (tempb(i)-sum)/(tempa(i,i));
        sum = 0;
    }
    return x;
};

/*
template <typename T>
inline void lu_decomp(matrix<T> a, matrix<T> &l, matrix<T> &u, matrix<T> &p){
    //assumptions about inputs l u and p are zeros
    p.set_diag(1);
    u = a;
    if(a.rows!=a.cols){
        std::cerr << "matrix not square" << '\n';
        return;
    }
    //partial pivoting for numerical stability
    int biggest_row;
    for(int i =0; i < u.cols;i++){
        biggest_row = i;
        for(int j = i+1; j < u.rows;j++){
            if(u(j,i) > u(biggest_row,i)){
                biggest_row = j;
            }
        }
        u.row_swap(biggest_row, i);
        p.row_swap(biggest_row, i);
    }
    for(int i = 0; i < u.rows;i++ ){
        
        for(int j = i+1; j <u.cols;j++){
            T m = (u(j,i)/u(i,i))*(-1);
            u.row_add(j,i,m);
            l(j,i) = m*(-1);
        }
    }
}*/

template <typename T>
inline T determinant(matrix<T> a){
    T det = 1;
    if(a.rows!=a.cols){
        std::cerr << "matrix not square" << '\n';
        return det;
    }
    matrix<T> u(a.rows,a.cols);
    u = a;
    int biggest_row =0;
    int swap_count=0;
    for(int i =0; i < u.cols;i++){
        biggest_row = i;
        for(int j = i+1; j < u.rows;j++){
            if(u(j,i) > u(biggest_row,i)){
                biggest_row = j;
            }
        }
        u.row_swap(biggest_row, i);
        if(biggest_row != i){
            swap_count++;
        }
    }
    for(int i = 0; i < u.rows;i++ ){
        for(int j = i+1; j <u.cols;j++){
            T m = (u(j,i)/u(i,i))*(-1);
            u.row_add(j,i,m);
        }
    }
    for(int i = 0; i< u.rows;i++){
        det *= u(i,i);
        if(det == 0){
            return det;
        }
    }
    if(swap_count%2 == 1){
        det *= -1;
    }
    return det;
};

template<typename T >
constexpr T norm2(vec<T> a){
    return sqrt(dot(a,a));
}
template<typename T>
inline T dot(vec<T>a, vec<T>b){
    T sum = 0;
    if(a.size != b.size){
        std::cerr << "dim do not match, dotProduct" << '\n';
        return sum;
    }
    for(int i = 0; i < a.size; i ++){
        sum +=(a(i))*(b(i));
    }
    return sum;
}


template <typename T>
inline matrix<T> alt_house_holder(vec<T> a, int c){
    matrix<T> A(a.size, a.size);
    vec<T> tempv(a.size);
    T alpha = 0;
    T beta = 0;
    tempv = a;
    for(int i =0; i < c; i++){
        tempv(i) = 0;
    }
    for(int i = c; i < A.cols; i++){
        for(int j = i; j < A.rows;j++){
            T entry = tempv(i) * tempv(j);
            if(i == j){
                A(i,i) = entry;
                alpha += entry;
            }else{
                A(i,j) = entry;
                A(j,i) = entry;
            }
        }
    }
    //now we have the uncorrected uuT and alpha
    tempv *= sqrt(alpha); // create correction vector
    tempv(c) *=2;
    //apply corrective matrix
    for(int i = c; i < tempv.size; i++){
        if(i == c){A(i,i)+= (alpha + tempv(i));}
        else{A(i,c)+=tempv(i);A(c,i)+=tempv(i);}
    }
    beta = 2*alpha + tempv(c);
    //A is corrected
    //Apply normalization factor of 1/alpha
    //Apply factor -2
    if(beta != 0){
        A/= beta;    
    }  
    A*=-2;
    for(int i =0; i<A.cols; i++){
        A(i,i) += 1;
    }
    return A;
}


template <typename T>
inline vec<T> get_col(matrix<T> a, int c){
    vec<T> temp(a.rows);
    for(int i =0; i < temp.size; i++){
        temp(i) = a(i,c);
    }
    return temp;
};

template <typename T>
inline void qr(matrix<T> a, matrix<T> &q, matrix<T> &r){
    matrix<T> h(a.rows,a.rows);
    h.set_diag(1);
    matrix<T> temp(a.rows, a.rows);
    r = a;
    for(int c = 0; c < a.cols; c++){
        vec<T> tempv = get_col(r,c);
        temp = alt_house_holder(tempv,c);
        r = temp * r;
        h = h*temp;
    }
    q = h;
};


template <typename T>
//all of these are of size x.size
inline vec<T> interpolate(vec<T> x, vec<T> fx)
{
    vec<T> poly(x.size);
    vec<T> temp(x.size);
    T prod;
    int i, j, k;
    for (i = 0; i < x.size; i++)
        poly(i) = 0;
    for (i = 0; i < x.size; i++) {
        prod = 1;
        for (j = 0; j < x.size; j++) {
            temp(j) = 0;
        }
        for (j = 0; j < x.size; j++) {
            if (i == j)
                continue;
            prod *= (x(i) - x(j));
        }
        prod = fx(i) / prod;
        temp(0) = prod;
        for (j = 0; j < x.size; j++) {
            if (i == j)
                continue;
            for (k = x.size - 1; k > 0; k--) {
                temp(k) += temp(k - 1);
                temp(k - 1) *= (-x(j));
            }
        }
        for (j = 0; j < x.size; j++) {
            poly(j) += temp(j);
        }
    }
    return poly;
};

// store each iteration in a vector; and then display the vector 
template <typename T>
inline T euler(T (*f)(T, T), T h, T a, T b, T inital){
    int its = int((b-a)/h);
    //vec<T> val(its+1);
    T val;
    val = inital;
    T t = a;
    for(int i =1; i < its; i++){
        val = (val + h*f(t, val));
        t+=h;
    }
    return val;
};

template <typename T>
inline T change_h(T diff, T h, T tol){
    T temp1 = 0.9*h*std::pow(tol/diff, 1./6.0);
    T temp2 = 10*h;
    return temp1<temp2 ? temp1 : temp2;
};

//i will want to do the h calc In this function
template <typename T>
inline T rkdps(T (*f)(T,T), T &h, T y, T t ,matrix<T> &a, vec<T> &b, vec<T> &b2 ,vec<T> &c, T &diff,T tol){
    vec<T> temp(2);
    vec<T> k(7);
    k(0) = h * f(t,y);
    k(1) = h * f(t + c(1)*h, y + h*a(1,0)*k(0));
    k(2) = h * f(t + c(2)*h, y + h*a(2,0)*k(0) + h*a(2,1)*k(1));
    k(3) = h * f(t + c(3)*h, y + h*a(3,0)*k(0) + h*a(3,1)*k(1) + h*a(3,2)*k(2));
    k(4) = h * f(t + c(4)*h, y + h*a(4,0)*k(0) + h*a(4,1)*k(1) + h*a(4,2)*k(2) + h*a(4,3)*k(3));
    k(5) = h * f(t + c(5)*h, y + h*a(5,0)*k(0) + h*a(5,1)*k(1) + h*a(5,2)*k(2) + h*a(5,3)*k(3) + h*a(5,4)*k(4));
    k(6) = h * f(t + c(6)*h, y + h*a(6,0)*k(0) + h*a(6,1)*k(1) + h*a(6,2)*k(2) + h*a(6,3)*k(3) + h*a(6,4)*k(4) + a(6,5)*k(5));
    temp(0) = y +  b(0)*k(0) + b(1)*k(1) + b(2)*k(2) + b(3)*k(3) + b(4)*k(4) + b(5)*k(5);
    temp(1) = y +  b2(0)*k(0) + b2(1)*k(1) + b2(2)*k(2) + b2(3)*k(3) + b2(4)*k(4) + b2(5)*k(5) + b2(6)*k(6);
    diff = norm2(temp);
    //diff = abs(temp(1)-temp(0));
    if(diff > h*tol){
        h=change_h(diff, h, tol);
        std::cout << "hit\n"; 
    }
    return temp(0);
};



//One of the matrices or vectors is either not doing a proper deep copy, or 
template <typename T>
inline T rkdp(T (*f)(T,T), T h, T a, T b, T inital, T tol){
    matrix<T>aii(7,7);
    vec<T>bi(7),bi2(7),ci(7);
    aii(1,0) = 1./5.;
    aii(2,0) = 3./40.;         aii(2,1) = 9./40.;
    aii(3,0) = 44./45.;        aii(3,1) = -56./15.;        aii(3,2) = 32./9.;
    aii(4,0) = 19372.0/6561.0; aii(4,1) = -25360./2187.;   aii(4,2) = 64448.0/6561.0; aii(4,3) = -212./729.;
    aii(5,0) = 9017./3168.;    aii(5,1) = -355./33.;       aii(5,2) = 46732./5247.;   aii(5,3) = 49./176.;    aii(5,4) = -5103./18656.;
    aii(6,0) = 35./384.;       aii(6,1) = 0.;              aii(6,2) = 500./1113.;     aii(6,3) = 125.0/192.;  aii(6,4) = -2187.0/6784.0;  aii(6,5) = 11./84.;
    bi(0) = 35./384.;          bi(1)=0.;                   bi(2)=500./1113.;          bi(3) = 125./192.;      bi(4)=-2187.0/6784.0;       bi(5)=11./84.;       bi(6)=0.;
    bi2(0) = 5179.0/57600.0;   bi2(1)=0.0;                 bi2(2)=7571.0/16695.0;     bi2(3) = 393.0/640.0;   bi2(4)=-92097.0/339200.0;   bi2(5)=187./2100.;   bi2(6)=1./40.;
    ci(0) = 0.0;               ci(1)=0.2;                  ci(2)=0.3;                 ci(3) = 0.8;            ci(4)=8.0/9.0;              ci(5)=1.0;           ci(6)=1.0;
    T y = inital; 
    T t = a;
    T diff;
    while(t < b){
        //std::cout << h << '\n';
        y = rkdps(f, h, y, t, aii,bi,bi2,ci,diff, tol);
        t+=h;
    }
    return y;
};



template <typename T>
inline T rk4(T (*f)(T,T), T h, T a, T b, T inital){
    int its = int((b-a)/h);
    T k0,k1,k2,k3,val;
    T t = a;
    T h2 = h*0.5;
    val = inital;
    for(int i =0; i < its; i++){
        k0 =  f(t, val);
        k1 = f(t+h2, val+h2*k0);
        k2 = f(t+h2, val+h2*k1);
        k3 = f(t+h, val+h*k2);
        val = val + h*(1.0/6.0)*(k0 + k1*2.f + k2*2.f + k3);
        t += h;
    }
    return val;
};
template <typename T>
inline void rk4(T (*f)(T,T), T h, T a, T b, T inital, vec<T> &val){
    int its = int((b-a)/h);
    //vec<T> vals(its);
    T k0,k1,k2,k3;
    T t = a;
    T h2 = h*0.5;
    val(0) = inital;
    for(int i =0; i < its; i++){
        k0 =  f(t, val(i));
        k1 = f(t+h2, val(i)+h2*k0);
        k2 = f(t+h2, val(i)+h2*k1);
        k3 = f(t+h, val(i)+h*k2);
        val(i+1) = val(i) + h*(1.0/6.0)*(k0 + k1*2.f + k2*2.f + k3);
        t += h;
    }
};

template <typename T>
inline T rk4(T (*f)(T,T), T h, T a, T val){
    T k0,k1,k2,k3;
    T h2 = h*0.5f;
    k1 = f(a+h2, val+(h2*k0));
    k2 = f(a+h2, val+(h2*k1));
    k3 = f(a+h, val+(h*k2));
    return val + h*(1.0/6.0)*(k0 + k2*2.f + k3*2.f + k3);
};

template <typename T>
inline vec<T> rk4(vec<T> (*f)(vec<T>,vec<T>), T h, vec<T> a, vec<T> val){
    vec<T> k0(a.size),k1(a.size),k2(a.size),k3(a.size);
    T h2 = h*0.5;
    k0 = f(a, val);
    k1 = f(a+h2, val+(h2*k0));
    k2 = f(a+h2, val+(h2*k1));
    k3 = f(a+h, val+(h*k2));
    return val + h*(1.0/6.0)*(k0 + k2*2.f + k3*2.f + k3);
};


template <typename T>
inline vec<T> AB4(T (*f)(T,T), T h, T a, T b, T inital){
    int its = int((b-a)/h);
    vec<T> val(its+1);
    vec<T> func_vals(4);
    val(0) = inital;
    rk4(f, h, a, a+3*h, inital,val);
    T t=a;
    func_vals(0) = f(t,val(0));t+=h;
    func_vals(1) = f(t,val(1));t+=h;
    func_vals(2) = f(t,val(2));t+=h;
    func_vals(3) = f(t,val(3));
    int i = 4;
    while(t<b){
        val(i) = val(i-1) + (h)*(55.*func_vals(3) - 59.*func_vals(2) + 37.*func_vals(1) -9.*func_vals(0))/(24.);
        t+=h;
        func_vals(0) = func_vals(1);func_vals(1) = func_vals(2);func_vals(2) = func_vals(3);
        func_vals(3) = f(t,val(i));
        i++;
    }
    return val;
};

template <typename T>
inline vec<T> pece(T (*f)(T,T), T h, T a, T b, T inital){
    int its = int((b-a)/h);
    vec<T> val(its+1);
    vec<T> func_vals(4);
    val(0) = inital;
    rk4(f, h, a, a+3*h, inital,val);
    T t=a;
    func_vals(0) = f(t,val(0));t+=h;
    func_vals(1) = f(t,val(1));t+=h;
    func_vals(2) = f(t,val(2));t+=h;
    func_vals(3) = f(t,val(3));
    int i = 4;
    T temp;
    while(t<b){
        val(i) = val(i-1) + (h)*(55.*func_vals(3) - 59.*func_vals(2) + 37.*func_vals(1) -9.*func_vals(0))/(24.);
        t+=h;
        temp = f(t,val(i));
        val(i) = val(i-1) + (h)*(251.0*temp + 646.0*func_vals(3) -264.0*func_vals(2)+106.0*func_vals(1)-19.0*func_vals(0))/(720.0);
        func_vals(0) = func_vals(1);func_vals(1) = func_vals(2);func_vals(2) = func_vals(3);
        func_vals(3) = f(t,val(i));
        i++;
    }
    return val;
};

template <typename T>
inline vec<T> AB4(T (*f)(T,T), T h, T a, T b, vec<T> inital){
    int its = int((b-a)/h);
    vec<T> val(its+1);
    vec<T> func_vals(4);
    val(0)= inital(0);
    T t=a;
    func_vals(0) = f(t,inital(0));t+=h;
    func_vals(1) = f(t,inital(1));t+=h;
    func_vals(2) = f(t,inital(2));t+=h;
    func_vals(3) = f(t,inital(3));
    int i = 4;
    while(t<b){
        val(i) = val(i-1) + (h)*(55.*func_vals(3) - 59.*func_vals(2) + 37.*func_vals(1) -9.*func_vals(0))/(24.);
        t+=h;
        func_vals(0) = func_vals(1);func_vals(1) = func_vals(2);func_vals(2) = func_vals(3);
        func_vals(3) = f(t,val(i));
        i++;
    }
    return val;
};
template <typename T>
inline vec<T> pece(T (*f)(T,T), T h, T a, T b, vec<T> inital){
    int its = int((b-a)/h);
    vec<T> val(its+1);
    vec<T> func_vals(4);
    val(0)= inital(0);
    T t=a;
    func_vals(0) = f(t,inital(0));t+=h;
    func_vals(1) = f(t,inital(1));t+=h;
    func_vals(2) = f(t,inital(2));t+=h;
    func_vals(3) = f(t,inital(3));
    int i = 4;
    T temp;
    while(t<b){
        val(i) = val(i-1) + (h)*(55.*func_vals(3) - 59.*func_vals(2) + 37.*func_vals(1) -9.*func_vals(0))/(24.);
        t+=h;
        temp = f(t,val(i));
        val(i) = val(i-1) + (h)*(251.0*temp + 646.0*func_vals(3) -264.0*func_vals(2)+106.0*func_vals(1)-19.0*func_vals(0))/(720.0);
        func_vals(0) = func_vals(1);func_vals(1) = func_vals(2);func_vals(2) = func_vals(3);
        func_vals(3) = f(t,val(i));
        i++;
    }
    return val;
};

// for the finite difference code, i will basically need to generate a matrix and then use my linear solver to get the results
//to generate this matrix, the function that is being used as an input must return a vector of size 3, that corresponds:
//p(xi), q(xi), r(xi)
template <typename T>
inline vec<T> finite_differences(vec<T> (*f)(T), T h, T a, T b, T y0, T yn){
    int its = (int)(((b-a)/h)+0.5f);//0.5f for rounding this caused an issue with the code from last week
    T x = a+h;
    T h2 = h*0.5;
    vec<T> b_vec(its);
    b_vec(its-1) = yn;
    b_vec(0) = y0;
    matrix<T>  A_mat(its,its);//default constructor allocates all the elements to be 0
    A_mat(0,0) = 1;A_mat(its-1,its-1) = 1;//setting the corners to 1;
    vec<T> temp(3);//this is a place holder so I don't allocate memory every iteration of the loop, and only need to do one function call;
    for(int i =1; i < its-2; i++){//bounded so it doesn't disturb the preset elements
        temp = f(x);
        b_vec(i) = h*h*temp(2);
        A_mat(i,i) = 1;
        A_mat(i,i-1) = h2*(temp(0)+1)/(-2.-temp(1));
        A_mat(i,i+1) = h2*(temp(0)-1)/(2.+temp(1)); 
        x+=h;
    }

    return lin_solver(A_mat,b_vec);

};




//matrix transpose
//projection
//orthagonalize vector
//1,2,inf norms vect
//1 and inf norm for mat
//normalize
//eigen value 
//cubic spline 
//numeric diff 
//numeric int 
//least squares