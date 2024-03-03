#pragma once
#include <iostream>
#include <cmath>
#include <string>
#include<cstdlib>
using namespace std;

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
        matrix(int n, int m);//will inizalize to all 0s
        matrix(int n, int m, _init_list_with_square_brackets<T> l);
        matrix(matrix<T> &a);   
        inline void printout();
        inline void set_diag(T c){for(int i = 0; i < rows; i++){data[(i*cols)+i]=c;}};
        inline void row_swap(int row1, int row2);
        inline void row_add(int target, int op, T fac);
        inline matrix<T> operator=(matrix<T> a) noexcept
        {
            for(int i = 0; i < a.rows*a.cols; i ++){
                this->data[i] = a.data[i];
            }
            return *this;
        };
        inline matrix<T> operator=(_init_list_with_square_brackets<T> l) noexcept{
            if(l.size() != rows*cols){
                cerr << "list wrong size" << endl;
                return *this;
            }
            for(int i = 0; i < l.size(); i ++){
                data[i] = l[i];
            }
            return *this;
        }
        inline friend matrix<T> operator*(matrix<T> a, matrix<T> b){
            matrix<T> temp(a.rows, b.cols);
            if((a.cols!=b.rows)){
                cerr << "dim mis match" << endl;
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
                cerr << "dim mis match" << endl;
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
                cerr << "dim do not match" << endl;
                return data[0];
            }
            return data[(a*cols)+b];
        };
        inline const T&operator()(int a, int b) const {
            if(a > rows-1 || b > cols-1){
                cerr << "dim do not match" << endl;
                return data[0];
            }
            return data[(a*cols)+b];
        };
        ~matrix(){
            free(data);
            rows =0;
            cols = 0;
        }
        int rows;
        int cols;
    private:
        T *data;

};
template <typename T>
inline matrix<T>::matrix(int n, int m)
{
    data = (T*)calloc(n*m,sizeof(T));
    rows = n;
    cols = m;
};
template <typename T>
inline matrix<T>::matrix(matrix<T> &a){
    rows = a.rows;
    cols = a.cols;
    data = (T*)calloc(rows * cols, sizeof(T));
    for(int i =0; i < a.rows*a.cols; i++){
        data[i] = a.data[i];
    }

};
template<class T> 
inline void matrix<T>::printout(){
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols; j++){
            cout << data[(cols*i)+j] << ' ';
        }
        cout << endl;
    }
    cout << endl;
};
template <typename T>
inline void matrix<T>::row_swap(int row1, int row2)
{
    if((row1>rows-1)||(row2>rows-1)){
        cerr << "out of bounds row" << " row1: " << row1 << " row2: "<< row2 <<endl;
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
        cerr << "out of bounds row" << " row1:" << target << " row2:"<<  op<<endl;
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
        inline void entry(T *v, int size);
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
            if(l.size() != size){
                cerr << "list wrong size" << endl;
                return *this;
            }
            for(int i = 0; i < l.size(); i ++){
                data[i] = l[i];
            }
            return *this;
        }
        inline friend vec<T> operator+(vec<T> a, vec<T> b){
            vec<T> temp(a.size);
            if(a.size!=b.size){
                cerr << "dim do not match" << endl;
                return temp;
            }
            for(int i = 0; i < a.size;i++){
                temp(i) = (a(i)+b(i));
            }
            return temp;
        }
        inline friend vec<T> operator-(vec<T> a, vec<T> b){
            vec<T> temp(a.size);
            if(a.size!=b.size){
                cerr << "dim do not match" << endl;
                return temp;
            }
            for(int i = 0; i < a.size;i++){
                temp(i) = (a(i)-b(i));
            }
            return temp;
        }
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
                cerr << "Dim do not match" << endl;
                return temp;
            }
            T sum  = 0;
            for(int i = 0; i<a.rows; i++){
                sum = 0;
                for(int j = 0; j < a.cols; j++){
                    sum += a(i,j) * b(j);
                    cout << "a("<< i << ',' << j << "):" << a(i,j) << " * b(" << j << "):" << b(j) << " = " <<(a(i,j)*b(j)) <<endl;
                }
                temp(i) = sum;
                cout << "sum:" <<sum << endl;
            }
            return temp;
        }
        inline vec<T> operator*=(T scalar){
            for(int i =0; i < this->size; i++){
                this->data[i] *= scalar;
            }
            return *this;
        }
        // x = A/b; system's solver;
        inline T& operator()(int a) {
            if(a > size-1){
                cerr << "dim do not match" << endl;
                return data[0];
            }else{
                return data[a];    
            }
            
        };
        inline const T&operator()(int a) const {
            if(a > size-1){
                cerr << "dim do not match" << endl;
                return data[0];
            }else{
                return data[a];
            }
        };
        ~vec(){
            size = 0;
            free(data);
        }
        int size;   
    private:
        T *data;
};

template <typename T>
inline vec<T>::vec(int n)
{
    data = (T*)calloc(n,sizeof(T));
    size = n;
}
template <typename T>
inline vec<T>::vec(vec<T> &a){
    size = a.size;
    data = (T*)malloc(size*sizeof(T));
    for(int i =0; i < size; i++){
        data[i] = a.data[i];
    }
}
template <typename T>
inline vec<T>::vec(_init_list_with_square_brackets<T> l) {
    size = l.size();
    data = (T*)malloc(size*sizeof(T));
    for(int i =0; i < size; i++){
        data[i] = l[i];
    }
};
template <typename T>
inline void vec<T>::printout(){
    for(int i=0; i < size; i++){
        cout << data[i] << endl;
    }
    cout << endl;
};


template <typename T>
inline T norm1(vec<T> a){
    T sum = 0;
    for(int i =0; i < a.size; i++){
        sum += a(i);
    }
    return sum;
};

template <typename T> 
inline vec<T> lin_solver(matrix<T> a, vec<T> b){
    matrix<T> tempa(a.rows, a.cols);
    vec<T> tempb(b.size), x(b.size); 
    //Two if else statements for size matching
    if(a.rows!=a.cols){
        cerr << "matrix not square" << endl;
        return x;
    }
    if(a.rows!=b.size){
        cerr << "dim of A and b do not match" << endl;
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


template <typename T>
inline vec<T> lin_solver_np(matrix<T> a, vec<T> b){
    matrix<T> tempa(a.rows, a.cols);
    vec<T> tempb(b.size), x(b.size); 
    //Two if else statements for size matching
    if(a.rows!=a.cols){
        cerr << "matrix not square" << endl;
        return x;
    }
    if(a.rows!=b.size){
        cerr << "dim of A and b do not match" << endl;
        return x;
    }
    //partial pivoting
    //iterate through each coloum
    tempa = a;
    tempb = b;
    if(a(0,0)==0){
        tempa.row_swap(0,1);
        tempb.row_swap(0,1);
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

template <typename T>
inline void lu_decomp(matrix<T> a, matrix<T> &l, matrix<T> &u, matrix<T> &p){
    //assumptions about inputs l u and p are zeros
    p.set_diag(1);
    u = a;
    if(a.rows!=a.cols){
        cerr << "matrix not square" << endl;
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
}

template <typename T>
inline T determinant(matrix<T> a){
    T det = 1;
    if(a.rows!=a.cols){
        cerr << "matrix not square" << endl;
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
inline T norm2(vec<T> a){
    return sqrt(dot(a,a));
}
template<typename T>
inline T dot(vec<T>a, vec<T>b){
    T sum = 0;
    if(a.size != b.size){
        cerr << "dim do not match, dotProduct" << endl;
        return sum;
    }
    for(int i = 0; i < a.size; i ++){
        sum +=(a(i))*(b(i));
    }
    return sum;
}

template <typename T>
inline matrix<T> house_holder(vec<T> a){
    int n = a.size;
    matrix<T> temp(n,n);
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            temp(i,j) = -2 *  a(i) * a(j);
        }
    }
    
    return temp;
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
}

/*
template <typename T>
inline T nInt(vec<T> (*f)(vec<T>), T a, T b){
    if(typeid(T) == int){

    }else if (typeid(T) == float){

    }else if(typeif(T) == double)
}
*/

// store each iteration in a vector; and then display the vector 
template <typename T>
inline vec<T> euler(T (*f)(T, T), T h, T a, T b, T inital){
    int its = int((b-a)/h);
    vec<T> val(its+1);
    val(0) = inital;
    for(int i =1; i < val.size; i++){
        val(i) = val(i-1) + h*f(h*i, val(i-1));
    }
    return val;
}

template <typename T>
inline vec<T> rk4(T (*f)(T,T), T h, T a, T b, T inital){
    int its = int((b-a)/h);
    vec<T> val(its+1);
    vec<T> k(4);
    vec<T> cons(4);
    T t = a;
    T h2 = h*0.5f;
    val(0) = inital;
    for(int i =1; i < val.size; i){
        k(0) =  f(t, val);
        k(1) = f(t+h2, val(i-1)+h2*k(0));
        k(2) = f(t+h2, val(i-1)+h2*k(1));
        k(3) = f(t+h2, val(i-1)+h*k(2));
        val(i) = val(i-1) + h*(k(0)*cons(0) + k(1)*cons(1) + k(2)*cons(3) + k(3)*cons(3));
        t += h;
    }
    return val;
}



template <typename T>//this one needs an rk function call
inline T AB4(T (*f)(T,T), T h, T a, T b, T inital){
    int its = int((b-a)/h);
    vec<T> val(its+1);
    vec<T> func_vals(4);
    val(0) = inital;
    func_vals = rk4(f, h, a, a+3*h, inital);
    for(int i = 0; i < its; i++){
        
    }
    //loop time
    for()
}

template <typename T>
inline T AB4(vec<T> (*f)(T,T), T h, T a, T b, vec<T> inital){

}




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