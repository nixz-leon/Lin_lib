#pragma once
#ifndef iostream
#include <iostream>
#endif
#ifndef string
#include <string>
#endif
#ifndef cstring
#include <cstring>
#endif


template <typename T> class matrix{
    public: //this are memory functions
        matrix();//default 
        matrix(int n, int m); //normal
        matrix(const matrix<T> &other);// copy
        matrix<T>& operator=(matrix<T>& other); //copy assigment
        //matrix<T>& operator=(matrix<T> other);
        matrix(matrix<T> &&other) noexcept;//move
        matrix<T>& operator=(matrix<T>&& other) noexcept;// move assignment 
        T& operator()(int row, int col);
        const T& operator()(int row, int col)const;
        void printout();
        ~matrix();
        int row;
        int col;
    private:
        T *data;
    public: //these will be the math operator functions and accessor functions
        
        inline void set_diag(T c);
        inline void row_swap(int a, int b);
        inline void row_add(int a, int b, T fac); 
        inline matrix<T> operator/=(T c);
        inline matrix<T> operator *=(T c);


        inline friend matrix<T> operator+(matrix<T> a, matrix<T> b){
            if((a.col!=b.col)||(a.row!=b.row)){std::cout << "Tried to add a " << a.row << ','<<a.col<< " matrix with a " << b.row << ',' << b.col << " matrix\n";exit(0);}
            matrix temp(&a);
            for(int i =0; i < a.row;i++){for(int j = 0; j < a.col;j++){temp(i,j)=a(i,j)+b(i,j);}}
            return temp;
        };
        inline friend matrix<T> operator-(matrix<T> a, matrix<T> b){
            if((a.col!=b.col)||(a.row!=b.row)){std::cout << "Tried to subtract a " << b.row << ','<<b.col<< " matrix from a " << a.row << ',' << a.col << " matrix\n";exit(0);}
            matrix temp(&a);
            for(int i =0; i < a.row;i++){for(int j = 0; j < a.col;j++){temp(i,j)=a(i,j)-b(i,j);}}
            return temp;
        };
        inline friend matrix<T> operator*(matrix<T> a, matrix<T> b){
            if((a.cols!=b.rows)){std::cout << "Tried to multiply a " << a.row << ','<<a.col<< " matrix by a " << b.row << ',' << b.col << " matrix\n";exit(0);}
            matrix<T> temp(&a);
            for(int i = 0; i < a.rows; i++){for(int j = 0; j < b.cols; j++){for(int k=0; k < b.rows; k++ ){temp(i,j) += a(i,k)*b(k,j);}}}
            return temp;
        };
        inline friend matrix<T> operator*(matrix<T> a, T b){matrix<T>temp(&a);for(int i =0; i < a.row;i++){for(int j =0; j<a.col;j++){temp = a(i,j)*b;}}return temp;};
        inline friend matrix<T> operator*(T b, matrix<T> a){matrix<T>temp(&a);for(int i =0; i < a.row;i++){for(int j =0; j<a.col;j++){temp =a(i,j)*b;}}return temp;};

};



template <typename T>
matrix<T>::matrix(){
    data = nullptr;
    row = 0;
    col = 0;
}
template <typename T>
matrix<T>::matrix(int n, int m){
    //data = new T[n*m];
    data = (T*)std::calloc(n*m,sizeof(T));
    row = n;
    col = m;
};
template <typename T>
matrix<T>::matrix(const matrix<T> &other){
    data = new T[other.row*other.col];//suposed deepcopy
    std::memcpy(data, other.data, sizeof(T)*other.col*other.row);
    row = other.row;
    col = other.col;
};
template <typename T>
inline matrix<T> &matrix<T>::operator=(matrix<T> &other){
    this->row = other.row;
    this->col = other.col;
    std::memcpy(this->data, other.data, sizeof(T)*other.col*other.row);
    return *this;
}
/*template <typename T>
inline matrix<T> &matrix<T>::operator=(matrix<T> other){
    std::memcpy(this->data, other.data, sizeof(T)*other.col*other.row);
    return *this;
};
*/

template <typename T>
matrix<T>::matrix(matrix<T> &&other) noexcept: data(nullptr), row(0),col(0){
    data = other.data; // reassigning ownership of pointer
    row = other.row; // redefining rows
    col = other.col; 
    other.data = nullptr; //releasing ownership of pointer from object passed in
    other.col = 0; // setting size to 0 to relfect size;
    other.row = 0;
};
template <typename  T>
inline matrix<T> &matrix<T>::operator=(matrix<T> &&other)noexcept{
    if(this != &other){
        delete[] data;
        data  = other.data;
        row = other.row;
        col = other.col;
        other.data = nullptr;
        other.row = 0;
        other.col = 0;
    }
    return *this;
}

template <typename T>
inline T &matrix<T>::operator()(int r, int c)
{
    if(row*col == 0){std::cout << "empty matrix\n";exit(0);}//checking if the object can be accessed in the first place 
    if((r > row-1)||(c > col-1)){std::cout << "Tried to access elm " << r << ',' << c << ", But size is " << row << ',' << col << '\n';exit(0);} // bounds checking
    if(r < 0){std::cout << "tried to access row " << r << " which is negative\n";}
    if(c < 0){std::cout << "tried to access row " << c << " which is negative\n";}
    return data[(r*col)+c];
}
template <typename T>
inline const T &matrix<T>::operator()(int r, int c) const
{
    if(row*col == 0){std::cout << "empty matrix\n";exit(0);}//checking if the object can be accessed in the first place 
    if((r > row-1)||(c > col-1)){std::cout << "Tried to access elm " << r << ',' << c << ", But size is " << row << ',' << col << '\n';exit(0);} // bounds checking
    if(r < 0){std::cout << "tried to access row " << r << " which is negative\n";}
    if(c < 0){std::cout << "tried to access row " << c << " which is negative\n";}
    return data[(r*col)+c];;
}

template <typename T>
inline void matrix<T>::printout()
{
    if(data !=nullptr){
        std::cout << '\n';
        for(int i = 0; i < row; i++){
            for(int j = 0; j < col; j++){
            std::cout << data[(col*i)+j] << ' ';
            }
        std::cout << '\n';
        }
    }
}
template <typename T>
inline matrix<T>::~matrix()
{
    delete[] data;
    row = 0;
    col = 0;
}

template <typename T>
inline void matrix<T>::set_diag(T c)
{
    for(int i = 0; i < row; i++){
        data[(i*col)+i]=c;
    }
}

template <typename T>
inline void matrix<T>::row_swap(int a, int b)
{
    if(a > row-1){std::cout << "Tried to access elm " << a << ", But size is " << row-1 << '\n';exit(0);}
    if(b > row-1){std::cout << "Tried to access elm " << b << ", But size is " << row-1 << '\n';exit(0);}
    if(a==b){return;}
    for(int i =0; i < col; i++){
        T temp = data[(a*col)+i];
        data[(a*col)+i] = data[(b*col)+i];
        data[(b*col)+i] = temp;
    } 
    
}

template <typename T>
inline void matrix<T>::row_add(int a, int b, T fac)
{
    if(a > row-1){std::cout << "Tried to access elm " << a << ", But size is " << row-1 << '\n';exit(0);}
    if(b > row-1){std::cout << "Tried to access elm " << b << ", But size is " << row-1 << '\n';exit(0);}
    if(a==b){return;}
    for(int i =0; i < col; i++){
        data[(a*col)+i] += fac*data[(b*col)+i];
    }
}

template <typename T>
inline matrix<T> matrix<T>::operator/=(T c)
{
    for(int i=0;i<this->row*this->col;i++){this->data[i]/=c;}
    return *this;
}

template <typename T>
inline matrix<T> matrix<T>::operator*=(T c)
{
    for(int i=0;i<this->row*this->col;i++){this->data[i]*=c;}
    return *this;
}

