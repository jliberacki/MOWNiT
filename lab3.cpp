#include <vector>
#include <iostream>
#include <math.h>

template <typename T> class AGHMatrix 
{
private:
    std::vector<std::vector<T>> matrix;
    unsigned rows;
    unsigned cols;

public:
    AGHMatrix(const std::vector<std::vector<T>>& matrix);
    AGHMatrix(unsigned _rows, unsigned _cols, const T& _initial);
    AGHMatrix(const AGHMatrix<T>& rhs);
    virtual ~AGHMatrix() = default;

    // Operator overloading, for "standard" mathematical matrix operations                                                                                                                                                          
    AGHMatrix<T>& operator=(const AGHMatrix<T>& rhs);

    // Matrix mathematical operations                                                                                                                                                                                               
    AGHMatrix<T> operator+(const AGHMatrix<T>& rhs);
    AGHMatrix<T> operator*(const AGHMatrix<T>& rhs);

    //helpful cut operation
    AGHMatrix<T> cut(const unsigned& row, const unsigned& col) const;

    // Access the individual elements                                                                                                                                                                                               
    T& operator()(const unsigned& row, const unsigned& col);
    const T& operator()(const unsigned& row, const unsigned& col) const;
    
    // Printing matrix
    std::ostream& operator<<(const AGHMatrix<T>& matrix);

    // Access the row and column sizes                                                                                                                                                                                              
    unsigned get_rows() const;
    unsigned get_cols() const;

    //homework functions
    bool if_symetric() const;
    double det() const;
    AGHMatrix<T> transpose() const;
    void LU_factorization() const;
    void Cholesky_factorization() const;
    std::vector<double> gaussianElimination() const;
};

template<typename T>
AGHMatrix<T>::AGHMatrix(const std::vector<std::vector<T>>& mat) 
{
  matrix.resize(mat.size());
  for (unsigned i = 0; i < mat.size(); i++) 
  {
    matrix[i].resize(mat[i].size());
    for(unsigned j = 0; j < mat[i].size(); j++)
    {
      matrix[i][j] = mat[i][j];
    }
  }
  rows = matrix.size();
  cols = matrix[1].size();
}

// Parameter Constructor                                                                                                                                                      
template<typename T>
AGHMatrix<T>::AGHMatrix(unsigned _rows, unsigned _cols, const T& _initial) 
{
  matrix.resize(_rows);
  for (unsigned i=0; i<matrix.size(); i++) 
  {
    matrix[i].resize(_cols, _initial);
  }
  rows = _rows;
  cols = _cols;
}

// Copy Constructor                                                                                                                                                           
template<typename T>
AGHMatrix<T>::AGHMatrix(const AGHMatrix<T>& rhs) 
{
  matrix = rhs.matrix;
  rows = rhs.get_rows();
  cols = rhs.get_cols();
}

// Get the number of rows of the matrix                                                                                                                                       
template<typename T>
unsigned AGHMatrix<T>::get_rows() const 
{
  return this->rows;
}

// Get the number of columns of the matrix                                                                                                                                    
template<typename T>
unsigned AGHMatrix<T>::get_cols() const 
{
  return this->cols;
}

// Assignment Operator                                                                                                                                                        
template<typename T>
AGHMatrix<T>& AGHMatrix<T>::operator=(const AGHMatrix<T>& rhs) 
{
  if (&rhs == this)
    return *this;

  unsigned new_rows = rhs.get_rows();
  unsigned new_cols = rhs.get_cols();

  matrix.resize(new_rows);
  for (unsigned i=0; i<matrix.size(); i++) 
  {
    matrix[i].resize(new_cols);
  }

  for (unsigned i=0; i<new_rows; i++) 
  {
    for (unsigned j=0; j<new_cols; j++) 
    {
      matrix[i][j] = rhs(i, j);
    }
  }
  rows = new_rows;
  cols = new_cols;

  return *this;
}

// Access the individual elements                                                                                                                                             
template<typename T>
T& AGHMatrix<T>::operator()(const unsigned& row, const unsigned& col) 
{
  return this->matrix[row][col];
}

// Access the individual elements (const)                                                                                                                                     
template<typename T>
const T& AGHMatrix<T>::operator()(const unsigned& row, const unsigned& col) const 
{
  return this->matrix[row][col];
}

// Addition of two matrices                                                                                                                                                   
template<typename T>
AGHMatrix<T> AGHMatrix<T>::operator+(const AGHMatrix<T>& rhs) 
{
  // Task 1 - implement addition of two matrices
  if(this->rows != rhs.get_rows() || this->cols != rhs.get_cols()) {
      throw "Matrixes must be same size";
  } else {
      for(int i=0; i<rhs.get_rows(); i++){
          for(int j=0; j<rhs.get_cols(); j++){
              matrix[i][j]+=rhs(i,j);
          }
      }
  }
  return *this;
}

// Left multiplication of this matrix and another                                                                                                                              
template<typename T>
AGHMatrix<T> AGHMatrix<T>::operator*(const AGHMatrix<T>& rhs) 
{
  // Task 1 - implement multiplication of two matrices
  if(this->cols != rhs.get_rows()) {
      throw "Matrix1 columns number must be same size as Matrix2 rows number";
  } else {
      AGHMatrix<double> res(this->rows, rhs.get_cols(), 0);
      for (int row = 0; row < this->rows; row++) {  
           for (int col = 0; col < rhs.get_cols(); col++) {  
               for (int inner = 0; inner < this->cols; inner++) {  
                   res(row,col) += matrix[row][inner] * rhs(inner,col);  
               }  
           }   
       } 
       return AGHMatrix<double>(res);
  }
}

// Printing matrix                                                                                                                        
template<typename T>
std::ostream& operator<<(std::ostream& stream, const AGHMatrix<T>& matrix) 
{
  for (int i=0; i<matrix.get_rows(); i++) 
  { 
    for (int j=0; j<matrix.get_cols(); j++) 
    {
        stream << matrix(i,j) << ", ";
    }
    stream << std::endl;
  }
    stream << std::endl;
}

template<typename T>
bool AGHMatrix<T>::if_symetric() const {
    if(this->cols != this->rows) return false;
    for(int row=0; row<this->rows; row++){
        for(int col=row; col<this->cols; col++){
            if(matrix[row][col]!=matrix[col][row]) return false;
        }
    }
    return true;
}

template<typename T>
AGHMatrix<T> AGHMatrix<T>::transpose() const {
    AGHMatrix<double> res(this->cols, this->rows, 0);
    for (int row=0; row < this->rows; row++) {
        for (int col=0; col < this->cols; col++) {
            res(col,row) = matrix[row][col];
        }
    }
    return AGHMatrix<double> (res);
}

template<typename T>
AGHMatrix<T> AGHMatrix<T>::cut(const unsigned& cut_row, const unsigned& cut_col) const {
    AGHMatrix<double> res(*this);
    res.matrix.erase(res.matrix.begin()+cut_row);
    for(int i=0; i<this->rows-1; i++){
        res.matrix[i].erase(res.matrix[i].begin()+cut_col);
    }
    res.cols -= 1;
    res.rows -= 1;
    return  AGHMatrix<double>(res);
}

template<typename T>
double AGHMatrix<T>::det() const {
    if (this->rows != this->cols) {
        throw "Matrix must be square to calculate det";
    }
    if(matrix.size() == 1) return matrix[0][0];
    else {
        double result=0;
        int sign = 1;

        for (int i=0; i<matrix.size(); i++) { 
            AGHMatrix<double> cofactor = this->cut(0,i);
            result += sign * matrix[0][i] * cofactor.det();
            sign *= -1;
        }
        return result; 
    }
}

template<typename T>
void AGHMatrix<T>::LU_factorization() const {

    if (this->rows != this->cols) {
        throw "Matrix must be square to factorize";
    }

    AGHMatrix<double> lower(this->rows, this->cols, 0);
    AGHMatrix<double> upper(this->rows, this->cols, 0);

    for (int row = 0; row < this->rows; row++) { 
        for (int col = row; col < this->cols; col++) { 
            double sum = 0; 
            for (int i = 0; i < row; i++) 
                sum += lower(row,i) * upper(i,col); 
            upper(row,col) = matrix[row][col] - sum; 
        } 
  
        for (int col = row; col < this->cols; col++) { 
            if (row == col) 
                lower(row,row) = 1;
            else { 
                double sum = 0; 
                for (int i = 0; i < row; i++) 
                    sum += lower(col,i) * upper(i,row); 
                lower(col,row) = (matrix[col][row] - sum) / upper(row,row); 
            } 
        } 
    }
    std::cout << "Upper:" << std::endl << upper << std::endl;
    std::cout << "Lower:" << std::endl << lower << std::endl;
}

template<typename T>
void AGHMatrix<T>::Cholesky_factorization() const {

    if (this->rows != this->cols) {
        throw "Matrix must be square to factorize";
    }

    AGHMatrix<double> lower(this->rows, this->cols, 0); 

    for (int row = 0; row < this->rows; row++) { 
        for (int col = 0; col <= row; col++) { 
            double sum = 0; 
            if (col == row)
            { 
                for (int i = 0; i < col; i++) 
                    sum += pow(lower(col,i), 2); 
                lower(col,col) = sqrt(matrix[col][col] - sum); 
            } else { 
                for (int i = 0; i < col; i++) 
                    sum += (lower(row,i) * lower(col,i)); 
                lower(row,col) = (matrix[row][col] - sum) / lower(col,col); 
            } 
        } 
    }
    AGHMatrix<double> upper = lower.transpose();
    std::cout << "Lower:" << std::endl << lower << std::endl;
    std::cout << "Upper:" << std::endl << upper << std::endl;
}

template<typename T>
std::vector<double> AGHMatrix<T>::gaussianElimination() const {

    AGHMatrix<double> tmp(*this);

    if(tmp.get_rows()+1 != tmp.get_cols()) throw "Matrix must be n x n+1 size";

    for (int i=0; i<tmp.get_rows(); i++) {
        // Szukamy elementu maksymalnego w danej kolumnie
        double max = abs(tmp(i,i));
        int index = i;
        for (int k=i+1; k<tmp.get_rows(); k++) {
            if (abs(tmp(k,i)) > max) {
                max = abs(tmp(k,i));
                index = k;
            }
        }

        // Zamieniamy aktualny wiersz z tym maksymalnym
        if (index != i) std::swap(tmp.matrix[i], tmp.matrix[index]);

        // Od każdego wiersza poniżej odejmujemy wielokrotność aktualnego aby je wyzerować
        for (int k=i+1; k<tmp.get_rows(); k++) {
            double subtractor = -tmp(k,i)/tmp(i,i);
            for (int j=i; j<tmp.get_cols(); j++) {
                if (i==j) {
                    tmp(k,j) = 0;
                } else {
                    tmp(k,j) += subtractor * tmp(i,j);
                }
            }
        }
    }

    // Po wyznaczeniu macierzą trójkątnej wyznaczamy poszczególne współczynniki
    std::vector<double> res(tmp.get_rows());
    for (int i=tmp.get_rows()-1; i>=0; i--) {
        res[i] = tmp(i,tmp.get_rows())/tmp(i,i);
        for (int k=i-1;k>=0; k--) {
            tmp(k,tmp.get_rows()) -= tmp(k,i) * res[i];
        }
    }
    return res;
}


int main() 
{
    try {
        // initialize matrices using init value
        AGHMatrix<double> mat1(5, 5, 1.2);
        AGHMatrix<double> mat2(5, 5, 2.8);

        // Uncomment when implemented
        AGHMatrix<double> mat3 = mat1 + mat2;
        std::cout << mat3;

        // initialize matrix using specified values
        std::vector<std::vector<double>> init { { 2.0, 1.0, 3.0 }, 
                                                { -1.0, 2.0, 4.0 } }; 

        std::vector<std::vector<double>> init2 { { 1.0, 3.0 }, 
                                                { 2.0, -2.0 },
                                                { -1.0, 4.0 } }; 
        
        std::vector<std::vector<double>> init3 { { 1.0, 2.0, 3.0 }, 
                                                { 6.0, 5.0, 4.0 },
                                                { 3.0, 7.0, 2.0 } }; 

        std::vector<std::vector<double>> init_LU {{ 5.0, 3.0, 2.0 }, 
                                            { 1.0, 2.0, 0.0 }, 
                                            { 3.0, 0.0, 4.0 }};

        std::vector<std::vector<double>> init_cholesky {{ 4.0, 12.0, -16.0 }, 
                                                    { 12.0, 37.0, -43.0 }, 
                                                    { -16.0, -43.0, 98.0 }};

        std::vector<std::vector<double>> init_gauss {{ 0.0001, -5.0300, 5.8090, 7.8320, 9.5740 }, 
                                                    { 2.2660, 1.9950,  1.2120, 8.0080, 7.2190 }, 
                                                    { 8.8500, 5.6810,  4.5520, 1.3020, 5.7300 },
                                                    { 6.7750, -2.253,  2.9080, 3.9700, 6.2910 }};

        // //zad1.2
        AGHMatrix<double> mat4(init);
        AGHMatrix<double> mat5(init2);
        AGHMatrix<double> mat_det(init3);
        std::cout << mat4;
        std::cout << mat5;
        AGHMatrix<double> mat0 = mat4 * mat5;
        std::cout << mat0;

        //zad2.1
        std::cout << mat1;
        std::cout << mat1.if_symetric() << std::endl;
        mat1(2,0)=0;
        std::cout << mat1;
        std::cout << mat1.if_symetric() << std::endl;

        //zad2.2 cut
        AGHMatrix<double> mat_cut = mat0.cut(0,0);
        std::cout << mat_cut << std::endl;

        //zad2.2 det
        std::cout << mat_det << std::endl;
        std::cout << mat_det.det() << std::endl;

        //zad3
        AGHMatrix<double> mat_lu(init_LU);
        std::cout << mat_lu << std::endl;
        mat_lu.LU_factorization();

        //zad4
        AGHMatrix<double> mat_cholesky(init_cholesky);
        std::cout << mat_cholesky << std::endl;
        mat_cholesky.Cholesky_factorization();

        //zad5
        AGHMatrix<double> mat_gauss(init_gauss);
        std::cout<< mat_gauss << std::endl;
        std::vector<double> res = mat_gauss.gaussianElimination();
        for (int i=0; i<res.size(); i++) {
            std::cout << res[i] << " ";
        }
    } catch (const char* msg) {
        std::cerr << msg << std::endl;
    }
    std::cout << std::endl;
    system ("PAUSE");
}