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
    AGHMatrix<T> jacoby(const AGHMatrix<T>& x, int iterations) const;
    AGHMatrix<T> gauss_seidel(const AGHMatrix<T>& x, int iterations) const;
    AGHMatrix<T> sor(const AGHMatrix<T>& x, int iterations) const;
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
AGHMatrix<T> AGHMatrix<T>::jacoby(const AGHMatrix<T>& b, int iterations) const {

    AGHMatrix<double> x(b);
    AGHMatrix<double> x_copy(b);

    if(x.get_cols() != 1) {
        throw "Input matrix must have one column";
    }
    if (this->rows != this->cols || matrix.size() != b.get_rows()) {
        throw "Matrix must have same amount of cols and rows, and must be same amount of rows as input matrix";
    }

    for(int i = 0; i < iterations; i++) {
        for(int row = 0; row < this->rows; row++) {
            double sum = 0;
            for(int col = 0; col < this->cols; col++) {
                if(col != row) {
                    sum += matrix[row][col] * x_copy(col,0);
                }
            }
            x(row,0) = (b(row,0) - sum) / matrix[row][row];
        }
        x_copy = x;
    }

    return x;
}

template<typename T>
AGHMatrix<T> AGHMatrix<T>::gauss_seidel(const AGHMatrix<T>& b, int iterations) const {

    AGHMatrix<double> x(b);
    AGHMatrix<double> x_copy(b);

    if(x.get_cols() != 1) {
        throw "Input matrix must have one column";
    }
    if (this->rows != this->cols || matrix.size() != b.get_rows()) {
        throw "Matrix must have same amount of cols and rows, and must be same amount of rows as input matrix";
    }

    for(int i = 0; i < iterations; i++) {
        for(int row = 0; row < this->rows; row++) {
            double sum = 0;
            for(int col = 0; col < this->cols; col++) {
                if(col > row) {
                    sum += matrix[row][col] * x_copy(col,0);
                }
                if(col < row) {
                    sum += matrix[row][col] * x(col,0);
                }
            }
            x(row,0) = (b(row,0) - sum) / matrix[row][row];
        }
        x_copy = x;
    }

    return x;
}

template<typename T>
AGHMatrix<T> AGHMatrix<T>::sor(const AGHMatrix<T>& b, int iterations) const {
    double omega = 1.25; //na podstawie: http://marta.certik.cz/NM/Asor.pdf
    AGHMatrix<double> x(b);
    AGHMatrix<double> x_copy(b);

    if(x.get_cols() != 1) {
        throw "Input matrix must have one column";
    }
    if (this->rows != this->cols || matrix.size() != b.get_rows()) {
        throw "Matrix must have same amount of cols and rows, and must be same amount of rows as input matrix";
    }

    for(int i = 0; i < iterations; i++) {
        for(int row = 0; row < this->rows; row++) {
            double sum = 0;
            for(int col = 0; col < this->cols; col++) {
                if(col > row) {
                    sum += matrix[row][col] * x_copy(col,0);
                }
                if(col < row) {
                    sum += matrix[row][col] * x(col,0);
                }
            }
            x(row,0) = (1 - omega) * x_copy(row,0) + omega * (b(row,0) - sum) / matrix[row][row];
        }
        x_copy = x;
    }

    return x;
}

int main() 
{
    try {
        // initialize matrices using init value
        AGHMatrix<double> mat1(5, 5, 1.2);

        // initialize matrix using specified values
        std::vector<std::vector<double>> init11 { { 3.0, -1.0, 1.0 }, 
                                                { -1.0, 3.0, -1.0 },
                                                { 1.0, -1.0, 3.0 } }; 

        std::vector<std::vector<double>> init12 { { -1.0 }, 
                                                { 7.0 },
                                                { -7.0 } }; 

        std::vector<std::vector<double>> solution1 { { 1.0 }, 
                                                { 2.0 },
                                                { -2.0 } };

        std::vector<std::vector<double>> init21 { { 5.0, -3.0 }, 
                                                { 1.0, -2.0 } }; 

        std::vector<std::vector<double>> init22 { { 21.0 }, 
                                                { 7.0 } }; 

        std::vector<std::vector<double>> solution2 { { 3.0 }, 
                                                { -2.0 } };  

        std::vector<std::vector<double>> init31 { { 4.0, -1.0, -0.2, 2.0 }, 
                                                { -1.0, 5.0, 0.0, -2.0 },
                                                { 0.2, 1.0, 10.0, -1.0 },
                                                { 0.0, -2.0, -1.0, 4.0 } }; 

        std::vector<std::vector<double>> init32 { { 21.6 }, 
                                                { 36.0 },
                                                { -20.6 },
                                                { -11.0 } };

        std::vector<std::vector<double>> solution3 { { 6.0 }, 
                                                { 9.0 },
                                                { -3.0 },
                                                { 1.0 } };  

        std::vector<std::vector<double>> init41 {{ 5.02, 2.01, -0.98 }, 
                                                    { 3.03, 6.95, 3.04 }, 
                                                    { 1.01, -3.99, 5.98 }};
        
        std::vector<std::vector<double>> init42 {{ 2.05 }, 
                                                    { -1.02 }, 
                                                    { 0.98 }};

        std::vector<std::vector<double>> solution4 {{ 0.50774 }, 
                                                    { -0.31141 }, 
                                                    { -0.12966 }};                                            

        std::vector<std::vector<double>> init51 {{ 1.0, 0.0, 0.0, 0.0, 0.0 }, 
                                                    { 0.0, 1.0,  0.0,0.0, 0.0 }, 
                                                    { 0.0, 0.0,  1.0, 0.0, 0.0 },
                                                    { 0.0, 0.0,  0.0, 1.0, 0.0 },
                                                    { 0.0, 0.0,  0.0, 0.0, 1.0 }};
        
        std::vector<std::vector<double>> init52 {{ 1.0 }, 
                                                    { 1.0 }, 
                                                    { 1.0 }, 
                                                    { 1.0 }, 
                                                    { 1.0 }};

        std::vector<std::vector<double>> solution5 {{ 1.0 }, 
                                                    { 1.0 }, 
                                                    { 1.0 }, 
                                                    { 1.0 }, 
                                                    { 1.0 }};
        
        AGHMatrix<double> mat11(init11);
        AGHMatrix<double> mat12(init12);
        AGHMatrix<double> sol1(solution1);

        AGHMatrix<double> mat21(init21);
        AGHMatrix<double> mat22(init22);
        AGHMatrix<double> sol2(solution2);
        
        AGHMatrix<double> mat31(init31);
        AGHMatrix<double> mat32(init32);
        AGHMatrix<double> sol3(solution3);

        AGHMatrix<double> mat41(init41);
        AGHMatrix<double> mat42(init42);
        AGHMatrix<double> sol4(solution4);

        AGHMatrix<double> mat51(init51);
        AGHMatrix<double> mat52(init52);
        AGHMatrix<double> sol5(solution5);

        AGHMatrix<double> res11 = mat11.jacoby(mat12, 10);
        AGHMatrix<double> res12 = mat11.gauss_seidel(mat12, 10);
        AGHMatrix<double> res13 = mat11.sor(mat12, 10);

        AGHMatrix<double> res21 = mat21.jacoby(mat22, 10);
        AGHMatrix<double> res22 = mat21.gauss_seidel(mat22, 10);
        AGHMatrix<double> res23 = mat21.sor(mat22, 10);

        AGHMatrix<double> res31 = mat31.jacoby(mat32, 10);
        AGHMatrix<double> res32 = mat31.gauss_seidel(mat32, 10);
        AGHMatrix<double> res33 = mat31.sor(mat32, 10);

        AGHMatrix<double> res41 = mat41.jacoby(mat42, 10);
        AGHMatrix<double> res42 = mat41.gauss_seidel(mat42, 10);
        AGHMatrix<double> res43 = mat41.sor(mat42, 10);

        AGHMatrix<double> res51 = mat51.jacoby(mat52, 10);
        AGHMatrix<double> res52 = mat51.gauss_seidel(mat52, 10);
        AGHMatrix<double> res53 = mat51.sor(mat52, 10);

        std::cout << "A:" << std::endl << mat51 << std::endl;
        std::cout << "B:" << std::endl <<mat52 << std::endl;
        std::cout << "jacoby:" << std::endl <<res51 << std::endl;
        std::cout << "gauss:" << std::endl <<res52 << std::endl;
        std::cout << "sor:" << std::endl <<res53 << std::endl;
        std::cout << "actual solution:" << std::endl <<sol5 << std::endl;
    } catch (const char* msg) {
        std::cerr << msg << std::endl;
    }
    std::cout << std::endl;
    system ("PAUSE");
}