
#ifndef PROJECTO_MATRIX_H
#define PROJECTO_MATRIX_H

class Matrix
{
public:
    Matrix(int fil, int col);
    Matrix(int fil, int col, double v[], int n);
    Matrix(const Matrix& m);
    ~Matrix();

    Matrix& operator=(const Matrix& matrix2);
    Matrix  operator+(const Matrix& matrix2);
    Matrix  operator-(const Matrix& matrix2);
    Matrix  operator*(const Matrix& matrix2);
    double& operator()(const int i, const int j) const;
    Matrix operator*(double producto);

    Matrix transpuesta();

    double dot(const Matrix &matrix2);
    double norma();

    int getNumCol();
    int getNumFil();

    Matrix getFil(int fil);
    Matrix getCol(int col);


    void print();

private:
    void initMatrix();

private:
    int fil;
    int col;
    double **matrix;
};

#endif //PROJECTO_MATRIX_H
