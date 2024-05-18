#include "../include/Matrix.h"
#include <iostream>
#include <iomanip>
#include <cmath>

Matrix::Matrix(int fil, int col) : fil(fil), col(col) {
    initMatrix();
}

Matrix::Matrix(int fil, int col, double v[], int n) : fil(fil), col(col) {
    initMatrix();

    int k = 0;

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++) {
            if (k < n)
                matrix[i][j] = v[k++];
            else
                matrix[i][j] = 0;
        }
}

Matrix::Matrix(const Matrix &m) {
    *this = m;
}

Matrix::~Matrix() {
    for (int i = 0; i < fil; i++)
        delete[] matrix[i];

    delete[] matrix;
}

void Matrix::initMatrix() {
    matrix = new double *[fil];
    for (int i = 0; i < fil; i++)
        matrix[i] = new double[col];

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            matrix[i][j] = 0.0;
}

Matrix &Matrix::operator=(const Matrix &matrix2) {
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            this->matrix[i][j] = matrix2.matrix[i][j];

    return *this;
}

Matrix Matrix::operator+(const Matrix &matrix2) {
    Matrix result(fil, col);

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] + matrix2.matrix[i][j];

    return result;
}

Matrix Matrix::operator-(const Matrix &matrix2) {
    Matrix result(fil, col);

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] - matrix2.matrix[i][j];

    return result;
}

Matrix Matrix::operator*(const Matrix &matrix2) {
    Matrix result(fil, col);

    for (int i = 0; i < this->fil; i++) {
        for (int j = 0; j < matrix2.col; j++) {
            result.matrix[i][j] = 0;
            for (int k = 0; k < this->col; k++) {
                result.matrix[i][j] = result.matrix[i][j] + this->matrix[i][k] * matrix2.matrix[k][j];
            }
        }
    }

    return result;
}

Matrix Matrix::operator*(double producto) {
    Matrix result(fil, col);

    for (int i = 0; i < this->fil; i++) {
        for (int j = 0; j < this->col; j++) {
            result.matrix[i][j] = matrix[i][j] * producto;
        }
    }
    return result;
}


double &Matrix::operator()(const int i, const int j) const {
    return matrix[i - 1][j - 1];
}

double Matrix::norma() {
    double n = 0.0;
    for (int i = 0; i < fil; i++) {
        for (int j = 0; j < col; j++) {
            n = n + matrix[i][j] * matrix[i][j];
        }
    }
    return sqrt(n);
}


double Matrix::dot(const Matrix &matrix2) {
    double n = 0.0;
    for (int i = 0; i < col; i++) {
        n += matrix[0][i] * matrix2.matrix[0][i];
    }
    return n;
}


void Matrix::print() {
    for (int i = 0; i < fil; i++) {
        for (int j = 0; j < col; j++) {
            std::cout << std::fixed << std::setprecision(14) << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}


int Matrix::getNumCol() {
    return this->col;
}

int Matrix::getNumFil() {
    return this->fil;

}

Matrix Matrix::getFil(int fil) {
    Matrix result(1, this->col);

    for (int i = 0; i < this->col; i++) {
        result.matrix[0][i] = matrix[fil - 1][i];

    }
    return result;

}

Matrix Matrix::getCol(int col) {
    Matrix result(1, this->fil);

    for (int i = 0; i < this->fil; i++) {
        result.matrix[0][i] = matrix[i][col -1];

    }
    return result;

}

Matrix Matrix::transpuesta() {
    Matrix result(col, fil);
    for (int i = 0; i < fil; ++i) {
        for (int j = 0; j < col; ++j) {
            result(j + 1, i + 1) = matrix[i][j];
        }
    }
    return result;
}

bool Matrix::equals(const Matrix &matrix2, int e) {
    if (fil != matrix2.fil || col != matrix2.col){
        return false;
    } else{
        for (int i = 0; i < fil; i++) {
            for (int j = 0; j < col; j++) {
                if (fabs(matrix[i][j] - matrix2.matrix[i][j]) > pow(10, -e)){
                    return false;
                }

            }
        }
        return true;
    }
}


Matrix Matrix::getPrimeraFil(int inic,int fin) {
    int lon = fin - inic +1;
    Matrix result(1, lon);

    for (int i = 0; i < lon; i++) {
        result.matrix[0][i] = matrix[i][0];

    }
    return result;

}

Matrix Matrix::concatenar(const Matrix &matrix2) {
    Matrix result(fil, col + matrix2.col);

    for (int i = 0; i < fil; ++i) {
        for (int j = 0; j < col; ++j) {
            result.matrix[i][j] = matrix[i][j];
        }
    }
    for (int i = 0; i < matrix2.fil; ++i) {
        for (int j = 0; j < matrix2.col; ++j) {
            result.matrix[i][col + j] = matrix2.matrix[i][j];
        }
    }

    return result;
}




