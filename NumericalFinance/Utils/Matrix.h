#include <vector>

class Matrix {
public:
    std::vector<std::vector<double>> data;
    // Constructors
    Matrix(); // Default constructor
    explicit Matrix(const std::vector<std::vector<double>>& vect); // Initializing constructor

    // Member Functions
    void setMatrixData(const Matrix& other_mat);
    void addRow(const std::vector<double>& newRow, int newRowIndex);
    Matrix getTranspose();
    bool luDecomposition(const std::vector<std::vector<double>>& matrix, std::vector<std::vector<double>>& L, std::vector<std::vector<double>>& U);
    double determinantLU();
    Matrix inverseLU();

    // Overloaded Operators
    Matrix operator-(const Matrix &other) const; // Difference of two matrix objects
    Matrix operator+(const Matrix &other) const; // Sums two matrix objects
    Matrix operator*(const Matrix &other) const; // Multiply two matrix objects

private:
};