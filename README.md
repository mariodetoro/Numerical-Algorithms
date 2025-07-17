# 🧮 Numerical Algorithms

This repository contains Python implementations of fundamental numerical methods used in scientific computing and engineering. Each script corresponds to a specific class of problems, such as solving differential equations, linear systems, finding eigenvalues, and performing interpolation or integration.

---

## 📂 Repository Contents

### 📌 Differential Equations
- **`Differential_equations_Euler_method.py`**  
  Implements Euler's method for solving first-order ordinary differential equations (ODEs). Suitable for initial value problems.

---

### 📌 Calculus and Integration
- **`Differentiation_integrals.py`**  
  Contains numerical differentiation and integration methods such as:
  - Finite difference methods for derivatives
  - Trapezoidal and Simpson’s Rule for numerical integration

---

### 📌 Linear Algebra and Systems of Equations
- **`Linear_equation_system.py`**  
  Solves systems of linear equations using direct or iterative methods (e.g., Gaussian elimination or matrix inversion).

- **`Jacobi_GaussSeidel_Matrix_Resolution.py`**  
  Implements the **Jacobi** and **Gauss-Seidel** iterative methods for solving linear systems, useful for large sparse matrices.

- **`QR_method.py`**  
  Uses the QR algorithm to compute eigenvalues of a matrix iteratively through QR decomposition.

- **`Eigenvalues_calculation.py`**  
  Provides basic numerical methods for estimating eigenvalues, possibly including the Power Method or other iterative techniques.

---

### 📌 Non-linear Equations
- **`Non_linear_equations.py`**  
  Implements numerical root-finding methods for nonlinear equations such as:
  - Bisection method  
  - Newton-Raphson method  
  - Secant method

---

### 📌 Interpolation
- **`Interpolation.py`**  
  Implements interpolation techniques including:
  - Lagrange polynomials  
  - Newton’s divided differences  
  - Linear and cubic spline interpolation

---

## 🛠️ How to Use

1. Clone the repository
2. Navigate to the folder
3. Run any script using Python
   - Modify the parameters in each script to match your input problem or dataset.
  
## Dependencies

Most scripts use only standard Python libraries like math and numpy. If necessary, install dependencies via:

pip install numpy matplotlib
