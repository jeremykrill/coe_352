import numpy as np
from scipy.linalg import svd

""" Three separate classes to deal with the nuances of each case """
class Fixed_Fixed:
    def __init__(self,masses_array,springs_array):
        self.masses = masses_array
        self.springs = springs_array
        self.At = self.At(springs)
        self.C = self.C(springs)
        self.A = self.A(self.At)
        self.K = np.dot(np.dot(self.At,self.C),self.A) # K = A^TCA
        self.K_inv = np.linalg.inv(self.K)
        self.u = self.displacements()
        self.e = self.elongations()
        self.w = self.internal_stresses()

    def At(self,springs):
        n_springs = len(springs)
        At = np.zeros([n_springs - 1, n_springs])
        for i in range(0,n_springs - 1):
            for j in range(0,n_springs):
                if(i == j):
                    At[i,j] = 1
                elif(j == i + 1):
                    At[i,j] = -1
        return At

    def C(self,springs):
        n_springs = len(springs)
        C = np.zeros([n_springs, n_springs])
        for i in range(0,n_springs):
            for j in range(0,n_springs):
                if(i == j):
                    C[i,j] = springs[i]
        return C

    def A(self,At):
        return np.transpose(At)

    def displacements(self):
        weights = self.masses * 9.81 # multiply masses by 9.81 to convert to weights
        u = np.dot(self.K_inv,np.transpose(weights))
        return u

    def elongations(self):
        e = np.dot(self.A,self.u)
        return e

    def internal_stresses(self):
        w = np.dot(self.C,self.e)
        return w

    def singular_values(self,matrix):
        U, s, Vt = np.linalg.svd(matrix)
        return s

    def eigenvalues(self,matrix):
        lambdas, v = np.linalg.eig(matrix)
        return lambdas

    def condition_number(self,matrix):
        """ Check if matrix is square """
        if(matrix.shape[0] == matrix.shape[1]):
            positive_eigenvalues = True
            for val in self.eigenvalues(matrix):
                if val <= 0:
                    positive_eigenvalues = False
            if (positive_eigenvalues):
                print("Its condition number is: ", np.max(self.eigenvalues(matrix))/np.min(self.eigenvalues(matrix)))
            else:
                print("Its condition number is: ", np.max(self.singular_values(matrix))/np.min(self.singular_values(matrix)))
        else:
            print("Its condition number is: ", np.max(self.singular_values(matrix))/np.min(self.singular_values(matrix)))

    def output(self):
        print("Matrix A is: ", self.A)
        print("Its singular values are:", self.singular_values(self.A))
        self.condition_number(self.A)
        print("Matrix A^T is: ", self.At)
        print("Its singular values are:", self.singular_values(self.At))
        self.condition_number(self.At)
        print("Matrix C is: ", self.C)
        print("Its singular values are:", self.singular_values(self.C))
        print("Its eigenvalues are: ", self.eigenvalues(self.C))
        self.condition_number(self.C)
        """ Solve for displacements using the force balance equation, u = K^-1 F """
        print("The force balance equation is: ", self.u, " = ", self.K_inv, self.masses*9.81)
        print("The elongation equation is: ", self.e, " = ", self.A, self.u)
        print("The internal stress equation is: ", self.w, " = ", self.C, self.e)

class Fixed_Free:
    def __init__(self,masses_array,springs_array):
        self.masses = masses_array
        self.springs = springs_array
        self.At = self.At(springs)
        self.C = self.C(springs)
        self.A = self.A(self.At)
        self.K = np.dot(np.dot(self.At,self.C),self.A) # K = A^TCA
        self.K_inv = np.linalg.inv(self.K)
        self.u = self.displacements()
        self.e = self.elongations()
        self.w = self.internal_stresses()

    def At(self,springs):
        n_springs = len(springs)
        At = np.zeros([n_springs, n_springs])
        for i in range(0,n_springs):
            for j in range(0,n_springs):
                if(i == j):
                    At[i,j] = 1
                elif(j == i + 1):
                    At[i,j] = -1
        return At

    def C(self,springs):
        n_springs = len(springs)
        C = np.zeros([n_springs, n_springs])
        for i in range(0,n_springs):
            for j in range(0,n_springs):
                if(i == j):
                    C[i,j] = springs[i]
        return C

    def A(self,At):
        return np.transpose(At)

    def displacements(self):
        weights = self.masses * 9.81 # multiply masses by 9.81 to convert to weights
        u = np.dot(self.K_inv,np.transpose(weights))
        return u

    def elongations(self):
        e = np.dot(self.A,self.u)
        return e

    def internal_stresses(self):
        w = np.dot(self.C,self.e)
        return w

    def singular_values(self,matrix):
        U, s, Vt = np.linalg.svd(matrix)
        return s

    def eigenvalues(self,matrix):
        lambdas, v = np.linalg.eig(matrix)
        return lambdas

    def condition_number(self,matrix):
        """ Check if matrix is square """
        if(matrix.shape[0] == matrix.shape[1]):
            positive_eigenvalues = True
            for val in self.eigenvalues(matrix):
                if val <= 0:
                    positive_eigenvalues = False
            if (positive_eigenvalues):
                print("Its condition number is: ", np.max(self.eigenvalues(matrix))/np.min(self.eigenvalues(matrix)))
            else:
                print("Its condition number is: ", np.max(self.singular_values(matrix))/np.min(self.singular_values(matrix)))
        else:
            print("Its condition number is: ", np.max(self.singular_values(matrix))/np.min(self.singular_values(matrix)))

    def output(self):
        print("Matrix A is: ", self.A)
        print("Its singular values are:", self.singular_values(self.A))
        print("Its eigenvalues are: ", self.eigenvalues(self.A))
        self.condition_number(self.A)
        print("Matrix A^T is: ", self.At)
        print("Its singular values are:", self.singular_values(self.At))
        print("Its eigenvalues are: ", self.eigenvalues(self.At))
        self.condition_number(self.At)
        print("Matrix C is: ", self.C)
        print("Its singular values are:", self.singular_values(self.C))
        print("Its eigenvalues are: ", self.eigenvalues(self.C))
        self.condition_number(self.C)
        """ Solve for displacements using the force balance equation, u = K^-1 F """
        print("The force balance equation is: ", self.u, " = ", self.K_inv, self.masses*9.81)
        print("The elongation equation is: ", self.e, " = ", self.A, self.u)
        print("The internal stress equation is: ", self.w, " = ", self.C, self.e)

class Free_Free:
    def __init__(self,masses_array,springs_array):
        self.masses = masses_array
        self.springs = springs_array
        self.At = self.At(springs)
        self.C = self.C(springs)
        self.A = self.A(self.At)
        self.K = np.dot(np.dot(self.At,self.C),self.A) # K = A^TCA
        self.K_inv = np.linalg.pinv(self.K)
        self.u = self.displacements()
        self.e = self.elongations()
        self.w = self.internal_stresses()

    def At(self,springs):
        n_springs = len(springs)
        At = np.zeros([n_springs + 1, n_springs])
        for i in range(0,n_springs + 1):
            for j in range(0,n_springs):
                if(i == j):
                    At[i,j] = -1
                elif(j == i - 1):
                    At[i,j] = 1
        return At

    def C(self,springs):
        n_springs = len(springs)
        C = np.zeros([n_springs, n_springs])
        for i in range(0,n_springs):
            for j in range(0,n_springs):
                if(i == j):
                    C[i,j] = springs[i]
        return C

    def A(self,At):
        return np.transpose(At)

    def displacements(self):
        weights = self.masses * 9.81 # multiply masses by 9.81 to convert to weights
        u = np.dot(self.K_inv,np.transpose(weights))
        return u

    def elongations(self):
        e = np.dot(self.A,self.u)
        return e

    def internal_stresses(self):
        w = np.dot(self.C,self.e)
        return w

    def singular_values(self,matrix):
        U, s, Vt = np.linalg.svd(matrix)
        return s

    def eigenvalues(self,matrix):
        lambdas, v = np.linalg.eig(matrix)
        return lambdas

    def condition_number(self,matrix):
        """ Check if matrix is square """
        if(matrix.shape[0] == matrix.shape[1]):
            positive_eigenvalues = True
            for val in self.eigenvalues(matrix):
                if val <= 0:
                    positive_eigenvalues = False
            if (positive_eigenvalues):
                print("Its condition number is: ", np.max(self.eigenvalues(matrix))/np.min(self.eigenvalues(matrix)))
            else:
                print("Its condition number is: ", np.max(self.singular_values(matrix))/np.min(self.singular_values(matrix)))
        else:
            print("Its condition number is: ", np.max(self.singular_values(matrix))/np.min(self.singular_values(matrix)))

    def output(self):
        print("Matrix A is: ", self.A)
        print("Its singular values are:", self.singular_values(self.A))
        self.condition_number(self.A)
        print("Matrix A^T is: ", self.At)
        print("Its singular values are:", self.singular_values(self.At))
        self.condition_number(self.At)
        print("Matrix C is: ", self.C)
        print("Its singular values are:", self.singular_values(self.C))
        print("Its eigenvalues are: ", self.eigenvalues(self.C))
        self.condition_number(self.C)
        """ Solve for displacements using the force balance equation, u = K^-1 F """
        print("The force balance equation is: ", self.u, " = ", self.K_inv, self.masses*9.81)
        print("The elongation equation is: ", self.e, " = ", self.A, self.u)
        print("The internal stress equation is: ", self.w, " = ", self.C, self.e)


if __name__ == '__main__':
    print("The program determines boundary conditions based on the masses and spring constants input.")
    print("If # masses = # spring constants - 1, the system is fixed-fixed.")
    print("If # masses = # spring constants, the system is fixed-free.")
    print("If # masses = # spring constants + 1, the system is free-free.")
    input("Press any key to begin: ")
    print("Enter masses. Once fully entered, an input of '-1' lets the program know you're done.")
    print("Constants must be positive values. The program will not consider the mass if it is not.")
    masses = np.array([])
    m = 0
    count = 1
    while(m >= 0):
        m = float(input("Enter mass " + str(count) + " : "))
        if m >= 0:
            masses = np.append(masses,[m])
            count += 1
    print("Enter spring constants. Once fully entered, an input of '-1' lets the program know you're done.")
    print("Constants must be positive values. The program will not consider the spring if it is not.")
    springs = np.array([])
    sp = 0
    count = 1
    while(float(sp) >= 0):
        sp = float(input("Enter spring " + str(count) + " constant: "))
        if sp >= 0:
            springs = np.append(springs,[sp])
        count += 1
    print("Inputs for masses and springs are:")
    print(masses)
    print(springs)
    if len(springs) == len(masses) + 1:
        system = Fixed_Fixed(masses,springs)
        print("Fixed-Fixed system")
        system.output()
    elif len(springs) == len(masses):
        system = Fixed_Free(masses,springs)
        print("Fixed-Free system")
        system.output()
    elif len(springs) == len(masses) - 1:
        system = Free_Free(masses,springs)
        print("Free-Free system")
        print("Note how the force balance always results in displacements that sum to zero, no matter the values input.")
        print("This is not the case for the other two systems.")
        system.output()
    else:
        print("Invalid input.")
