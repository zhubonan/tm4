# -*- coding: utf-8 -*-
"""
Created on Wed Dec 30 13:53:35 2015
Temp file
@author: Bonan
"""

#Disable transfer matrix calculation method as the class now represent an segment of the whole
    def update_T(self):
        """
        Calcualte transfer matrix for each interface and get the overall one
        Then calculate coupled refelctivity and transmisivity
        """
        if len(self.D) == self.N + 2 and len(self.P) == self.N + 1:
            self.T = [np.linalg.solve(self.D[i], self.D[i+1].dot(self.P[i])) for i in range(self.N +1)]
        else:
            print("Mismatched D and P stack")
        self.T_total = stack_dot(self.T)
        self.coeff = calc_coeff(self.T_total)
        self.coeff_modulus = self.coeff.copy()
        for i in self.coeff_modulus:
            self.coeff_modulus[i] = np.abs(self.coeff_modulus[i])**2
    def LR_basis(self):
        """
        Caculate T_total and coefficients for LR polarised light
        """
        a = 1/math.sqrt(2)
        b = - 1j / math.sqrt(2)
        M = np.array([[a,0,a,0],[0,a,0,a],[b,0,-b,0],[0,b,0,-b]])
        N = np.array([[a,0,a,0],[0,0,0,0],[b,0,-b,0],[0,0,0,0]])
        self.T_total_LR = np.linalg.solve(M, self.T_total.dot(N))
        self.coeff_LR = calc_coeff(self.T_total_LR)
        self.coeff_modulus_LR = self.coeff_LR.copy()
        for i in self.coeff_modulus_LR:
            self.coeff_modulus_LR[i] = np.abs(self.coeff_modulus_LR[i])**2       