class Material:
    def __init__(self, rho, cp, k):
        """
        Create a material with thermal properties.

        Arguments:
        rho - density of the material in kg/m^3
        cp - heat capacity of the material in J/kg
        k - heat conduction of the material in W/mK
        """
        self.rho = rho
        self.cp = cp
        self.k = k
