class Conversion:
    """Converts the inputs of an electrochemical simulation into dimensionless equivalents"""
    
    def __init__(self, E, tmax, d):
        self.E = E
        self.tmax = tmax
        #alpha, K0, K1, K2
        
class Deconversion:
    """Converts the outputs of an electrochemical simulation into quantities with dimensionality"""
    
    def __init__(self, theta, flux):
        self.theta = theta
        self.flux = flux