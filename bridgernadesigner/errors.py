import warnings

class TargetNotDNAError(Exception):
    def __init__(self, message="The target sequence must be DNA."):
        self.message = message
        super().__init__(self.message)

class DonorNotDNAError(Exception):
    def __init__(self, message="The donor sequence must be DNA."):
        self.message = message
        super().__init__(self.message)

class TargetLengthError(Exception):
    def __init__(self, message="The target sequence must be exactly 14 nucleotides long."):
        self.message = message
        super().__init__(self.message)

class DonorLengthError(Exception):
    def __init__(self, message="The donor sequence must be exactly 14 nucleotides long."):
        self.message = message
        super().__init__(self.message)

def CoreMismatchWarning():
    warnings.warn("The target and donor cores should be identical.")

def NotNTCoreWarning():
    warnings.warn("Core does not follow the expected NT format, recombination may not work.")

def P6P7Warning():
    warnings.warn("Target position 6-7 and donor position 6-7 match, efficiency may be dramatically reduced.")

def P7Warning():
    warnings.warn("Target position 7 and donor position 7 match, efficiency may be significantly reduced.")

def P6Warning():
    warnings.warn("Target position 6 and donor position 6 match, efficiency may be reduced.")