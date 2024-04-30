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

class CoreMismatchError(Exception):
    def __init__(self, message="That target and donor cores must be the same."):
        self.message = message
        super().__init__(self.message)

def DonorP7CWarning():
    warnings.warn("Donor P7 is C, this was found to be very inefficient in a screen.")

def MatchingP6P7Warning():
    warnings.warn("Target P6-P7 and Donor P6-P7 match, efficiency is unclear. DBL HSG forced to be anti-complementary.")

def NotCTorGTCoreWarning():
    warnings.warn("Core is not CT or GT, efficiency is unclear.")

def NotNTCoreWarning():
    warnings.warn("Core does not follow the expected NT format, likely inefficient.")