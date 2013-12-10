class HGVSError(Exception):
    pass

class HGVSParseError(HGVSError):
    pass

class InvalidIntervalError(HGVSError):
    pass

###
#class UTAError(Exception):
#    pass
#
#class DatabaseError(UTAError):
#    pass
#
#class InvalidTranscriptError(UTAError):
#    pass
#
#
#
class InvalidHGVSVariantError(HGVSError):
   pass
