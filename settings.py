from models import *

COMPONENTS = {
    "pointSource": {
        "id": 1,
        "variables": [ 
            "RA", "Dec", "amp", "offset"
        ],
        "analytical"  :   True,
        "spectrum"    :  False, 
        "var_to_print":[
            'RA (deg)','Dec (deg)','Amplitude (Jy)','Offset'
        ]
    },
   "gaussSource": {
        "id": 2,
        "variables": [ 
            "RA", "Dec", "amp", "Major", "e", "Angle", "offset"
        ],
        "analytical"  :   True,
        "spectrum"    :  False, 
        "var_to_print":[
            'RA (deg)','Dec (deg)','Amplitude (Jy)','Major (deg)','e','Angle (deg)','Offset'
        ]
    },
    "gaussSurface": {
        "id": 3,
        "variables": [ 
            "RA", "Dec", "amp", "Major", "e", "Angle", "offset", "Temp"
        ],
        "analytical"  :   True,
        "spectrum"    :  False, 
        "var_to_print":[
            'RA (deg)','Dec (deg)','Amplitude (Jy)','Major (deg)','e','Angle (deg)','Offset','Temperature (keV)'
        ]
    },
    "betaPressure": {
        "id": 4,
        "function": betaProfile,
        "variables": [ 
            "RA", "Dec", "amp", "Major", "e", "Angle", "offset", "Temp", "beta", 'z'
        ],
        "analytical"  :  False,
        "spectrum"    :  False, 
        "var_to_print":[
            'RA (deg)','Dec (deg)','Amplitude (keV/cm3)','Major (deg)','e','Angle (deg)','Offset','Temperature (keV)','Beta','z'
        ]
    },
    "gnfwPressure": {
        "id": 5,
        "function": gnfwProfile,
        "variables": [ 
            "RA", "Dec", "amp", "Major", "e", "Angle", "offset", "Temp", "beta", "gamma", "z"
        ],
        "analytical"  :  False,
        "spectrum"    :  False, 
        "var_to_print":[
            'RA (deg)','Dec (deg)','Amplitude (keV/cm3)','Major (deg)','e','Angle (deg)','Offset','Temperature (keV)','Beta', 'Gamma', 'z'
        ]
    },
    "A10Pressure": {
        "id": 6,
        "function": a10Profile,
        "variables": [  #major = c500?, also 'ap' and 'major' can't be found back 
            "RA", "Dec", "mass", "c500", "e", "Angle", "offset", "Temp", "alpha", "beta", "gamma", "P0", "Alpha_p", "z", "bias"
        ],
        "analytical"  :  False,
        "spectrum"    :  False, 
        "var_to_print":[
            'RA (deg)','Dec (deg)','log10(M500/Msun)','c500','e','Angle (deg)','Offset','Temperature (keV)','Alpha','Beta','Gamma','P0','Alpha_p','z','bias'
        ]
    },
    "powerLaw": {
        "id": 7,
        "variables": [ 
            "spec_index"
        ],
        "analytical"  :  True,
        "spectrum"    :  True, 
        "var_to_print":[
            'SpecIndex'
        ]
    },
    "powerLawMod": {
        "id": 8,
        "variables": [ 
            "spec_index", "spec_curv"
        ],
        "analytical"  :  True,
        "spectrum"    :  True, 
        "var_to_print":[
            'SpecIndex','SpecCurv'
        ]
    },
    "powerDust": {
        "id": 9,
        "variables": [ 
            "spec_index", "spec_curv", "Temp", "beta", "z", "kappa0", "nu0"
        ],
        "analytical"  :  True,
        "spectrum"    :  True, 
        "var_to_print":[
            'SpecIndex','log(Mass/M_sun)','Temperature (keV)','Beta','z','kappa0','nu0'
        ]
    },
    "tSZ": {
        'id':  10, 
        "variables": [ 
        ],
        "analytical"  :  True,
        "spectrum"    :  True, 
        "var_to_print":[

        ]
    }
}
