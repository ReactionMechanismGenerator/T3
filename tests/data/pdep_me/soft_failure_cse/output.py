#Thermo used:  
#H(34)  SMILES: [H]
#Thermo library: Klippenstein_Glarborg2016
#CH2(T)(48)  SMILES: [CH2]
#Thermo library: thermo_DFT_CCSDTF12_BAC
#C1rad(2)  SMILES: [CH3]
#Thermo library: thermo_DFT_CCSDTF12_BAC

#Path Reactions used:  
#H(34) + CH2(T)(48) <=> C1rad(2)
#Kinetics taken from the arrheniusHigh attribute of a Troe/Lindemann exprssion. Originally from reaction library FFCM1(-)

#   =========== =========== =========== =========== =========== =========== 
#         T / P   1.184e-01   4.153e-01   3.162e+00   2.408e+01   8.445e+01
#   =========== =========== =========== =========== =========== =========== 
#       302.491  -0.000e+00  -0.000e+00   4.043e-65   5.828e-65   6.302e-65
#       323.355   2.531e-54  -0.000e+00   5.153e-60   8.242e-60   8.993e-60
#       370.585  -0.000e+00   2.396e-50   1.593e-50   2.693e-50   3.004e-50
#       457.988   2.985e-38   1.050e-38   3.752e-38   7.103e-38   8.293e-38
#       614.983   6.430e-26   1.965e-25   8.091e-25   1.833e-24   2.341e-24
#       900.017   8.174e-14   2.594e-13   1.248e-12   3.621e-12   5.423e-12
#        1394.8   8.168e-05   2.690e-04   1.504e-03   5.739e-03   1.063e-02
#       1985.54   3.305e+00   1.113e+01   6.853e+01   3.174e+02   6.942e+02
#   =========== =========== =========== =========== =========== =========== 
pdepreaction(
    reactants = ['C1rad(2)'],
    products = ['H(34)', 'CH2(T)(48)'],
    kinetics = Chebyshev(
        coeffs = [
            [
                None,
                None,
                None,
                None,
            ],
            [
                None,
                None,
                None,
                None,
            ],
            [
                None,
                None,
                None,
                None,
            ],
            [
                None,
                None,
                None,
                None,
            ],
            [
                None,
                None,
                None,
                None,
            ],
            [
                None,
                None,
                None,
                None,
            ],
        ],
        kunits = 's^-1',
        Tmin = (300, 'K'),
        Tmax = (2100, 'K'),
        Pmin = (0.1, 'bar'),
        Pmax = (100, 'bar'),
    ),
)

