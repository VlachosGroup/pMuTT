# -*- coding: utf-8 -*-
"""
pMuTT.models.empirical.shomate

Operations related to Shomate polynomials
"""

from pMuTT.models.empirical import BaseThermo

class Shomate(BaseThermo):
    """Stores the information for an individual Shomate specie
    Inherits from pMuTT.models.empirical.BaseThermo

    The thermodynamic properties are calculated using the following form:

    :math:`\\frac{c_P}{R}=\\frac{1}{R}\\bigg(A+Bt+Ct^2+Dt^3+\\frac{E}{t^2}
    \\bigg)`
    :math:`\\frac{H}{RT}=\\frac{1}{RT}\\bigg(At+B\\frac{t^2}{2}+C\\frac{t^3}{3}
    +D\\frac{t^4}{4}-\\frac{E}{t}+F\\bigg)`
    :math:`\\frac{S}{R}=\\frac{1}{R}\\bigg(A\\ln(t)+Bt+C\\frac{t^2}{2}+D
    \\frac{t^3}{3}-\\frac{E}{2t^2}+G\\bigg)`
    
    where :math:`t=\\frac{T}{1000}` in K

    Attributes
    ----------
        shomate_phase
    """
    def __init__(self, name, shomate_, **kwargs):
        super().__init__(name=name, **kwargs)

        pass