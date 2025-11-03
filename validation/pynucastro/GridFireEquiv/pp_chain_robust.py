import numba
import numpy as np
from scipy import constants
from numba.experimental import jitclass

from pynucastro.rates import TableIndex, TableInterpolator, TabularRate, Tfactors
from pynucastro.screening import PlasmaState, ScreenFactors

jn = 0
jp = 1
jd = 2
jt = 3
jhe3 = 4
jhe4 = 5
jhe6 = 6
jli6 = 7
jli7 = 8
jbe7 = 9
jbe9 = 10
jb8 = 11
jc11 = 12
jc12 = 13
jc13 = 14
jc14 = 15
jn13 = 16
jn14 = 17
jn15 = 18
jo14 = 19
jo15 = 20
jo16 = 21
jo17 = 22
jo18 = 23
jf15 = 24
jf16 = 25
jf17 = 26
jf18 = 27
jf19 = 28
jne18 = 29
jne19 = 30
jne20 = 31
jne21 = 32
jna21 = 33
jna22 = 34
jna23 = 35
jna24 = 36
jmg22 = 37
jmg23 = 38
jmg24 = 39
jmg25 = 40
jmg26 = 41
jal25 = 42
jal26 = 43
jal27 = 44
jal28 = 45
jsi26 = 46
jsi28 = 47
jsi29 = 48
jsi30 = 49
jp29 = 50
jp30 = 51
jp31 = 52
jp32 = 53
js31 = 54
js32 = 55
js34 = 56
jcl32 = 57
jcl34 = 58
jcl35 = 59
jar35 = 60
nnuc = 61

A = np.zeros((nnuc), dtype=np.int32)

A[jn] = 1
A[jp] = 1
A[jd] = 2
A[jt] = 3
A[jhe3] = 3
A[jhe4] = 4
A[jhe6] = 6
A[jli6] = 6
A[jli7] = 7
A[jbe7] = 7
A[jbe9] = 9
A[jb8] = 8
A[jc11] = 11
A[jc12] = 12
A[jc13] = 13
A[jc14] = 14
A[jn13] = 13
A[jn14] = 14
A[jn15] = 15
A[jo14] = 14
A[jo15] = 15
A[jo16] = 16
A[jo17] = 17
A[jo18] = 18
A[jf15] = 15
A[jf16] = 16
A[jf17] = 17
A[jf18] = 18
A[jf19] = 19
A[jne18] = 18
A[jne19] = 19
A[jne20] = 20
A[jne21] = 21
A[jna21] = 21
A[jna22] = 22
A[jna23] = 23
A[jna24] = 24
A[jmg22] = 22
A[jmg23] = 23
A[jmg24] = 24
A[jmg25] = 25
A[jmg26] = 26
A[jal25] = 25
A[jal26] = 26
A[jal27] = 27
A[jal28] = 28
A[jsi26] = 26
A[jsi28] = 28
A[jsi29] = 29
A[jsi30] = 30
A[jp29] = 29
A[jp30] = 30
A[jp31] = 31
A[jp32] = 32
A[js31] = 31
A[js32] = 32
A[js34] = 34
A[jcl32] = 32
A[jcl34] = 34
A[jcl35] = 35
A[jar35] = 35

Z = np.zeros((nnuc), dtype=np.int32)

Z[jn] = 0
Z[jp] = 1
Z[jd] = 1
Z[jt] = 1
Z[jhe3] = 2
Z[jhe4] = 2
Z[jhe6] = 2
Z[jli6] = 3
Z[jli7] = 3
Z[jbe7] = 4
Z[jbe9] = 4
Z[jb8] = 5
Z[jc11] = 6
Z[jc12] = 6
Z[jc13] = 6
Z[jc14] = 6
Z[jn13] = 7
Z[jn14] = 7
Z[jn15] = 7
Z[jo14] = 8
Z[jo15] = 8
Z[jo16] = 8
Z[jo17] = 8
Z[jo18] = 8
Z[jf15] = 9
Z[jf16] = 9
Z[jf17] = 9
Z[jf18] = 9
Z[jf19] = 9
Z[jne18] = 10
Z[jne19] = 10
Z[jne20] = 10
Z[jne21] = 10
Z[jna21] = 11
Z[jna22] = 11
Z[jna23] = 11
Z[jna24] = 11
Z[jmg22] = 12
Z[jmg23] = 12
Z[jmg24] = 12
Z[jmg25] = 12
Z[jmg26] = 12
Z[jal25] = 13
Z[jal26] = 13
Z[jal27] = 13
Z[jal28] = 13
Z[jsi26] = 14
Z[jsi28] = 14
Z[jsi29] = 14
Z[jsi30] = 14
Z[jp29] = 15
Z[jp30] = 15
Z[jp31] = 15
Z[jp32] = 15
Z[js31] = 16
Z[js32] = 16
Z[js34] = 16
Z[jcl32] = 17
Z[jcl34] = 17
Z[jcl35] = 17
Z[jar35] = 18

# masses in ergs
mass = np.zeros((nnuc), dtype=np.float64)

mass[jn] = 0.001505349762871528
mass[jp] = 0.0015040963047307696
mass[jd] = 0.0030058819195053215
mass[jt] = 0.004501206494525079
mass[jhe3] = 0.004501176706825056
mass[jhe4] = 0.0059735574859708365
mass[jhe6] = 0.008982694165199753
mass[jli6] = 0.008977078184259593
mass[jli7] = 0.010470810414554471
mass[jbe7] = 0.010472191322584432
mass[jbe9] = 0.013449944991877256
mass[jb8] = 0.011976069136782909
mass[jc11] = 0.016433661161513517
mass[jc12] = 0.017909017027273523
mass[jc13] = 0.01940644192976114
mass[jc14] = 0.020898691600487556
mass[jn13] = 0.01940999951603316
mass[jn14] = 0.020898440897976135
mass[jn15] = 0.022386433805845516
mass[jo14] = 0.020906683078094165
mass[jo15] = 0.02239084645968795
mass[jo16] = 0.023871099855618198
mass[jo17] = 0.025369811672200496
mass[jo18] = 0.026862271330925704
mass[jf15] = 0.02241281454438738
mass[jf16] = 0.023895792605265975
mass[jf17] = 0.025374234423440733
mass[jf18] = 0.026864924401329426
mass[jf19] = 0.028353560476732827
mass[jne18] = 0.026872045275379234
mass[jne19] = 0.02835875072008801
mass[jne20] = 0.029837079292893483
mass[jne21] = 0.03133159647374143
mass[jna21] = 0.0313372792660881
mass[jna22] = 0.03282489638134515
mass[jna23] = 0.03431034746033777
mass[jna24] = 0.03580454709025751
mass[jmg22] = 0.032832557028702955
mass[jmg23] = 0.03431684618276469
mass[jmg24] = 0.03579571000144743
mass[jmg25] = 0.03728931497629966
mass[jmg26] = 0.03877689166223153
mass[jal25] = 0.03729616718134972
mass[jal26] = 0.03878330742475342
mass[jal27] = 0.040267735904095256
mass[jal28] = 0.041760708567121475
mass[jsi26] = 0.0387914290824159
mass[jsi28] = 0.041753271139658626
mass[jsi29] = 0.04324504469459885
mass[jsi30] = 0.044733396647367994
mass[jp29] = 0.0432529631025368
mass[jp30] = 0.0447401772223051
mass[jp31] = 0.04622580258463001
mass[jp32] = 0.04771843802856788
mass[js31] = 0.04623445115717168
mass[js32] = 0.04771569725373347
mass[js34] = 0.05069425905627418
mass[jcl32] = 0.047736014194634656
mass[jcl34] = 0.05070305756947744
mass[jcl35] = 0.05218814819637857
mass[jar35] = 0.05219770715067765

names = []
names.append("n")
names.append("H1")
names.append("H2")
names.append("H3")
names.append("He3")
names.append("He4")
names.append("He6")
names.append("Li6")
names.append("Li7")
names.append("Be7")
names.append("Be9")
names.append("B8")
names.append("C11")
names.append("C12")
names.append("C13")
names.append("C14")
names.append("N13")
names.append("N14")
names.append("N15")
names.append("O14")
names.append("O15")
names.append("O16")
names.append("O17")
names.append("O18")
names.append("F15")
names.append("F16")
names.append("F17")
names.append("F18")
names.append("F19")
names.append("Ne18")
names.append("Ne19")
names.append("Ne20")
names.append("Ne21")
names.append("Na21")
names.append("Na22")
names.append("Na23")
names.append("Na24")
names.append("Mg22")
names.append("Mg23")
names.append("Mg24")
names.append("Mg25")
names.append("Mg26")
names.append("Al25")
names.append("Al26")
names.append("Al27")
names.append("Al28")
names.append("Si26")
names.append("Si28")
names.append("Si29")
names.append("Si30")
names.append("P29")
names.append("P30")
names.append("P31")
names.append("P32")
names.append("S31")
names.append("S32")
names.append("S34")
names.append("Cl32")
names.append("Cl34")
names.append("Cl35")
names.append("Ar35")

def to_composition(Y):
    """Convert an array of molar fractions to a Composition object."""
    from pynucastro import Composition, Nucleus
    nuclei = [Nucleus.from_cache(name) for name in names]
    comp = Composition(nuclei)
    for i, nuc in enumerate(nuclei):
        comp.X[nuc] = Y[i] * A[i]
    return comp


def energy_release(dY):
    """return the energy release in erg/g (/s if dY is actually dY/dt)"""
    enuc = 0.0
    for i, y in enumerate(dY):
        enuc += y * mass[i]
    enuc *= -1*constants.Avogadro
    return enuc

@jitclass([
    ("n__p__weak__wc12", numba.float64),
    ("t__He3__weak__wc12", numba.float64),
    ("He3__t__weak__electron_capture", numba.float64),
    ("He6__Li6__weak__wc12", numba.float64),
    ("Be7__Li7__weak__electron_capture", numba.float64),
    ("C14__N14__weak__wc12", numba.float64),
    ("N13__C13__weak__wc12", numba.float64),
    ("O14__N14__weak__wc12", numba.float64),
    ("O15__N15__weak__wc12", numba.float64),
    ("F17__O17__weak__wc12", numba.float64),
    ("F18__O18__weak__wc12", numba.float64),
    ("Ne18__F18__weak__wc12", numba.float64),
    ("Ne19__F19__weak__wc12", numba.float64),
    ("Na21__Ne21__weak__wc12", numba.float64),
    ("Na24__Mg24__weak__wc12", numba.float64),
    ("Mg22__Na22__weak__wc12", numba.float64),
    ("Mg23__Na23__weak__wc12", numba.float64),
    ("Al25__Mg25__weak__wc12", numba.float64),
    ("Al26__Mg26__weak__wc12", numba.float64),
    ("Al28__Si28__weak__wc12", numba.float64),
    ("Si26__Al26__weak__wc12", numba.float64),
    ("P29__Si29__weak__wc12", numba.float64),
    ("P30__Si30__weak__wc12", numba.float64),
    ("P32__S32__weak__wc12", numba.float64),
    ("S31__P31__weak__wc12", numba.float64),
    ("Cl32__S32__weak__wc12", numba.float64),
    ("Cl34__S34__weak__wc12", numba.float64),
    ("Ar35__Cl35__weak__wc12", numba.float64),
    ("d__n_p", numba.float64),
    ("t__n_d", numba.float64),
    ("He3__p_d", numba.float64),
    ("He4__n_He3", numba.float64),
    ("He4__p_t", numba.float64),
    ("He4__d_d", numba.float64),
    ("Li6__He4_d", numba.float64),
    ("Li7__n_Li6", numba.float64),
    ("Li7__He4_t", numba.float64),
    ("Be7__p_Li6", numba.float64),
    ("Be7__He4_He3", numba.float64),
    ("B8__p_Be7", numba.float64),
    ("B8__He4_He4__weak__wc12", numba.float64),
    ("C11__He4_Be7", numba.float64),
    ("C12__n_C11", numba.float64),
    ("C13__n_C12", numba.float64),
    ("C14__n_C13", numba.float64),
    ("N13__p_C12", numba.float64),
    ("N14__n_N13", numba.float64),
    ("N14__p_C13", numba.float64),
    ("N15__n_N14", numba.float64),
    ("N15__p_C14", numba.float64),
    ("O14__p_N13", numba.float64),
    ("O15__n_O14", numba.float64),
    ("O15__p_N14", numba.float64),
    ("O16__n_O15", numba.float64),
    ("O16__p_N15", numba.float64),
    ("O16__He4_C12", numba.float64),
    ("O17__n_O16", numba.float64),
    ("O18__n_O17", numba.float64),
    ("O18__He4_C14", numba.float64),
    ("F15__p_O14", numba.float64),
    ("F16__n_F15", numba.float64),
    ("F16__p_O15", numba.float64),
    ("F17__n_F16", numba.float64),
    ("F17__p_O16", numba.float64),
    ("F18__n_F17", numba.float64),
    ("F18__p_O17", numba.float64),
    ("F18__He4_N14", numba.float64),
    ("F19__n_F18", numba.float64),
    ("F19__p_O18", numba.float64),
    ("F19__He4_N15", numba.float64),
    ("Ne18__p_F17", numba.float64),
    ("Ne18__He4_O14", numba.float64),
    ("Ne19__n_Ne18", numba.float64),
    ("Ne19__p_F18", numba.float64),
    ("Ne19__He4_O15", numba.float64),
    ("Ne20__n_Ne19", numba.float64),
    ("Ne20__p_F19", numba.float64),
    ("Ne20__He4_O16", numba.float64),
    ("Ne21__n_Ne20", numba.float64),
    ("Ne21__He4_O17", numba.float64),
    ("Na21__p_Ne20", numba.float64),
    ("Na21__He4_F17", numba.float64),
    ("Na22__n_Na21", numba.float64),
    ("Na22__p_Ne21", numba.float64),
    ("Na22__He4_F18", numba.float64),
    ("Na23__n_Na22", numba.float64),
    ("Na23__He4_F19", numba.float64),
    ("Na24__n_Na23", numba.float64),
    ("Mg22__p_Na21", numba.float64),
    ("Mg22__He4_Ne18", numba.float64),
    ("Mg23__n_Mg22", numba.float64),
    ("Mg23__p_Na22", numba.float64),
    ("Mg23__He4_Ne19", numba.float64),
    ("Mg24__n_Mg23", numba.float64),
    ("Mg24__p_Na23", numba.float64),
    ("Mg24__He4_Ne20", numba.float64),
    ("Mg25__n_Mg24", numba.float64),
    ("Mg25__p_Na24", numba.float64),
    ("Mg25__He4_Ne21", numba.float64),
    ("Mg26__n_Mg25", numba.float64),
    ("Al25__p_Mg24", numba.float64),
    ("Al25__He4_Na21", numba.float64),
    ("Al26__n_Al25", numba.float64),
    ("Al26__p_Mg25", numba.float64),
    ("Al26__He4_Na22", numba.float64),
    ("Al27__n_Al26", numba.float64),
    ("Al27__p_Mg26", numba.float64),
    ("Al27__He4_Na23", numba.float64),
    ("Al28__n_Al27", numba.float64),
    ("Al28__He4_Na24", numba.float64),
    ("Si26__p_Al25", numba.float64),
    ("Si26__He4_Mg22", numba.float64),
    ("Si28__p_Al27", numba.float64),
    ("Si28__He4_Mg24", numba.float64),
    ("Si29__n_Si28", numba.float64),
    ("Si29__p_Al28", numba.float64),
    ("Si29__He4_Mg25", numba.float64),
    ("Si30__n_Si29", numba.float64),
    ("Si30__He4_Mg26", numba.float64),
    ("P29__p_Si28", numba.float64),
    ("P29__He4_Al25", numba.float64),
    ("P30__n_P29", numba.float64),
    ("P30__p_Si29", numba.float64),
    ("P30__He4_Al26", numba.float64),
    ("P31__n_P30", numba.float64),
    ("P31__p_Si30", numba.float64),
    ("P31__He4_Al27", numba.float64),
    ("P32__n_P31", numba.float64),
    ("P32__He4_Al28", numba.float64),
    ("S31__p_P30", numba.float64),
    ("S32__n_S31", numba.float64),
    ("S32__p_P31", numba.float64),
    ("S32__He4_Si28", numba.float64),
    ("S34__He4_Si30", numba.float64),
    ("Cl32__p_S31", numba.float64),
    ("Cl32__p_P31__weak__wc12", numba.float64),
    ("Cl32__He4_Si28__weak__wc12", numba.float64),
    ("Cl34__He4_P30", numba.float64),
    ("Cl35__n_Cl34", numba.float64),
    ("Cl35__p_S34", numba.float64),
    ("Cl35__He4_P31", numba.float64),
    ("Ar35__p_Cl34", numba.float64),
    ("Ar35__He4_S31", numba.float64),
    ("He6__n_n_He4", numba.float64),
    ("Li6__n_p_He4", numba.float64),
    ("Be9__n_He4_He4", numba.float64),
    ("C12__He4_He4_He4", numba.float64),
    ("n_p__d", numba.float64),
    ("p_p__d__weak__bet_pos_", numba.float64),
    ("p_p__d__weak__electron_capture", numba.float64),
    ("n_d__t", numba.float64),
    ("p_d__He3", numba.float64),
    ("d_d__He4", numba.float64),
    ("He4_d__Li6", numba.float64),
    ("p_t__He4", numba.float64),
    ("He4_t__Li7", numba.float64),
    ("n_He3__He4", numba.float64),
    ("p_He3__He4__weak__bet_pos_", numba.float64),
    ("He4_He3__Be7", numba.float64),
    ("n_Li6__Li7", numba.float64),
    ("p_Li6__Be7", numba.float64),
    ("p_Be7__B8", numba.float64),
    ("He4_Be7__C11", numba.float64),
    ("n_C11__C12", numba.float64),
    ("n_C12__C13", numba.float64),
    ("p_C12__N13", numba.float64),
    ("He4_C12__O16", numba.float64),
    ("n_C13__C14", numba.float64),
    ("p_C13__N14", numba.float64),
    ("p_C14__N15", numba.float64),
    ("He4_C14__O18", numba.float64),
    ("n_N13__N14", numba.float64),
    ("p_N13__O14", numba.float64),
    ("n_N14__N15", numba.float64),
    ("p_N14__O15", numba.float64),
    ("He4_N14__F18", numba.float64),
    ("p_N15__O16", numba.float64),
    ("He4_N15__F19", numba.float64),
    ("n_O14__O15", numba.float64),
    ("p_O14__F15", numba.float64),
    ("He4_O14__Ne18", numba.float64),
    ("n_O15__O16", numba.float64),
    ("p_O15__F16", numba.float64),
    ("He4_O15__Ne19", numba.float64),
    ("n_O16__O17", numba.float64),
    ("p_O16__F17", numba.float64),
    ("He4_O16__Ne20", numba.float64),
    ("n_O17__O18", numba.float64),
    ("p_O17__F18", numba.float64),
    ("He4_O17__Ne21", numba.float64),
    ("p_O18__F19", numba.float64),
    ("n_F15__F16", numba.float64),
    ("n_F16__F17", numba.float64),
    ("n_F17__F18", numba.float64),
    ("p_F17__Ne18", numba.float64),
    ("He4_F17__Na21", numba.float64),
    ("n_F18__F19", numba.float64),
    ("p_F18__Ne19", numba.float64),
    ("He4_F18__Na22", numba.float64),
    ("p_F19__Ne20", numba.float64),
    ("He4_F19__Na23", numba.float64),
    ("n_Ne18__Ne19", numba.float64),
    ("He4_Ne18__Mg22", numba.float64),
    ("n_Ne19__Ne20", numba.float64),
    ("He4_Ne19__Mg23", numba.float64),
    ("n_Ne20__Ne21", numba.float64),
    ("p_Ne20__Na21", numba.float64),
    ("He4_Ne20__Mg24", numba.float64),
    ("p_Ne21__Na22", numba.float64),
    ("He4_Ne21__Mg25", numba.float64),
    ("n_Na21__Na22", numba.float64),
    ("p_Na21__Mg22", numba.float64),
    ("He4_Na21__Al25", numba.float64),
    ("n_Na22__Na23", numba.float64),
    ("p_Na22__Mg23", numba.float64),
    ("He4_Na22__Al26", numba.float64),
    ("n_Na23__Na24", numba.float64),
    ("p_Na23__Mg24", numba.float64),
    ("He4_Na23__Al27", numba.float64),
    ("p_Na24__Mg25", numba.float64),
    ("He4_Na24__Al28", numba.float64),
    ("n_Mg22__Mg23", numba.float64),
    ("He4_Mg22__Si26", numba.float64),
    ("n_Mg23__Mg24", numba.float64),
    ("n_Mg24__Mg25", numba.float64),
    ("p_Mg24__Al25", numba.float64),
    ("He4_Mg24__Si28", numba.float64),
    ("n_Mg25__Mg26", numba.float64),
    ("p_Mg25__Al26", numba.float64),
    ("He4_Mg25__Si29", numba.float64),
    ("p_Mg26__Al27", numba.float64),
    ("He4_Mg26__Si30", numba.float64),
    ("n_Al25__Al26", numba.float64),
    ("p_Al25__Si26", numba.float64),
    ("He4_Al25__P29", numba.float64),
    ("n_Al26__Al27", numba.float64),
    ("He4_Al26__P30", numba.float64),
    ("n_Al27__Al28", numba.float64),
    ("p_Al27__Si28", numba.float64),
    ("He4_Al27__P31", numba.float64),
    ("p_Al28__Si29", numba.float64),
    ("He4_Al28__P32", numba.float64),
    ("n_Si28__Si29", numba.float64),
    ("p_Si28__P29", numba.float64),
    ("He4_Si28__S32", numba.float64),
    ("n_Si29__Si30", numba.float64),
    ("p_Si29__P30", numba.float64),
    ("p_Si30__P31", numba.float64),
    ("He4_Si30__S34", numba.float64),
    ("n_P29__P30", numba.float64),
    ("n_P30__P31", numba.float64),
    ("p_P30__S31", numba.float64),
    ("He4_P30__Cl34", numba.float64),
    ("n_P31__P32", numba.float64),
    ("p_P31__S32", numba.float64),
    ("He4_P31__Cl35", numba.float64),
    ("n_S31__S32", numba.float64),
    ("p_S31__Cl32", numba.float64),
    ("He4_S31__Ar35", numba.float64),
    ("p_S34__Cl35", numba.float64),
    ("n_Cl34__Cl35", numba.float64),
    ("p_Cl34__Ar35", numba.float64),
    ("d_d__n_He3", numba.float64),
    ("d_d__p_t", numba.float64),
    ("p_t__n_He3", numba.float64),
    ("p_t__d_d", numba.float64),
    ("d_t__n_He4", numba.float64),
    ("He4_t__n_Li6", numba.float64),
    ("n_He3__p_t", numba.float64),
    ("n_He3__d_d", numba.float64),
    ("d_He3__p_He4", numba.float64),
    ("t_He3__d_He4", numba.float64),
    ("He4_He3__p_Li6", numba.float64),
    ("n_He4__d_t", numba.float64),
    ("p_He4__d_He3", numba.float64),
    ("d_He4__t_He3", numba.float64),
    ("He4_He4__n_Be7", numba.float64),
    ("He4_He4__p_Li7", numba.float64),
    ("n_Li6__He4_t", numba.float64),
    ("p_Li6__He4_He3", numba.float64),
    ("d_Li6__n_Be7", numba.float64),
    ("d_Li6__p_Li7", numba.float64),
    ("He4_Li6__p_Be9", numba.float64),
    ("p_Li7__n_Be7", numba.float64),
    ("p_Li7__d_Li6", numba.float64),
    ("p_Li7__He4_He4", numba.float64),
    ("t_Li7__n_Be9", numba.float64),
    ("n_Be7__p_Li7", numba.float64),
    ("n_Be7__d_Li6", numba.float64),
    ("n_Be7__He4_He4", numba.float64),
    ("n_Be9__t_Li7", numba.float64),
    ("p_Be9__He4_Li6", numba.float64),
    ("He4_Be9__n_C12", numba.float64),
    ("He4_B8__p_C11", numba.float64),
    ("p_C11__He4_B8", numba.float64),
    ("He4_C11__n_O14", numba.float64),
    ("He4_C11__p_N14", numba.float64),
    ("n_C12__He4_Be9", numba.float64),
    ("He4_C12__n_O15", numba.float64),
    ("He4_C12__p_N15", numba.float64),
    ("C12_C12__n_Mg23", numba.float64),
    ("C12_C12__p_Na23", numba.float64),
    ("C12_C12__He4_Ne20", numba.float64),
    ("p_C13__n_N13", numba.float64),
    ("d_C13__n_N14", numba.float64),
    ("He4_C13__n_O16", numba.float64),
    ("p_C14__n_N14", numba.float64),
    ("d_C14__n_N15", numba.float64),
    ("He4_C14__n_O17", numba.float64),
    ("n_N13__p_C13", numba.float64),
    ("He4_N13__n_F16", numba.float64),
    ("He4_N13__p_O16", numba.float64),
    ("n_N14__p_C14", numba.float64),
    ("n_N14__d_C13", numba.float64),
    ("p_N14__n_O14", numba.float64),
    ("p_N14__He4_C11", numba.float64),
    ("He4_N14__n_F17", numba.float64),
    ("He4_N14__p_O17", numba.float64),
    ("n_N15__d_C14", numba.float64),
    ("p_N15__n_O15", numba.float64),
    ("p_N15__He4_C12", numba.float64),
    ("He4_N15__n_F18", numba.float64),
    ("He4_N15__p_O18", numba.float64),
    ("n_O14__p_N14", numba.float64),
    ("n_O14__He4_C11", numba.float64),
    ("He4_O14__p_F17", numba.float64),
    ("n_O15__p_N15", numba.float64),
    ("n_O15__He4_C12", numba.float64),
    ("p_O15__n_F15", numba.float64),
    ("He4_O15__n_Ne18", numba.float64),
    ("He4_O15__p_F18", numba.float64),
    ("n_O16__He4_C13", numba.float64),
    ("p_O16__n_F16", numba.float64),
    ("p_O16__He4_N13", numba.float64),
    ("He4_O16__n_Ne19", numba.float64),
    ("He4_O16__p_F19", numba.float64),
    ("C12_O16__p_Al27", numba.float64),
    ("C12_O16__He4_Mg24", numba.float64),
    ("O16_O16__n_S31", numba.float64),
    ("O16_O16__p_P31", numba.float64),
    ("O16_O16__He4_Si28", numba.float64),
    ("n_O17__He4_C14", numba.float64),
    ("p_O17__n_F17", numba.float64),
    ("p_O17__He4_N14", numba.float64),
    ("He4_O17__n_Ne20", numba.float64),
    ("p_O18__n_F18", numba.float64),
    ("p_O18__He4_N15", numba.float64),
    ("He4_O18__n_Ne21", numba.float64),
    ("n_F15__p_O15", numba.float64),
    ("He4_F15__p_Ne18", numba.float64),
    ("n_F16__p_O16", numba.float64),
    ("n_F16__He4_N13", numba.float64),
    ("He4_F16__p_Ne19", numba.float64),
    ("n_F17__p_O17", numba.float64),
    ("n_F17__He4_N14", numba.float64),
    ("p_F17__He4_O14", numba.float64),
    ("He4_F17__p_Ne20", numba.float64),
    ("n_F18__p_O18", numba.float64),
    ("n_F18__He4_N15", numba.float64),
    ("p_F18__n_Ne18", numba.float64),
    ("p_F18__He4_O15", numba.float64),
    ("He4_F18__n_Na21", numba.float64),
    ("He4_F18__p_Ne21", numba.float64),
    ("p_F19__n_Ne19", numba.float64),
    ("p_F19__He4_O16", numba.float64),
    ("He4_F19__n_Na22", numba.float64),
    ("n_Ne18__p_F18", numba.float64),
    ("n_Ne18__He4_O15", numba.float64),
    ("p_Ne18__He4_F15", numba.float64),
    ("He4_Ne18__p_Na21", numba.float64),
    ("n_Ne19__p_F19", numba.float64),
    ("n_Ne19__He4_O16", numba.float64),
    ("p_Ne19__He4_F16", numba.float64),
    ("He4_Ne19__n_Mg22", numba.float64),
    ("He4_Ne19__p_Na22", numba.float64),
    ("n_Ne20__He4_O17", numba.float64),
    ("p_Ne20__He4_F17", numba.float64),
    ("He4_Ne20__n_Mg23", numba.float64),
    ("He4_Ne20__p_Na23", numba.float64),
    ("He4_Ne20__C12_C12", numba.float64),
    ("C12_Ne20__n_S31", numba.float64),
    ("C12_Ne20__p_P31", numba.float64),
    ("C12_Ne20__He4_Si28", numba.float64),
    ("n_Ne21__He4_O18", numba.float64),
    ("p_Ne21__n_Na21", numba.float64),
    ("p_Ne21__He4_F18", numba.float64),
    ("He4_Ne21__n_Mg24", numba.float64),
    ("He4_Ne21__p_Na24", numba.float64),
    ("n_Na21__p_Ne21", numba.float64),
    ("n_Na21__He4_F18", numba.float64),
    ("p_Na21__He4_Ne18", numba.float64),
    ("He4_Na21__p_Mg24", numba.float64),
    ("n_Na22__He4_F19", numba.float64),
    ("p_Na22__n_Mg22", numba.float64),
    ("p_Na22__He4_Ne19", numba.float64),
    ("He4_Na22__n_Al25", numba.float64),
    ("He4_Na22__p_Mg25", numba.float64),
    ("p_Na23__n_Mg23", numba.float64),
    ("p_Na23__He4_Ne20", numba.float64),
    ("p_Na23__C12_C12", numba.float64),
    ("He4_Na23__n_Al26", numba.float64),
    ("He4_Na23__p_Mg26", numba.float64),
    ("p_Na24__n_Mg24", numba.float64),
    ("p_Na24__He4_Ne21", numba.float64),
    ("He4_Na24__n_Al27", numba.float64),
    ("n_Mg22__p_Na22", numba.float64),
    ("n_Mg22__He4_Ne19", numba.float64),
    ("He4_Mg22__p_Al25", numba.float64),
    ("n_Mg23__p_Na23", numba.float64),
    ("n_Mg23__He4_Ne20", numba.float64),
    ("n_Mg23__C12_C12", numba.float64),
    ("He4_Mg23__n_Si26", numba.float64),
    ("He4_Mg23__p_Al26", numba.float64),
    ("n_Mg24__p_Na24", numba.float64),
    ("n_Mg24__He4_Ne21", numba.float64),
    ("p_Mg24__He4_Na21", numba.float64),
    ("He4_Mg24__p_Al27", numba.float64),
    ("He4_Mg24__C12_O16", numba.float64),
    ("p_Mg25__n_Al25", numba.float64),
    ("p_Mg25__He4_Na22", numba.float64),
    ("He4_Mg25__n_Si28", numba.float64),
    ("He4_Mg25__p_Al28", numba.float64),
    ("p_Mg26__n_Al26", numba.float64),
    ("p_Mg26__He4_Na23", numba.float64),
    ("He4_Mg26__n_Si29", numba.float64),
    ("n_Al25__p_Mg25", numba.float64),
    ("n_Al25__He4_Na22", numba.float64),
    ("p_Al25__He4_Mg22", numba.float64),
    ("He4_Al25__p_Si28", numba.float64),
    ("n_Al26__p_Mg26", numba.float64),
    ("n_Al26__He4_Na23", numba.float64),
    ("p_Al26__n_Si26", numba.float64),
    ("p_Al26__He4_Mg23", numba.float64),
    ("He4_Al26__n_P29", numba.float64),
    ("He4_Al26__p_Si29", numba.float64),
    ("n_Al27__He4_Na24", numba.float64),
    ("p_Al27__He4_Mg24", numba.float64),
    ("p_Al27__C12_O16", numba.float64),
    ("He4_Al27__n_P30", numba.float64),
    ("He4_Al27__p_Si30", numba.float64),
    ("p_Al28__n_Si28", numba.float64),
    ("p_Al28__He4_Mg25", numba.float64),
    ("He4_Al28__n_P31", numba.float64),
    ("n_Si26__p_Al26", numba.float64),
    ("n_Si26__He4_Mg23", numba.float64),
    ("He4_Si26__p_P29", numba.float64),
    ("n_Si28__p_Al28", numba.float64),
    ("n_Si28__He4_Mg25", numba.float64),
    ("p_Si28__He4_Al25", numba.float64),
    ("He4_Si28__n_S31", numba.float64),
    ("He4_Si28__p_P31", numba.float64),
    ("He4_Si28__C12_Ne20", numba.float64),
    ("He4_Si28__O16_O16", numba.float64),
    ("n_Si29__He4_Mg26", numba.float64),
    ("p_Si29__n_P29", numba.float64),
    ("p_Si29__He4_Al26", numba.float64),
    ("He4_Si29__n_S32", numba.float64),
    ("He4_Si29__p_P32", numba.float64),
    ("p_Si30__n_P30", numba.float64),
    ("p_Si30__He4_Al27", numba.float64),
    ("n_P29__p_Si29", numba.float64),
    ("n_P29__He4_Al26", numba.float64),
    ("p_P29__He4_Si26", numba.float64),
    ("He4_P29__n_Cl32", numba.float64),
    ("He4_P29__p_S32", numba.float64),
    ("n_P30__p_Si30", numba.float64),
    ("n_P30__He4_Al27", numba.float64),
    ("n_P31__He4_Al28", numba.float64),
    ("p_P31__n_S31", numba.float64),
    ("p_P31__He4_Si28", numba.float64),
    ("p_P31__C12_Ne20", numba.float64),
    ("p_P31__O16_O16", numba.float64),
    ("He4_P31__n_Cl34", numba.float64),
    ("He4_P31__p_S34", numba.float64),
    ("p_P32__n_S32", numba.float64),
    ("p_P32__He4_Si29", numba.float64),
    ("He4_P32__n_Cl35", numba.float64),
    ("n_S31__p_P31", numba.float64),
    ("n_S31__He4_Si28", numba.float64),
    ("n_S31__C12_Ne20", numba.float64),
    ("n_S31__O16_O16", numba.float64),
    ("He4_S31__p_Cl34", numba.float64),
    ("n_S32__p_P32", numba.float64),
    ("n_S32__He4_Si29", numba.float64),
    ("p_S32__n_Cl32", numba.float64),
    ("p_S32__He4_P29", numba.float64),
    ("He4_S32__n_Ar35", numba.float64),
    ("He4_S32__p_Cl35", numba.float64),
    ("p_S34__n_Cl34", numba.float64),
    ("p_S34__He4_P31", numba.float64),
    ("n_Cl32__p_S32", numba.float64),
    ("n_Cl32__He4_P29", numba.float64),
    ("He4_Cl32__p_Ar35", numba.float64),
    ("n_Cl34__p_S34", numba.float64),
    ("n_Cl34__He4_P31", numba.float64),
    ("p_Cl34__He4_S31", numba.float64),
    ("n_Cl35__He4_P32", numba.float64),
    ("p_Cl35__n_Ar35", numba.float64),
    ("p_Cl35__He4_S32", numba.float64),
    ("n_Ar35__p_Cl35", numba.float64),
    ("n_Ar35__He4_S32", numba.float64),
    ("p_Ar35__He4_Cl32", numba.float64),
    ("p_d__n_p_p", numba.float64),
    ("t_t__n_n_He4", numba.float64),
    ("t_He3__n_p_He4", numba.float64),
    ("He3_He3__p_p_He4", numba.float64),
    ("d_Li7__n_He4_He4", numba.float64),
    ("d_Be7__p_He4_He4", numba.float64),
    ("p_Be9__d_He4_He4", numba.float64),
    ("n_B8__p_He4_He4", numba.float64),
    ("n_C11__He4_He4_He4", numba.float64),
    ("t_Li7__n_n_He4_He4", numba.float64),
    ("He3_Li7__n_p_He4_He4", numba.float64),
    ("t_Be7__n_p_He4_He4", numba.float64),
    ("He3_Be7__p_p_He4_He4", numba.float64),
    ("p_Be9__n_p_He4_He4", numba.float64),
    ("n_n_He4__He6", numba.float64),
    ("n_p_He4__Li6", numba.float64),
    ("n_He4_He4__Be9", numba.float64),
    ("He4_He4_He4__C12", numba.float64),
    ("n_p_p__p_d", numba.float64),
    ("n_n_He4__t_t", numba.float64),
    ("n_p_He4__t_He3", numba.float64),
    ("p_p_He4__He3_He3", numba.float64),
    ("n_He4_He4__d_Li7", numba.float64),
    ("p_He4_He4__n_B8", numba.float64),
    ("p_He4_He4__d_Be7", numba.float64),
    ("d_He4_He4__p_Be9", numba.float64),
    ("He4_He4_He4__n_C11", numba.float64),
    ("n_n_He4_He4__t_Li7", numba.float64),
    ("n_p_He4_He4__He3_Li7", numba.float64),
    ("n_p_He4_He4__t_Be7", numba.float64),
    ("n_p_He4_He4__p_Be9", numba.float64),
    ("p_p_He4_He4__He3_Be7", numba.float64),
])
class RateEval:
    def __init__(self):
        self.n__p__weak__wc12 = np.nan
        self.t__He3__weak__wc12 = np.nan
        self.He3__t__weak__electron_capture = np.nan
        self.He6__Li6__weak__wc12 = np.nan
        self.Be7__Li7__weak__electron_capture = np.nan
        self.C14__N14__weak__wc12 = np.nan
        self.N13__C13__weak__wc12 = np.nan
        self.O14__N14__weak__wc12 = np.nan
        self.O15__N15__weak__wc12 = np.nan
        self.F17__O17__weak__wc12 = np.nan
        self.F18__O18__weak__wc12 = np.nan
        self.Ne18__F18__weak__wc12 = np.nan
        self.Ne19__F19__weak__wc12 = np.nan
        self.Na21__Ne21__weak__wc12 = np.nan
        self.Na24__Mg24__weak__wc12 = np.nan
        self.Mg22__Na22__weak__wc12 = np.nan
        self.Mg23__Na23__weak__wc12 = np.nan
        self.Al25__Mg25__weak__wc12 = np.nan
        self.Al26__Mg26__weak__wc12 = np.nan
        self.Al28__Si28__weak__wc12 = np.nan
        self.Si26__Al26__weak__wc12 = np.nan
        self.P29__Si29__weak__wc12 = np.nan
        self.P30__Si30__weak__wc12 = np.nan
        self.P32__S32__weak__wc12 = np.nan
        self.S31__P31__weak__wc12 = np.nan
        self.Cl32__S32__weak__wc12 = np.nan
        self.Cl34__S34__weak__wc12 = np.nan
        self.Ar35__Cl35__weak__wc12 = np.nan
        self.d__n_p = np.nan
        self.t__n_d = np.nan
        self.He3__p_d = np.nan
        self.He4__n_He3 = np.nan
        self.He4__p_t = np.nan
        self.He4__d_d = np.nan
        self.Li6__He4_d = np.nan
        self.Li7__n_Li6 = np.nan
        self.Li7__He4_t = np.nan
        self.Be7__p_Li6 = np.nan
        self.Be7__He4_He3 = np.nan
        self.B8__p_Be7 = np.nan
        self.B8__He4_He4__weak__wc12 = np.nan
        self.C11__He4_Be7 = np.nan
        self.C12__n_C11 = np.nan
        self.C13__n_C12 = np.nan
        self.C14__n_C13 = np.nan
        self.N13__p_C12 = np.nan
        self.N14__n_N13 = np.nan
        self.N14__p_C13 = np.nan
        self.N15__n_N14 = np.nan
        self.N15__p_C14 = np.nan
        self.O14__p_N13 = np.nan
        self.O15__n_O14 = np.nan
        self.O15__p_N14 = np.nan
        self.O16__n_O15 = np.nan
        self.O16__p_N15 = np.nan
        self.O16__He4_C12 = np.nan
        self.O17__n_O16 = np.nan
        self.O18__n_O17 = np.nan
        self.O18__He4_C14 = np.nan
        self.F15__p_O14 = np.nan
        self.F16__n_F15 = np.nan
        self.F16__p_O15 = np.nan
        self.F17__n_F16 = np.nan
        self.F17__p_O16 = np.nan
        self.F18__n_F17 = np.nan
        self.F18__p_O17 = np.nan
        self.F18__He4_N14 = np.nan
        self.F19__n_F18 = np.nan
        self.F19__p_O18 = np.nan
        self.F19__He4_N15 = np.nan
        self.Ne18__p_F17 = np.nan
        self.Ne18__He4_O14 = np.nan
        self.Ne19__n_Ne18 = np.nan
        self.Ne19__p_F18 = np.nan
        self.Ne19__He4_O15 = np.nan
        self.Ne20__n_Ne19 = np.nan
        self.Ne20__p_F19 = np.nan
        self.Ne20__He4_O16 = np.nan
        self.Ne21__n_Ne20 = np.nan
        self.Ne21__He4_O17 = np.nan
        self.Na21__p_Ne20 = np.nan
        self.Na21__He4_F17 = np.nan
        self.Na22__n_Na21 = np.nan
        self.Na22__p_Ne21 = np.nan
        self.Na22__He4_F18 = np.nan
        self.Na23__n_Na22 = np.nan
        self.Na23__He4_F19 = np.nan
        self.Na24__n_Na23 = np.nan
        self.Mg22__p_Na21 = np.nan
        self.Mg22__He4_Ne18 = np.nan
        self.Mg23__n_Mg22 = np.nan
        self.Mg23__p_Na22 = np.nan
        self.Mg23__He4_Ne19 = np.nan
        self.Mg24__n_Mg23 = np.nan
        self.Mg24__p_Na23 = np.nan
        self.Mg24__He4_Ne20 = np.nan
        self.Mg25__n_Mg24 = np.nan
        self.Mg25__p_Na24 = np.nan
        self.Mg25__He4_Ne21 = np.nan
        self.Mg26__n_Mg25 = np.nan
        self.Al25__p_Mg24 = np.nan
        self.Al25__He4_Na21 = np.nan
        self.Al26__n_Al25 = np.nan
        self.Al26__p_Mg25 = np.nan
        self.Al26__He4_Na22 = np.nan
        self.Al27__n_Al26 = np.nan
        self.Al27__p_Mg26 = np.nan
        self.Al27__He4_Na23 = np.nan
        self.Al28__n_Al27 = np.nan
        self.Al28__He4_Na24 = np.nan
        self.Si26__p_Al25 = np.nan
        self.Si26__He4_Mg22 = np.nan
        self.Si28__p_Al27 = np.nan
        self.Si28__He4_Mg24 = np.nan
        self.Si29__n_Si28 = np.nan
        self.Si29__p_Al28 = np.nan
        self.Si29__He4_Mg25 = np.nan
        self.Si30__n_Si29 = np.nan
        self.Si30__He4_Mg26 = np.nan
        self.P29__p_Si28 = np.nan
        self.P29__He4_Al25 = np.nan
        self.P30__n_P29 = np.nan
        self.P30__p_Si29 = np.nan
        self.P30__He4_Al26 = np.nan
        self.P31__n_P30 = np.nan
        self.P31__p_Si30 = np.nan
        self.P31__He4_Al27 = np.nan
        self.P32__n_P31 = np.nan
        self.P32__He4_Al28 = np.nan
        self.S31__p_P30 = np.nan
        self.S32__n_S31 = np.nan
        self.S32__p_P31 = np.nan
        self.S32__He4_Si28 = np.nan
        self.S34__He4_Si30 = np.nan
        self.Cl32__p_S31 = np.nan
        self.Cl32__p_P31__weak__wc12 = np.nan
        self.Cl32__He4_Si28__weak__wc12 = np.nan
        self.Cl34__He4_P30 = np.nan
        self.Cl35__n_Cl34 = np.nan
        self.Cl35__p_S34 = np.nan
        self.Cl35__He4_P31 = np.nan
        self.Ar35__p_Cl34 = np.nan
        self.Ar35__He4_S31 = np.nan
        self.He6__n_n_He4 = np.nan
        self.Li6__n_p_He4 = np.nan
        self.Be9__n_He4_He4 = np.nan
        self.C12__He4_He4_He4 = np.nan
        self.n_p__d = np.nan
        self.p_p__d__weak__bet_pos_ = np.nan
        self.p_p__d__weak__electron_capture = np.nan
        self.n_d__t = np.nan
        self.p_d__He3 = np.nan
        self.d_d__He4 = np.nan
        self.He4_d__Li6 = np.nan
        self.p_t__He4 = np.nan
        self.He4_t__Li7 = np.nan
        self.n_He3__He4 = np.nan
        self.p_He3__He4__weak__bet_pos_ = np.nan
        self.He4_He3__Be7 = np.nan
        self.n_Li6__Li7 = np.nan
        self.p_Li6__Be7 = np.nan
        self.p_Be7__B8 = np.nan
        self.He4_Be7__C11 = np.nan
        self.n_C11__C12 = np.nan
        self.n_C12__C13 = np.nan
        self.p_C12__N13 = np.nan
        self.He4_C12__O16 = np.nan
        self.n_C13__C14 = np.nan
        self.p_C13__N14 = np.nan
        self.p_C14__N15 = np.nan
        self.He4_C14__O18 = np.nan
        self.n_N13__N14 = np.nan
        self.p_N13__O14 = np.nan
        self.n_N14__N15 = np.nan
        self.p_N14__O15 = np.nan
        self.He4_N14__F18 = np.nan
        self.p_N15__O16 = np.nan
        self.He4_N15__F19 = np.nan
        self.n_O14__O15 = np.nan
        self.p_O14__F15 = np.nan
        self.He4_O14__Ne18 = np.nan
        self.n_O15__O16 = np.nan
        self.p_O15__F16 = np.nan
        self.He4_O15__Ne19 = np.nan
        self.n_O16__O17 = np.nan
        self.p_O16__F17 = np.nan
        self.He4_O16__Ne20 = np.nan
        self.n_O17__O18 = np.nan
        self.p_O17__F18 = np.nan
        self.He4_O17__Ne21 = np.nan
        self.p_O18__F19 = np.nan
        self.n_F15__F16 = np.nan
        self.n_F16__F17 = np.nan
        self.n_F17__F18 = np.nan
        self.p_F17__Ne18 = np.nan
        self.He4_F17__Na21 = np.nan
        self.n_F18__F19 = np.nan
        self.p_F18__Ne19 = np.nan
        self.He4_F18__Na22 = np.nan
        self.p_F19__Ne20 = np.nan
        self.He4_F19__Na23 = np.nan
        self.n_Ne18__Ne19 = np.nan
        self.He4_Ne18__Mg22 = np.nan
        self.n_Ne19__Ne20 = np.nan
        self.He4_Ne19__Mg23 = np.nan
        self.n_Ne20__Ne21 = np.nan
        self.p_Ne20__Na21 = np.nan
        self.He4_Ne20__Mg24 = np.nan
        self.p_Ne21__Na22 = np.nan
        self.He4_Ne21__Mg25 = np.nan
        self.n_Na21__Na22 = np.nan
        self.p_Na21__Mg22 = np.nan
        self.He4_Na21__Al25 = np.nan
        self.n_Na22__Na23 = np.nan
        self.p_Na22__Mg23 = np.nan
        self.He4_Na22__Al26 = np.nan
        self.n_Na23__Na24 = np.nan
        self.p_Na23__Mg24 = np.nan
        self.He4_Na23__Al27 = np.nan
        self.p_Na24__Mg25 = np.nan
        self.He4_Na24__Al28 = np.nan
        self.n_Mg22__Mg23 = np.nan
        self.He4_Mg22__Si26 = np.nan
        self.n_Mg23__Mg24 = np.nan
        self.n_Mg24__Mg25 = np.nan
        self.p_Mg24__Al25 = np.nan
        self.He4_Mg24__Si28 = np.nan
        self.n_Mg25__Mg26 = np.nan
        self.p_Mg25__Al26 = np.nan
        self.He4_Mg25__Si29 = np.nan
        self.p_Mg26__Al27 = np.nan
        self.He4_Mg26__Si30 = np.nan
        self.n_Al25__Al26 = np.nan
        self.p_Al25__Si26 = np.nan
        self.He4_Al25__P29 = np.nan
        self.n_Al26__Al27 = np.nan
        self.He4_Al26__P30 = np.nan
        self.n_Al27__Al28 = np.nan
        self.p_Al27__Si28 = np.nan
        self.He4_Al27__P31 = np.nan
        self.p_Al28__Si29 = np.nan
        self.He4_Al28__P32 = np.nan
        self.n_Si28__Si29 = np.nan
        self.p_Si28__P29 = np.nan
        self.He4_Si28__S32 = np.nan
        self.n_Si29__Si30 = np.nan
        self.p_Si29__P30 = np.nan
        self.p_Si30__P31 = np.nan
        self.He4_Si30__S34 = np.nan
        self.n_P29__P30 = np.nan
        self.n_P30__P31 = np.nan
        self.p_P30__S31 = np.nan
        self.He4_P30__Cl34 = np.nan
        self.n_P31__P32 = np.nan
        self.p_P31__S32 = np.nan
        self.He4_P31__Cl35 = np.nan
        self.n_S31__S32 = np.nan
        self.p_S31__Cl32 = np.nan
        self.He4_S31__Ar35 = np.nan
        self.p_S34__Cl35 = np.nan
        self.n_Cl34__Cl35 = np.nan
        self.p_Cl34__Ar35 = np.nan
        self.d_d__n_He3 = np.nan
        self.d_d__p_t = np.nan
        self.p_t__n_He3 = np.nan
        self.p_t__d_d = np.nan
        self.d_t__n_He4 = np.nan
        self.He4_t__n_Li6 = np.nan
        self.n_He3__p_t = np.nan
        self.n_He3__d_d = np.nan
        self.d_He3__p_He4 = np.nan
        self.t_He3__d_He4 = np.nan
        self.He4_He3__p_Li6 = np.nan
        self.n_He4__d_t = np.nan
        self.p_He4__d_He3 = np.nan
        self.d_He4__t_He3 = np.nan
        self.He4_He4__n_Be7 = np.nan
        self.He4_He4__p_Li7 = np.nan
        self.n_Li6__He4_t = np.nan
        self.p_Li6__He4_He3 = np.nan
        self.d_Li6__n_Be7 = np.nan
        self.d_Li6__p_Li7 = np.nan
        self.He4_Li6__p_Be9 = np.nan
        self.p_Li7__n_Be7 = np.nan
        self.p_Li7__d_Li6 = np.nan
        self.p_Li7__He4_He4 = np.nan
        self.t_Li7__n_Be9 = np.nan
        self.n_Be7__p_Li7 = np.nan
        self.n_Be7__d_Li6 = np.nan
        self.n_Be7__He4_He4 = np.nan
        self.n_Be9__t_Li7 = np.nan
        self.p_Be9__He4_Li6 = np.nan
        self.He4_Be9__n_C12 = np.nan
        self.He4_B8__p_C11 = np.nan
        self.p_C11__He4_B8 = np.nan
        self.He4_C11__n_O14 = np.nan
        self.He4_C11__p_N14 = np.nan
        self.n_C12__He4_Be9 = np.nan
        self.He4_C12__n_O15 = np.nan
        self.He4_C12__p_N15 = np.nan
        self.C12_C12__n_Mg23 = np.nan
        self.C12_C12__p_Na23 = np.nan
        self.C12_C12__He4_Ne20 = np.nan
        self.p_C13__n_N13 = np.nan
        self.d_C13__n_N14 = np.nan
        self.He4_C13__n_O16 = np.nan
        self.p_C14__n_N14 = np.nan
        self.d_C14__n_N15 = np.nan
        self.He4_C14__n_O17 = np.nan
        self.n_N13__p_C13 = np.nan
        self.He4_N13__n_F16 = np.nan
        self.He4_N13__p_O16 = np.nan
        self.n_N14__p_C14 = np.nan
        self.n_N14__d_C13 = np.nan
        self.p_N14__n_O14 = np.nan
        self.p_N14__He4_C11 = np.nan
        self.He4_N14__n_F17 = np.nan
        self.He4_N14__p_O17 = np.nan
        self.n_N15__d_C14 = np.nan
        self.p_N15__n_O15 = np.nan
        self.p_N15__He4_C12 = np.nan
        self.He4_N15__n_F18 = np.nan
        self.He4_N15__p_O18 = np.nan
        self.n_O14__p_N14 = np.nan
        self.n_O14__He4_C11 = np.nan
        self.He4_O14__p_F17 = np.nan
        self.n_O15__p_N15 = np.nan
        self.n_O15__He4_C12 = np.nan
        self.p_O15__n_F15 = np.nan
        self.He4_O15__n_Ne18 = np.nan
        self.He4_O15__p_F18 = np.nan
        self.n_O16__He4_C13 = np.nan
        self.p_O16__n_F16 = np.nan
        self.p_O16__He4_N13 = np.nan
        self.He4_O16__n_Ne19 = np.nan
        self.He4_O16__p_F19 = np.nan
        self.C12_O16__p_Al27 = np.nan
        self.C12_O16__He4_Mg24 = np.nan
        self.O16_O16__n_S31 = np.nan
        self.O16_O16__p_P31 = np.nan
        self.O16_O16__He4_Si28 = np.nan
        self.n_O17__He4_C14 = np.nan
        self.p_O17__n_F17 = np.nan
        self.p_O17__He4_N14 = np.nan
        self.He4_O17__n_Ne20 = np.nan
        self.p_O18__n_F18 = np.nan
        self.p_O18__He4_N15 = np.nan
        self.He4_O18__n_Ne21 = np.nan
        self.n_F15__p_O15 = np.nan
        self.He4_F15__p_Ne18 = np.nan
        self.n_F16__p_O16 = np.nan
        self.n_F16__He4_N13 = np.nan
        self.He4_F16__p_Ne19 = np.nan
        self.n_F17__p_O17 = np.nan
        self.n_F17__He4_N14 = np.nan
        self.p_F17__He4_O14 = np.nan
        self.He4_F17__p_Ne20 = np.nan
        self.n_F18__p_O18 = np.nan
        self.n_F18__He4_N15 = np.nan
        self.p_F18__n_Ne18 = np.nan
        self.p_F18__He4_O15 = np.nan
        self.He4_F18__n_Na21 = np.nan
        self.He4_F18__p_Ne21 = np.nan
        self.p_F19__n_Ne19 = np.nan
        self.p_F19__He4_O16 = np.nan
        self.He4_F19__n_Na22 = np.nan
        self.n_Ne18__p_F18 = np.nan
        self.n_Ne18__He4_O15 = np.nan
        self.p_Ne18__He4_F15 = np.nan
        self.He4_Ne18__p_Na21 = np.nan
        self.n_Ne19__p_F19 = np.nan
        self.n_Ne19__He4_O16 = np.nan
        self.p_Ne19__He4_F16 = np.nan
        self.He4_Ne19__n_Mg22 = np.nan
        self.He4_Ne19__p_Na22 = np.nan
        self.n_Ne20__He4_O17 = np.nan
        self.p_Ne20__He4_F17 = np.nan
        self.He4_Ne20__n_Mg23 = np.nan
        self.He4_Ne20__p_Na23 = np.nan
        self.He4_Ne20__C12_C12 = np.nan
        self.C12_Ne20__n_S31 = np.nan
        self.C12_Ne20__p_P31 = np.nan
        self.C12_Ne20__He4_Si28 = np.nan
        self.n_Ne21__He4_O18 = np.nan
        self.p_Ne21__n_Na21 = np.nan
        self.p_Ne21__He4_F18 = np.nan
        self.He4_Ne21__n_Mg24 = np.nan
        self.He4_Ne21__p_Na24 = np.nan
        self.n_Na21__p_Ne21 = np.nan
        self.n_Na21__He4_F18 = np.nan
        self.p_Na21__He4_Ne18 = np.nan
        self.He4_Na21__p_Mg24 = np.nan
        self.n_Na22__He4_F19 = np.nan
        self.p_Na22__n_Mg22 = np.nan
        self.p_Na22__He4_Ne19 = np.nan
        self.He4_Na22__n_Al25 = np.nan
        self.He4_Na22__p_Mg25 = np.nan
        self.p_Na23__n_Mg23 = np.nan
        self.p_Na23__He4_Ne20 = np.nan
        self.p_Na23__C12_C12 = np.nan
        self.He4_Na23__n_Al26 = np.nan
        self.He4_Na23__p_Mg26 = np.nan
        self.p_Na24__n_Mg24 = np.nan
        self.p_Na24__He4_Ne21 = np.nan
        self.He4_Na24__n_Al27 = np.nan
        self.n_Mg22__p_Na22 = np.nan
        self.n_Mg22__He4_Ne19 = np.nan
        self.He4_Mg22__p_Al25 = np.nan
        self.n_Mg23__p_Na23 = np.nan
        self.n_Mg23__He4_Ne20 = np.nan
        self.n_Mg23__C12_C12 = np.nan
        self.He4_Mg23__n_Si26 = np.nan
        self.He4_Mg23__p_Al26 = np.nan
        self.n_Mg24__p_Na24 = np.nan
        self.n_Mg24__He4_Ne21 = np.nan
        self.p_Mg24__He4_Na21 = np.nan
        self.He4_Mg24__p_Al27 = np.nan
        self.He4_Mg24__C12_O16 = np.nan
        self.p_Mg25__n_Al25 = np.nan
        self.p_Mg25__He4_Na22 = np.nan
        self.He4_Mg25__n_Si28 = np.nan
        self.He4_Mg25__p_Al28 = np.nan
        self.p_Mg26__n_Al26 = np.nan
        self.p_Mg26__He4_Na23 = np.nan
        self.He4_Mg26__n_Si29 = np.nan
        self.n_Al25__p_Mg25 = np.nan
        self.n_Al25__He4_Na22 = np.nan
        self.p_Al25__He4_Mg22 = np.nan
        self.He4_Al25__p_Si28 = np.nan
        self.n_Al26__p_Mg26 = np.nan
        self.n_Al26__He4_Na23 = np.nan
        self.p_Al26__n_Si26 = np.nan
        self.p_Al26__He4_Mg23 = np.nan
        self.He4_Al26__n_P29 = np.nan
        self.He4_Al26__p_Si29 = np.nan
        self.n_Al27__He4_Na24 = np.nan
        self.p_Al27__He4_Mg24 = np.nan
        self.p_Al27__C12_O16 = np.nan
        self.He4_Al27__n_P30 = np.nan
        self.He4_Al27__p_Si30 = np.nan
        self.p_Al28__n_Si28 = np.nan
        self.p_Al28__He4_Mg25 = np.nan
        self.He4_Al28__n_P31 = np.nan
        self.n_Si26__p_Al26 = np.nan
        self.n_Si26__He4_Mg23 = np.nan
        self.He4_Si26__p_P29 = np.nan
        self.n_Si28__p_Al28 = np.nan
        self.n_Si28__He4_Mg25 = np.nan
        self.p_Si28__He4_Al25 = np.nan
        self.He4_Si28__n_S31 = np.nan
        self.He4_Si28__p_P31 = np.nan
        self.He4_Si28__C12_Ne20 = np.nan
        self.He4_Si28__O16_O16 = np.nan
        self.n_Si29__He4_Mg26 = np.nan
        self.p_Si29__n_P29 = np.nan
        self.p_Si29__He4_Al26 = np.nan
        self.He4_Si29__n_S32 = np.nan
        self.He4_Si29__p_P32 = np.nan
        self.p_Si30__n_P30 = np.nan
        self.p_Si30__He4_Al27 = np.nan
        self.n_P29__p_Si29 = np.nan
        self.n_P29__He4_Al26 = np.nan
        self.p_P29__He4_Si26 = np.nan
        self.He4_P29__n_Cl32 = np.nan
        self.He4_P29__p_S32 = np.nan
        self.n_P30__p_Si30 = np.nan
        self.n_P30__He4_Al27 = np.nan
        self.n_P31__He4_Al28 = np.nan
        self.p_P31__n_S31 = np.nan
        self.p_P31__He4_Si28 = np.nan
        self.p_P31__C12_Ne20 = np.nan
        self.p_P31__O16_O16 = np.nan
        self.He4_P31__n_Cl34 = np.nan
        self.He4_P31__p_S34 = np.nan
        self.p_P32__n_S32 = np.nan
        self.p_P32__He4_Si29 = np.nan
        self.He4_P32__n_Cl35 = np.nan
        self.n_S31__p_P31 = np.nan
        self.n_S31__He4_Si28 = np.nan
        self.n_S31__C12_Ne20 = np.nan
        self.n_S31__O16_O16 = np.nan
        self.He4_S31__p_Cl34 = np.nan
        self.n_S32__p_P32 = np.nan
        self.n_S32__He4_Si29 = np.nan
        self.p_S32__n_Cl32 = np.nan
        self.p_S32__He4_P29 = np.nan
        self.He4_S32__n_Ar35 = np.nan
        self.He4_S32__p_Cl35 = np.nan
        self.p_S34__n_Cl34 = np.nan
        self.p_S34__He4_P31 = np.nan
        self.n_Cl32__p_S32 = np.nan
        self.n_Cl32__He4_P29 = np.nan
        self.He4_Cl32__p_Ar35 = np.nan
        self.n_Cl34__p_S34 = np.nan
        self.n_Cl34__He4_P31 = np.nan
        self.p_Cl34__He4_S31 = np.nan
        self.n_Cl35__He4_P32 = np.nan
        self.p_Cl35__n_Ar35 = np.nan
        self.p_Cl35__He4_S32 = np.nan
        self.n_Ar35__p_Cl35 = np.nan
        self.n_Ar35__He4_S32 = np.nan
        self.p_Ar35__He4_Cl32 = np.nan
        self.p_d__n_p_p = np.nan
        self.t_t__n_n_He4 = np.nan
        self.t_He3__n_p_He4 = np.nan
        self.He3_He3__p_p_He4 = np.nan
        self.d_Li7__n_He4_He4 = np.nan
        self.d_Be7__p_He4_He4 = np.nan
        self.p_Be9__d_He4_He4 = np.nan
        self.n_B8__p_He4_He4 = np.nan
        self.n_C11__He4_He4_He4 = np.nan
        self.t_Li7__n_n_He4_He4 = np.nan
        self.He3_Li7__n_p_He4_He4 = np.nan
        self.t_Be7__n_p_He4_He4 = np.nan
        self.He3_Be7__p_p_He4_He4 = np.nan
        self.p_Be9__n_p_He4_He4 = np.nan
        self.n_n_He4__He6 = np.nan
        self.n_p_He4__Li6 = np.nan
        self.n_He4_He4__Be9 = np.nan
        self.He4_He4_He4__C12 = np.nan
        self.n_p_p__p_d = np.nan
        self.n_n_He4__t_t = np.nan
        self.n_p_He4__t_He3 = np.nan
        self.p_p_He4__He3_He3 = np.nan
        self.n_He4_He4__d_Li7 = np.nan
        self.p_He4_He4__n_B8 = np.nan
        self.p_He4_He4__d_Be7 = np.nan
        self.d_He4_He4__p_Be9 = np.nan
        self.He4_He4_He4__n_C11 = np.nan
        self.n_n_He4_He4__t_Li7 = np.nan
        self.n_p_He4_He4__He3_Li7 = np.nan
        self.n_p_He4_He4__t_Be7 = np.nan
        self.n_p_He4_He4__p_Be9 = np.nan
        self.p_p_He4_He4__He3_Be7 = np.nan

@numba.njit()
def ye(Y):
    return np.sum(Z * Y)/np.sum(A * Y)

@numba.njit()
def n__p__weak__wc12(rate_eval, tf):
    # n --> p
    rate = 0.0

    # wc12w
    rate += np.exp(  -6.78161)

    rate_eval.n__p__weak__wc12 = rate

@numba.njit()
def t__He3__weak__wc12(rate_eval, tf):
    # t --> He3
    rate = 0.0

    # wc12w
    rate += np.exp(  -20.1456)

    rate_eval.t__He3__weak__wc12 = rate

@numba.njit()
def He3__t__weak__electron_capture(rate_eval, tf):
    # He3 --> t
    rate = 0.0

    #   ecw
    rate += np.exp(  -32.462 + -0.21338*tf.T9i + -0.821581*tf.T913i + 11.1241*tf.T913
                  + -0.577338*tf.T9 + 0.0290471*tf.T953 + -0.262705*tf.lnT9)

    rate_eval.He3__t__weak__electron_capture = rate

@numba.njit()
def He6__Li6__weak__wc12(rate_eval, tf):
    # He6 --> Li6
    rate = 0.0

    # wc12w
    rate += np.exp(  -0.144619)

    rate_eval.He6__Li6__weak__wc12 = rate

@numba.njit()
def Be7__Li7__weak__electron_capture(rate_eval, tf):
    # Be7 --> Li7
    rate = 0.0

    #   ecw
    rate += np.exp(  -23.8328 + 3.02033*tf.T913
                  + -0.0742132*tf.T9 + -0.00792386*tf.T953 + -0.650113*tf.lnT9)

    rate_eval.Be7__Li7__weak__electron_capture = rate

@numba.njit()
def C14__N14__weak__wc12(rate_eval, tf):
    # C14 --> N14
    rate = 0.0

    # wc12w
    rate += np.exp(  -26.2827)

    rate_eval.C14__N14__weak__wc12 = rate

@numba.njit()
def N13__C13__weak__wc12(rate_eval, tf):
    # N13 --> C13
    rate = 0.0

    # wc12w
    rate += np.exp(  -6.7601)

    rate_eval.N13__C13__weak__wc12 = rate

@numba.njit()
def O14__N14__weak__wc12(rate_eval, tf):
    # O14 --> N14
    rate = 0.0

    # wc12w
    rate += np.exp(  -4.62354)

    rate_eval.O14__N14__weak__wc12 = rate

@numba.njit()
def O15__N15__weak__wc12(rate_eval, tf):
    # O15 --> N15
    rate = 0.0

    # wc12w
    rate += np.exp(  -5.17053)

    rate_eval.O15__N15__weak__wc12 = rate

@numba.njit()
def F17__O17__weak__wc12(rate_eval, tf):
    # F17 --> O17
    rate = 0.0

    # wc12w
    rate += np.exp(  -4.53318)

    rate_eval.F17__O17__weak__wc12 = rate

@numba.njit()
def F18__O18__weak__wc12(rate_eval, tf):
    # F18 --> O18
    rate = 0.0

    # wc12w
    rate += np.exp(  -9.15982)

    rate_eval.F18__O18__weak__wc12 = rate

@numba.njit()
def Ne18__F18__weak__wc12(rate_eval, tf):
    # Ne18 --> F18
    rate = 0.0

    # wc12w
    rate += np.exp(  -0.879336)

    rate_eval.Ne18__F18__weak__wc12 = rate

@numba.njit()
def Ne19__F19__weak__wc12(rate_eval, tf):
    # Ne19 --> F19
    rate = 0.0

    # wc12w
    rate += np.exp(  -3.21142)

    rate_eval.Ne19__F19__weak__wc12 = rate

@numba.njit()
def Na21__Ne21__weak__wc12(rate_eval, tf):
    # Na21 --> Ne21
    rate = 0.0

    # wc12w
    rate += np.exp(  -3.48003)

    rate_eval.Na21__Ne21__weak__wc12 = rate

@numba.njit()
def Na24__Mg24__weak__wc12(rate_eval, tf):
    # Na24 --> Mg24
    rate = 0.0

    # wc12w
    rate += np.exp(  -11.2633)

    rate_eval.Na24__Mg24__weak__wc12 = rate

@numba.njit()
def Mg22__Na22__weak__wc12(rate_eval, tf):
    # Mg22 --> Na22
    rate = 0.0

    # wc12w
    rate += np.exp(  -1.72235)

    rate_eval.Mg22__Na22__weak__wc12 = rate

@numba.njit()
def Mg23__Na23__weak__wc12(rate_eval, tf):
    # Mg23 --> Na23
    rate = 0.0

    # wc12w
    rate += np.exp(  -2.79132)

    rate_eval.Mg23__Na23__weak__wc12 = rate

@numba.njit()
def Al25__Mg25__weak__wc12(rate_eval, tf):
    # Al25 --> Mg25
    rate = 0.0

    # wc12w
    rate += np.exp(  -2.33781)

    rate_eval.Al25__Mg25__weak__wc12 = rate

@numba.njit()
def Al26__Mg26__weak__wc12(rate_eval, tf):
    # Al26 --> Mg26
    rate = 0.0

    # wc12w
    rate += np.exp(  -4.62175 + -2.64931*tf.T9i + -0.025978*tf.T913
                  + -0.0291284*tf.T9 + 0.00389774*tf.T953)

    rate_eval.Al26__Mg26__weak__wc12 = rate

@numba.njit()
def Al28__Si28__weak__wc12(rate_eval, tf):
    # Al28 --> Si28
    rate = 0.0

    # wc12w
    rate += np.exp(  -5.26435)

    rate_eval.Al28__Si28__weak__wc12 = rate

@numba.njit()
def Si26__Al26__weak__wc12(rate_eval, tf):
    # Si26 --> Al26
    rate = 0.0

    # wc12w
    rate += np.exp(  -1.16851)

    rate_eval.Si26__Al26__weak__wc12 = rate

@numba.njit()
def P29__Si29__weak__wc12(rate_eval, tf):
    # P29 --> Si29
    rate = 0.0

    # wc12w
    rate += np.exp(  -1.78721)

    rate_eval.P29__Si29__weak__wc12 = rate

@numba.njit()
def P30__Si30__weak__wc12(rate_eval, tf):
    # P30 --> Si30
    rate = 0.0

    # wc12w
    rate += np.exp(  -5.37715)

    rate_eval.P30__Si30__weak__wc12 = rate

@numba.njit()
def P32__S32__weak__wc12(rate_eval, tf):
    # P32 --> S32
    rate = 0.0

    # wc12w
    rate += np.exp(  -14.389)

    rate_eval.P32__S32__weak__wc12 = rate

@numba.njit()
def S31__P31__weak__wc12(rate_eval, tf):
    # S31 --> P31
    rate = 0.0

    # wc12w
    rate += np.exp(  -1.31042)

    rate_eval.S31__P31__weak__wc12 = rate

@numba.njit()
def Cl32__S32__weak__wc12(rate_eval, tf):
    # Cl32 --> S32
    rate = 0.0

    # wc12w
    rate += np.exp(  0.843349)

    rate_eval.Cl32__S32__weak__wc12 = rate

@numba.njit()
def Cl34__S34__weak__wc12(rate_eval, tf):
    # Cl34 --> S34
    rate = 0.0

    # wc12w
    rate += np.exp(  -0.791781)

    rate_eval.Cl34__S34__weak__wc12 = rate

@numba.njit()
def Ar35__Cl35__weak__wc12(rate_eval, tf):
    # Ar35 --> Cl35
    rate = 0.0

    # wc12w
    rate += np.exp(  -0.943126)

    rate_eval.Ar35__Cl35__weak__wc12 = rate

@numba.njit()
def d__n_p(rate_eval, tf):
    # d --> n + p
    rate = 0.0

    # an06n
    rate += np.exp(  34.6293 + -25.815*tf.T9i + -2.70618*tf.T913
                  + 0.11718*tf.T9 + -0.00312788*tf.T953 + 1.96913*tf.lnT9)
    # an06n
    rate += np.exp(  31.1075 + -25.815*tf.T9i + -0.0102082*tf.T913
                  + -0.0893959*tf.T9 + 0.00696704*tf.T953 + 2.5*tf.lnT9)
    # an06n
    rate += np.exp(  33.0154 + -25.815*tf.T9i + -2.30472*tf.T913
                  + -0.887862*tf.T9 + 0.137663*tf.T953 + 1.5*tf.lnT9)

    rate_eval.d__n_p = rate

@numba.njit()
def t__n_d(rate_eval, tf):
    # t --> n + d
    rate = 0.0

    # nk06n
    rate += np.exp(  30.1124 + -72.6136*tf.T9i
                  + 2.5*tf.lnT9)
    # nk06n
    rate += np.exp(  28.869 + -72.6136*tf.T9i
                  + 1.575*tf.lnT9)

    rate_eval.t__n_d = rate

@numba.njit()
def He3__p_d(rate_eval, tf):
    # He3 --> p + d
    rate = 0.0

    # de04 
    rate += np.exp(  32.4383 + -63.7435*tf.T9i + -3.7208*tf.T913i + 0.198654*tf.T913
                  + 1.83333*tf.lnT9)
    # de04n
    rate += np.exp(  31.032 + -63.7435*tf.T9i + -3.7208*tf.T913i + 0.871782*tf.T913
                  + 0.833333*tf.lnT9)

    rate_eval.He3__p_d = rate

@numba.njit()
def He4__n_He3(rate_eval, tf):
    # He4 --> n + He3
    rate = 0.0

    # ka02n
    rate += np.exp(  33.0131 + -238.79*tf.T9i + -1.50147*tf.T913
                  + 2.5*tf.lnT9)
    # ka02n
    rate += np.exp(  29.4845 + -238.79*tf.T9i
                  + 1.5*tf.lnT9)

    rate_eval.He4__n_He3 = rate

@numba.njit()
def He4__p_t(rate_eval, tf):
    # He4 --> p + t
    rate = 0.0

    # cf88n
    rate += np.exp(  33.7327 + -229.932*tf.T9i + -3.869*tf.T913i + 1.45482*tf.T913
                  + 0.577246*tf.T9 + -0.112199*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.He4__p_t = rate

@numba.njit()
def He4__d_d(rate_eval, tf):
    # He4 --> d + d
    rate = 0.0

    # nacrn
    rate += np.exp(  28.2984 + -276.744*tf.T9i + -4.26166*tf.T913i + -0.119233*tf.T913
                  + 0.778829*tf.T9 + -0.0925203*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.He4__d_d = rate

@numba.njit()
def Li6__He4_d(rate_eval, tf):
    # Li6 --> He4 + d
    rate = 0.0

    # tu19r
    rate += np.exp(  27.5672 + -24.9919*tf.T9i)
    # tu19n
    rate += np.exp(  22.7676 + -17.1028*tf.T9i + -7.55198*tf.T913i + 5.77546*tf.T913
                  + -0.487854*tf.T9 + 0.032833*tf.T953 + 0.376948*tf.lnT9)

    rate_eval.Li6__He4_d = rate

@numba.njit()
def Li7__n_Li6(rate_eval, tf):
    # Li7 --> n + Li6
    rate = 0.0

    # jz10n
    rate += np.exp(  32.2347 + -84.1369*tf.T9i
                  + 1.5*tf.lnT9)

    rate_eval.Li7__n_Li6 = rate

@numba.njit()
def Li7__He4_t(rate_eval, tf):
    # Li7 --> He4 + t
    rate = 0.0

    # de04 
    rate += np.exp(  36.7442 + -28.6283*tf.T9i + -8.0805*tf.T913i + -0.217514*tf.T913
                  + -0.114859*tf.T9 + 0.0470043*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Li7__He4_t = rate

@numba.njit()
def Be7__p_Li6(rate_eval, tf):
    # Be7 --> p + Li6
    rate = 0.0

    # nacrn
    rate += np.exp(  37.4661 + -65.0548*tf.T9i + -8.4372*tf.T913i + -0.515473*tf.T913
                  + 0.0285578*tf.T9 + 0.00879731*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Be7__p_Li6 = rate

@numba.njit()
def Be7__He4_He3(rate_eval, tf):
    # Be7 --> He4 + He3
    rate = 0.0

    # cd08n
    rate += np.exp(  38.7379 + -18.4059*tf.T9i + -12.8271*tf.T913i + -0.0308225*tf.T913
                  + -0.654685*tf.T9 + 0.0896331*tf.T953 + 0.833333*tf.lnT9)
    # cd08n
    rate += np.exp(  40.8355 + -18.4059*tf.T9i + -12.8271*tf.T913i + -3.8126*tf.T913
                  + 0.0942285*tf.T9 + -0.00301018*tf.T953 + 2.83333*tf.lnT9)

    rate_eval.Be7__He4_He3 = rate

@numba.njit()
def B8__p_Be7(rate_eval, tf):
    # B8 --> p + Be7
    rate = 0.0

    # nacrr
    rate += np.exp(  31.0163 + -8.93482*tf.T9i)
    # nacrn
    rate += np.exp(  35.8138 + -1.58982*tf.T9i + -10.264*tf.T913i + -0.203472*tf.T913
                  + 0.121083*tf.T9 + -0.00700063*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.B8__p_Be7 = rate

@numba.njit()
def B8__He4_He4__weak__wc12(rate_eval, tf):
    # B8 --> He4 + He4
    rate = 0.0

    # wc12w
    rate += np.exp(  -0.105148)

    rate_eval.B8__He4_He4__weak__wc12 = rate

@numba.njit()
def C11__He4_Be7(rate_eval, tf):
    # C11 --> He4 + Be7
    rate = 0.0

    # nacrr
    rate += np.exp(  33.8476 + -94.0424*tf.T9i)
    # nacrn
    rate += np.exp(  48.6414 + -87.5443*tf.T9i + -23.214*tf.T913i + -3.74943*tf.T913
                  + 1.23242*tf.T9 + -0.195665*tf.T953 + 0.833333*tf.lnT9)
    # nacrr
    rate += np.exp(  36.2917 + -97.7214*tf.T9i + 0.685881*tf.T913
                  + -0.697071*tf.T9 + 0.13274*tf.T953)

    rate_eval.C11__He4_Be7 = rate

@numba.njit()
def C12__n_C11(rate_eval, tf):
    # C12 --> n + C11
    rate = 0.0

    # bb92r
    rate += np.exp(  33.0631 + -218.179*tf.T9i + -5.21741e-12*tf.T913i + 2.32794e-11*tf.T913
                  + -2.03119e-12*tf.T9 + 7.09452e-14*tf.T953)
    # bb92n
    rate += np.exp(  35.3287 + -217.262*tf.T9i + 1.18875e-10*tf.T913i + -3.10977e-10*tf.T913
                  + 2.0423e-11*tf.T9 + -6.59293e-13*tf.T953 + 1.5*tf.lnT9)
    # bb92r
    rate += np.exp(  38.8258 + -222.832*tf.T9i + 1.38808e-09*tf.T913i + -3.66292e-09*tf.T913
                  + 2.41904e-10*tf.T9 + -7.82223e-12*tf.T953)

    rate_eval.C12__n_C11 = rate

@numba.njit()
def C13__n_C12(rate_eval, tf):
    # C13 --> n + C12
    rate = 0.0

    # ks03 
    rate += np.exp(  30.8808 + -57.4077*tf.T9i + 1.49573*tf.T913i + -0.841102*tf.T913
                  + 0.0340543*tf.T9 + -0.0026392*tf.T953 + 3.1662*tf.lnT9)

    rate_eval.C13__n_C12 = rate

@numba.njit()
def C14__n_C13(rate_eval, tf):
    # C14 --> n + C13
    rate = 0.0

    # ks03 
    rate += np.exp(  59.5926 + -95.0156*tf.T9i + 18.3578*tf.T913i + -46.5786*tf.T913
                  + 2.58472*tf.T9 + -0.118622*tf.T953 + 21.4142*tf.lnT9)

    rate_eval.C14__n_C13 = rate

@numba.njit()
def N13__p_C12(rate_eval, tf):
    # N13 --> p + C12
    rate = 0.0

    # ls09r
    rate += np.exp(  40.4354 + -26.326*tf.T9i + -5.10735*tf.T913i + -2.24111*tf.T913
                  + 0.148883*tf.T9)
    # ls09n
    rate += np.exp(  40.0408 + -22.5475*tf.T9i + -13.692*tf.T913i + -0.230881*tf.T913
                  + 4.44362*tf.T9 + -3.15898*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.N13__p_C12 = rate

@numba.njit()
def N14__n_N13(rate_eval, tf):
    # N14 --> n + N13
    rate = 0.0

    # wiesr
    rate += np.exp(  19.5584 + -125.474*tf.T9i + 9.44873e-10*tf.T913i + -2.33713e-09*tf.T913
                  + 1.97507e-10*tf.T9 + -1.49747e-11*tf.T953)
    # wiesn
    rate += np.exp(  37.1268 + -122.484*tf.T9i + 1.72241e-10*tf.T913i + -5.62522e-10*tf.T913
                  + 5.59212e-11*tf.T9 + -4.6549e-12*tf.T953 + 1.5*tf.lnT9)

    rate_eval.N14__n_N13 = rate

@numba.njit()
def N14__p_C13(rate_eval, tf):
    # N14 --> p + C13
    rate = 0.0

    # nacrr
    rate += np.exp(  37.1528 + -93.4071*tf.T9i + -0.196703*tf.T913
                  + 0.142126*tf.T9 + -0.0238912*tf.T953)
    # nacrr
    rate += np.exp(  38.3716 + -101.18*tf.T9i)
    # nacrn
    rate += np.exp(  41.7046 + -87.6256*tf.T9i + -13.72*tf.T913i + -0.450018*tf.T913
                  + 3.70823*tf.T9 + -1.70545*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.N14__p_C13 = rate

@numba.njit()
def N15__n_N14(rate_eval, tf):
    # N15 --> n + N14
    rate = 0.0

    # ks03 
    rate += np.exp(  34.1728 + -125.726*tf.T9i + 1.396*tf.T913i + -3.47552*tf.T913
                  + 0.351773*tf.T9 + -0.0229344*tf.T953 + 2.52161*tf.lnT9)

    rate_eval.N15__n_N14 = rate

@numba.njit()
def N15__p_C14(rate_eval, tf):
    # N15 --> p + C14
    rate = 0.0

    # il10r
    rate += np.exp(  40.0115 + -119.975*tf.T9i + -10.658*tf.T913i + 1.73644*tf.T913
                  + -0.350498*tf.T9 + 0.0279902*tf.T953)
    # il10n
    rate += np.exp(  43.0281 + -118.452*tf.T9i + -13.9619*tf.T913i + -4.34315*tf.T913
                  + 6.64922*tf.T9 + -3.22592*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.N15__p_C14 = rate

@numba.njit()
def O14__p_N13(rate_eval, tf):
    # O14 --> p + N13
    rate = 0.0

    # lg06r
    rate += np.exp(  35.2849 + -59.8313*tf.T9i + 1.57122*tf.T913i)
    # lg06n
    rate += np.exp(  42.4234 + -53.7053*tf.T9i + -15.1676*tf.T913i + 0.0955166*tf.T913
                  + 3.0659*tf.T9 + -0.507339*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.O14__p_N13 = rate

@numba.njit()
def O15__n_O14(rate_eval, tf):
    # O15 --> n + O14
    rate = 0.0

    # rpsmr
    rate += np.exp(  32.7811 + -153.419*tf.T9i + -1.38986*tf.T913i + 1.74662*tf.T913
                  + -0.0276897*tf.T9 + 0.00321014*tf.T953 + 0.438778*tf.lnT9)

    rate_eval.O15__n_O14 = rate

@numba.njit()
def O15__p_N14(rate_eval, tf):
    # O15 --> p + N14
    rate = 0.0

    # im05r
    rate += np.exp(  30.7435 + -89.5667*tf.T9i
                  + 1.5682*tf.lnT9)
    # im05r
    rate += np.exp(  31.6622 + -87.6737*tf.T9i)
    # im05n
    rate += np.exp(  44.1246 + -84.6757*tf.T9i + -15.193*tf.T913i + -4.63975*tf.T913
                  + 9.73458*tf.T9 + -9.55051*tf.T953 + 1.83333*tf.lnT9)
    # im05n
    rate += np.exp(  41.0177 + -84.6757*tf.T9i + -15.193*tf.T913i + -0.161954*tf.T913
                  + -7.52123*tf.T9 + -0.987565*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.O15__p_N14 = rate

@numba.njit()
def O16__n_O15(rate_eval, tf):
    # O16 --> n + O15
    rate = 0.0

    # rpsmr
    rate += np.exp(  32.3869 + -181.759*tf.T9i + -1.11761*tf.T913i + 1.0167*tf.T913
                  + 0.0449976*tf.T9 + -0.00204682*tf.T953 + 0.710783*tf.lnT9)

    rate_eval.O16__n_O15 = rate

@numba.njit()
def O16__p_N15(rate_eval, tf):
    # O16 --> p + N15
    rate = 0.0

    # li10r
    rate += np.exp(  38.8465 + -150.962*tf.T9i
                  + 0.0459037*tf.T9)
    # li10r
    rate += np.exp(  30.8927 + -143.656*tf.T9i)
    # li10n
    rate += np.exp(  44.3197 + -140.732*tf.T9i + -15.24*tf.T913i + 0.334926*tf.T913
                  + 4.59088*tf.T9 + -4.78468*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.O16__p_N15 = rate

@numba.njit()
def O16__He4_C12(rate_eval, tf):
    # O16 --> He4 + C12
    rate = 0.0

    # nac2 
    rate += np.exp(  94.3131 + -84.503*tf.T9i + 58.9128*tf.T913i + -148.273*tf.T913
                  + 9.08324*tf.T9 + -0.541041*tf.T953 + 71.8554*tf.lnT9)
    # nac2 
    rate += np.exp(  279.295 + -84.9515*tf.T9i + 103.411*tf.T913i + -420.567*tf.T913
                  + 64.0874*tf.T9 + -12.4624*tf.T953 + 138.803*tf.lnT9)

    rate_eval.O16__He4_C12 = rate

@numba.njit()
def O17__n_O16(rate_eval, tf):
    # O17 --> n + O16
    rate = 0.0

    # ks03 
    rate += np.exp(  29.0385 + -48.0574*tf.T9i + -2.11246*tf.T913i + 4.87742*tf.T913
                  + -0.314426*tf.T9 + 0.0169515*tf.T953 + 0.515216*tf.lnT9)

    rate_eval.O17__n_O16 = rate

@numba.njit()
def O18__n_O17(rate_eval, tf):
    # O18 --> n + O17
    rate = 0.0

    # bb92n
    rate += np.exp(  29.682 + -93.3549*tf.T9i + 7.48144e-11*tf.T913i + -2.47239e-10*tf.T913
                  + 2.48052e-11*tf.T9 + -2.07736e-12*tf.T953 + 1.5*tf.lnT9)
    # bb92r
    rate += np.exp(  45.5915 + -95.4845*tf.T9i + 16.8052*tf.T913i + -30.138*tf.T913
                  + 1.14711*tf.T9 + -0.0220312*tf.T953 + 15.2452*tf.lnT9)

    rate_eval.O18__n_O17 = rate

@numba.njit()
def O18__He4_C14(rate_eval, tf):
    # O18 --> He4 + C14
    rate = 0.0

    # il10r
    rate += np.exp(  36.546 + -82.6514*tf.T9i + -3.83188*tf.T913
                  + 1.64358*tf.T9 + -0.177785*tf.T953)
    # il10r
    rate += np.exp(  0.910093 + -74.3219*tf.T9i)
    # il10n
    rate += np.exp(  43.2028 + -72.2531*tf.T9i + -31.7222*tf.T913i + 11.3923*tf.T913
                  + -9.92249*tf.T9 + -2.0*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.O18__He4_C14 = rate

@numba.njit()
def F15__p_O14(rate_eval, tf):
    # F15 --> p + O14
    rate = 0.0

    # rpsmr
    rate += np.exp(  26.6052 + -2.37708*tf.T9i + 71.6134*tf.T913i + -76.9654*tf.T913
                  + 3.54377*tf.T9 + -0.159804*tf.T953 + 47.8922*tf.lnT9)

    rate_eval.F15__p_O14 = rate

@numba.njit()
def F16__n_F15(rate_eval, tf):
    # F16 --> n + F15
    rate = 0.0

    # rpsmr
    rate += np.exp(  32.8445 + -164.391*tf.T9i + -1.8316*tf.T913i + 2.8333*tf.T913
                  + -0.0925718*tf.T9 + 0.00805423*tf.T953 + 0.0590139*tf.lnT9)

    rate_eval.F16__n_F15 = rate

@numba.njit()
def F16__p_O15(rate_eval, tf):
    # F16 --> p + O15
    rate = 0.0

    # rpsmr
    rate += np.exp(  21.2899 + -2.29538*tf.T9i + 74.5955*tf.T913i + -73.794*tf.T913
                  + 2.87936*tf.T9 + -0.117988*tf.T953 + 49.0366*tf.lnT9)

    rate_eval.F16__p_O15 = rate

@numba.njit()
def F17__n_F16(rate_eval, tf):
    # F17 --> n + F16
    rate = 0.0

    # rpsmr
    rate += np.exp(  31.1876 + -194.947*tf.T9i + -0.79526*tf.T913i + 1.84363*tf.T913
                  + -0.123684*tf.T9 + 0.0149635*tf.T953 + 0.795305*tf.lnT9)

    rate_eval.F17__n_F16 = rate

@numba.njit()
def F17__p_O16(rate_eval, tf):
    # F17 --> p + O16
    rate = 0.0

    # ia08n
    rate += np.exp(  40.9135 + -6.96583*tf.T9i + -16.696*tf.T913i + -1.16252*tf.T913
                  + 0.267703*tf.T9 + -0.0338411*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.F17__p_O16 = rate

@numba.njit()
def F18__n_F17(rate_eval, tf):
    # F18 --> n + F17
    rate = 0.0

    # rpsmr
    rate += np.exp(  27.7411 + -106.154*tf.T9i + -2.71353*tf.T913i + 6.19871*tf.T913
                  + -0.26836*tf.T9 + 0.0115489*tf.T953 + -1.02002*tf.lnT9)

    rate_eval.F18__n_F17 = rate

@numba.njit()
def F18__p_O17(rate_eval, tf):
    # F18 --> p + O17
    rate = 0.0

    # il10r
    rate += np.exp(  33.7037 + -71.2889*tf.T9i + 2.31435*tf.T913
                  + -0.302835*tf.T9 + 0.020133*tf.T953)
    # il10r
    rate += np.exp(  11.2362 + -65.8069*tf.T9i)
    # il10n
    rate += np.exp(  40.2061 + -65.0606*tf.T9i + -16.4035*tf.T913i + 4.31885*tf.T913
                  + -0.709921*tf.T9 + -2.0*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.F18__p_O17 = rate

@numba.njit()
def F18__He4_N14(rate_eval, tf):
    # F18 --> He4 + N14
    rate = 0.0

    # il10r
    rate += np.exp(  38.6146 + -62.1948*tf.T9i + -5.6227*tf.T913i)
    # il10r
    rate += np.exp(  24.9119 + -56.3896*tf.T9i)
    # il10n
    rate += np.exp(  46.249 + -51.2292*tf.T9i + -36.2504*tf.T913i
                  + -5.0*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.F18__He4_N14 = rate

@numba.njit()
def F19__n_F18(rate_eval, tf):
    # F19 --> n + F18
    rate = 0.0

    # rpsmr
    rate += np.exp(  24.8451 + -120.998*tf.T9i + -6.65988*tf.T913i + 18.3243*tf.T913
                  + -1.47263*tf.T9 + 0.0955082*tf.T953 + -5.25505*tf.lnT9)

    rate_eval.F19__n_F18 = rate

@numba.njit()
def F19__p_O18(rate_eval, tf):
    # F19 --> p + O18
    rate = 0.0

    # il10r
    rate += np.exp(  30.2003 + -99.501*tf.T9i + 3.99059*tf.T913
                  + -0.593127*tf.T9 + 0.0877534*tf.T953)
    # il10r
    rate += np.exp(  28.008 + -94.4325*tf.T9i)
    # il10r
    rate += np.exp(  -12.0764 + -93.0204*tf.T9i)
    # il10n
    rate += np.exp(  42.8485 + -92.7757*tf.T9i + -16.7246*tf.T913i
                  + -3.0*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.F19__p_O18 = rate

@numba.njit()
def F19__He4_N15(rate_eval, tf):
    # F19 --> He4 + N15
    rate = 0.0

    # il10r
    rate += np.exp(  -4.06142 + -50.7773*tf.T9i + 35.4292*tf.T913
                  + -5.5767*tf.T9 + 0.441293*tf.T953)
    # il10r
    rate += np.exp(  28.2717 + -53.5621*tf.T9i)
    # il10r
    rate += np.exp(  15.3186 + -50.7554*tf.T9i)
    # il10n
    rate += np.exp(  50.1291 + -46.5774*tf.T9i + -36.2324*tf.T913i
                  + -2.0*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.F19__He4_N15 = rate

@numba.njit()
def Ne18__p_F17(rate_eval, tf):
    # Ne18 --> p + F17
    rate = 0.0

    # cb09 
    rate += np.exp(  52.9895 + -50.492*tf.T9i + -21.3249*tf.T913i + -0.230774*tf.T913
                  + 0.917931*tf.T9 + -0.0440377*tf.T953 + -5.86014*tf.lnT9)
    # cb09 
    rate += np.exp(  17.5646 + -45.5647*tf.T9i + -14.2191*tf.T913i + 34.0647*tf.T913
                  + -16.5698*tf.T9 + 2.48116*tf.T953 + -0.63376*tf.lnT9)

    rate_eval.Ne18__p_F17 = rate

@numba.njit()
def Ne18__He4_O14(rate_eval, tf):
    # Ne18 --> He4 + O14
    rate = 0.0

    # wh87r
    rate += np.exp(  20.0156 + -71.5044*tf.T9i
                  + 6.5*tf.lnT9)
    # wh87r
    rate += np.exp(  28.2415 + -81.9554*tf.T9i)
    # wh87r
    rate += np.exp(  22.5609 + -71.0754*tf.T9i)
    # wh87n
    rate += np.exp(  51.158 + -59.3454*tf.T9i + -39.38*tf.T913i + -0.0772187*tf.T913
                  + -0.635361*tf.T9 + 0.106236*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Ne18__He4_O14 = rate

@numba.njit()
def Ne19__n_Ne18(rate_eval, tf):
    # Ne19 --> n + Ne18
    rate = 0.0

    # ths8r
    rate += np.exp(  32.1805 + -135.042*tf.T9i + 0.937162*tf.T913
                  + -0.0221952*tf.T9 + -0.00101206*tf.T953 + 1.5*tf.lnT9)

    rate_eval.Ne19__n_Ne18 = rate

@numba.njit()
def Ne19__p_F18(rate_eval, tf):
    # Ne19 --> p + F18
    rate = 0.0

    # il10n
    rate += np.exp(  81.4385 + -74.3988*tf.T9i + -21.4023*tf.T913i + -93.766*tf.T913
                  + 179.258*tf.T9 + -202.561*tf.T953 + 0.833333*tf.lnT9)
    # il10r
    rate += np.exp(  18.1729 + -77.2902*tf.T9i + 13.1683*tf.T913
                  + -1.92023*tf.T9 + 0.16901*tf.T953)
    # il10r
    rate += np.exp(  -5.41887 + -74.7977*tf.T9i + 22.4903*tf.T913
                  + 0.307872*tf.T9 + -0.296226*tf.T953)

    rate_eval.Ne19__p_F18 = rate

@numba.njit()
def Ne19__He4_O15(rate_eval, tf):
    # Ne19 --> He4 + O15
    rate = 0.0

    # dc11r
    rate += np.exp(  -7.51212 + -45.1578*tf.T9i + -3.24609*tf.T913i + 44.4647*tf.T913
                  + -9.79962*tf.T9 + 0.841782*tf.T953)
    # dc11r
    rate += np.exp(  24.6922 + -46.8378*tf.T9i)
    # dc11n
    rate += np.exp(  51.0289 + -40.9534*tf.T9i + -39.578*tf.T913i
                  + -3.0*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Ne19__He4_O15 = rate

@numba.njit()
def Ne20__n_Ne19(rate_eval, tf):
    # Ne20 --> n + Ne19
    rate = 0.0

    # ths8r
    rate += np.exp(  30.7283 + -195.706*tf.T9i + 1.57592*tf.T913
                  + -0.11175*tf.T9 + 0.00226473*tf.T953 + 1.5*tf.lnT9)

    rate_eval.Ne20__n_Ne19 = rate

@numba.njit()
def Ne20__p_F19(rate_eval, tf):
    # Ne20 --> p + F19
    rate = 0.0

    # nacrr
    rate += np.exp(  18.691 + -156.781*tf.T9i + 31.6442*tf.T913i + -58.6563*tf.T913
                  + 67.7365*tf.T9 + -22.9721*tf.T953)
    # nacrr
    rate += np.exp(  36.7036 + -150.75*tf.T9i + -11.3832*tf.T913i + 5.47872*tf.T913
                  + -1.07203*tf.T9 + 0.11196*tf.T953)
    # nacrn
    rate += np.exp(  42.6027 + -149.037*tf.T9i + -18.116*tf.T913i + -1.4622*tf.T913
                  + 6.95113*tf.T9 + -2.90366*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Ne20__p_F19 = rate

@numba.njit()
def Ne20__He4_O16(rate_eval, tf):
    # Ne20 --> He4 + O16
    rate = 0.0

    # co10r
    rate += np.exp(  28.6431 + -65.246*tf.T9i)
    # co10n
    rate += np.exp(  48.6604 + -54.8875*tf.T9i + -39.7262*tf.T913i + -0.210799*tf.T913
                  + 0.442879*tf.T9 + -0.0797753*tf.T953 + 0.833333*tf.lnT9)
    # co10r
    rate += np.exp(  34.2658 + -67.6518*tf.T9i + -3.65925*tf.T913
                  + 0.714224*tf.T9 + -0.00107508*tf.T953)

    rate_eval.Ne20__He4_O16 = rate

@numba.njit()
def Ne21__n_Ne20(rate_eval, tf):
    # Ne21 --> n + Ne20
    rate = 0.0

    # ka02n
    rate += np.exp(  30.8228 + -78.458*tf.T9i
                  + 1.5*tf.lnT9)
    # ka02r
    rate += np.exp(  34.9807 + -80.162*tf.T9i)

    rate_eval.Ne21__n_Ne20 = rate

@numba.njit()
def Ne21__He4_O17(rate_eval, tf):
    # Ne21 --> He4 + O17
    rate = 0.0

    # be13r
    rate += np.exp(  0.0906657 + -90.782*tf.T9i + 123.363*tf.T913i + -87.4351*tf.T913
                  + -3.40974e-06*tf.T9 + -57.0469*tf.T953 + 83.7218*tf.lnT9)
    # be13r
    rate += np.exp(  -91.954 + -98.9487*tf.T9i + 3.31162e-08*tf.T913i + 130.258*tf.T913
                  + -7.92551e-05*tf.T9 + -4.13772*tf.T953 + -41.2753*tf.lnT9)
    # be13r
    rate += np.exp(  27.3205 + -91.2722*tf.T9i + 2.87641*tf.T913i + -3.54489*tf.T913
                  + -2.11222e-08*tf.T9 + -3.90649e-09*tf.T953 + 6.25778*tf.lnT9)

    rate_eval.Ne21__He4_O17 = rate

@numba.njit()
def Na21__p_Ne20(rate_eval, tf):
    # Na21 --> p + Ne20
    rate = 0.0

    # ly18 
    rate += np.exp(  230.123 + -28.3722*tf.T9i + 15.325*tf.T913i + -294.859*tf.T913
                  + 107.692*tf.T9 + -46.2072*tf.T953 + 59.3398*tf.lnT9)
    # ly18 
    rate += np.exp(  28.0772 + -37.0575*tf.T9i + 20.5893*tf.T913i + -17.5841*tf.T913
                  + 0.243226*tf.T9 + -0.000231418*tf.T953 + 14.3398*tf.lnT9)
    # ly18 
    rate += np.exp(  252.265 + -32.6731*tf.T9i + 258.57*tf.T913i + -506.387*tf.T913
                  + 22.1576*tf.T9 + -0.721182*tf.T953 + 231.788*tf.lnT9)
    # ly18 
    rate += np.exp(  195320.0 + -89.3596*tf.T9i + 21894.7*tf.T913i + -319153.0*tf.T913
                  + 224369.0*tf.T9 + -188049.0*tf.T953 + 48704.9*tf.lnT9)

    rate_eval.Na21__p_Ne20 = rate

@numba.njit()
def Na21__He4_F17(rate_eval, tf):
    # Na21 --> He4 + F17
    rate = 0.0

    # rpsmr
    rate += np.exp(  66.3334 + -77.8653*tf.T9i + 15.559*tf.T913i + -68.3231*tf.T913
                  + 2.54275*tf.T9 + -0.0989207*tf.T953 + 38.3877*tf.lnT9)

    rate_eval.Na21__He4_F17 = rate

@numba.njit()
def Na22__n_Na21(rate_eval, tf):
    # Na22 --> n + Na21
    rate = 0.0

    # ths8r
    rate += np.exp(  30.6596 + -128.457*tf.T9i + 1.31089*tf.T913
                  + -0.164931*tf.T9 + 0.00903374*tf.T953 + 1.5*tf.lnT9)

    rate_eval.Na22__n_Na21 = rate

@numba.njit()
def Na22__p_Ne21(rate_eval, tf):
    # Na22 --> p + Ne21
    rate = 0.0

    # il10r
    rate += np.exp(  24.8334 + -79.6093*tf.T9i)
    # il10r
    rate += np.exp(  -24.579 + -78.4059*tf.T9i)
    # il10n
    rate += np.exp(  42.146 + -78.2097*tf.T9i + -19.2096*tf.T913i
                  + -1.0*tf.T953 + 0.833333*tf.lnT9)
    # il10r
    rate += np.exp(  -16.4098 + -82.4235*tf.T9i + 21.1176*tf.T913i + 34.0411*tf.T913
                  + -4.45593*tf.T9 + 0.328613*tf.T953)

    rate_eval.Na22__p_Ne21 = rate

@numba.njit()
def Na22__He4_F18(rate_eval, tf):
    # Na22 --> He4 + F18
    rate = 0.0

    # rpsmr
    rate += np.exp(  59.3224 + -100.236*tf.T9i + 18.8956*tf.T913i + -65.6134*tf.T913
                  + 1.71114*tf.T9 + -0.0260999*tf.T953 + 39.3396*tf.lnT9)

    rate_eval.Na22__He4_F18 = rate

@numba.njit()
def Na23__n_Na22(rate_eval, tf):
    # Na23 --> n + Na22
    rate = 0.0

    # ths8r
    rate += np.exp(  37.0665 + -144.113*tf.T9i + 1.02148*tf.T913
                  + -0.334638*tf.T9 + 0.0258708*tf.T953 + 1.5*tf.lnT9)

    rate_eval.Na23__n_Na22 = rate

@numba.njit()
def Na23__He4_F19(rate_eval, tf):
    # Na23 --> He4 + F19
    rate = 0.0

    # rpsmr
    rate += np.exp(  76.8979 + -123.578*tf.T9i + 39.7219*tf.T913i + -100.401*tf.T913
                  + 3.15808*tf.T9 + -0.0629822*tf.T953 + 55.9823*tf.lnT9)

    rate_eval.Na23__He4_F19 = rate

@numba.njit()
def Na24__n_Na23(rate_eval, tf):
    # Na24 --> n + Na23
    rate = 0.0

    # ks03 
    rate += np.exp(  44.327 + -80.9446*tf.T9i + 16.0285*tf.T913i + -25.7668*tf.T913
                  + 1.07385*tf.T9 + -0.0293163*tf.T953 + 14.1777*tf.lnT9)

    rate_eval.Na24__n_Na23 = rate

@numba.njit()
def Mg22__p_Na21(rate_eval, tf):
    # Mg22 --> p + Na21
    rate = 0.0

    # il10r
    rate += np.exp(  31.3024 + -71.3515*tf.T9i + 4.75873*tf.T913
                  + -0.708067*tf.T9 + 0.0523665*tf.T953)
    # il10r
    rate += np.exp(  30.1015 + -66.2607*tf.T9i)
    # il10n
    rate += np.exp(  43.6781 + -63.8733*tf.T9i + -20.7329*tf.T913i
                  + -1.5*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Mg22__p_Na21 = rate

@numba.njit()
def Mg22__He4_Ne18(rate_eval, tf):
    # Mg22 --> He4 + Ne18
    rate = 0.0

    # ths8r
    rate += np.exp(  57.6776 + -94.4496*tf.T9i + -46.4859*tf.T913i + 0.956741*tf.T913
                  + -0.914402*tf.T9 + 0.0722478*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Mg22__He4_Ne18 = rate

@numba.njit()
def Mg23__n_Mg22(rate_eval, tf):
    # Mg23 --> n + Mg22
    rate = 0.0

    # ths8r
    rate += np.exp(  31.5736 + -152.577*tf.T9i + 2.72435*tf.T913
                  + -0.470498*tf.T9 + 0.0357362*tf.T953 + 1.5*tf.lnT9)

    rate_eval.Mg23__n_Mg22 = rate

@numba.njit()
def Mg23__p_Na22(rate_eval, tf):
    # Mg23 --> p + Na22
    rate = 0.0

    # il10r
    rate += np.exp(  7.95641 + -88.7434*tf.T9i)
    # il10r
    rate += np.exp(  -1.07519 + -88.4655*tf.T9i)
    # il10r
    rate += np.exp(  12.9256 + -90.3923*tf.T9i + 4.86658*tf.T913i + 16.4592*tf.T913
                  + -1.95129*tf.T9 + 0.132972*tf.T953)

    rate_eval.Mg23__p_Na22 = rate

@numba.njit()
def Mg23__He4_Ne19(rate_eval, tf):
    # Mg23 --> He4 + Ne19
    rate = 0.0

    # ths8r
    rate += np.exp(  61.3121 + -111.985*tf.T9i + -46.6346*tf.T913i + -1.1007*tf.T913
                  + -0.794097*tf.T9 + 0.0813036*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Mg23__He4_Ne19 = rate

@numba.njit()
def Mg24__n_Mg23(rate_eval, tf):
    # Mg24 --> n + Mg23
    rate = 0.0

    # ths8r
    rate += np.exp(  32.0344 + -191.835*tf.T9i + 2.66964*tf.T913
                  + -0.448904*tf.T9 + 0.0326505*tf.T953 + 1.5*tf.lnT9)

    rate_eval.Mg24__n_Mg23 = rate

@numba.njit()
def Mg24__p_Na23(rate_eval, tf):
    # Mg24 --> p + Na23
    rate = 0.0

    # il10r
    rate += np.exp(  34.0876 + -138.968*tf.T9i + -0.360588*tf.T913
                  + 1.4187*tf.T9 + -0.184061*tf.T953)
    # il10r
    rate += np.exp(  20.0024 + -137.3*tf.T9i)
    # il10n
    rate += np.exp(  43.9357 + -135.688*tf.T9i + -20.6428*tf.T913i + 1.52954*tf.T913
                  + 2.7487*tf.T9 + -1.0*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Mg24__p_Na23 = rate

@numba.njit()
def Mg24__He4_Ne20(rate_eval, tf):
    # Mg24 --> He4 + Ne20
    rate = 0.0

    # il10r
    rate += np.exp(  16.0203 + -120.895*tf.T9i + 16.9229*tf.T913
                  + -2.57325*tf.T9 + 0.208997*tf.T953)
    # il10r
    rate += np.exp(  26.8017 + -117.334*tf.T9i)
    # il10r
    rate += np.exp(  -13.8869 + -110.62*tf.T9i)
    # il10n
    rate += np.exp(  49.3244 + -108.114*tf.T9i + -46.2525*tf.T913i + 5.58901*tf.T913
                  + 7.61843*tf.T9 + -3.683*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Mg24__He4_Ne20 = rate

@numba.njit()
def Mg25__n_Mg24(rate_eval, tf):
    # Mg25 --> n + Mg24
    rate = 0.0

    # ks03 
    rate += np.exp(  86.4748 + -84.9032*tf.T9i + -0.142939*tf.T913i + -57.7499*tf.T913
                  + 7.01981*tf.T9 + -0.582057*tf.T953 + 14.3133*tf.lnT9)

    rate_eval.Mg25__n_Mg24 = rate

@numba.njit()
def Mg25__p_Na24(rate_eval, tf):
    # Mg25 --> p + Na24
    rate = 0.0

    # rath 
    rate += np.exp(  87.7599 + -142.543*tf.T9i + 88.4533*tf.T913i + -156.066*tf.T913
                  + 8.77114*tf.T9 + -0.500891*tf.T953 + 78.7574*tf.lnT9)

    rate_eval.Mg25__p_Na24 = rate

@numba.njit()
def Mg25__He4_Ne21(rate_eval, tf):
    # Mg25 --> He4 + Ne21
    rate = 0.0

    # cf88r
    rate += np.exp(  50.668 + -136.725*tf.T9i + -29.4583*tf.T913
                  + 14.6328*tf.T9 + -3.47392*tf.T953)
    # cf88n
    rate += np.exp(  61.1178 + -114.676*tf.T9i + -46.89*tf.T913i + -0.72642*tf.T913
                  + -0.76406*tf.T9 + 0.0797483*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Mg25__He4_Ne21 = rate

@numba.njit()
def Mg26__n_Mg25(rate_eval, tf):
    # Mg26 --> n + Mg25
    rate = 0.0

    # ks03 
    rate += np.exp(  63.7787 + -128.778*tf.T9i + 9.392*tf.T913i + -36.6784*tf.T913
                  + 3.09567*tf.T9 + -0.223882*tf.T953 + 13.8852*tf.lnT9)

    rate_eval.Mg26__n_Mg25 = rate

@numba.njit()
def Al25__p_Mg24(rate_eval, tf):
    # Al25 --> p + Mg24
    rate = 0.0

    # il10r
    rate += np.exp(  30.093 + -28.8453*tf.T9i + -1.57811*tf.T913
                  + 1.52232*tf.T9 + -0.183001*tf.T953)
    # il10n
    rate += np.exp(  41.7494 + -26.3608*tf.T9i + -22.0227*tf.T913i + 0.361297*tf.T913
                  + 2.61292*tf.T9 + -1.0*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Al25__p_Mg24 = rate

@numba.njit()
def Al25__He4_Na21(rate_eval, tf):
    # Al25 --> He4 + Na21
    rate = 0.0

    # ths8r
    rate += np.exp(  59.7257 + -106.262*tf.T9i + -49.9709*tf.T913i + 1.63835*tf.T913
                  + -1.18562*tf.T9 + 0.101965*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Al25__He4_Na21 = rate

@numba.njit()
def Al26__n_Al25(rate_eval, tf):
    # Al26 --> n + Al25
    rate = 0.0

    # ths8r
    rate += np.exp(  30.9667 + -131.891*tf.T9i + 1.17141*tf.T913
                  + -0.162515*tf.T9 + 0.0126275*tf.T953 + 1.5*tf.lnT9)

    rate_eval.Al26__n_Al25 = rate

@numba.njit()
def Al26__p_Mg25(rate_eval, tf):
    # Al26 --> p + Mg25
    rate = 0.0

    # il10r
    rate += np.exp(  25.2686 + -76.4067*tf.T9i + 8.46334*tf.T913
                  + -0.907024*tf.T9 + 0.0642981*tf.T953)
    # il10r
    rate += np.exp(  27.2591 + -73.903*tf.T9i + -88.9297*tf.T913
                  + 302.948*tf.T9 + -346.461*tf.T953)
    # il10r
    rate += np.exp(  -14.1555 + -73.6126*tf.T9i)

    rate_eval.Al26__p_Mg25 = rate

@numba.njit()
def Al26__He4_Na22(rate_eval, tf):
    # Al26 --> He4 + Na22
    rate = 0.0

    # ths8r
    rate += np.exp(  60.7692 + -109.695*tf.T9i + -50.0924*tf.T913i + -0.390826*tf.T913
                  + -0.99531*tf.T9 + 0.101354*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Al26__He4_Na22 = rate

@numba.njit()
def Al27__n_Al26(rate_eval, tf):
    # Al27 --> n + Al26
    rate = 0.0

    # ks03 
    rate += np.exp(  39.0178 + -151.532*tf.T9i + -0.171158*tf.T913i + -1.77283*tf.T913
                  + 0.206192*tf.T9 + -0.0191705*tf.T953 + 1.63961*tf.lnT9)

    rate_eval.Al27__n_Al26 = rate

@numba.njit()
def Al27__p_Mg26(rate_eval, tf):
    # Al27 --> p + Mg26
    rate = 0.0

    # il10r
    rate += np.exp(  27.118 + -99.3406*tf.T9i + 6.78105*tf.T913
                  + -1.25771*tf.T9 + 0.140754*tf.T953)
    # il10r
    rate += np.exp(  -5.3594 + -96.8701*tf.T9i + 35.6312*tf.T913
                  + -5.27265*tf.T9 + 0.392932*tf.T953)
    # il10r
    rate += np.exp(  -62.6356 + -96.4509*tf.T9i + 251.281*tf.T913
                  + -730.009*tf.T9 + -224.016*tf.T953)

    rate_eval.Al27__p_Mg26 = rate

@numba.njit()
def Al27__He4_Na23(rate_eval, tf):
    # Al27 --> He4 + Na23
    rate = 0.0

    # ths8r
    rate += np.exp(  69.2185 + -117.109*tf.T9i + -50.2042*tf.T913i + -1.64239*tf.T913
                  + -1.59995*tf.T9 + 0.184933*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Al27__He4_Na23 = rate

@numba.njit()
def Al28__n_Al27(rate_eval, tf):
    # Al28 --> n + Al27
    rate = 0.0

    # ks03 
    rate += np.exp(  43.2942 + -89.7267*tf.T9i + 7.49524*tf.T913i + -14.6211*tf.T913
                  + 0.721304*tf.T9 + -0.0288351*tf.T953 + 8.08346*tf.lnT9)

    rate_eval.Al28__n_Al27 = rate

@numba.njit()
def Al28__He4_Na24(rate_eval, tf):
    # Al28 --> He4 + Na24
    rate = 0.0

    # rath 
    rate += np.exp(  26.7427 + -128.003*tf.T9i + -22.0373*tf.T913i + 12.2654*tf.T913
                  + -3.20073*tf.T9 + 0.258909*tf.T953 + 6.02728*tf.lnT9)

    rate_eval.Al28__He4_Na24 = rate

@numba.njit()
def Si26__p_Al25(rate_eval, tf):
    # Si26 --> p + Al25
    rate = 0.0

    # li20r
    rate += np.exp(  34.1845 + -68.7787*tf.T9i
                  + 0.01455*tf.lnT9)
    # li20r
    rate += np.exp(  30.8265 + -76.4488*tf.T9i
                  + 2.98721*tf.lnT9)
    # li20r
    rate += np.exp(  19.2308 + -65.7211*tf.T9i
                  + 1.25571*tf.lnT9)

    rate_eval.Si26__p_Al25 = rate

@numba.njit()
def Si26__He4_Mg22(rate_eval, tf):
    # Si26 --> He4 + Mg22
    rate = 0.0

    # ths8r
    rate += np.exp(  60.748 + -106.443*tf.T9i + -53.0847*tf.T913i + 0.625592*tf.T913
                  + -0.893632*tf.T9 + 0.0691552*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Si26__He4_Mg22 = rate

@numba.njit()
def Si28__p_Al27(rate_eval, tf):
    # Si28 --> p + Al27
    rate = 0.0

    # il10n
    rate += np.exp(  46.5494 + -134.445*tf.T9i + -23.2205*tf.T913i
                  + -2.0*tf.T953 + 0.833333*tf.lnT9)
    # il10r
    rate += np.exp(  111.466 + -134.832*tf.T9i + -26.8327*tf.T913i + -116.137*tf.T913
                  + 0.00950567*tf.T9 + 0.00999755*tf.T953)
    # il10r
    rate += np.exp(  11.7765 + -136.349*tf.T9i + 23.8634*tf.T913
                  + -3.70135*tf.T9 + 0.28964*tf.T953)

    rate_eval.Si28__p_Al27 = rate

@numba.njit()
def Si28__He4_Mg24(rate_eval, tf):
    # Si28 --> He4 + Mg24
    rate = 0.0

    # st08r
    rate += np.exp(  32.9006 + -131.488*tf.T9i)
    # st08r
    rate += np.exp(  -25.6886 + -128.693*tf.T9i + 21.3721*tf.T913i + 37.7649*tf.T913
                  + -4.10635*tf.T9 + 0.249618*tf.T953)

    rate_eval.Si28__He4_Mg24 = rate

@numba.njit()
def Si29__n_Si28(rate_eval, tf):
    # Si29 --> n + Si28
    rate = 0.0

    # ka02 
    rate += np.exp(  31.7355 + -98.3365*tf.T9i
                  + 1.5*tf.lnT9)
    # ka02r
    rate += np.exp(  29.8758 + -98.7165*tf.T9i + 7.68863*tf.T913
                  + -1.7991*tf.T9)

    rate_eval.Si29__n_Si28 = rate

@numba.njit()
def Si29__p_Al28(rate_eval, tf):
    # Si29 --> p + Al28
    rate = 0.0

    # rath 
    rate += np.exp(  107.692 + -146.264*tf.T9i + 114.903*tf.T913i + -203.401*tf.T913
                  + 11.5793*tf.T9 + -0.672783*tf.T953 + 100.918*tf.lnT9)

    rate_eval.Si29__p_Al28 = rate

@numba.njit()
def Si29__He4_Mg25(rate_eval, tf):
    # Si29 --> He4 + Mg25
    rate = 0.0

    # cf88n
    rate += np.exp(  66.3395 + -129.123*tf.T9i + -53.41*tf.T913i + -1.83266*tf.T913
                  + -0.573073*tf.T9 + 0.0462678*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Si29__He4_Mg25 = rate

@numba.njit()
def Si30__n_Si29(rate_eval, tf):
    # Si30 --> n + Si29
    rate = 0.0

    # ka02n
    rate += np.exp(  36.1504 + -123.112*tf.T9i + 0.650904*tf.T913
                  + 1.5*tf.lnT9)
    # ka02r
    rate += np.exp(  33.9492 + -123.292*tf.T9i + 5.50678*tf.T913
                  + -2.85656*tf.T9)

    rate_eval.Si30__n_Si29 = rate

@numba.njit()
def Si30__He4_Mg26(rate_eval, tf):
    # Si30 --> He4 + Mg26
    rate = 0.0

    # cf88r
    rate += np.exp(  26.2068 + -142.235*tf.T9i + -1.87411*tf.T913
                  + 3.41299*tf.T9 + -0.43226*tf.T953)
    # cf88n
    rate += np.exp(  70.7561 + -123.518*tf.T9i + -53.7518*tf.T913i + -4.8647*tf.T913
                  + -1.51467*tf.T9 + 0.833333*tf.lnT9)

    rate_eval.Si30__He4_Mg26 = rate

@numba.njit()
def P29__p_Si28(rate_eval, tf):
    # P29 --> p + Si28
    rate = 0.0

    # il10n
    rate += np.exp(  39.1379 + -31.8984*tf.T9i + -23.8173*tf.T913i + 7.08203*tf.T913
                  + -1.44753*tf.T9 + 0.0804296*tf.T953 + 0.833333*tf.lnT9)
    # il10r
    rate += np.exp(  28.6997 + -36.0408*tf.T9i)

    rate_eval.P29__p_Si28 = rate

@numba.njit()
def P29__He4_Al25(rate_eval, tf):
    # P29 --> He4 + Al25
    rate = 0.0

    # ths8r
    rate += np.exp(  63.8779 + -121.399*tf.T9i + -56.3424*tf.T913i + 0.542998*tf.T913
                  + -0.721716*tf.T9 + 0.0469712*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.P29__He4_Al25 = rate

@numba.njit()
def P30__n_P29(rate_eval, tf):
    # P30 --> n + P29
    rate = 0.0

    # ths8r
    rate += np.exp(  32.0379 + -131.355*tf.T9i + 0.15555*tf.T913
                  + 0.155359*tf.T9 + -0.0208019*tf.T953 + 1.5*tf.lnT9)

    rate_eval.P30__n_P29 = rate

@numba.njit()
def P30__p_Si29(rate_eval, tf):
    # P30 --> p + Si29
    rate = 0.0

    # il10r
    rate += np.exp(  22.0015 + -68.2607*tf.T9i + 14.0921*tf.T913
                  + -3.92096*tf.T9 + 0.447706*tf.T953)
    # il10r
    rate += np.exp(  9.77935 + -66.1716*tf.T9i)
    # il10n
    rate += np.exp(  39.7677 + -64.9214*tf.T9i + -23.9101*tf.T913i + 10.7796*tf.T913
                  + -3.04181*tf.T9 + 0.274565*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.P30__p_Si29 = rate

@numba.njit()
def P30__He4_Al26(rate_eval, tf):
    # P30 --> He4 + Al26
    rate = 0.0

    # ths8r
    rate += np.exp(  69.1545 + -120.863*tf.T9i + -56.4422*tf.T913i + -2.44848*tf.T913
                  + -1.17578*tf.T9 + 0.150757*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.P30__He4_Al26 = rate

@numba.njit()
def P31__n_P30(rate_eval, tf):
    # P31 --> n + P30
    rate = 0.0

    # ths8r
    rate += np.exp(  36.8808 + -142.87*tf.T9i + 0.909911*tf.T913
                  + -0.162367*tf.T9 + 0.00668293*tf.T953 + 1.5*tf.lnT9)

    rate_eval.P31__n_P30 = rate

@numba.njit()
def P31__p_Si30(rate_eval, tf):
    # P31 --> p + Si30
    rate = 0.0

    # de20r
    rate += np.exp(  18.7213 + -85.9282*tf.T9i
                  + 6.49034*tf.lnT9)
    # de20r
    rate += np.exp(  4.04359 + -85.6217*tf.T9i
                  + 2.80331*tf.lnT9)
    # de20r
    rate += np.exp(  -1115.38 + -180.553*tf.T9i
                  + -895.258*tf.lnT9)
    # de20r
    rate += np.exp(  32.9288 + -90.2661*tf.T9i
                  + -0.070816*tf.lnT9)
    # de20r
    rate += np.exp(  -11.8961 + -85.2694*tf.T9i
                  + -0.128387*tf.lnT9)
    # de20r
    rate += np.exp(  35.3329 + -91.3175*tf.T9i
                  + 0.3809*tf.lnT9)
    # de20r
    rate += np.exp(  31.761 + -89.8588*tf.T9i
                  + 2.7883*tf.lnT9)
    # de20r
    rate += np.exp(  -311.303 + -85.8097*tf.T9i
                  + -77.047*tf.lnT9)

    rate_eval.P31__p_Si30 = rate

@numba.njit()
def P31__He4_Al27(rate_eval, tf):
    # P31 --> He4 + Al27
    rate = 0.0

    # ths8r
    rate += np.exp(  73.2168 + -112.206*tf.T9i + -56.5351*tf.T913i + -0.896208*tf.T913
                  + -1.72024*tf.T9 + 0.185409*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.P31__He4_Al27 = rate

@numba.njit()
def P32__n_P31(rate_eval, tf):
    # P32 --> n + P31
    rate = 0.0

    # ka02r
    rate += np.exp(  32.7384 + -92.3401*tf.T9i + 2.13185*tf.T913)
    # ka02n
    rate += np.exp(  35.7859 + -92.0933*tf.T9i + 0.327031*tf.T913
                  + 2.5*tf.lnT9)

    rate_eval.P32__n_P31 = rate

@numba.njit()
def P32__He4_Al28(rate_eval, tf):
    # P32 --> He4 + Al28
    rate = 0.0

    # rath 
    rate += np.exp(  25.2167 + -117.232*tf.T9i + -14.2992*tf.T913i + 6.40835*tf.T913
                  + -3.59742*tf.T9 + 0.308804*tf.T953 + 11.594*tf.lnT9)

    rate_eval.P32__He4_Al28 = rate

@numba.njit()
def S31__p_P30(rate_eval, tf):
    # S31 --> p + P30
    rate = 0.0

    # mb07 
    rate += np.exp(  -7732.57 + -60.6616*tf.T9i + -1999.51*tf.T913i + 11886.5*tf.T913
                  + -2668.72*tf.T9 + 354.294*tf.T953 + -2898.95*tf.lnT9)
    # mb07 
    rate += np.exp(  48.4487 + -80.154*tf.T9i + 156.029*tf.T913i + -174.377*tf.T913
                  + 7.4644*tf.T9 + -0.342232*tf.T953 + 100.758*tf.lnT9)

    rate_eval.S31__p_P30 = rate

@numba.njit()
def S32__n_S31(rate_eval, tf):
    # S32 --> n + S31
    rate = 0.0

    # ths8r
    rate += np.exp(  31.9171 + -174.56*tf.T9i + 1.71463*tf.T913
                  + -0.221804*tf.T9 + 0.00880104*tf.T953 + 1.5*tf.lnT9)

    rate_eval.S32__n_S31 = rate

@numba.njit()
def S32__p_P31(rate_eval, tf):
    # S32 --> p + P31
    rate = 0.0

    # il10r
    rate += np.exp(  25.1729 + -106.637*tf.T9i + 8.09341*tf.T913
                  + -0.615971*tf.T9 + 0.031159*tf.T953)
    # il10r
    rate += np.exp(  21.6829 + -105.119*tf.T9i)
    # il10n
    rate += np.exp(  43.6109 + -102.86*tf.T9i + -25.3278*tf.T913i + 6.4931*tf.T913
                  + -9.27513*tf.T9 + -0.610439*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.S32__p_P31 = rate

@numba.njit()
def S32__He4_Si28(rate_eval, tf):
    # S32 --> He4 + Si28
    rate = 0.0

    # ths8r
    rate += np.exp(  72.813 + -80.626*tf.T9i + -59.4896*tf.T913i + 4.47205*tf.T913
                  + -4.78989*tf.T9 + 0.557201*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.S32__He4_Si28 = rate

@numba.njit()
def S34__He4_Si30(rate_eval, tf):
    # S34 --> He4 + Si30
    rate = 0.0

    # ths8r
    rate += np.exp(  72.4051 + -91.9515*tf.T9i + -59.6559*tf.T913i + 3.70141*tf.T913
                  + -3.12537*tf.T9 + 0.307626*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.S34__He4_Si30 = rate

@numba.njit()
def Cl32__p_S31(rate_eval, tf):
    # Cl32 --> p + S31
    rate = 0.0

    # mm11r
    rate += np.exp(  11.4998 + -20.041*tf.T9i)
    # mm11n
    rate += np.exp(  29.5568 + -18.2655*tf.T9i + -25.3278*tf.T913i + 42.7226*tf.T913
                  + -102.23*tf.T9 + -2.0*tf.T953 + 0.833333*tf.lnT9)
    # mm11r
    rate += np.exp(  26.5247 + -24.6248*tf.T9i + 3.202*tf.T913
                  + 0.101684*tf.T9 + -0.0334894*tf.T953)

    rate_eval.Cl32__p_S31 = rate

@numba.njit()
def Cl32__p_P31__weak__wc12(rate_eval, tf):
    # Cl32 --> p + P31
    rate = 0.0

    # wc12w
    rate += np.exp(  -7.26758)

    rate_eval.Cl32__p_P31__weak__wc12 = rate

@numba.njit()
def Cl32__He4_Si28__weak__wc12(rate_eval, tf):
    # Cl32 --> He4 + Si28
    rate = 0.0

    # wc12w
    rate += np.exp(  -6.75675)

    rate_eval.Cl32__He4_Si28__weak__wc12 = rate

@numba.njit()
def Cl34__He4_P30(rate_eval, tf):
    # Cl34 --> He4 + P30
    rate = 0.0

    # ths8r
    rate += np.exp(  71.335 + -77.3338*tf.T9i + -62.4643*tf.T913i + -3.19028*tf.T913
                  + -0.832633*tf.T9 + 0.0987525*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Cl34__He4_P30 = rate

@numba.njit()
def Cl35__n_Cl34(rate_eval, tf):
    # Cl35 --> n + Cl34
    rate = 0.0

    # ths8r
    rate += np.exp(  34.9299 + -146.74*tf.T9i + 0.990222*tf.T913
                  + -0.146686*tf.T9 + 0.00560251*tf.T953 + 1.5*tf.lnT9)

    rate_eval.Cl35__n_Cl34 = rate

@numba.njit()
def Cl35__p_S34(rate_eval, tf):
    # Cl35 --> p + S34
    rate = 0.0

    # se19r
    rate += np.exp(  22.7104 + -77.1013*tf.T9i)
    # se19r
    rate += np.exp(  -37.7488 + -101.622*tf.T9i + 454.53*tf.T913i + -375.557*tf.T913
                  + 12.6533*tf.T9 + -0.408677*tf.T953 + 246.618*tf.lnT9)
    # se19r
    rate += np.exp(  14.5754 + -75.9594*tf.T9i)
    # se19r
    rate += np.exp(  30.5303 + -79.6891*tf.T9i)
    # se19r
    rate += np.exp(  7.39781 + -75.3374*tf.T9i)
    # se19r
    rate += np.exp(  29.2773 + -79.0392*tf.T9i)
    # se19r
    rate += np.exp(  4.33927 + -75.1053*tf.T9i)
    # se19r
    rate += np.exp(  28.0819 + -78.5751*tf.T9i)

    rate_eval.Cl35__p_S34 = rate

@numba.njit()
def Cl35__He4_P31(rate_eval, tf):
    # Cl35 --> He4 + P31
    rate = 0.0

    # ths8r
    rate += np.exp(  74.6679 + -81.2033*tf.T9i + -62.5433*tf.T913i + -2.95026*tf.T913
                  + -0.89652*tf.T9 + 0.0774126*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Cl35__He4_P31 = rate

@numba.njit()
def Ar35__p_Cl34(rate_eval, tf):
    # Ar35 --> p + Cl34
    rate = 0.0

    # ths8r
    rate += np.exp(  57.222 + -68.4272*tf.T9i + -27.8896*tf.T913i + -0.491144*tf.T913
                  + -2.47595*tf.T9 + 0.34035*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Ar35__p_Cl34 = rate

@numba.njit()
def Ar35__He4_S31(rate_eval, tf):
    # Ar35 --> He4 + S31
    rate = 0.0

    # ths8r
    rate += np.exp(  80.3587 + -74.59*tf.T9i + -65.2934*tf.T913i + -13.3991*tf.T913
                  + 1.97512*tf.T9 + -0.208606*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Ar35__He4_S31 = rate

@numba.njit()
def He6__n_n_He4(rate_eval, tf):
    # He6 --> n + n + He4
    rate = 0.0

    # cf88r
    rate += np.exp(  22.178 + -20.8994*tf.T9i + 0.694279*tf.T913i + -3.33326*tf.T913
                  + 0.507932*tf.T9 + -0.0427342*tf.T953 + 2.0*tf.lnT9)

    rate_eval.He6__n_n_He4 = rate

@numba.njit()
def Li6__n_p_He4(rate_eval, tf):
    # Li6 --> n + p + He4
    rate = 0.0

    # cf88r
    rate += np.exp(  33.4196 + -62.2896*tf.T9i + 1.44987*tf.T913i + -1.42759*tf.T913
                  + 0.0454035*tf.T9 + 0.00471161*tf.T953 + 2.0*tf.lnT9)

    rate_eval.Li6__n_p_He4 = rate

@numba.njit()
def Be9__n_He4_He4(rate_eval, tf):
    # Be9 --> n + He4 + He4
    rate = 0.0

    # ac12n
    rate += np.exp(  37.273 + -18.2597*tf.T9i + -13.3317*tf.T913i + 13.2237*tf.T913
                  + -9.06339*tf.T9 + 2.33333*tf.lnT9)
    # ac12r
    rate += np.exp(  38.6902 + -19.2792*tf.T9i + -1.56673*tf.T913i + -5.43497*tf.T913
                  + 0.673807*tf.T9 + -0.041014*tf.T953 + 1.5*tf.lnT9)

    rate_eval.Be9__n_He4_He4 = rate

@numba.njit()
def C12__He4_He4_He4(rate_eval, tf):
    # C12 --> 3 He4
    rate = 0.0

    # fy05r
    rate += np.exp(  22.394 + -88.5493*tf.T9i + -13.49*tf.T913i + 21.4259*tf.T913
                  + -1.34769*tf.T9 + 0.0879816*tf.T953 + -10.1653*tf.lnT9)
    # fy05r
    rate += np.exp(  34.9561 + -85.4472*tf.T9i + -23.57*tf.T913i + 20.4886*tf.T913
                  + -12.9882*tf.T9 + -20.0*tf.T953 + 0.83333*tf.lnT9)
    # fy05n
    rate += np.exp(  45.7734 + -84.4227*tf.T9i + -37.06*tf.T913i + 29.3493*tf.T913
                  + -115.507*tf.T9 + -10.0*tf.T953 + 1.66667*tf.lnT9)

    rate_eval.C12__He4_He4_He4 = rate

@numba.njit()
def n_p__d(rate_eval, tf):
    # n + p --> d
    rate = 0.0

    # an06n
    rate += np.exp(  12.3687 + -2.70618*tf.T913
                  + 0.11718*tf.T9 + -0.00312788*tf.T953 + 0.469127*tf.lnT9)
    # an06n
    rate += np.exp(  10.7548 + -2.30472*tf.T913
                  + -0.887862*tf.T9 + 0.137663*tf.T953)
    # an06n
    rate += np.exp(  8.84688 + -0.0102082*tf.T913
                  + -0.0893959*tf.T9 + 0.00696704*tf.T953 + 1.0*tf.lnT9)

    rate_eval.n_p__d = rate

@numba.njit()
def p_p__d__weak__bet_pos_(rate_eval, tf):
    # p + p --> d
    rate = 0.0

    # bet+w
    rate += np.exp(  -34.7863 + -3.51193*tf.T913i + 3.10086*tf.T913
                  + -0.198314*tf.T9 + 0.0126251*tf.T953 + -1.02517*tf.lnT9)

    rate_eval.p_p__d__weak__bet_pos_ = rate

@numba.njit()
def p_p__d__weak__electron_capture(rate_eval, tf):
    # p + p --> d
    rate = 0.0

    #   ecw
    rate += np.exp(  -43.6499 + -0.00246064*tf.T9i + -2.7507*tf.T913i + -0.424877*tf.T913
                  + 0.015987*tf.T9 + -0.000690875*tf.T953 + -0.207625*tf.lnT9)

    rate_eval.p_p__d__weak__electron_capture = rate

@numba.njit()
def n_d__t(rate_eval, tf):
    # d + n --> t
    rate = 0.0

    # nk06n
    rate += np.exp(  6.60935
                  + 1.0*tf.lnT9)
    # nk06n
    rate += np.exp(  5.36598
                  + 0.075*tf.lnT9)

    rate_eval.n_d__t = rate

@numba.njit()
def p_d__He3(rate_eval, tf):
    # d + p --> He3
    rate = 0.0

    # de04 
    rate += np.exp(  8.93525 + -3.7208*tf.T913i + 0.198654*tf.T913
                  + 0.333333*tf.lnT9)
    # de04n
    rate += np.exp(  7.52898 + -3.7208*tf.T913i + 0.871782*tf.T913
                  + -0.666667*tf.lnT9)

    rate_eval.p_d__He3 = rate

@numba.njit()
def d_d__He4(rate_eval, tf):
    # d + d --> He4
    rate = 0.0

    # nacrn
    rate += np.exp(  3.78177 + -4.26166*tf.T913i + -0.119233*tf.T913
                  + 0.778829*tf.T9 + -0.0925203*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.d_d__He4 = rate

@numba.njit()
def He4_d__Li6(rate_eval, tf):
    # d + He4 --> Li6
    rate = 0.0

    # tu19r
    rate += np.exp(  4.12313 + -7.889*tf.T9i
                  + -1.5*tf.lnT9)
    # tu19n
    rate += np.exp(  -0.676485 + 6.3911e-05*tf.T9i + -7.55198*tf.T913i + 5.77546*tf.T913
                  + -0.487854*tf.T9 + 0.032833*tf.T953 + -1.12305*tf.lnT9)

    rate_eval.He4_d__Li6 = rate

@numba.njit()
def p_t__He4(rate_eval, tf):
    # t + p --> He4
    rate = 0.0

    # cf88n
    rate += np.exp(  9.76526 + -3.869*tf.T913i + 1.45482*tf.T913
                  + 0.577246*tf.T9 + -0.112199*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_t__He4 = rate

@numba.njit()
def He4_t__Li7(rate_eval, tf):
    # t + He4 --> Li7
    rate = 0.0

    # de04 
    rate += np.exp(  13.6162 + -8.0805*tf.T913i + -0.217514*tf.T913
                  + -0.114859*tf.T9 + 0.0470043*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_t__Li7 = rate

@numba.njit()
def n_He3__He4(rate_eval, tf):
    # He3 + n --> He4
    rate = 0.0

    # ka02n
    rate += np.exp(  9.04572 + -1.50147*tf.T913
                  + 1.0*tf.lnT9)
    # ka02n
    rate += np.exp(  5.51711)

    rate_eval.n_He3__He4 = rate

@numba.njit()
def p_He3__He4__weak__bet_pos_(rate_eval, tf):
    # He3 + p --> He4
    rate = 0.0

    # bet+w
    rate += np.exp(  -27.7611 + -4.30107e-12*tf.T9i + -6.141*tf.T913i + -1.93473e-09*tf.T913
                  + 2.04145e-10*tf.T9 + -1.80372e-11*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_He3__He4__weak__bet_pos_ = rate

@numba.njit()
def He4_He3__Be7(rate_eval, tf):
    # He3 + He4 --> Be7
    rate = 0.0

    # cd08n
    rate += np.exp(  17.7075 + -12.8271*tf.T913i + -3.8126*tf.T913
                  + 0.0942285*tf.T9 + -0.00301018*tf.T953 + 1.33333*tf.lnT9)
    # cd08n
    rate += np.exp(  15.6099 + -12.8271*tf.T913i + -0.0308225*tf.T913
                  + -0.654685*tf.T9 + 0.0896331*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_He3__Be7 = rate

@numba.njit()
def n_Li6__Li7(rate_eval, tf):
    # Li6 + n --> Li7
    rate = 0.0

    # jz10n
    rate += np.exp(  9.04782)

    rate_eval.n_Li6__Li7 = rate

@numba.njit()
def p_Li6__Be7(rate_eval, tf):
    # Li6 + p --> Be7
    rate = 0.0

    # nacrn
    rate += np.exp(  14.2792 + -8.4372*tf.T913i + -0.515473*tf.T913
                  + 0.0285578*tf.T9 + 0.00879731*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Li6__Be7 = rate

@numba.njit()
def p_Be7__B8(rate_eval, tf):
    # Be7 + p --> B8
    rate = 0.0

    # nacrr
    rate += np.exp(  7.73399 + -7.345*tf.T9i
                  + -1.5*tf.lnT9)
    # nacrn
    rate += np.exp(  12.5315 + -10.264*tf.T913i + -0.203472*tf.T913
                  + 0.121083*tf.T9 + -0.00700063*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Be7__B8 = rate

@numba.njit()
def He4_Be7__C11(rate_eval, tf):
    # Be7 + He4 --> C11
    rate = 0.0

    # nacrr
    rate += np.exp(  11.8776 + -10.177*tf.T9i + 0.685881*tf.T913
                  + -0.697071*tf.T9 + 0.13274*tf.T953 + -1.5*tf.lnT9)
    # nacrr
    rate += np.exp(  9.43348 + -6.498*tf.T9i
                  + -1.5*tf.lnT9)
    # nacrn
    rate += np.exp(  24.2273 + -23.214*tf.T913i + -3.74943*tf.T913
                  + 1.23242*tf.T9 + -0.195665*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Be7__C11 = rate

@numba.njit()
def n_C11__C12(rate_eval, tf):
    # C11 + n --> C12
    rate = 0.0

    # bb92n
    rate += np.exp(  10.3672 + -0.0095137*tf.T9i + 1.18875e-10*tf.T913i + -3.10977e-10*tf.T913
                  + 2.0423e-11*tf.T9 + -6.59293e-13*tf.T953)
    # bb92r
    rate += np.exp(  13.8643 + -5.57951*tf.T9i + 1.38808e-09*tf.T913i + -3.66292e-09*tf.T913
                  + 2.41904e-10*tf.T9 + -7.82223e-12*tf.T953 + -1.5*tf.lnT9)
    # bb92r
    rate += np.exp(  8.10155 + -0.926514*tf.T9i + -5.21741e-12*tf.T913i + 2.32794e-11*tf.T913
                  + -2.03119e-12*tf.T9 + 7.09452e-14*tf.T953 + -1.5*tf.lnT9)

    rate_eval.n_C11__C12 = rate

@numba.njit()
def n_C12__C13(rate_eval, tf):
    # C12 + n --> C13
    rate = 0.0

    # ks03 
    rate += np.exp(  7.98821 + -0.00836732*tf.T9i + 1.49573*tf.T913i + -0.841102*tf.T913
                  + 0.0340543*tf.T9 + -0.0026392*tf.T953 + 1.6662*tf.lnT9)

    rate_eval.n_C12__C13 = rate

@numba.njit()
def p_C12__N13(rate_eval, tf):
    # C12 + p --> N13
    rate = 0.0

    # ls09n
    rate += np.exp(  17.1482 + -13.692*tf.T913i + -0.230881*tf.T913
                  + 4.44362*tf.T9 + -3.15898*tf.T953 + -0.666667*tf.lnT9)
    # ls09r
    rate += np.exp(  17.5428 + -3.77849*tf.T9i + -5.10735*tf.T913i + -2.24111*tf.T913
                  + 0.148883*tf.T9 + -1.5*tf.lnT9)

    rate_eval.p_C12__N13 = rate

@numba.njit()
def He4_C12__O16(rate_eval, tf):
    # C12 + He4 --> O16
    rate = 0.0

    # nac2 
    rate += np.exp(  254.634 + -1.84097*tf.T9i + 103.411*tf.T913i + -420.567*tf.T913
                  + 64.0874*tf.T9 + -12.4624*tf.T953 + 137.303*tf.lnT9)
    # nac2 
    rate += np.exp(  69.6526 + -1.39254*tf.T9i + 58.9128*tf.T913i + -148.273*tf.T913
                  + 9.08324*tf.T9 + -0.541041*tf.T953 + 70.3554*tf.lnT9)

    rate_eval.He4_C12__O16 = rate

@numba.njit()
def n_C13__C14(rate_eval, tf):
    # C13 + n --> C14
    rate = 0.0

    # ks03 
    rate += np.exp(  35.3048 + -0.133687*tf.T9i + 18.3578*tf.T913i + -46.5786*tf.T913
                  + 2.58472*tf.T9 + -0.118622*tf.T953 + 19.9142*tf.lnT9)

    rate_eval.n_C13__C14 = rate

@numba.njit()
def p_C13__N14(rate_eval, tf):
    # C13 + p --> N14
    rate = 0.0

    # nacrn
    rate += np.exp(  18.5155 + -13.72*tf.T913i + -0.450018*tf.T913
                  + 3.70823*tf.T9 + -1.70545*tf.T953 + -0.666667*tf.lnT9)
    # nacrr
    rate += np.exp(  13.9637 + -5.78147*tf.T9i + -0.196703*tf.T913
                  + 0.142126*tf.T9 + -0.0238912*tf.T953 + -1.5*tf.lnT9)
    # nacrr
    rate += np.exp(  15.1825 + -13.5543*tf.T9i
                  + -1.5*tf.lnT9)

    rate_eval.p_C13__N14 = rate

@numba.njit()
def p_C14__N15(rate_eval, tf):
    # C14 + p --> N15
    rate = 0.0

    # il10r
    rate += np.exp(  17.1024 + -1.52341*tf.T9i + -10.658*tf.T913i + 1.73644*tf.T913
                  + -0.350498*tf.T9 + 0.0279902*tf.T953 + -1.5*tf.lnT9)
    # il10n
    rate += np.exp(  20.119 + -13.9619*tf.T913i + -4.34315*tf.T913
                  + 6.64922*tf.T9 + -3.22592*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_C14__N15 = rate

@numba.njit()
def He4_C14__O18(rate_eval, tf):
    # C14 + He4 --> O18
    rate = 0.0

    # il10r
    rate += np.exp(  11.8309 + -10.3983*tf.T9i + -3.83188*tf.T913
                  + 1.64358*tf.T9 + -0.177785*tf.T953 + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  -23.805 + -2.06876*tf.T9i
                  + -1.5*tf.lnT9)
    # il10n
    rate += np.exp(  18.4877 + -31.7222*tf.T913i + 11.3923*tf.T913
                  + -9.92249*tf.T9 + -2.0*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_C14__O18 = rate

@numba.njit()
def n_N13__N14(rate_eval, tf):
    # N13 + n --> N14
    rate = 0.0

    # wiesn
    rate += np.exp(  13.9377 + -0.0054652*tf.T9i + 1.72241e-10*tf.T913i + -5.62522e-10*tf.T913
                  + 5.59212e-11*tf.T9 + -4.6549e-12*tf.T953)
    # wiesr
    rate += np.exp(  -3.63074 + -2.99547*tf.T9i + 9.44873e-10*tf.T913i + -2.33713e-09*tf.T913
                  + 1.97507e-10*tf.T9 + -1.49747e-11*tf.T953 + -1.5*tf.lnT9)

    rate_eval.n_N13__N14 = rate

@numba.njit()
def p_N13__O14(rate_eval, tf):
    # N13 + p --> O14
    rate = 0.0

    # lg06r
    rate += np.exp(  10.9971 + -6.12602*tf.T9i + 1.57122*tf.T913i
                  + -1.5*tf.lnT9)
    # lg06n
    rate += np.exp(  18.1356 + -15.1676*tf.T913i + 0.0955166*tf.T913
                  + 3.0659*tf.T9 + -0.507339*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_N13__O14 = rate

@numba.njit()
def n_N14__N15(rate_eval, tf):
    # N14 + n --> N15
    rate = 0.0

    # ks03 
    rate += np.exp(  10.1651 + -0.0114078*tf.T9i + 1.396*tf.T913i + -3.47552*tf.T913
                  + 0.351773*tf.T9 + -0.0229344*tf.T953 + 1.02161*tf.lnT9)

    rate_eval.n_N14__N15 = rate

@numba.njit()
def p_N14__O15(rate_eval, tf):
    # N14 + p --> O15
    rate = 0.0

    # im05r
    rate += np.exp(  6.73578 + -4.891*tf.T9i
                  + 0.0682*tf.lnT9)
    # im05r
    rate += np.exp(  7.65444 + -2.998*tf.T9i
                  + -1.5*tf.lnT9)
    # im05n
    rate += np.exp(  20.1169 + -15.193*tf.T913i + -4.63975*tf.T913
                  + 9.73458*tf.T9 + -9.55051*tf.T953 + 0.333333*tf.lnT9)
    # im05n
    rate += np.exp(  17.01 + -15.193*tf.T913i + -0.161954*tf.T913
                  + -7.52123*tf.T9 + -0.987565*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_N14__O15 = rate

@numba.njit()
def He4_N14__F18(rate_eval, tf):
    # N14 + He4 --> F18
    rate = 0.0

    # il10r
    rate += np.exp(  13.8995 + -10.9656*tf.T9i + -5.6227*tf.T913i
                  + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  0.196838 + -5.16034*tf.T9i
                  + -1.5*tf.lnT9)
    # il10n
    rate += np.exp(  21.5339 + -36.2504*tf.T913i
                  + -5.0*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_N14__F18 = rate

@numba.njit()
def p_N15__O16(rate_eval, tf):
    # N15 + p --> O16
    rate = 0.0

    # li10r
    rate += np.exp(  14.5444 + -10.2295*tf.T9i
                  + 0.0459037*tf.T9 + -1.5*tf.lnT9)
    # li10r
    rate += np.exp(  6.59056 + -2.92315*tf.T9i
                  + -1.5*tf.lnT9)
    # li10n
    rate += np.exp(  20.0176 + -15.24*tf.T913i + 0.334926*tf.T913
                  + 4.59088*tf.T9 + -4.78468*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_N15__O16 = rate

@numba.njit()
def He4_N15__F19(rate_eval, tf):
    # N15 + He4 --> F19
    rate = 0.0

    # il10r
    rate += np.exp(  -28.7989 + -4.19986*tf.T9i + 35.4292*tf.T913
                  + -5.5767*tf.T9 + 0.441293*tf.T953 + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  3.5342 + -6.98462*tf.T9i
                  + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  -9.41892 + -4.17795*tf.T9i
                  + -1.5*tf.lnT9)
    # il10n
    rate += np.exp(  25.3916 + -36.2324*tf.T913i
                  + -2.0*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_N15__F19 = rate

@numba.njit()
def n_O14__O15(rate_eval, tf):
    # O14 + n --> O15
    rate = 0.0

    # rpsmr
    rate += np.exp(  9.87196 + 0.0160481*tf.T9i + -1.38986*tf.T913i + 1.74662*tf.T913
                  + -0.0276897*tf.T9 + 0.00321014*tf.T953 + -1.06122*tf.lnT9)

    rate_eval.n_O14__O15 = rate

@numba.njit()
def p_O14__F15(rate_eval, tf):
    # O14 + p --> F15
    rate = 0.0

    # rpsmr
    rate += np.exp(  3.69607 + -19.5633*tf.T9i + 71.6134*tf.T913i + -76.9654*tf.T913
                  + 3.54377*tf.T9 + -0.159804*tf.T953 + 46.3922*tf.lnT9)

    rate_eval.p_O14__F15 = rate

@numba.njit()
def He4_O14__Ne18(rate_eval, tf):
    # O14 + He4 --> Ne18
    rate = 0.0

    # wh87n
    rate += np.exp(  26.4429 + -39.38*tf.T913i + -0.0772187*tf.T913
                  + -0.635361*tf.T9 + 0.106236*tf.T953 + -0.666667*tf.lnT9)
    # wh87r
    rate += np.exp(  -4.69948 + -12.159*tf.T9i
                  + 5.0*tf.lnT9)
    # wh87r
    rate += np.exp(  3.52636 + -22.61*tf.T9i
                  + -1.5*tf.lnT9)
    # wh87r
    rate += np.exp(  -2.15417 + -11.73*tf.T9i
                  + -1.5*tf.lnT9)

    rate_eval.He4_O14__Ne18 = rate

@numba.njit()
def n_O15__O16(rate_eval, tf):
    # O15 + n --> O16
    rate = 0.0

    # rpsmr
    rate += np.exp(  8.08476 + 0.0135346*tf.T9i + -1.11761*tf.T913i + 1.0167*tf.T913
                  + 0.0449976*tf.T9 + -0.00204682*tf.T953 + -0.789217*tf.lnT9)

    rate_eval.n_O15__O16 = rate

@numba.njit()
def p_O15__F16(rate_eval, tf):
    # O15 + p --> F16
    rate = 0.0

    # rpsmr
    rate += np.exp(  -3.01222 + -8.51539*tf.T9i + 74.5955*tf.T913i + -73.794*tf.T913
                  + 2.87936*tf.T9 + -0.117988*tf.T953 + 47.5366*tf.lnT9)

    rate_eval.p_O15__F16 = rate

@numba.njit()
def He4_O15__Ne19(rate_eval, tf):
    # O15 + He4 --> Ne19
    rate = 0.0

    # dc11r
    rate += np.exp(  -32.2496 + -4.20439*tf.T9i + -3.24609*tf.T913i + 44.4647*tf.T913
                  + -9.79962*tf.T9 + 0.841782*tf.T953 + -1.5*tf.lnT9)
    # dc11r
    rate += np.exp(  -0.0452465 + -5.88439*tf.T9i
                  + -1.5*tf.lnT9)
    # dc11n
    rate += np.exp(  26.2914 + -39.578*tf.T913i
                  + -3.0*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_O15__Ne19 = rate

@numba.njit()
def n_O16__O17(rate_eval, tf):
    # O16 + n --> O17
    rate = 0.0

    # ks03 
    rate += np.exp(  7.21546 + 0.0235015*tf.T9i + -2.11246*tf.T913i + 4.87742*tf.T913
                  + -0.314426*tf.T9 + 0.0169515*tf.T953 + -0.984784*tf.lnT9)

    rate_eval.n_O16__O17 = rate

@numba.njit()
def p_O16__F17(rate_eval, tf):
    # O16 + p --> F17
    rate = 0.0

    # ia08n
    rate += np.exp(  19.0904 + -16.696*tf.T913i + -1.16252*tf.T913
                  + 0.267703*tf.T9 + -0.0338411*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_O16__F17 = rate

@numba.njit()
def He4_O16__Ne20(rate_eval, tf):
    # O16 + He4 --> Ne20
    rate = 0.0

    # co10r
    rate += np.exp(  3.88571 + -10.3585*tf.T9i
                  + -1.5*tf.lnT9)
    # co10n
    rate += np.exp(  23.903 + -39.7262*tf.T913i + -0.210799*tf.T913
                  + 0.442879*tf.T9 + -0.0797753*tf.T953 + -0.666667*tf.lnT9)
    # co10r
    rate += np.exp(  9.50848 + -12.7643*tf.T9i + -3.65925*tf.T913
                  + 0.714224*tf.T9 + -0.00107508*tf.T953 + -1.5*tf.lnT9)

    rate_eval.He4_O16__Ne20 = rate

@numba.njit()
def n_O17__O18(rate_eval, tf):
    # O17 + n --> O18
    rate = 0.0

    # bb92n
    rate += np.exp(  4.27026 + -5.69067e-13*tf.T9i + 7.48144e-11*tf.T913i + -2.47239e-10*tf.T913
                  + 2.48052e-11*tf.T9 + -2.07736e-12*tf.T953 + 8.43048e-11*tf.lnT9)
    # bb92r
    rate += np.exp(  20.1798 + -2.12961*tf.T9i + 16.8052*tf.T913i + -30.138*tf.T913
                  + 1.14711*tf.T9 + -0.0220312*tf.T953 + 13.7452*tf.lnT9)

    rate_eval.n_O17__O18 = rate

@numba.njit()
def p_O17__F18(rate_eval, tf):
    # O17 + p --> F18
    rate = 0.0

    # il10r
    rate += np.exp(  9.39048 + -6.22828*tf.T9i + 2.31435*tf.T913
                  + -0.302835*tf.T9 + 0.020133*tf.T953 + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  -13.077 + -0.746296*tf.T9i
                  + -1.5*tf.lnT9)
    # il10n
    rate += np.exp(  15.8929 + -16.4035*tf.T913i + 4.31885*tf.T913
                  + -0.709921*tf.T9 + -2.0*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_O17__F18 = rate

@numba.njit()
def He4_O17__Ne21(rate_eval, tf):
    # O17 + He4 --> Ne21
    rate = 0.0

    # be13r
    rate += np.exp(  -25.0898 + -5.50926*tf.T9i + 123.363*tf.T913i + -87.4351*tf.T913
                  + -3.40974e-06*tf.T9 + -57.0469*tf.T953 + 82.2218*tf.lnT9)
    # be13r
    rate += np.exp(  -117.134 + -13.6759*tf.T9i + 3.31162e-08*tf.T913i + 130.258*tf.T913
                  + -7.92551e-05*tf.T9 + -4.13772*tf.T953 + -42.7753*tf.lnT9)
    # be13r
    rate += np.exp(  2.14 + -5.99952*tf.T9i + 2.87641*tf.T913i + -3.54489*tf.T913
                  + -2.11222e-08*tf.T9 + -3.90649e-09*tf.T953 + 4.75778*tf.lnT9)

    rate_eval.He4_O17__Ne21 = rate

@numba.njit()
def p_O18__F19(rate_eval, tf):
    # O18 + p --> F19
    rate = 0.0

    # il10n
    rate += np.exp(  19.917 + -16.7246*tf.T913i
                  + -3.0*tf.T953 + -0.666667*tf.lnT9)
    # il10r
    rate += np.exp(  7.26876 + -6.7253*tf.T9i + 3.99059*tf.T913
                  + -0.593127*tf.T9 + 0.0877534*tf.T953 + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  5.07648 + -1.65681*tf.T9i
                  + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  -35.0079 + -0.244743*tf.T9i
                  + -1.5*tf.lnT9)

    rate_eval.p_O18__F19 = rate

@numba.njit()
def n_F15__F16(rate_eval, tf):
    # F15 + n --> F16
    rate = 0.0

    # rpsmr
    rate += np.exp(  8.54244 + 0.0214707*tf.T9i + -1.8316*tf.T913i + 2.8333*tf.T913
                  + -0.0925718*tf.T9 + 0.00805423*tf.T953 + -1.44099*tf.lnT9)

    rate_eval.n_F15__F16 = rate

@numba.njit()
def n_F16__F17(rate_eval, tf):
    # F16 + n --> F17
    rate = 0.0

    # rpsmr
    rate += np.exp(  9.36454 + 0.00885357*tf.T9i + -0.79526*tf.T913i + 1.84363*tf.T913
                  + -0.123684*tf.T9 + 0.0149635*tf.T953 + -0.704695*tf.lnT9)

    rate_eval.n_F16__F17 = rate

@numba.njit()
def n_F17__F18(rate_eval, tf):
    # F17 + n --> F18
    rate = 0.0

    # rpsmr
    rate += np.exp(  3.42798 + 0.0273799*tf.T9i + -2.71353*tf.T913i + 6.19871*tf.T913
                  + -0.26836*tf.T9 + 0.0115489*tf.T953 + -2.52002*tf.lnT9)

    rate_eval.n_F17__F18 = rate

@numba.njit()
def p_F17__Ne18(rate_eval, tf):
    # F17 + p --> Ne18
    rate = 0.0

    # cb09 
    rate += np.exp(  -7.84708 + -0.0323504*tf.T9i + -14.2191*tf.T913i + 34.0647*tf.T913
                  + -16.5698*tf.T9 + 2.48116*tf.T953 + -2.13376*tf.lnT9)
    # cb09 
    rate += np.exp(  27.5778 + -4.95969*tf.T9i + -21.3249*tf.T913i + -0.230774*tf.T913
                  + 0.917931*tf.T9 + -0.0440377*tf.T953 + -7.36014*tf.lnT9)

    rate_eval.p_F17__Ne18 = rate

@numba.njit()
def He4_F17__Na21(rate_eval, tf):
    # F17 + He4 --> Na21
    rate = 0.0

    # rpsmr
    rate += np.exp(  41.1529 + -1.72817*tf.T9i + 15.559*tf.T913i + -68.3231*tf.T913
                  + 2.54275*tf.T9 + -0.0989207*tf.T953 + 36.8877*tf.lnT9)

    rate_eval.He4_F17__Na21 = rate

@numba.njit()
def n_F18__F19(rate_eval, tf):
    # F18 + n --> F19
    rate = 0.0

    # rpsmr
    rate += np.exp(  0.814964 + 0.0604483*tf.T9i + -6.65988*tf.T913i + 18.3243*tf.T913
                  + -1.47263*tf.T9 + 0.0955082*tf.T953 + -6.75505*tf.lnT9)

    rate_eval.n_F18__F19 = rate

@numba.njit()
def p_F18__Ne19(rate_eval, tf):
    # F18 + p --> Ne19
    rate = 0.0

    # il10r
    rate += np.exp(  -5.85727 + -2.89147*tf.T9i + 13.1683*tf.T913
                  + -1.92023*tf.T9 + 0.16901*tf.T953 + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  -29.449 + -0.39895*tf.T9i + 22.4903*tf.T913
                  + 0.307872*tf.T9 + -0.296226*tf.T953 + -1.5*tf.lnT9)
    # il10n
    rate += np.exp(  57.4084 + -21.4023*tf.T913i + -93.766*tf.T913
                  + 179.258*tf.T9 + -202.561*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_F18__Ne19 = rate

@numba.njit()
def He4_F18__Na22(rate_eval, tf):
    # F18 + He4 --> Na22
    rate = 0.0

    # rpsmr
    rate += np.exp(  35.3786 + -1.82957*tf.T9i + 18.8956*tf.T913i + -65.6134*tf.T913
                  + 1.71114*tf.T9 + -0.0260999*tf.T953 + 37.8396*tf.lnT9)

    rate_eval.He4_F18__Na22 = rate

@numba.njit()
def p_F19__Ne20(rate_eval, tf):
    # F19 + p --> Ne20
    rate = 0.0

    # nacrr
    rate += np.exp(  -5.63093 + -7.74414*tf.T9i + 31.6442*tf.T913i + -58.6563*tf.T913
                  + 67.7365*tf.T9 + -22.9721*tf.T953 + -1.5*tf.lnT9)
    # nacrr
    rate += np.exp(  12.3816 + -1.71383*tf.T9i + -11.3832*tf.T913i + 5.47872*tf.T913
                  + -1.07203*tf.T9 + 0.11196*tf.T953 + -1.5*tf.lnT9)
    # nacrn
    rate += np.exp(  18.2807 + -18.116*tf.T913i + -1.4622*tf.T913
                  + 6.95113*tf.T9 + -2.90366*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_F19__Ne20 = rate

@numba.njit()
def He4_F19__Na23(rate_eval, tf):
    # F19 + He4 --> Na23
    rate = 0.0

    # rpsmr
    rate += np.exp(  52.7856 + -2.11408*tf.T9i + 39.7219*tf.T913i + -100.401*tf.T913
                  + 3.15808*tf.T9 + -0.0629822*tf.T953 + 54.4823*tf.lnT9)

    rate_eval.He4_F19__Na23 = rate

@numba.njit()
def n_Ne18__Ne19(rate_eval, tf):
    # Ne18 + n --> Ne19
    rate = 0.0

    # ths8r
    rate += np.exp(  9.249 + 0.937162*tf.T913
                  + -0.0221952*tf.T9 + -0.00101206*tf.T953)

    rate_eval.n_Ne18__Ne19 = rate

@numba.njit()
def He4_Ne18__Mg22(rate_eval, tf):
    # Ne18 + He4 --> Mg22
    rate = 0.0

    # ths8r
    rate += np.exp(  32.8865 + -46.4859*tf.T913i + 0.956741*tf.T913
                  + -0.914402*tf.T9 + 0.0722478*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Ne18__Mg22 = rate

@numba.njit()
def n_Ne19__Ne20(rate_eval, tf):
    # Ne19 + n --> Ne20
    rate = 0.0

    # ths8r
    rate += np.exp(  6.40633 + 1.57592*tf.T913
                  + -0.11175*tf.T9 + 0.00226473*tf.T953)

    rate_eval.n_Ne19__Ne20 = rate

@numba.njit()
def He4_Ne19__Mg23(rate_eval, tf):
    # Ne19 + He4 --> Mg23
    rate = 0.0

    # ths8r
    rate += np.exp(  37.1998 + -46.6346*tf.T913i + -1.1007*tf.T913
                  + -0.794097*tf.T9 + 0.0813036*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Ne19__Mg23 = rate

@numba.njit()
def n_Ne20__Ne21(rate_eval, tf):
    # Ne20 + n --> Ne21
    rate = 0.0

    # ka02r
    rate += np.exp(  12.7344 + -1.70393*tf.T9i
                  + -1.5*tf.lnT9)
    # ka02n
    rate += np.exp(  8.57651)

    rate_eval.n_Ne20__Ne21 = rate

@numba.njit()
def p_Ne20__Na21(rate_eval, tf):
    # Ne20 + p --> Na21
    rate = 0.0

    # ly18 
    rate += np.exp(  195297.0 + -61.14*tf.T9i + 21894.7*tf.T913i + -319153.0*tf.T913
                  + 224369.0*tf.T9 + -188049.0*tf.T953 + 48703.4*tf.lnT9)
    # ly18 
    rate += np.exp(  207.877 + -0.152711*tf.T9i + 15.325*tf.T913i + -294.859*tf.T913
                  + 107.692*tf.T9 + -46.2072*tf.T953 + 57.8398*tf.lnT9)
    # ly18 
    rate += np.exp(  5.83103 + -8.838*tf.T9i + 20.5893*tf.T913i + -17.5841*tf.T913
                  + 0.243226*tf.T9 + -0.000231418*tf.T953 + 12.8398*tf.lnT9)
    # ly18 
    rate += np.exp(  230.019 + -4.45358*tf.T9i + 258.57*tf.T913i + -506.387*tf.T913
                  + 22.1576*tf.T9 + -0.721182*tf.T953 + 230.288*tf.lnT9)

    rate_eval.p_Ne20__Na21 = rate

@numba.njit()
def He4_Ne20__Mg24(rate_eval, tf):
    # Ne20 + He4 --> Mg24
    rate = 0.0

    # il10r
    rate += np.exp(  -8.79827 + -12.7809*tf.T9i + 16.9229*tf.T913
                  + -2.57325*tf.T9 + 0.208997*tf.T953 + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  1.98307 + -9.22026*tf.T9i
                  + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  -38.7055 + -2.50605*tf.T9i
                  + -1.5*tf.lnT9)
    # il10n
    rate += np.exp(  24.5058 + -46.2525*tf.T913i + 5.58901*tf.T913
                  + 7.61843*tf.T9 + -3.683*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Ne20__Mg24 = rate

@numba.njit()
def p_Ne21__Na22(rate_eval, tf):
    # Ne21 + p --> Na22
    rate = 0.0

    # il10n
    rate += np.exp(  19.0696 + -19.2096*tf.T913i
                  + -1.0*tf.T953 + -0.666667*tf.lnT9)
    # il10r
    rate += np.exp(  -39.4862 + -4.21385*tf.T9i + 21.1176*tf.T913i + 34.0411*tf.T913
                  + -4.45593*tf.T9 + 0.328613*tf.T953 + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  1.75704 + -1.39957*tf.T9i
                  + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  -47.6554 + -0.19618*tf.T9i
                  + -1.5*tf.lnT9)

    rate_eval.p_Ne21__Na22 = rate

@numba.njit()
def He4_Ne21__Mg25(rate_eval, tf):
    # Ne21 + He4 --> Mg25
    rate = 0.0

    # cf88r
    rate += np.exp(  26.2429 + -22.049*tf.T9i + -29.4583*tf.T913
                  + 14.6328*tf.T9 + -3.47392*tf.T953 + -1.5*tf.lnT9)
    # cf88n
    rate += np.exp(  36.6927 + -46.89*tf.T913i + -0.72642*tf.T913
                  + -0.76406*tf.T9 + 0.0797483*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Ne21__Mg25 = rate

@numba.njit()
def n_Na21__Na22(rate_eval, tf):
    # Na21 + n --> Na22
    rate = 0.0

    # ths8r
    rate += np.exp(  7.58321 + 1.31089*tf.T913
                  + -0.164931*tf.T9 + 0.00903374*tf.T953)

    rate_eval.n_Na21__Na22 = rate

@numba.njit()
def p_Na21__Mg22(rate_eval, tf):
    # Na21 + p --> Mg22
    rate = 0.0

    # il10r
    rate += np.exp(  6.28009 + -7.47826*tf.T9i + 4.75873*tf.T913
                  + -0.708067*tf.T9 + 0.0523665*tf.T953 + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  5.07923 + -2.38741*tf.T9i
                  + -1.5*tf.lnT9)
    # il10n
    rate += np.exp(  18.6558 + -20.7329*tf.T913i
                  + -1.5*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Na21__Mg22 = rate

@numba.njit()
def He4_Na21__Al25(rate_eval, tf):
    # Na21 + He4 --> Al25
    rate = 0.0

    # ths8r
    rate += np.exp(  35.3006 + -49.9709*tf.T913i + 1.63835*tf.T913
                  + -1.18562*tf.T9 + 0.101965*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Na21__Al25 = rate

@numba.njit()
def n_Na22__Na23(rate_eval, tf):
    # Na22 + n --> Na23
    rate = 0.0

    # ths8r
    rate += np.exp(  12.8678 + 1.02148*tf.T913
                  + -0.334638*tf.T9 + 0.0258708*tf.T953)

    rate_eval.n_Na22__Na23 = rate

@numba.njit()
def p_Na22__Mg23(rate_eval, tf):
    # Na22 + p --> Mg23
    rate = 0.0

    # il10r
    rate += np.exp(  -11.2731 + -2.42669*tf.T9i + 4.86658*tf.T913i + 16.4592*tf.T913
                  + -1.95129*tf.T9 + 0.132972*tf.T953 + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  -16.2423 + -0.777841*tf.T9i
                  + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  -25.2739 + -0.499888*tf.T9i
                  + -1.5*tf.lnT9)

    rate_eval.p_Na22__Mg23 = rate

@numba.njit()
def He4_Na22__Al26(rate_eval, tf):
    # Na22 + He4 --> Al26
    rate = 0.0

    # ths8r
    rate += np.exp(  36.3797 + -50.0924*tf.T913i + -0.390826*tf.T913
                  + -0.99531*tf.T9 + 0.101354*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Na22__Al26 = rate

@numba.njit()
def n_Na23__Na24(rate_eval, tf):
    # Na23 + n --> Na24
    rate = 0.0

    # ks03 
    rate += np.exp(  21.496 + -0.185436*tf.T9i + 16.0285*tf.T913i + -25.7668*tf.T913
                  + 1.07385*tf.T9 + -0.0293163*tf.T953 + 12.6777*tf.lnT9)

    rate_eval.n_Na23__Na24 = rate

@numba.njit()
def p_Na23__Mg24(rate_eval, tf):
    # Na23 + p --> Mg24
    rate = 0.0

    # il10r
    rate += np.exp(  9.0594 + -3.28029*tf.T9i + -0.360588*tf.T913
                  + 1.4187*tf.T9 + -0.184061*tf.T953 + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  -5.02585 + -1.61219*tf.T9i
                  + -1.5*tf.lnT9)
    # il10n
    rate += np.exp(  18.9075 + -20.6428*tf.T913i + 1.52954*tf.T913
                  + 2.7487*tf.T9 + -1.0*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Na23__Mg24 = rate

@numba.njit()
def He4_Na23__Al27(rate_eval, tf):
    # Na23 + He4 --> Al27
    rate = 0.0

    # ths8r
    rate += np.exp(  44.7724 + -50.2042*tf.T913i + -1.64239*tf.T913
                  + -1.59995*tf.T9 + 0.184933*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Na23__Al27 = rate

@numba.njit()
def p_Na24__Mg25(rate_eval, tf):
    # Na24 + p --> Mg25
    rate = 0.0

    # rath 
    rate += np.exp(  63.7099 + -2.54619*tf.T9i + 88.4533*tf.T913i + -156.066*tf.T913
                  + 8.77114*tf.T9 + -0.500891*tf.T953 + 77.2574*tf.lnT9)

    rate_eval.p_Na24__Mg25 = rate

@numba.njit()
def He4_Na24__Al28(rate_eval, tf):
    # Na24 + He4 --> Al28
    rate = 0.0

    # rath 
    rate += np.exp(  1.63055 + -2.00108*tf.T9i + -22.0373*tf.T913i + 12.2654*tf.T913
                  + -3.20073*tf.T9 + 0.258909*tf.T953 + 4.52728*tf.lnT9)

    rate_eval.He4_Na24__Al28 = rate

@numba.njit()
def n_Mg22__Mg23(rate_eval, tf):
    # Mg22 + n --> Mg23
    rate = 0.0

    # ths8r
    rate += np.exp(  9.32079 + 2.72435*tf.T913
                  + -0.470498*tf.T9 + 0.0357362*tf.T953)

    rate_eval.n_Mg22__Mg23 = rate

@numba.njit()
def He4_Mg22__Si26(rate_eval, tf):
    # Mg22 + He4 --> Si26
    rate = 0.0

    # ths8r
    rate += np.exp(  35.9065 + -53.0847*tf.T913i + 0.625592*tf.T913
                  + -0.893632*tf.T9 + 0.0691552*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Mg22__Si26 = rate

@numba.njit()
def n_Mg23__Mg24(rate_eval, tf):
    # Mg23 + n --> Mg24
    rate = 0.0

    # ths8r
    rate += np.exp(  7.00613 + 2.66964*tf.T913
                  + -0.448904*tf.T9 + 0.0326505*tf.T953)

    rate_eval.n_Mg23__Mg24 = rate

@numba.njit()
def n_Mg24__Mg25(rate_eval, tf):
    # Mg24 + n --> Mg25
    rate = 0.0

    # ks03 
    rate += np.exp(  64.622 + 0.161296*tf.T9i + -0.142939*tf.T913i + -57.7499*tf.T913
                  + 7.01981*tf.T9 + -0.582057*tf.T953 + 12.8133*tf.lnT9)

    rate_eval.n_Mg24__Mg25 = rate

@numba.njit()
def p_Mg24__Al25(rate_eval, tf):
    # Mg24 + p --> Al25
    rate = 0.0

    # il10r
    rate += np.exp(  8.24021 + -2.48451*tf.T9i + -1.57811*tf.T913
                  + 1.52232*tf.T9 + -0.183001*tf.T953 + -1.5*tf.lnT9)
    # il10n
    rate += np.exp(  19.8966 + -22.0227*tf.T913i + 0.361297*tf.T913
                  + 2.61292*tf.T9 + -1.0*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Mg24__Al25 = rate

@numba.njit()
def He4_Mg24__Si28(rate_eval, tf):
    # Mg24 + He4 --> Si28
    rate = 0.0

    # st08r
    rate += np.exp(  8.03977 + -15.629*tf.T9i
                  + -1.5*tf.lnT9)
    # st08r
    rate += np.exp(  -50.5494 + -12.8332*tf.T9i + 21.3721*tf.T913i + 37.7649*tf.T913
                  + -4.10635*tf.T9 + 0.249618*tf.T953 + -1.5*tf.lnT9)

    rate_eval.He4_Mg24__Si28 = rate

@numba.njit()
def n_Mg25__Mg26(rate_eval, tf):
    # Mg25 + n --> Mg26
    rate = 0.0

    # ks03 
    rate += np.exp(  38.34 + -0.0457591*tf.T9i + 9.392*tf.T913i + -36.6784*tf.T913
                  + 3.09567*tf.T9 + -0.223882*tf.T953 + 12.3852*tf.lnT9)

    rate_eval.n_Mg25__Mg26 = rate

@numba.njit()
def p_Mg25__Al26(rate_eval, tf):
    # Mg25 + p --> Al26
    rate = 0.0

    # il10r
    rate += np.exp(  4.21826 + -0.71983*tf.T9i + -88.9297*tf.T913
                  + 302.948*tf.T9 + -346.461*tf.T953 + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  -37.1963 + -0.429366*tf.T9i
                  + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  2.22778 + -3.22353*tf.T9i + 8.46334*tf.T913
                  + -0.907024*tf.T9 + 0.0642981*tf.T953 + -1.5*tf.lnT9)

    rate_eval.p_Mg25__Al26 = rate

@numba.njit()
def He4_Mg25__Si29(rate_eval, tf):
    # Mg25 + He4 --> Si29
    rate = 0.0

    # cf88n
    rate += np.exp(  40.3715 + -53.41*tf.T913i + -1.83266*tf.T913
                  + -0.573073*tf.T9 + 0.0462678*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Mg25__Si29 = rate

@numba.njit()
def p_Mg26__Al27(rate_eval, tf):
    # Mg26 + p --> Al27
    rate = 0.0

    # il10r
    rate += np.exp(  -27.2168 + -0.888689*tf.T9i + 35.6312*tf.T913
                  + -5.27265*tf.T9 + 0.392932*tf.T953 + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  -84.493 + -0.469464*tf.T9i + 251.281*tf.T913
                  + -730.009*tf.T9 + -224.016*tf.T953 + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  5.26056 + -3.35921*tf.T9i + 6.78105*tf.T913
                  + -1.25771*tf.T9 + 0.140754*tf.T953 + -1.5*tf.lnT9)

    rate_eval.p_Mg26__Al27 = rate

@numba.njit()
def He4_Mg26__Si30(rate_eval, tf):
    # Mg26 + He4 --> Si30
    rate = 0.0

    # cf88r
    rate += np.exp(  1.32941 + -18.7164*tf.T9i + -1.87411*tf.T913
                  + 3.41299*tf.T9 + -0.43226*tf.T953 + -1.5*tf.lnT9)
    # cf88n
    rate += np.exp(  45.8787 + -53.7518*tf.T913i + -4.8647*tf.T913
                  + -1.51467*tf.T9 + -0.666667*tf.lnT9)

    rate_eval.He4_Mg26__Si30 = rate

@numba.njit()
def n_Al25__Al26(rate_eval, tf):
    # Al25 + n --> Al26
    rate = 0.0

    # ths8r
    rate += np.exp(  7.92587 + 1.17141*tf.T913
                  + -0.162515*tf.T9 + 0.0126275*tf.T953)

    rate_eval.n_Al25__Al26 = rate

@numba.njit()
def p_Al25__Si26(rate_eval, tf):
    # Al25 + p --> Si26
    rate = 0.0

    # li20r
    rate += np.exp(  8.74592 + -4.78862*tf.T9i
                  + -1.48545*tf.lnT9)
    # li20r
    rate += np.exp(  5.38793 + -12.4587*tf.T9i
                  + 1.48721*tf.lnT9)
    # li20r
    rate += np.exp(  -6.20781 + -1.73102*tf.T9i
                  + -0.244294*tf.lnT9)

    rate_eval.p_Al25__Si26 = rate

@numba.njit()
def He4_Al25__P29(rate_eval, tf):
    # Al25 + He4 --> P29
    rate = 0.0

    # ths8r
    rate += np.exp(  37.9099 + -56.3424*tf.T913i + 0.542998*tf.T913
                  + -0.721716*tf.T9 + 0.0469712*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Al25__P29 = rate

@numba.njit()
def n_Al26__Al27(rate_eval, tf):
    # Al26 + n --> Al27
    rate = 0.0

    # ks03 
    rate += np.exp(  14.7625 + 0.00350938*tf.T9i + -0.171158*tf.T913i + -1.77283*tf.T913
                  + 0.206192*tf.T9 + -0.0191705*tf.T953 + 0.139609*tf.lnT9)

    rate_eval.n_Al26__Al27 = rate

@numba.njit()
def He4_Al26__P30(rate_eval, tf):
    # Al26 + He4 --> P30
    rate = 0.0

    # ths8r
    rate += np.exp(  42.9778 + -56.4422*tf.T913i + -2.44848*tf.T913
                  + -1.17578*tf.T9 + 0.150757*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Al26__P30 = rate

@numba.njit()
def n_Al27__Al28(rate_eval, tf):
    # Al27 + n --> Al28
    rate = 0.0

    # ks03 
    rate += np.exp(  19.7971 + -0.0784751*tf.T9i + 7.49524*tf.T913i + -14.6211*tf.T913
                  + 0.721304*tf.T9 + -0.0288351*tf.T953 + 6.58346*tf.lnT9)

    rate_eval.n_Al27__Al28 = rate

@numba.njit()
def p_Al27__Si28(rate_eval, tf):
    # Al27 + p --> Si28
    rate = 0.0

    # il10r
    rate += np.exp(  -13.6664 + -1.90396*tf.T9i + 23.8634*tf.T913
                  + -3.70135*tf.T9 + 0.28964*tf.T953 + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  86.0234 + -0.387313*tf.T9i + -26.8327*tf.T913i + -116.137*tf.T913
                  + 0.00950567*tf.T9 + 0.00999755*tf.T953 + -1.5*tf.lnT9)
    # il10n
    rate += np.exp(  21.1065 + -23.2205*tf.T913i
                  + -2.0*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Al27__Si28 = rate

@numba.njit()
def He4_Al27__P31(rate_eval, tf):
    # Al27 + He4 --> P31
    rate = 0.0

    # ths8r
    rate += np.exp(  47.2333 + -56.5351*tf.T913i + -0.896208*tf.T913
                  + -1.72024*tf.T9 + 0.185409*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Al27__P31 = rate

@numba.njit()
def p_Al28__Si29(rate_eval, tf):
    # Al28 + p --> Si29
    rate = 0.0

    # rath 
    rate += np.exp(  82.7866 + -3.13386*tf.T9i + 114.903*tf.T913i + -203.401*tf.T913
                  + 11.5793*tf.T9 + -0.672783*tf.T953 + 99.4177*tf.lnT9)

    rate_eval.p_Al28__Si29 = rate

@numba.njit()
def He4_Al28__P32(rate_eval, tf):
    # Al28 + He4 --> P32
    rate = 0.0

    # rath 
    rate += np.exp(  -0.522392 + -2.59086*tf.T9i + -14.2992*tf.T913i + 6.40835*tf.T913
                  + -3.59742*tf.T9 + 0.308804*tf.T953 + 10.094*tf.lnT9)

    rate_eval.He4_Al28__P32 = rate

@numba.njit()
def n_Si28__Si29(rate_eval, tf):
    # Si28 + n --> Si29
    rate = 0.0

    # ka02r
    rate += np.exp(  6.9158 + -0.38*tf.T9i + 7.68863*tf.T913
                  + -1.7991*tf.T9 + -1.5*tf.lnT9)
    # ka02 
    rate += np.exp(  8.77553)

    rate_eval.n_Si28__Si29 = rate

@numba.njit()
def p_Si28__P29(rate_eval, tf):
    # Si28 + p --> P29
    rate = 0.0

    # il10n
    rate += np.exp(  16.1779 + -23.8173*tf.T913i + 7.08203*tf.T913
                  + -1.44753*tf.T9 + 0.0804296*tf.T953 + -0.666667*tf.lnT9)
    # il10r
    rate += np.exp(  5.73975 + -4.14232*tf.T9i
                  + -1.5*tf.lnT9)

    rate_eval.p_Si28__P29 = rate

@numba.njit()
def He4_Si28__S32(rate_eval, tf):
    # Si28 + He4 --> S32
    rate = 0.0

    # ths8r
    rate += np.exp(  47.9212 + -59.4896*tf.T913i + 4.47205*tf.T913
                  + -4.78989*tf.T9 + 0.557201*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Si28__S32 = rate

@numba.njit()
def n_Si29__Si30(rate_eval, tf):
    # Si29 + n --> Si30
    rate = 0.0

    # ka02r
    rate += np.exp(  9.60115 + -0.179366*tf.T9i + 5.50678*tf.T913
                  + -2.85656*tf.T9 + -1.5*tf.lnT9)
    # ka02n
    rate += np.exp(  11.8023 + 0.650904*tf.T913)

    rate_eval.n_Si29__Si30 = rate

@numba.njit()
def p_Si29__P30(rate_eval, tf):
    # Si29 + p --> P30
    rate = 0.0

    # il10r
    rate += np.exp(  -1.24791 + -3.33929*tf.T9i + 14.0921*tf.T913
                  + -3.92096*tf.T9 + 0.447706*tf.T953 + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  -13.4701 + -1.25026*tf.T9i
                  + -1.5*tf.lnT9)
    # il10n
    rate += np.exp(  16.5182 + -23.9101*tf.T913i + 10.7796*tf.T913
                  + -3.04181*tf.T9 + 0.274565*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Si29__P30 = rate

@numba.njit()
def p_Si30__P31(rate_eval, tf):
    # Si30 + p --> P31
    rate = 0.0

    # de20r
    rate += np.exp(  12.3695 + -6.64105*tf.T9i
                  + -1.1191*tf.lnT9)
    # de20r
    rate += np.exp(  8.79766 + -5.18231*tf.T9i
                  + 1.2883*tf.lnT9)
    # de20r
    rate += np.exp(  -334.266 + -1.13327*tf.T9i
                  + -78.547*tf.lnT9)
    # de20r
    rate += np.exp(  -4.24208 + -1.25174*tf.T9i
                  + 4.99034*tf.lnT9)
    # de20r
    rate += np.exp(  -18.9198 + -0.945261*tf.T9i
                  + 1.30331*tf.lnT9)
    # de20r
    rate += np.exp(  -1138.34 + -95.8769*tf.T9i
                  + -896.758*tf.lnT9)
    # de20r
    rate += np.exp(  9.96544 + -5.58963*tf.T9i
                  + -1.57082*tf.lnT9)
    # de20r
    rate += np.exp(  -34.8594 + -0.592934*tf.T9i
                  + -1.62839*tf.lnT9)

    rate_eval.p_Si30__P31 = rate

@numba.njit()
def He4_Si30__S34(rate_eval, tf):
    # Si30 + He4 --> S34
    rate = 0.0

    # ths8r
    rate += np.exp(  47.5008 + -59.6559*tf.T913i + 3.70141*tf.T913
                  + -3.12537*tf.T9 + 0.307626*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Si30__S34 = rate

@numba.njit()
def n_P29__P30(rate_eval, tf):
    # P29 + n --> P30
    rate = 0.0

    # ths8r
    rate += np.exp(  8.78841 + 0.15555*tf.T913
                  + 0.155359*tf.T9 + -0.0208019*tf.T953)

    rate_eval.n_P29__P30 = rate

@numba.njit()
def n_P30__P31(rate_eval, tf):
    # P30 + n --> P31
    rate = 0.0

    # ths8r
    rate += np.exp(  12.8187 + 0.909911*tf.T913
                  + -0.162367*tf.T9 + 0.00668293*tf.T953)

    rate_eval.n_P30__P31 = rate

@numba.njit()
def p_P30__S31(rate_eval, tf):
    # P30 + p --> S31
    rate = 0.0

    # mb07 
    rate += np.exp(  -7756.63 + 10.5093*tf.T9i + -1999.51*tf.T913i + 11886.5*tf.T913
                  + -2668.72*tf.T9 + 354.294*tf.T953 + -2900.45*tf.lnT9)
    # mb07 
    rate += np.exp(  24.3866 + -8.98316*tf.T9i + 156.029*tf.T913i + -174.377*tf.T913
                  + 7.4644*tf.T9 + -0.342232*tf.T953 + 99.2579*tf.lnT9)

    rate_eval.p_P30__S31 = rate

@numba.njit()
def He4_P30__Cl34(rate_eval, tf):
    # P30 + He4 --> Cl34
    rate = 0.0

    # ths8r
    rate += np.exp(  45.3321 + -62.4643*tf.T913i + -3.19028*tf.T913
                  + -0.832633*tf.T9 + 0.0987525*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_P30__Cl34 = rate

@numba.njit()
def n_P31__P32(rate_eval, tf):
    # P31 + n --> P32
    rate = 0.0

    # ka02r
    rate += np.exp(  9.48567 + -0.246803*tf.T9i + 2.13185*tf.T913
                  + -1.5*tf.lnT9)
    # ka02n
    rate += np.exp(  12.5332 + 0.327031*tf.T913
                  + 1.0*tf.lnT9)

    rate_eval.n_P31__P32 = rate

@numba.njit()
def p_P31__S32(rate_eval, tf):
    # P31 + p --> S32
    rate = 0.0

    # il10r
    rate += np.exp(  0.821556 + -3.77704*tf.T9i + 8.09341*tf.T913
                  + -0.615971*tf.T9 + 0.031159*tf.T953 + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  -2.66839 + -2.25958*tf.T9i
                  + -1.5*tf.lnT9)
    # il10n
    rate += np.exp(  19.2596 + -25.3278*tf.T913i + 6.4931*tf.T913
                  + -9.27513*tf.T9 + -0.610439*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_P31__S32 = rate

@numba.njit()
def He4_P31__Cl35(rate_eval, tf):
    # P31 + He4 --> Cl35
    rate = 0.0

    # ths8r
    rate += np.exp(  50.451 + -62.5433*tf.T913i + -2.95026*tf.T913
                  + -0.89652*tf.T9 + 0.0774126*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_P31__Cl35 = rate

@numba.njit()
def n_S31__S32(rate_eval, tf):
    # S31 + n --> S32
    rate = 0.0

    # ths8r
    rate += np.exp(  7.56582 + 1.71463*tf.T913
                  + -0.221804*tf.T9 + 0.00880104*tf.T953)

    rate_eval.n_S31__S32 = rate

@numba.njit()
def p_S31__Cl32(rate_eval, tf):
    # S31 + p --> Cl32
    rate = 0.0

    # mm11r
    rate += np.exp(  -11.7529 + -1.77549*tf.T9i
                  + -1.5*tf.lnT9)
    # mm11n
    rate += np.exp(  6.30412 + -25.3278*tf.T913i + 42.7226*tf.T913
                  + -102.23*tf.T9 + -2.0*tf.T953 + -0.666667*tf.lnT9)
    # mm11r
    rate += np.exp(  3.272 + -6.35927*tf.T9i + 3.202*tf.T913
                  + 0.101684*tf.T9 + -0.0334894*tf.T953 + -1.5*tf.lnT9)

    rate_eval.p_S31__Cl32 = rate

@numba.njit()
def He4_S31__Ar35(rate_eval, tf):
    # S31 + He4 --> Ar35
    rate = 0.0

    # ths8r
    rate += np.exp(  56.1418 + -65.2934*tf.T913i + -13.3991*tf.T913
                  + 1.97512*tf.T9 + -0.208606*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_S31__Ar35 = rate

@numba.njit()
def p_S34__Cl35(rate_eval, tf):
    # S34 + p --> Cl35
    rate = 0.0

    # se19r
    rate += np.exp(  -7.70046 + -2.02615*tf.T9i
                  + -1.5*tf.lnT9)
    # se19r
    rate += np.exp(  8.25443 + -5.75583*tf.T9i
                  + -1.5*tf.lnT9)
    # se19r
    rate += np.exp(  -14.8781 + -1.40414*tf.T9i
                  + -1.5*tf.lnT9)
    # se19r
    rate += np.exp(  7.00142 + -5.10598*tf.T9i
                  + -1.5*tf.lnT9)
    # se19r
    rate += np.exp(  -17.9366 + -1.17205*tf.T9i
                  + -1.5*tf.lnT9)
    # se19r
    rate += np.exp(  5.80604 + -4.6418*tf.T9i
                  + -1.5*tf.lnT9)
    # se19r
    rate += np.exp(  0.43449 + -3.16803*tf.T9i
                  + -1.5*tf.lnT9)
    # se19r
    rate += np.exp(  -60.0247 + -27.689*tf.T9i + 454.53*tf.T913i + -375.557*tf.T913
                  + 12.6533*tf.T9 + -0.408677*tf.T953 + 245.118*tf.lnT9)

    rate_eval.p_S34__Cl35 = rate

@numba.njit()
def n_Cl34__Cl35(rate_eval, tf):
    # Cl34 + n --> Cl35
    rate = 0.0

    # ths8r
    rate += np.exp(  12.6539 + 0.990222*tf.T913
                  + -0.146686*tf.T9 + 0.00560251*tf.T953)

    rate_eval.n_Cl34__Cl35 = rate

@numba.njit()
def p_Cl34__Ar35(rate_eval, tf):
    # Cl34 + p --> Ar35
    rate = 0.0

    # ths8r
    rate += np.exp(  34.946 + -27.8896*tf.T913i + -0.491144*tf.T913
                  + -2.47595*tf.T9 + 0.34035*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Cl34__Ar35 = rate

@numba.njit()
def d_d__n_He3(rate_eval, tf):
    # d + d --> n + He3
    rate = 0.0

    # gi17n
    rate += np.exp(  19.0876 + -0.00019002*tf.T9i + -4.2292*tf.T913i + 1.6932*tf.T913
                  + -0.0855529*tf.T9 + -1.35709e-25*tf.T953 + -0.734513*tf.lnT9)

    rate_eval.d_d__n_He3 = rate

@numba.njit()
def d_d__p_t(rate_eval, tf):
    # d + d --> p + t
    rate = 0.0

    # go17n
    rate += np.exp(  18.8052 + 4.36209e-05*tf.T9i + -4.32296*tf.T913i + 1.91572*tf.T913
                  + -0.081562*tf.T9 + -3.28804e-22*tf.T953 + -0.879518*tf.lnT9)

    rate_eval.d_d__p_t = rate

@numba.njit()
def p_t__n_He3(rate_eval, tf):
    # t + p --> n + He3
    rate = 0.0

    # de04 
    rate += np.exp(  20.3787 + -8.86352*tf.T9i + -0.332788*tf.T913
                  + -0.700485*tf.T9 + 0.0976521*tf.T953)
    # de04 
    rate += np.exp(  19.2762 + -8.86352*tf.T9i + 0.0438557*tf.T913
                  + -0.201527*tf.T9 + 0.0153433*tf.T953 + 1.0*tf.lnT9)

    rate_eval.p_t__n_He3 = rate

@numba.njit()
def p_t__d_d(rate_eval, tf):
    # t + p --> d + d
    rate = 0.0

    # go17n
    rate += np.exp(  19.3545 + -46.799*tf.T9i + -4.32296*tf.T913i + 1.91572*tf.T913
                  + -0.081562*tf.T9 + -3.28804e-22*tf.T953 + -0.879518*tf.lnT9)

    rate_eval.p_t__d_d = rate

@numba.njit()
def d_t__n_He4(rate_eval, tf):
    # t + d --> n + He4
    rate = 0.0

    # de04 
    rate += np.exp(  39.3457 + -4.5244*tf.T913i + -16.4028*tf.T913
                  + 1.73103*tf.T9 + -0.122966*tf.T953 + 2.31304*tf.lnT9)
    # de04 
    rate += np.exp(  25.1794 + -4.5244*tf.T913i + 0.350337*tf.T913
                  + 0.58747*tf.T9 + -8.84909*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.d_t__n_He4 = rate

@numba.njit()
def He4_t__n_Li6(rate_eval, tf):
    # t + He4 --> n + Li6
    rate = 0.0

    # cf88r
    rate += np.exp(  21.7239 + -57.884*tf.T9i
                  + -1.5*tf.lnT9)
    # cf88n
    rate += np.exp(  19.0085 + -55.494*tf.T9i)

    rate_eval.He4_t__n_Li6 = rate

@numba.njit()
def n_He3__p_t(rate_eval, tf):
    # He3 + n --> p + t
    rate = 0.0

    # de04 
    rate += np.exp(  19.2762 + 0.0438557*tf.T913
                  + -0.201527*tf.T9 + 0.0153433*tf.T953 + 1.0*tf.lnT9)
    # de04 
    rate += np.exp(  20.3787 + -0.332788*tf.T913
                  + -0.700485*tf.T9 + 0.0976521*tf.T953)

    rate_eval.n_He3__p_t = rate

@numba.njit()
def n_He3__d_d(rate_eval, tf):
    # He3 + n --> d + d
    rate = 0.0

    # gi17n
    rate += np.exp(  19.6369 + -37.9358*tf.T9i + -4.2292*tf.T913i + 1.6932*tf.T913
                  + -0.0855529*tf.T9 + -1.35709e-25*tf.T953 + -0.734513*tf.lnT9)

    rate_eval.n_He3__d_d = rate

@numba.njit()
def d_He3__p_He4(rate_eval, tf):
    # He3 + d --> p + He4
    rate = 0.0

    # de04 
    rate += np.exp(  41.2969 + -7.182*tf.T913i + -17.1349*tf.T913
                  + 1.36908*tf.T9 + -0.0814423*tf.T953 + 3.35395*tf.lnT9)
    # de04 
    rate += np.exp(  24.6839 + -7.182*tf.T913i + 0.473288*tf.T913
                  + 1.46847*tf.T9 + -27.9603*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.d_He3__p_He4 = rate

@numba.njit()
def t_He3__d_He4(rate_eval, tf):
    # He3 + t --> d + He4
    rate = 0.0

    # cf88n
    rate += np.exp(  22.4207 + -7.733*tf.T913i + -0.133473*tf.T913
                  + -0.294412*tf.T9 + 0.0310968*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.t_He3__d_He4 = rate

@numba.njit()
def He4_He3__p_Li6(rate_eval, tf):
    # He3 + He4 --> p + Li6
    rate = 0.0

    # pt05n
    rate += np.exp(  24.4064 + -46.6405*tf.T9i + -8.39481*tf.T913i + -0.165254*tf.T913
                  + -0.16936*tf.T9 + 0.0533676*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_He3__p_Li6 = rate

@numba.njit()
def n_He4__d_t(rate_eval, tf):
    # n + He4 --> d + t
    rate = 0.0

    # de04 
    rate += np.exp(  41.0525 + -204.112*tf.T9i + -4.5244*tf.T913i + -16.4028*tf.T913
                  + 1.73103*tf.T9 + -0.122966*tf.T953 + 2.31304*tf.lnT9)
    # de04 
    rate += np.exp(  26.8862 + -204.112*tf.T9i + -4.5244*tf.T913i + 0.350337*tf.T913
                  + 0.58747*tf.T9 + -8.84909*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.n_He4__d_t = rate

@numba.njit()
def p_He4__d_He3(rate_eval, tf):
    # p + He4 --> d + He3
    rate = 0.0

    # de04 
    rate += np.exp(  26.3907 + -212.977*tf.T9i + -7.182*tf.T913i + 0.473288*tf.T913
                  + 1.46847*tf.T9 + -27.9603*tf.T953 + -0.666667*tf.lnT9)
    # de04 
    rate += np.exp(  43.0037 + -212.977*tf.T9i + -7.182*tf.T913i + -17.1349*tf.T913
                  + 1.36908*tf.T9 + -0.0814423*tf.T953 + 3.35395*tf.lnT9)

    rate_eval.p_He4__d_He3 = rate

@numba.njit()
def d_He4__t_He3(rate_eval, tf):
    # d + He4 --> t + He3
    rate = 0.0

    # cf88n
    rate += np.exp(  22.8851 + -166.176*tf.T9i + -7.733*tf.T913i + -0.133473*tf.T913
                  + -0.294412*tf.T9 + 0.0310968*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.d_He4__t_He3 = rate

@numba.njit()
def He4_He4__n_Be7(rate_eval, tf):
    # He4 + He4 --> n + Be7
    rate = 0.0

    #  wagn
    rate += np.exp(  19.694 + -220.375*tf.T9i + -0.00210045*tf.T913
                  + 0.000176541*tf.T9 + -1.36797e-05*tf.T953 + 1.00083*tf.lnT9)

    rate_eval.He4_He4__n_Be7 = rate

@numba.njit()
def He4_He4__p_Li7(rate_eval, tf):
    # He4 + He4 --> p + Li7
    rate = 0.0

    # de04r
    rate += np.exp(  15.7864 + -205.79*tf.T9i
                  + -1.5*tf.lnT9)
    # de04 
    rate += np.exp(  13.4902 + -201.312*tf.T9i + -8.4727*tf.T913i + 0.417943*tf.T913
                  + 5.34565*tf.T9 + -4.8684*tf.T953 + -0.666667*tf.lnT9)
    # de04r
    rate += np.exp(  23.4325 + -227.465*tf.T9i
                  + -1.5*tf.lnT9)
    # de04 
    rate += np.exp(  21.9764 + -201.312*tf.T9i + -8.4727*tf.T913i + 0.297934*tf.T913
                  + 0.0582335*tf.T9 + -0.00413383*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_He4__p_Li7 = rate

@numba.njit()
def n_Li6__He4_t(rate_eval, tf):
    # Li6 + n --> He4 + t
    rate = 0.0

    # cf88r
    rate += np.exp(  21.665 + -2.39128*tf.T9i
                  + -1.5*tf.lnT9)
    # cf88n
    rate += np.exp(  18.9496 + -0.001281*tf.T9i)

    rate_eval.n_Li6__He4_t = rate

@numba.njit()
def p_Li6__He4_He3(rate_eval, tf):
    # Li6 + p --> He4 + He3
    rate = 0.0

    # pt05n
    rate += np.exp(  24.3475 + -8.39481*tf.T913i + -0.165254*tf.T913
                  + -0.16936*tf.T9 + 0.0533676*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Li6__He4_He3 = rate

@numba.njit()
def d_Li6__n_Be7(rate_eval, tf):
    # Li6 + d --> n + Be7
    rate = 0.0

    # mafon
    rate += np.exp(  28.0095 + -4.77456e-12*tf.T9i + -10.259*tf.T913i + -2.01559e-09*tf.T913
                  + 1.99542e-10*tf.T9 + -1.65595e-11*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.d_Li6__n_Be7 = rate

@numba.njit()
def d_Li6__p_Li7(rate_eval, tf):
    # Li6 + d --> p + Li7
    rate = 0.0

    # mafon
    rate += np.exp(  28.0231 + -10.135*tf.T913i
                  + -0.666667*tf.lnT9)

    rate_eval.d_Li6__p_Li7 = rate

@numba.njit()
def He4_Li6__p_Be9(rate_eval, tf):
    # Li6 + He4 --> p + Be9
    rate = 0.0

    # cf88n
    rate += np.exp(  25.5847 + -24.6712*tf.T9i + -10.359*tf.T913i + 0.102577*tf.T913
                  + 4.43544*tf.T9 + -5.97105*tf.T953 + -0.666667*tf.lnT9)
    # cf88r
    rate += np.exp(  19.8324 + -29.8312*tf.T9i
                  + -0.75*tf.lnT9)
    # cf88r
    rate += np.exp(  19.4366 + -27.7172*tf.T9i
                  + -1.0*tf.lnT9)

    rate_eval.He4_Li6__p_Be9 = rate

@numba.njit()
def p_Li7__n_Be7(rate_eval, tf):
    # Li7 + p --> n + Be7
    rate = 0.0

    # db18 
    rate += np.exp(  21.7899 + -19.0779*tf.T9i + -0.30254*tf.T913i + -0.3602*tf.T913
                  + 0.17472*tf.T9 + -0.0223*tf.T953 + -0.4581*tf.lnT9)

    rate_eval.p_Li7__n_Be7 = rate

@numba.njit()
def p_Li7__d_Li6(rate_eval, tf):
    # Li7 + p --> d + Li6
    rate = 0.0

    # mafon
    rate += np.exp(  28.9494 + -58.3239*tf.T9i + -10.135*tf.T913i
                  + -0.666667*tf.lnT9)

    rate_eval.p_Li7__d_Li6 = rate

@numba.njit()
def p_Li7__He4_He4(rate_eval, tf):
    # Li7 + p --> He4 + He4
    rate = 0.0

    # de04r
    rate += np.exp(  21.8999 + -26.1527*tf.T9i
                  + -1.5*tf.lnT9)
    # de04 
    rate += np.exp(  20.4438 + -8.4727*tf.T913i + 0.297934*tf.T913
                  + 0.0582335*tf.T9 + -0.00413383*tf.T953 + -0.666667*tf.lnT9)
    # de04r
    rate += np.exp(  14.2538 + -4.478*tf.T9i
                  + -1.5*tf.lnT9)
    # de04 
    rate += np.exp(  11.9576 + -8.4727*tf.T913i + 0.417943*tf.T913
                  + 5.34565*tf.T9 + -4.8684*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Li7__He4_He4 = rate

@numba.njit()
def t_Li7__n_Be9(rate_eval, tf):
    # Li7 + t --> n + Be9
    rate = 0.0

    # bk91 
    rate += np.exp(  30.2619 + -11.333*tf.T913i + -0.170964*tf.T913
                  + -6.30572*tf.T9 + 1.2248*tf.T953 + -0.666667*tf.lnT9)
    # bk91 
    rate += np.exp(  34.634 + -11.333*tf.T913i + -7.3964*tf.T913
                  + 0.947759*tf.T9 + -0.0839008*tf.T953 + 0.333333*tf.lnT9)

    rate_eval.t_Li7__n_Be9 = rate

@numba.njit()
def n_Be7__p_Li7(rate_eval, tf):
    # Be7 + n --> p + Li7
    rate = 0.0

    # db18 
    rate += np.exp(  21.7899 + 0.000728098*tf.T9i + -0.30254*tf.T913i + -0.3602*tf.T913
                  + 0.17472*tf.T9 + -0.0223*tf.T953 + -0.4581*tf.lnT9)

    rate_eval.n_Be7__p_Li7 = rate

@numba.njit()
def n_Be7__d_Li6(rate_eval, tf):
    # Be7 + n --> d + Li6
    rate = 0.0

    # mafon
    rate += np.exp(  28.9358 + -39.2438*tf.T9i + -10.259*tf.T913i + -2.01559e-09*tf.T913
                  + 1.99542e-10*tf.T9 + -1.65595e-11*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.n_Be7__d_Li6 = rate

@numba.njit()
def n_Be7__He4_He4(rate_eval, tf):
    # Be7 + n --> He4 + He4
    rate = 0.0

    #  wagn
    rate += np.exp(  18.1614 + -0.00210045*tf.T913
                  + 0.000176541*tf.T9 + -1.36797e-05*tf.T953 + 1.00083*tf.lnT9)

    rate_eval.n_Be7__He4_He4 = rate

@numba.njit()
def n_Be9__t_Li7(rate_eval, tf):
    # Be9 + n --> t + Li7
    rate = 0.0

    # bk91 
    rate += np.exp(  35.9049 + -121.123*tf.T9i + -11.333*tf.T913i + -7.3964*tf.T913
                  + 0.947759*tf.T9 + -0.0839008*tf.T953 + 0.333333*tf.lnT9)
    # bk91 
    rate += np.exp(  31.5328 + -121.123*tf.T9i + -11.333*tf.T913i + -0.170964*tf.T913
                  + -6.30572*tf.T9 + 1.2248*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.n_Be9__t_Li7 = rate

@numba.njit()
def p_Be9__He4_Li6(rate_eval, tf):
    # Be9 + p --> He4 + Li6
    rate = 0.0

    # cf88n
    rate += np.exp(  26.0751 + -10.359*tf.T913i + 0.102577*tf.T913
                  + 4.43544*tf.T9 + -5.97105*tf.T953 + -0.666667*tf.lnT9)
    # cf88r
    rate += np.exp(  20.3228 + -5.16*tf.T9i
                  + -0.75*tf.lnT9)
    # cf88r
    rate += np.exp(  19.927 + -3.046*tf.T9i
                  + -1.0*tf.lnT9)

    rate_eval.p_Be9__He4_Li6 = rate

@numba.njit()
def He4_Be9__n_C12(rate_eval, tf):
    # Be9 + He4 --> n + C12
    rate = 0.0

    # cf88r
    rate += np.exp(  11.744 + -4.179*tf.T9i
                  + -1.5*tf.lnT9)
    # cf88r
    rate += np.exp(  -1.48281 + -1.834*tf.T9i
                  + -1.5*tf.lnT9)
    # cf88r
    rate += np.exp(  -9.51959 + -1.184*tf.T9i
                  + -1.5*tf.lnT9)
    # cf88n
    rate += np.exp(  31.464 + -23.87*tf.T913i + 0.566698*tf.T913
                  + 44.0957*tf.T9 + -314.232*tf.T953 + -0.666667*tf.lnT9)
    # cf88r
    rate += np.exp(  19.2962 + -12.732*tf.T9i)

    rate_eval.He4_Be9__n_C12 = rate

@numba.njit()
def He4_B8__p_C11(rate_eval, tf):
    # B8 + He4 --> p + C11
    rate = 0.0

    # wiesr
    rate += np.exp(  207.835 + -9.09757*tf.T9i + 381.949*tf.T913i + -608.493*tf.T913
                  + 34.7229*tf.T9 + -2.02465*tf.T953 + 297.922*tf.lnT9)
    # wiesr
    rate += np.exp(  13.8524 + -1.0296*tf.T9i + 3.19293*tf.T913i + -37.5025*tf.T913
                  + 4.4781*tf.T9 + -0.401587*tf.T953 + 8.35133*tf.lnT9)

    rate_eval.He4_B8__p_C11 = rate

@numba.njit()
def p_C11__He4_B8(rate_eval, tf):
    # C11 + p --> He4 + B8
    rate = 0.0

    # wiesr
    rate += np.exp(  208.967 + -95.0429*tf.T9i + 381.949*tf.T913i + -608.493*tf.T913
                  + 34.7229*tf.T9 + -2.02465*tf.T953 + 297.922*tf.lnT9)
    # wiesr
    rate += np.exp(  14.9842 + -86.975*tf.T9i + 3.19293*tf.T913i + -37.5025*tf.T913
                  + 4.4781*tf.T9 + -0.401587*tf.T953 + 8.35133*tf.lnT9)

    rate_eval.p_C11__He4_B8 = rate

@numba.njit()
def He4_C11__n_O14(rate_eval, tf):
    # C11 + He4 --> n + O14
    rate = 0.0

    # rpsmr
    rate += np.exp(  15.1629 + -34.8289*tf.T9i + -1.31632*tf.T913i + 2.06431*tf.T913
                  + 0.0585225*tf.T9 + -0.00948426*tf.T953 + -1.11933*tf.lnT9)

    rate_eval.He4_C11__n_O14 = rate

@numba.njit()
def He4_C11__p_N14(rate_eval, tf):
    # C11 + He4 --> p + N14
    rate = 0.0

    # cf88n
    rate += np.exp(  36.613 + -31.883*tf.T913i + -0.361593*tf.T913
                  + -0.394216*tf.T9 + 0.0239162*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_C11__p_N14 = rate

@numba.njit()
def n_C12__He4_Be9(rate_eval, tf):
    # C12 + n --> He4 + Be9
    rate = 0.0

    # cf88r
    rate += np.exp(  0.858256 + -67.9913*tf.T9i
                  + -1.5*tf.lnT9)
    # cf88r
    rate += np.exp(  -7.17852 + -67.3413*tf.T9i
                  + -1.5*tf.lnT9)
    # cf88n
    rate += np.exp(  33.8051 + -66.1573*tf.T9i + -23.87*tf.T913i + 0.566698*tf.T913
                  + 44.0957*tf.T9 + -314.232*tf.T953 + -0.666667*tf.lnT9)
    # cf88r
    rate += np.exp(  21.6373 + -78.8893*tf.T9i)
    # cf88r
    rate += np.exp(  14.0851 + -70.3363*tf.T9i
                  + -1.5*tf.lnT9)

    rate_eval.n_C12__He4_Be9 = rate

@numba.njit()
def He4_C12__n_O15(rate_eval, tf):
    # C12 + He4 --> n + O15
    rate = 0.0

    # cf88n
    rate += np.exp(  17.0115 + -98.6615*tf.T9i + 0.124787*tf.T913
                  + 0.0588937*tf.T9 + -0.00679206*tf.T953)

    rate_eval.He4_C12__n_O15 = rate

@numba.njit()
def He4_C12__p_N15(rate_eval, tf):
    # C12 + He4 --> p + N15
    rate = 0.0

    # nacrn
    rate += np.exp(  27.118 + -57.6279*tf.T9i + -15.253*tf.T913i + 1.59318*tf.T913
                  + 2.4479*tf.T9 + -2.19708*tf.T953 + -0.666667*tf.lnT9)
    # nacrr
    rate += np.exp(  -6.93365 + -58.7917*tf.T9i + 22.7105*tf.T913
                  + -2.90707*tf.T9 + 0.205754*tf.T953 + -1.5*tf.lnT9)
    # nacrr
    rate += np.exp(  20.5388 + -65.034*tf.T9i
                  + -1.5*tf.lnT9)
    # nacrr
    rate += np.exp(  -5.2319 + -59.6491*tf.T9i + 30.8497*tf.T913
                  + -8.50433*tf.T9 + -1.54426*tf.T953 + -1.5*tf.lnT9)

    rate_eval.He4_C12__p_N15 = rate

@numba.njit()
def C12_C12__n_Mg23(rate_eval, tf):
    # C12 + C12 --> n + Mg23
    rate = 0.0

    # cf88r
    rate += np.exp(  -12.8056 + -30.1498*tf.T9i + 11.4826*tf.T913
                  + 1.82849*tf.T9 + -0.34844*tf.T953)

    rate_eval.C12_C12__n_Mg23 = rate

@numba.njit()
def C12_C12__p_Na23(rate_eval, tf):
    # C12 + C12 --> p + Na23
    rate = 0.0

    # cf88r
    rate += np.exp(  60.9649 + -84.165*tf.T913i + -1.4191*tf.T913
                  + -0.114619*tf.T9 + -0.070307*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.C12_C12__p_Na23 = rate

@numba.njit()
def C12_C12__He4_Ne20(rate_eval, tf):
    # C12 + C12 --> He4 + Ne20
    rate = 0.0

    # cf88r
    rate += np.exp(  61.2863 + -84.165*tf.T913i + -1.56627*tf.T913
                  + -0.0736084*tf.T9 + -0.072797*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.C12_C12__He4_Ne20 = rate

@numba.njit()
def p_C13__n_N13(rate_eval, tf):
    # C13 + p --> n + N13
    rate = 0.0

    # nacrn
    rate += np.exp(  17.7625 + -34.8483*tf.T9i + 1.26126*tf.T913
                  + -0.204952*tf.T9 + 0.0310523*tf.T953)

    rate_eval.p_C13__n_N13 = rate

@numba.njit()
def d_C13__n_N14(rate_eval, tf):
    # C13 + d --> n + N14
    rate = 0.0

    # bb92n
    rate += np.exp(  27.1993 + -0.00261944*tf.T9i + -16.8935*tf.T913i + 4.06445*tf.T913
                  + -1.1715*tf.T9 + 0.118556*tf.T953 + -1.13937*tf.lnT9)

    rate_eval.d_C13__n_N14 = rate

@numba.njit()
def He4_C13__n_O16(rate_eval, tf):
    # C13 + He4 --> n + O16
    rate = 0.0

    # gl12 
    rate += np.exp(  62.5775 + -0.0277331*tf.T9i + -32.3917*tf.T913i + -48.934*tf.T913
                  + 44.1843*tf.T9 + -20.8743*tf.T953 + 2.02494*tf.lnT9)
    # gl12 
    rate += np.exp(  79.3008 + -0.30489*tf.T9i + 7.43132*tf.T913i + -84.8689*tf.T913
                  + 3.65083*tf.T9 + -0.148015*tf.T953 + 37.6008*tf.lnT9)

    rate_eval.He4_C13__n_O16 = rate

@numba.njit()
def p_C14__n_N14(rate_eval, tf):
    # C14 + p --> n + N14
    rate = 0.0

    # cf88 
    rate += np.exp(  5.23589 + -7.26442*tf.T9i + 12.3428*tf.T913
                  + -2.70025*tf.T9 + 0.236625*tf.T953 + 1.0*tf.lnT9)
    # cf88n
    rate += np.exp(  14.7608 + -7.26442*tf.T9i + -4.33989*tf.T913
                  + 11.4311*tf.T9 + -11.7764*tf.T953)

    rate_eval.p_C14__n_N14 = rate

@numba.njit()
def d_C14__n_N15(rate_eval, tf):
    # C14 + d --> n + N15
    rate = 0.0

    # bk92 
    rate += np.exp(  33.5637 + -16.939*tf.T913i + -4.14392*tf.T913
                  + 0.438623*tf.T9 + -0.0354193*tf.T953 + 0.333333*tf.lnT9)
    # bk92 
    rate += np.exp(  30.6841 + -16.939*tf.T913i + -0.582342*tf.T913
                  + -8.17066*tf.T9 + 1.70865*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.d_C14__n_N15 = rate

@numba.njit()
def He4_C14__n_O17(rate_eval, tf):
    # C14 + He4 --> n + O17
    rate = 0.0

    # kg91r
    rate += np.exp(  18.7889 + -23.8496*tf.T9i
                  + -1.5*tf.lnT9)
    # kg91r
    rate += np.exp(  14.4267 + -23.0516*tf.T9i
                  + -1.5*tf.lnT9)
    # kg91n
    rate += np.exp(  11.0417 + -21.0906*tf.T9i)

    rate_eval.He4_C14__n_O17 = rate

@numba.njit()
def n_N13__p_C13(rate_eval, tf):
    # N13 + n --> p + C13
    rate = 0.0

    # nacrn
    rate += np.exp(  17.7625 + 1.26126*tf.T913
                  + -0.204952*tf.T9 + 0.0310523*tf.T953)

    rate_eval.n_N13__p_C13 = rate

@numba.njit()
def He4_N13__n_F16(rate_eval, tf):
    # N13 + He4 --> n + F16
    rate = 0.0

    # rpsmr
    rate += np.exp(  14.9908 + -127.425*tf.T9i + -1.70176*tf.T913i + 3.84453*tf.T913
                  + -0.121351*tf.T9 + 0.00264846*tf.T953 + -1.63623*tf.lnT9)

    rate_eval.He4_N13__n_F16 = rate

@numba.njit()
def He4_N13__p_O16(rate_eval, tf):
    # N13 + He4 --> p + O16
    rate = 0.0

    # cf88n
    rate += np.exp(  40.4644 + -35.829*tf.T913i + -0.530275*tf.T913
                  + -0.982462*tf.T9 + 0.0808059*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_N13__p_O16 = rate

@numba.njit()
def n_N14__p_C14(rate_eval, tf):
    # N14 + n --> p + C14
    rate = 0.0

    # cf88 
    rate += np.exp(  4.13728 + 12.3428*tf.T913
                  + -2.70025*tf.T9 + 0.236625*tf.T953 + 1.0*tf.lnT9)
    # cf88n
    rate += np.exp(  13.6622 + -4.33989*tf.T913
                  + 11.4311*tf.T9 + -11.7764*tf.T953)

    rate_eval.n_N14__p_C14 = rate

@numba.njit()
def n_N14__d_C13(rate_eval, tf):
    # N14 + n --> d + C13
    rate = 0.0

    # bb92n
    rate += np.exp(  28.1279 + -61.8182*tf.T9i + -16.8935*tf.T913i + 4.06445*tf.T913
                  + -1.1715*tf.T9 + 0.118556*tf.T953 + -1.13937*tf.lnT9)

    rate_eval.n_N14__d_C13 = rate

@numba.njit()
def p_N14__n_O14(rate_eval, tf):
    # N14 + p --> n + O14
    rate = 0.0

    # nacr 
    rate += np.exp(  11.3432 + -68.7567*tf.T9i + 5.48024*tf.T913
                  + -0.764072*tf.T9 + 0.0587804*tf.T953)

    rate_eval.p_N14__n_O14 = rate

@numba.njit()
def p_N14__He4_C11(rate_eval, tf):
    # N14 + p --> He4 + C11
    rate = 0.0

    # cf88n
    rate += np.exp(  37.9252 + -33.92*tf.T9i + -31.883*tf.T913i + -0.361593*tf.T913
                  + -0.394216*tf.T9 + 0.0239162*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_N14__He4_C11 = rate

@numba.njit()
def He4_N14__n_F17(rate_eval, tf):
    # N14 + He4 --> n + F17
    rate = 0.0

    # nacr 
    rate += np.exp(  19.2094 + -54.9473*tf.T9i + -0.712285*tf.T913
                  + 0.240317*tf.T9 + -0.0201705*tf.T953)

    rate_eval.He4_N14__n_F17 = rate

@numba.njit()
def He4_N14__p_O17(rate_eval, tf):
    # N14 + He4 --> p + O17
    rate = 0.0

    # il10r
    rate += np.exp(  5.13169 + -15.9452*tf.T9i
                  + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  -7.60954 + -14.5839*tf.T9i
                  + -1.5*tf.lnT9)
    # il10n
    rate += np.exp(  19.1771 + -13.8305*tf.T9i + -16.9078*tf.T913i
                  + -2.0*tf.T953 + -0.666667*tf.lnT9)
    # il10r
    rate += np.exp(  9.77209 + -18.7891*tf.T9i + 5.10182*tf.T913
                  + 0.379373*tf.T9 + -0.0672515*tf.T953 + -1.5*tf.lnT9)

    rate_eval.He4_N14__p_O17 = rate

@numba.njit()
def n_N15__d_C14(rate_eval, tf):
    # N15 + n --> d + C14
    rate = 0.0

    # bk92 
    rate += np.exp(  34.2122 + -92.6344*tf.T9i + -16.939*tf.T913i + -4.14392*tf.T913
                  + 0.438623*tf.T9 + -0.0354193*tf.T953 + 0.333333*tf.lnT9)
    # bk92 
    rate += np.exp(  31.3326 + -92.6344*tf.T9i + -16.939*tf.T913i + -0.582342*tf.T913
                  + -8.17066*tf.T9 + 1.70865*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.n_N15__d_C14 = rate

@numba.njit()
def p_N15__n_O15(rate_eval, tf):
    # N15 + p --> n + O15
    rate = 0.0

    # nacrn
    rate += np.exp(  18.3942 + -41.0335*tf.T9i + 0.331392*tf.T913
                  + 0.0171473*tf.T9)

    rate_eval.p_N15__n_O15 = rate

@numba.njit()
def p_N15__He4_C12(rate_eval, tf):
    # N15 + p --> He4 + C12
    rate = 0.0

    # nacrn
    rate += np.exp(  27.4764 + -15.253*tf.T913i + 1.59318*tf.T913
                  + 2.4479*tf.T9 + -2.19708*tf.T953 + -0.666667*tf.lnT9)
    # nacrr
    rate += np.exp(  -6.57522 + -1.1638*tf.T9i + 22.7105*tf.T913
                  + -2.90707*tf.T9 + 0.205754*tf.T953 + -1.5*tf.lnT9)
    # nacrr
    rate += np.exp(  20.8972 + -7.406*tf.T9i
                  + -1.5*tf.lnT9)
    # nacrr
    rate += np.exp(  -4.87347 + -2.02117*tf.T9i + 30.8497*tf.T913
                  + -8.50433*tf.T9 + -1.54426*tf.T953 + -1.5*tf.lnT9)

    rate_eval.p_N15__He4_C12 = rate

@numba.njit()
def He4_N15__n_F18(rate_eval, tf):
    # N15 + He4 --> n + F18
    rate = 0.0

    # cf88n
    rate += np.exp(  18.0938 + -74.4713*tf.T9i + 1.74308*tf.T913
                  + -1.15123*tf.T9 + 0.135196*tf.T953)

    rate_eval.He4_N15__n_F18 = rate

@numba.njit()
def He4_N15__p_O18(rate_eval, tf):
    # N15 + He4 --> p + O18
    rate = 0.0

    # il10r
    rate += np.exp(  8.46654 + -47.8616*tf.T9i
                  + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  -29.7104 + -46.4444*tf.T9i
                  + -1.5*tf.lnT9)
    # il10n
    rate += np.exp(  25.1611 + -46.1986*tf.T9i + -16.6979*tf.T913i
                  + -3.0*tf.T953 + -0.666667*tf.lnT9)
    # il10r
    rate += np.exp(  7.13756 + -51.5219*tf.T9i + 11.6568*tf.T913
                  + -2.16303*tf.T9 + 0.209965*tf.T953 + -1.5*tf.lnT9)

    rate_eval.He4_N15__p_O18 = rate

@numba.njit()
def n_O14__p_N14(rate_eval, tf):
    # O14 + n --> p + N14
    rate = 0.0

    # nacr 
    rate += np.exp(  12.4418 + 5.48024*tf.T913
                  + -0.764072*tf.T9 + 0.0587804*tf.T953)

    rate_eval.n_O14__p_N14 = rate

@numba.njit()
def n_O14__He4_C11(rate_eval, tf):
    # O14 + n --> He4 + C11
    rate = 0.0

    # rpsmr
    rate += np.exp(  17.5737 + 0.00781138*tf.T9i + -1.31632*tf.T913i + 2.06431*tf.T913
                  + 0.0585225*tf.T9 + -0.00948426*tf.T953 + -1.11933*tf.lnT9)

    rate_eval.n_O14__He4_C11 = rate

@numba.njit()
def He4_O14__p_F17(rate_eval, tf):
    # O14 + He4 --> p + F17
    rate = 0.0

    # Ha96r
    rate += np.exp(  12.1289 + -12.0223*tf.T9i
                  + -1.5*tf.lnT9)
    # Ha96r
    rate += np.exp(  18.6518 + -26.0*tf.T9i
                  + -1.5*tf.lnT9)
    # Ha96n
    rate += np.exp(  40.8358 + -39.388*tf.T913i + -17.4673*tf.T913
                  + 35.3029*tf.T9 + -24.8162*tf.T953 + -0.666667*tf.lnT9)
    # Ha96r
    rate += np.exp(  16.3087 + -22.51*tf.T9i
                  + -1.5*tf.lnT9)
    # Ha96r
    rate += np.exp(  11.1184 + -13.6*tf.T9i
                  + -1.5*tf.lnT9)
    # Ha96r
    rate += np.exp(  -106.091 + -0.453036*tf.T9i
                  + -1.5*tf.lnT9)

    rate_eval.He4_O14__p_F17 = rate

@numba.njit()
def n_O15__p_N15(rate_eval, tf):
    # O15 + n --> p + N15
    rate = 0.0

    # nacrn
    rate += np.exp(  18.3942 + 0.331392*tf.T913
                  + 0.0171473*tf.T9)

    rate_eval.n_O15__p_N15 = rate

@numba.njit()
def n_O15__He4_C12(rate_eval, tf):
    # O15 + n --> He4 + C12
    rate = 0.0

    # cf88n
    rate += np.exp(  17.3699 + 0.124787*tf.T913
                  + 0.0588937*tf.T9 + -0.00679206*tf.T953)

    rate_eval.n_O15__He4_C12 = rate

@numba.njit()
def p_O15__n_F15(rate_eval, tf):
    # O15 + p --> n + F15
    rate = 0.0

    # rpsmr
    rate += np.exp(  18.6182 + -170.574*tf.T9i + -3.25485*tf.T913i + 4.5547*tf.T913
                  + -0.151514*tf.T9 + 0.0064577*tf.T953 + -2.38438*tf.lnT9)

    rate_eval.p_O15__n_F15 = rate

@numba.njit()
def He4_O15__n_Ne18(rate_eval, tf):
    # O15 + He4 --> n + Ne18
    rate = 0.0

    # ths8r
    rate += np.exp(  16.2469 + -94.0888*tf.T9i + 0.377325*tf.T913
                  + 0.0968945*tf.T9 + -0.007769*tf.T953)

    rate_eval.He4_O15__n_Ne18 = rate

@numba.njit()
def He4_O15__p_F18(rate_eval, tf):
    # O15 + He4 --> p + F18
    rate = 0.0

    # il10n
    rate += np.exp(  61.2985 + -33.4459*tf.T9i + -21.4023*tf.T913i + -80.8891*tf.T913
                  + 134.6*tf.T9 + -126.504*tf.T953 + -0.666667*tf.lnT9)
    # il10r
    rate += np.exp(  1.04969 + -36.4627*tf.T9i + 13.3223*tf.T913
                  + -1.36696*tf.T9 + 0.0757363*tf.T953 + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  -32.4461 + -33.8223*tf.T9i + 61.738*tf.T913
                  + -108.29*tf.T9 + -34.2365*tf.T953 + -1.5*tf.lnT9)

    rate_eval.He4_O15__p_F18 = rate

@numba.njit()
def n_O16__He4_C13(rate_eval, tf):
    # O16 + n --> He4 + C13
    rate = 0.0

    # gl12 
    rate += np.exp(  64.3455 + -25.7388*tf.T9i + -32.3917*tf.T913i + -48.934*tf.T913
                  + 44.1843*tf.T9 + -20.8743*tf.T953 + 2.02494*tf.lnT9)
    # gl12 
    rate += np.exp(  81.0688 + -26.0159*tf.T9i + 7.43132*tf.T913i + -84.8689*tf.T913
                  + 3.65083*tf.T9 + -0.148015*tf.T953 + 37.6008*tf.lnT9)

    rate_eval.n_O16__He4_C13 = rate

@numba.njit()
def p_O16__n_F16(rate_eval, tf):
    # O16 + p --> n + F16
    rate = 0.0

    # rpsmr
    rate += np.exp(  18.9074 + -187.94*tf.T9i + -2.69385*tf.T913i + 3.18271*tf.T913
                  + -0.0491204*tf.T9 + 0.000327601*tf.T953 + -1.83412*tf.lnT9)

    rate_eval.p_O16__n_F16 = rate

@numba.njit()
def p_O16__He4_N13(rate_eval, tf):
    # O16 + p --> He4 + N13
    rate = 0.0

    # cf88n
    rate += np.exp(  42.2324 + -60.5523*tf.T9i + -35.829*tf.T913i + -0.530275*tf.T913
                  + -0.982462*tf.T9 + 0.0808059*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_O16__He4_N13 = rate

@numba.njit()
def He4_O16__n_Ne19(rate_eval, tf):
    # O16 + He4 --> n + Ne19
    rate = 0.0

    # ths8r
    rate += np.exp(  17.2055 + -140.818*tf.T9i + 1.70736*tf.T913
                  + -0.132579*tf.T9 + 0.00454218*tf.T953)

    rate_eval.He4_O16__n_Ne19 = rate

@numba.njit()
def He4_O16__p_F19(rate_eval, tf):
    # O16 + He4 --> p + F19
    rate = 0.0

    # nacr 
    rate += np.exp(  7.80363 + -96.6272*tf.T9i
                  + -1.5*tf.lnT9)
    # nacr 
    rate += np.exp(  -53.1397 + -94.2866*tf.T9i
                  + -1.5*tf.lnT9)
    # nacr 
    rate += np.exp(  25.8562 + -94.1589*tf.T9i + -18.116*tf.T913i
                  + 1.86674*tf.T9 + -7.5666*tf.T953 + -0.666667*tf.lnT9)
    # nacrr
    rate += np.exp(  13.9232 + -97.4449*tf.T9i
                  + -0.21103*tf.T9 + 2.87702*tf.lnT9)
    # nacr 
    rate += np.exp(  14.7601 + -97.9108*tf.T9i
                  + -1.5*tf.lnT9)

    rate_eval.He4_O16__p_F19 = rate

@numba.njit()
def C12_O16__p_Al27(rate_eval, tf):
    # O16 + C12 --> p + Al27
    rate = 0.0

    # cf88r
    rate += np.exp(  68.5253 + 0.205134*tf.T9i + -119.242*tf.T913i + 13.3667*tf.T913
                  + 0.295425*tf.T9 + -0.267288*tf.T953 + -9.91729*tf.lnT9)

    rate_eval.C12_O16__p_Al27 = rate

@numba.njit()
def C12_O16__He4_Mg24(rate_eval, tf):
    # O16 + C12 --> He4 + Mg24
    rate = 0.0

    # cf88r
    rate += np.exp(  48.5341 + 0.37204*tf.T9i + -133.413*tf.T913i + 50.1572*tf.T913
                  + -3.15987*tf.T9 + 0.0178251*tf.T953 + -23.7027*tf.lnT9)

    rate_eval.C12_O16__He4_Mg24 = rate

@numba.njit()
def O16_O16__n_S31(rate_eval, tf):
    # O16 + O16 --> n + S31
    rate = 0.0

    # cf88r
    rate += np.exp(  77.5491 + -0.373641*tf.T9i + -120.83*tf.T913i + -7.72334*tf.T913
                  + -2.27939*tf.T9 + 0.167655*tf.T953 + 7.62001*tf.lnT9)

    rate_eval.O16_O16__n_S31 = rate

@numba.njit()
def O16_O16__p_P31(rate_eval, tf):
    # O16 + O16 --> p + P31
    rate = 0.0

    # cf88r
    rate += np.exp(  85.2628 + 0.223453*tf.T9i + -145.844*tf.T913i + 8.72612*tf.T913
                  + -0.554035*tf.T9 + -0.137562*tf.T953 + -6.88807*tf.lnT9)

    rate_eval.O16_O16__p_P31 = rate

@numba.njit()
def O16_O16__He4_Si28(rate_eval, tf):
    # O16 + O16 --> He4 + Si28
    rate = 0.0

    # cf88r
    rate += np.exp(  97.2435 + -0.268514*tf.T9i + -119.324*tf.T913i + -32.2497*tf.T913
                  + 1.46214*tf.T9 + -0.200893*tf.T953 + 13.2148*tf.lnT9)

    rate_eval.O16_O16__He4_Si28 = rate

@numba.njit()
def n_O17__He4_C14(rate_eval, tf):
    # O17 + n --> He4 + C14
    rate = 0.0

    # kg91r
    rate += np.exp(  13.73 + -1.961*tf.T9i
                  + -1.5*tf.lnT9)
    # kg91n
    rate += np.exp(  10.345)
    # kg91r
    rate += np.exp(  18.0922 + -2.759*tf.T9i
                  + -1.5*tf.lnT9)

    rate_eval.n_O17__He4_C14 = rate

@numba.njit()
def p_O17__n_F17(rate_eval, tf):
    # O17 + p --> n + F17
    rate = 0.0

    #  wagn
    rate += np.exp(  19.0085 + -41.1187*tf.T9i)

    rate_eval.p_O17__n_F17 = rate

@numba.njit()
def p_O17__He4_N14(rate_eval, tf):
    # O17 + p --> He4 + N14
    rate = 0.0

    # il10r
    rate += np.exp(  -7.20763 + -0.753395*tf.T9i
                  + -1.5*tf.lnT9)
    # il10n
    rate += np.exp(  19.579 + -16.9078*tf.T913i
                  + -2.0*tf.T953 + -0.666667*tf.lnT9)
    # il10r
    rate += np.exp(  10.174 + -4.95865*tf.T9i + 5.10182*tf.T913
                  + 0.379373*tf.T9 + -0.0672515*tf.T953 + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  5.5336 + -2.11477*tf.T9i
                  + -1.5*tf.lnT9)

    rate_eval.p_O17__He4_N14 = rate

@numba.njit()
def He4_O17__n_Ne20(rate_eval, tf):
    # O17 + He4 --> n + Ne20
    rate = 0.0

    # nacrr
    rate += np.exp(  1.80342 + -13.8*tf.T9i + 12.6501*tf.T913
                  + -1.10938*tf.T9 + 0.0696232*tf.T953 + -1.5*tf.lnT9)
    # nacrr
    rate += np.exp(  7.45588 + -8.55*tf.T9i
                  + -1.5*tf.lnT9)
    # nacrn
    rate += np.exp(  40.621 + -39.918*tf.T913i
                  + 0.227017*tf.T9 + -0.900234*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_O17__n_Ne20 = rate

@numba.njit()
def p_O18__n_F18(rate_eval, tf):
    # O18 + p --> n + F18
    rate = 0.0

    #  wagn
    rate += np.exp(  20.1071 + -28.2841*tf.T9i)

    rate_eval.p_O18__n_F18 = rate

@numba.njit()
def p_O18__He4_N15(rate_eval, tf):
    # O18 + p --> He4 + N15
    rate = 0.0

    # il10r
    rate += np.exp(  -27.9044 + -0.245884*tf.T9i
                  + -1.5*tf.lnT9)
    # il10n
    rate += np.exp(  26.9671 + -16.6979*tf.T913i
                  + -3.0*tf.T953 + -0.666667*tf.lnT9)
    # il10r
    rate += np.exp(  8.94352 + -5.32335*tf.T9i + 11.6568*tf.T913
                  + -2.16303*tf.T9 + 0.209965*tf.T953 + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  10.2725 + -1.663*tf.T9i
                  + -1.5*tf.lnT9)

    rate_eval.p_O18__He4_N15 = rate

@numba.njit()
def He4_O18__n_Ne21(rate_eval, tf):
    # O18 + He4 --> n + Ne21
    rate = 0.0

    # nacr 
    rate += np.exp(  2.51575 + -8.41969*tf.T9i + 6.51299*tf.T913
                  + -0.334277*tf.T9 + 2.0*tf.lnT9)
    # nacr 
    rate += np.exp(  1.3415 + -15.444*tf.T9i + 11.656*tf.T913
                  + -1.24869*tf.T9 + 1.0*tf.lnT9)
    # nacrn
    rate += np.exp(  2.82796 + -8.085*tf.T9i)

    rate_eval.He4_O18__n_Ne21 = rate

@numba.njit()
def n_F15__p_O15(rate_eval, tf):
    # F15 + n --> p + O15
    rate = 0.0

    # rpsmr
    rate += np.exp(  18.6182 + 0.0583936*tf.T9i + -3.25485*tf.T913i + 4.5547*tf.T913
                  + -0.151514*tf.T9 + 0.0064577*tf.T953 + -2.38438*tf.lnT9)

    rate_eval.n_F15__p_O15 = rate

@numba.njit()
def He4_F15__p_Ne18(rate_eval, tf):
    # F15 + He4 --> p + Ne18
    rate = 0.0

    # rpsmr
    rate += np.exp(  46.8138 + -1.58356*tf.T9i + 10.8535*tf.T913i + -60.1083*tf.T913
                  + 2.19998*tf.T9 + -0.0925197*tf.T953 + 33.282*tf.lnT9)

    rate_eval.He4_F15__p_Ne18 = rate

@numba.njit()
def n_F16__p_O16(rate_eval, tf):
    # F16 + n --> p + O16
    rate = 0.0

    # rpsmr
    rate += np.exp(  18.9074 + 0.0532742*tf.T9i + -2.69385*tf.T913i + 3.18271*tf.T913
                  + -0.0491204*tf.T9 + 0.000327601*tf.T953 + -1.83412*tf.lnT9)

    rate_eval.n_F16__p_O16 = rate

@numba.njit()
def n_F16__He4_N13(rate_eval, tf):
    # F16 + n --> He4 + N13
    rate = 0.0

    # rpsmr
    rate += np.exp(  16.7588 + 0.00426483*tf.T9i + -1.70176*tf.T913i + 3.84453*tf.T913
                  + -0.121351*tf.T9 + 0.00264846*tf.T953 + -1.63623*tf.lnT9)

    rate_eval.n_F16__He4_N13 = rate

@numba.njit()
def He4_F16__p_Ne19(rate_eval, tf):
    # F16 + He4 --> p + Ne19
    rate = 0.0

    # rpsmr
    rate += np.exp(  48.4419 + -1.67449*tf.T9i + 14.6635*tf.T913i + -65.6394*tf.T913
                  + 2.4479*tf.T9 + -0.102505*tf.T953 + 36.1187*tf.lnT9)

    rate_eval.He4_F16__p_Ne19 = rate

@numba.njit()
def n_F17__p_O17(rate_eval, tf):
    # F17 + n --> p + O17
    rate = 0.0

    #  wagn
    rate += np.exp(  19.0085)

    rate_eval.n_F17__p_O17 = rate

@numba.njit()
def n_F17__He4_N14(rate_eval, tf):
    # F17 + n --> He4 + N14
    rate = 0.0

    # nacr 
    rate += np.exp(  19.6113 + -0.712285*tf.T913
                  + 0.240317*tf.T9 + -0.0201705*tf.T953)

    rate_eval.n_F17__He4_N14 = rate

@numba.njit()
def p_F17__He4_O14(rate_eval, tf):
    # F17 + p --> He4 + O14
    rate = 0.0

    # Ha96r
    rate += np.exp(  10.4217 + -27.4326*tf.T9i
                  + -1.5*tf.lnT9)
    # Ha96r
    rate += np.exp(  -106.788 + -14.2856*tf.T9i
                  + -1.5*tf.lnT9)
    # Ha96r
    rate += np.exp(  11.4322 + -25.8549*tf.T9i
                  + -1.5*tf.lnT9)
    # Ha96r
    rate += np.exp(  17.9551 + -39.8326*tf.T9i
                  + -1.5*tf.lnT9)
    # Ha96n
    rate += np.exp(  40.1391 + -13.8326*tf.T9i + -39.388*tf.T913i + -17.4673*tf.T913
                  + 35.3029*tf.T9 + -24.8162*tf.T953 + -0.666667*tf.lnT9)
    # Ha96r
    rate += np.exp(  15.612 + -36.3426*tf.T9i
                  + -1.5*tf.lnT9)

    rate_eval.p_F17__He4_O14 = rate

@numba.njit()
def He4_F17__p_Ne20(rate_eval, tf):
    # F17 + He4 --> p + Ne20
    rate = 0.0

    # nacr 
    rate += np.exp(  38.6287 + -43.18*tf.T913i + 4.46827*tf.T913
                  + -1.63915*tf.T9 + 0.123483*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_F17__p_Ne20 = rate

@numba.njit()
def n_F18__p_O18(rate_eval, tf):
    # F18 + n --> p + O18
    rate = 0.0

    #  wagn
    rate += np.exp(  19.0085)

    rate_eval.n_F18__p_O18 = rate

@numba.njit()
def n_F18__He4_N15(rate_eval, tf):
    # F18 + n --> He4 + N15
    rate = 0.0

    # cf88n
    rate += np.exp(  18.8011 + 1.74308*tf.T913
                  + -1.15123*tf.T9 + 0.135196*tf.T953)

    rate_eval.n_F18__He4_N15 = rate

@numba.njit()
def p_F18__n_Ne18(rate_eval, tf):
    # F18 + p --> n + Ne18
    rate = 0.0

    # ths8r
    rate += np.exp(  17.4443 + -60.6429*tf.T9i + 0.590026*tf.T913
                  + 0.0742242*tf.T9 + -0.0116856*tf.T953)

    rate_eval.p_F18__n_Ne18 = rate

@numba.njit()
def p_F18__He4_O15(rate_eval, tf):
    # F18 + p --> He4 + O15
    rate = 0.0

    # il10r
    rate += np.exp(  1.75704 + -3.01675*tf.T9i + 13.3223*tf.T913
                  + -1.36696*tf.T9 + 0.0757363*tf.T953 + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  -31.7388 + -0.376432*tf.T9i + 61.738*tf.T913
                  + -108.29*tf.T9 + -34.2365*tf.T953 + -1.5*tf.lnT9)
    # il10n
    rate += np.exp(  62.0058 + -21.4023*tf.T913i + -80.8891*tf.T913
                  + 134.6*tf.T9 + -126.504*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_F18__He4_O15 = rate

@numba.njit()
def He4_F18__n_Na21(rate_eval, tf):
    # F18 + He4 --> n + Na21
    rate = 0.0

    # ths8r
    rate += np.exp(  14.8808 + -30.0387*tf.T9i + 1.06229*tf.T913
                  + 0.212448*tf.T9 + -0.020685*tf.T953)

    rate_eval.He4_F18__n_Na21 = rate

@numba.njit()
def He4_F18__p_Ne21(rate_eval, tf):
    # F18 + He4 --> p + Ne21
    rate = 0.0

    # rpsmr
    rate += np.exp(  49.7863 + -1.84559*tf.T9i + 21.4461*tf.T913i + -73.252*tf.T913
                  + 2.42329*tf.T9 + -0.077278*tf.T953 + 40.7604*tf.lnT9)

    rate_eval.He4_F18__p_Ne21 = rate

@numba.njit()
def p_F19__n_Ne19(rate_eval, tf):
    # F19 + p --> n + Ne19
    rate = 0.0

    # nacr 
    rate += np.exp(  15.6732 + -46.6617*tf.T9i + 4.34441*tf.T913
                  + -1.71268*tf.T9)
    # nacr 
    rate += np.exp(  20.353 + -46.6617*tf.T9i + -1.58171*tf.T913
                  + 0.0324472*tf.T9 + 1.0*tf.lnT9)

    rate_eval.p_F19__n_Ne19 = rate

@numba.njit()
def p_F19__He4_O16(rate_eval, tf):
    # F19 + p --> He4 + O16
    rate = 0.0

    # nacr 
    rate += np.exp(  -52.7043 + -0.12765*tf.T9i
                  + -1.5*tf.lnT9)
    # nacr 
    rate += np.exp(  26.2916 + -18.116*tf.T913i
                  + 1.86674*tf.T9 + -7.5666*tf.T953 + -0.666667*tf.lnT9)
    # nacrr
    rate += np.exp(  14.3586 + -3.286*tf.T9i
                  + -0.21103*tf.T9 + 2.87702*tf.lnT9)
    # nacr 
    rate += np.exp(  15.1955 + -3.75185*tf.T9i
                  + -1.5*tf.lnT9)
    # nacr 
    rate += np.exp(  8.239 + -2.46828*tf.T9i
                  + -1.5*tf.lnT9)

    rate_eval.p_F19__He4_O16 = rate

@numba.njit()
def He4_F19__n_Na22(rate_eval, tf):
    # F19 + He4 --> n + Na22
    rate = 0.0

    # cf88 
    rate += np.exp(  6.67675 + -22.6172*tf.T9i + 6.22036*tf.T913
                  + -0.0209891*tf.T9 + -0.0689843*tf.T953 + 1.0*tf.lnT9)
    # cf88n
    rate += np.exp(  15.0925 + -22.6172*tf.T9i + -0.0532251*tf.T913
                  + 1.00044*tf.T9 + -0.13238*tf.T953)

    rate_eval.He4_F19__n_Na22 = rate

@numba.njit()
def n_Ne18__p_F18(rate_eval, tf):
    # Ne18 + n --> p + F18
    rate = 0.0

    # ths8r
    rate += np.exp(  18.5429 + 0.590026*tf.T913
                  + 0.0742242*tf.T9 + -0.0116856*tf.T953)

    rate_eval.n_Ne18__p_F18 = rate

@numba.njit()
def n_Ne18__He4_O15(rate_eval, tf):
    # Ne18 + n --> He4 + O15
    rate = 0.0

    # ths8r
    rate += np.exp(  18.0529 + 0.377325*tf.T913
                  + 0.0968945*tf.T9 + -0.007769*tf.T953)

    rate_eval.n_Ne18__He4_O15 = rate

@numba.njit()
def p_Ne18__He4_F15(rate_eval, tf):
    # Ne18 + p --> He4 + F15
    rate = 0.0

    # rpsmr
    rate += np.exp(  48.6197 + -78.2429*tf.T9i + 10.8535*tf.T913i + -60.1083*tf.T913
                  + 2.19998*tf.T9 + -0.0925197*tf.T953 + 33.282*tf.lnT9)

    rate_eval.p_Ne18__He4_F15 = rate

@numba.njit()
def He4_Ne18__p_Na21(rate_eval, tf):
    # Ne18 + He4 --> p + Na21
    rate = 0.0

    # mo14r
    rate += np.exp(  -23.5811 + -5.76874*tf.T9i + 8.5632e-08*tf.T913i + 24.8579*tf.T913
                  + 0.0823845*tf.T9 + -0.365374*tf.T953 + -2.21415e-06*tf.lnT9)
    # mo14r
    rate += np.exp(  -22.4389 + -10.2598*tf.T9i + 4.73034e-07*tf.T913i + 36.29*tf.T913
                  + -6.56565*tf.T9 + -2.96287e-06*tf.T953 + -9.00373e-07*tf.lnT9)
    # mo14r
    rate += np.exp(  -37.3979 + -3.8681*tf.T9i + 7.09521e-06*tf.T913i + 34.1789*tf.T913
                  + -1.72974*tf.T9 + -0.0395081*tf.T953 + -7.82759e-07*tf.lnT9)

    rate_eval.He4_Ne18__p_Na21 = rate

@numba.njit()
def n_Ne19__p_F19(rate_eval, tf):
    # Ne19 + n --> p + F19
    rate = 0.0

    # nacr 
    rate += np.exp(  15.6732 + 4.34441*tf.T913
                  + -1.71268*tf.T9)
    # nacr 
    rate += np.exp(  20.353 + -1.58171*tf.T913
                  + 0.0324472*tf.T9 + 1.0*tf.lnT9)

    rate_eval.n_Ne19__p_F19 = rate

@numba.njit()
def n_Ne19__He4_O16(rate_eval, tf):
    # Ne19 + n --> He4 + O16
    rate = 0.0

    # ths8r
    rate += np.exp(  17.6409 + 1.70736*tf.T913
                  + -0.132579*tf.T9 + 0.00454218*tf.T953)

    rate_eval.n_Ne19__He4_O16 = rate

@numba.njit()
def p_Ne19__He4_F16(rate_eval, tf):
    # Ne19 + p --> He4 + F16
    rate = 0.0

    # rpsmr
    rate += np.exp(  48.8773 + -48.8468*tf.T9i + 14.6635*tf.T913i + -65.6394*tf.T913
                  + 2.4479*tf.T9 + -0.102505*tf.T953 + 36.1187*tf.lnT9)

    rate_eval.p_Ne19__He4_F16 = rate

@numba.njit()
def He4_Ne19__n_Mg22(rate_eval, tf):
    # Ne19 + He4 --> n + Mg22
    rate = 0.0

    # ths8r
    rate += np.exp(  14.3989 + -40.5925*tf.T9i + 2.08592*tf.T913
                  + -0.188654*tf.T9 + 0.0101103*tf.T953)

    rate_eval.He4_Ne19__n_Mg22 = rate

@numba.njit()
def He4_Ne19__p_Na22(rate_eval, tf):
    # Ne19 + He4 --> p + Na22
    rate = 0.0

    # ths8r
    rate += np.exp(  43.1874 + -46.6346*tf.T913i + 0.866532*tf.T913
                  + -0.893541*tf.T9 + 0.0747971*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Ne19__p_Na22 = rate

@numba.njit()
def n_Ne20__He4_O17(rate_eval, tf):
    # Ne20 + n --> He4 + O17
    rate = 0.0

    # nacrr
    rate += np.exp(  4.7377 + -20.6002*tf.T9i + 12.6501*tf.T913
                  + -1.10938*tf.T9 + 0.0696232*tf.T953 + -1.5*tf.lnT9)
    # nacrr
    rate += np.exp(  10.3902 + -15.3502*tf.T9i
                  + -1.5*tf.lnT9)
    # nacrn
    rate += np.exp(  43.5553 + -6.80024*tf.T9i + -39.918*tf.T913i
                  + 0.227017*tf.T9 + -0.900234*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.n_Ne20__He4_O17 = rate

@numba.njit()
def p_Ne20__He4_F17(rate_eval, tf):
    # Ne20 + p --> He4 + F17
    rate = 0.0

    # nacr 
    rate += np.exp(  41.563 + -47.9266*tf.T9i + -43.18*tf.T913i + 4.46827*tf.T913
                  + -1.63915*tf.T9 + 0.123483*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Ne20__He4_F17 = rate

@numba.njit()
def He4_Ne20__n_Mg23(rate_eval, tf):
    # Ne20 + He4 --> n + Mg23
    rate = 0.0

    # ths8r
    rate += np.exp(  17.9544 + -83.7215*tf.T9i + 1.83199*tf.T913
                  + -0.290485*tf.T9 + 0.0242929*tf.T953)

    rate_eval.He4_Ne20__n_Mg23 = rate

@numba.njit()
def He4_Ne20__p_Na23(rate_eval, tf):
    # Ne20 + He4 --> p + Na23
    rate = 0.0

    # il10r
    rate += np.exp(  -6.37772 + -29.8896*tf.T9i + 19.7297*tf.T913
                  + -2.20987*tf.T9 + 0.153374*tf.T953 + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  0.227472 + -29.4348*tf.T9i
                  + -1.5*tf.lnT9)
    # il10n
    rate += np.exp(  19.1852 + -27.5738*tf.T9i + -20.0024*tf.T913i + 11.5988*tf.T913
                  + -1.37398*tf.T9 + -1.0*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Ne20__p_Na23 = rate

@numba.njit()
def He4_Ne20__C12_C12(rate_eval, tf):
    # Ne20 + He4 --> C12 + C12
    rate = 0.0

    # cf88r
    rate += np.exp(  61.4748 + -53.6267*tf.T9i + -84.165*tf.T913i + -1.56627*tf.T913
                  + -0.0736084*tf.T9 + -0.072797*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Ne20__C12_C12 = rate

@numba.njit()
def C12_Ne20__n_S31(rate_eval, tf):
    # Ne20 + C12 --> n + S31
    rate = 0.0

    # rolfr
    rate += np.exp(  -342.129 + -54.118*tf.T9i + 638.957*tf.T913i + -289.04*tf.T913
                  + -3.53144*tf.T9 + 0.648991*tf.T953 + 297.721*tf.lnT9)

    rate_eval.C12_Ne20__n_S31 = rate

@numba.njit()
def C12_Ne20__p_P31(rate_eval, tf):
    # Ne20 + C12 --> p + P31
    rate = 0.0

    # rolfr
    rate += np.exp(  -268.136 + -38.7624*tf.T9i + 361.154*tf.T913i + -92.643*tf.T913
                  + -9.98738*tf.T9 + 0.892737*tf.T953 + 161.042*tf.lnT9)

    rate_eval.C12_Ne20__p_P31 = rate

@numba.njit()
def C12_Ne20__He4_Si28(rate_eval, tf):
    # Ne20 + C12 --> He4 + Si28
    rate = 0.0

    # rolfr
    rate += np.exp(  -308.905 + -47.2175*tf.T9i + 514.197*tf.T913i + -200.896*tf.T913
                  + -6.42713*tf.T9 + 0.758256*tf.T953 + 236.359*tf.lnT9)

    rate_eval.C12_Ne20__He4_Si28 = rate

@numba.njit()
def n_Ne21__He4_O18(rate_eval, tf):
    # Ne21 + n --> He4 + O18
    rate = 0.0

    # nacr 
    rate += np.exp(  1.11027 + -7.35566*tf.T9i + 11.656*tf.T913
                  + -1.24869*tf.T9 + 1.0*tf.lnT9)
    # nacrn
    rate += np.exp(  2.59673 + 0.0033365*tf.T9i)
    # nacr 
    rate += np.exp(  2.28452 + -0.331353*tf.T9i + 6.51299*tf.T913
                  + -0.334277*tf.T9 + 2.0*tf.lnT9)

    rate_eval.n_Ne21__He4_O18 = rate

@numba.njit()
def p_Ne21__n_Na21(rate_eval, tf):
    # Ne21 + p --> n + Na21
    rate = 0.0

    # ths8r
    rate += np.exp(  17.9167 + -50.247*tf.T9i + 2.17917*tf.T913
                  + -0.303135*tf.T9 + 0.019321*tf.T953)

    rate_eval.p_Ne21__n_Na21 = rate

@numba.njit()
def p_Ne21__He4_F18(rate_eval, tf):
    # Ne21 + p --> He4 + F18
    rate = 0.0

    # rpsmr
    rate += np.exp(  50.6536 + -22.049*tf.T9i + 21.4461*tf.T913i + -73.252*tf.T913
                  + 2.42329*tf.T9 + -0.077278*tf.T953 + 40.7604*tf.lnT9)

    rate_eval.p_Ne21__He4_F18 = rate

@numba.njit()
def He4_Ne21__n_Mg24(rate_eval, tf):
    # Ne21 + He4 --> n + Mg24
    rate = 0.0

    # nacrr
    rate += np.exp(  -7.26831 + -13.2638*tf.T9i + 18.0748*tf.T913
                  + -0.981883*tf.T9 + -1.5*tf.lnT9)
    # nacrn
    rate += np.exp(  43.9762 + -46.88*tf.T913i + -0.536629*tf.T913
                  + 0.144715*tf.T9 + -0.197624*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Ne21__n_Mg24 = rate

@numba.njit()
def He4_Ne21__p_Na24(rate_eval, tf):
    # Ne21 + He4 --> p + Na24
    rate = 0.0

    # rath 
    rate += np.exp(  46.9899 + -27.7425*tf.T9i + 79.8849*tf.T913i + -124.834*tf.T913
                  + 6.1655*tf.T9 + -0.339742*tf.T953 + 67.7826*tf.lnT9)

    rate_eval.He4_Ne21__p_Na24 = rate

@numba.njit()
def n_Na21__p_Ne21(rate_eval, tf):
    # Na21 + n --> p + Ne21
    rate = 0.0

    # ths8r
    rate += np.exp(  17.9167 + 2.17917*tf.T913
                  + -0.303135*tf.T9 + 0.019321*tf.T953)

    rate_eval.n_Na21__p_Ne21 = rate

@numba.njit()
def n_Na21__He4_F18(rate_eval, tf):
    # Na21 + n --> He4 + F18
    rate = 0.0

    # ths8r
    rate += np.exp(  15.7482 + 1.06229*tf.T913
                  + 0.212448*tf.T9 + -0.020685*tf.T953)

    rate_eval.n_Na21__He4_F18 = rate

@numba.njit()
def p_Na21__He4_Ne18(rate_eval, tf):
    # Na21 + p --> He4 + Ne18
    rate = 0.0

    # mo14r
    rate += np.exp(  -23.8123 + -36.3851*tf.T9i + 8.5632e-08*tf.T913i + 24.8579*tf.T913
                  + 0.0823845*tf.T9 + -0.365374*tf.T953 + -2.21415e-06*tf.lnT9)
    # mo14r
    rate += np.exp(  -22.6701 + -40.8761*tf.T9i + 4.73034e-07*tf.T913i + 36.29*tf.T913
                  + -6.56565*tf.T9 + -2.96287e-06*tf.T953 + -9.00373e-07*tf.lnT9)
    # mo14r
    rate += np.exp(  -37.6291 + -34.4844*tf.T9i + 7.09521e-06*tf.T913i + 34.1789*tf.T913
                  + -1.72974*tf.T9 + -0.0395081*tf.T953 + -7.82759e-07*tf.lnT9)

    rate_eval.p_Na21__He4_Ne18 = rate

@numba.njit()
def He4_Na21__p_Mg24(rate_eval, tf):
    # Na21 + He4 --> p + Mg24
    rate = 0.0

    # nacr 
    rate += np.exp(  39.8144 + -49.9621*tf.T913i + 5.90498*tf.T913
                  + -1.6598*tf.T9 + 0.117817*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Na21__p_Mg24 = rate

@numba.njit()
def n_Na22__He4_F19(rate_eval, tf):
    # Na22 + n --> He4 + F19
    rate = 0.0

    # cf88 
    rate += np.exp(  6.59038 + 6.22036*tf.T913
                  + -0.0209891*tf.T9 + -0.0689843*tf.T953 + 1.0*tf.lnT9)
    # cf88n
    rate += np.exp(  15.0061 + -0.0532251*tf.T913
                  + 1.00044*tf.T9 + -0.13238*tf.T953)

    rate_eval.n_Na22__He4_F19 = rate

@numba.njit()
def p_Na22__n_Mg22(rate_eval, tf):
    # Na22 + p --> n + Mg22
    rate = 0.0

    # ths8r
    rate += np.exp(  15.9053 + -64.6118*tf.T9i + 2.67784*tf.T913
                  + -0.418718*tf.T9 + 0.031411*tf.T953)

    rate_eval.p_Na22__n_Mg22 = rate

@numba.njit()
def p_Na22__He4_Ne19(rate_eval, tf):
    # Na22 + p --> He4 + Ne19
    rate = 0.0

    # ths8r
    rate += np.exp(  43.101 + -24.0192*tf.T9i + -46.6346*tf.T913i + 0.866532*tf.T913
                  + -0.893541*tf.T9 + 0.0747971*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Na22__He4_Ne19 = rate

@numba.njit()
def He4_Na22__n_Al25(rate_eval, tf):
    # Na22 + He4 --> n + Al25
    rate = 0.0

    # ths8r
    rate += np.exp(  7.59058 + -22.1956*tf.T9i + 2.92382*tf.T913
                  + 0.706669*tf.T9 + -0.0950292*tf.T953)

    rate_eval.He4_Na22__n_Al25 = rate

@numba.njit()
def He4_Na22__p_Mg25(rate_eval, tf):
    # Na22 + He4 --> p + Mg25
    rate = 0.0

    # ths8r
    rate += np.exp(  44.973 + -50.0924*tf.T913i + 0.807739*tf.T913
                  + -0.956029*tf.T9 + 0.0793321*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Na22__p_Mg25 = rate

@numba.njit()
def p_Na23__n_Mg23(rate_eval, tf):
    # Na23 + p --> n + Mg23
    rate = 0.0

    # nacr 
    rate += np.exp(  19.4638 + -56.1542*tf.T9i + 0.993488*tf.T913
                  + -0.257094*tf.T9 + 0.0284334*tf.T953)

    rate_eval.p_Na23__n_Mg23 = rate

@numba.njit()
def p_Na23__He4_Ne20(rate_eval, tf):
    # Na23 + p --> He4 + Ne20
    rate = 0.0

    # il10r
    rate += np.exp(  -6.58736 + -2.31577*tf.T9i + 19.7297*tf.T913
                  + -2.20987*tf.T9 + 0.153374*tf.T953 + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  0.0178295 + -1.86103*tf.T9i
                  + -1.5*tf.lnT9)
    # il10n
    rate += np.exp(  18.9756 + -20.0024*tf.T913i + 11.5988*tf.T913
                  + -1.37398*tf.T9 + -1.0*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Na23__He4_Ne20 = rate

@numba.njit()
def p_Na23__C12_C12(rate_eval, tf):
    # Na23 + p --> C12 + C12
    rate = 0.0

    # cf88r
    rate += np.exp(  60.9438 + -26.0184*tf.T9i + -84.165*tf.T913i + -1.4191*tf.T913
                  + -0.114619*tf.T9 + -0.070307*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Na23__C12_C12 = rate

@numba.njit()
def He4_Na23__n_Al26(rate_eval, tf):
    # Na23 + He4 --> n + Al26
    rate = 0.0

    # ol11r
    rate += np.exp(  14.5219 + -34.8285*tf.T9i
                  + -1.5*tf.lnT9)
    # ol11r
    rate += np.exp(  13.4292 + -34.4845*tf.T9i
                  + -1.5*tf.lnT9)
    # ol11n
    rate += np.exp(  11.4506 + -34.4184*tf.T9i + 5.07134*tf.T913
                  + -0.557537*tf.T9 + 0.0451737*tf.T953)

    rate_eval.He4_Na23__n_Al26 = rate

@numba.njit()
def He4_Na23__p_Mg26(rate_eval, tf):
    # Na23 + He4 --> p + Mg26
    rate = 0.0

    # ths8r
    rate += np.exp(  44.527 + -50.2042*tf.T913i + 1.76141*tf.T913
                  + -1.36813*tf.T9 + 0.123087*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Na23__p_Mg26 = rate

@numba.njit()
def p_Na24__n_Mg24(rate_eval, tf):
    # Na24 + p --> n + Mg24
    rate = 0.0

    # rath 
    rate += np.exp(  64.9407 + -2.51894*tf.T9i + 85.1374*tf.T913i + -145.038*tf.T913
                  + 7.70922*tf.T9 + -0.427731*tf.T953 + 73.5924*tf.lnT9)

    rate_eval.p_Na24__n_Mg24 = rate

@numba.njit()
def p_Na24__He4_Ne21(rate_eval, tf):
    # Na24 + p --> He4 + Ne21
    rate = 0.0

    # rath 
    rate += np.exp(  47.365 + -2.4679*tf.T9i + 79.8849*tf.T913i + -124.834*tf.T913
                  + 6.1655*tf.T9 + -0.339742*tf.T953 + 67.7826*tf.lnT9)

    rate_eval.p_Na24__He4_Ne21 = rate

@numba.njit()
def He4_Na24__n_Al27(rate_eval, tf):
    # Na24 + He4 --> n + Al27
    rate = 0.0

    # rath 
    rate += np.exp(  67.5936 + -2.6585*tf.T9i + 49.3636*tf.T913i + -124.715*tf.T913
                  + 5.57965*tf.T9 + -0.275543*tf.T953 + 65.8717*tf.lnT9)

    rate_eval.He4_Na24__n_Al27 = rate

@numba.njit()
def n_Mg22__p_Na22(rate_eval, tf):
    # Mg22 + n --> p + Na22
    rate = 0.0

    # ths8r
    rate += np.exp(  17.8512 + 2.67784*tf.T913
                  + -0.418718*tf.T9 + 0.031411*tf.T953)

    rate_eval.n_Mg22__p_Na22 = rate

@numba.njit()
def n_Mg22__He4_Ne19(rate_eval, tf):
    # Mg22 + n --> He4 + Ne19
    rate = 0.0

    # ths8r
    rate += np.exp(  16.2584 + 2.08592*tf.T913
                  + -0.188654*tf.T9 + 0.0101103*tf.T953)

    rate_eval.n_Mg22__He4_Ne19 = rate

@numba.njit()
def He4_Mg22__p_Al25(rate_eval, tf):
    # Mg22 + He4 --> p + Al25
    rate = 0.0

    # ths8r
    rate += np.exp(  46.0189 + -53.0847*tf.T913i + 0.634581*tf.T913
                  + -0.829647*tf.T9 + 0.0650313*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Mg22__p_Al25 = rate

@numba.njit()
def n_Mg23__p_Na23(rate_eval, tf):
    # Mg23 + n --> p + Na23
    rate = 0.0

    # nacr 
    rate += np.exp(  19.4638 + 0.993488*tf.T913
                  + -0.257094*tf.T9 + 0.0284334*tf.T953)

    rate_eval.n_Mg23__p_Na23 = rate

@numba.njit()
def n_Mg23__He4_Ne20(rate_eval, tf):
    # Mg23 + n --> He4 + Ne20
    rate = 0.0

    # ths8r
    rate += np.exp(  17.7448 + 1.83199*tf.T913
                  + -0.290485*tf.T9 + 0.0242929*tf.T953)

    rate_eval.n_Mg23__He4_Ne20 = rate

@numba.njit()
def n_Mg23__C12_C12(rate_eval, tf):
    # Mg23 + n --> C12 + C12
    rate = 0.0

    # cf88r
    rate += np.exp(  -12.8267 + 11.4826*tf.T913
                  + 1.82849*tf.T9 + -0.34844*tf.T953)

    rate_eval.n_Mg23__C12_C12 = rate

@numba.njit()
def He4_Mg23__n_Si26(rate_eval, tf):
    # Mg23 + He4 --> n + Si26
    rate = 0.0

    # ths8r
    rate += np.exp(  13.5183 + -46.1342*tf.T9i + 1.08923*tf.T913
                  + -0.0248723*tf.T9 + 0.00450822*tf.T953)

    rate_eval.He4_Mg23__n_Si26 = rate

@numba.njit()
def He4_Mg23__p_Al26(rate_eval, tf):
    # Mg23 + He4 --> p + Al26
    rate = 0.0

    # ths8r
    rate += np.exp(  46.215 + -53.203*tf.T913i + 0.71292*tf.T913
                  + -0.892548*tf.T9 + 0.0709813*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Mg23__p_Al26 = rate

@numba.njit()
def n_Mg24__p_Na24(rate_eval, tf):
    # Mg24 + n --> p + Na24
    rate = 0.0

    # rath 
    rate += np.exp(  67.1379 + -57.443*tf.T9i + 85.1374*tf.T913i + -145.038*tf.T913
                  + 7.70922*tf.T9 + -0.427731*tf.T953 + 73.5924*tf.lnT9)

    rate_eval.n_Mg24__p_Na24 = rate

@numba.njit()
def n_Mg24__He4_Ne21(rate_eval, tf):
    # Mg24 + n --> He4 + Ne21
    rate = 0.0

    # nacrr
    rate += np.exp(  -4.69602 + -42.9133*tf.T9i + 18.0748*tf.T913
                  + -0.981883*tf.T9 + -1.5*tf.lnT9)
    # nacrn
    rate += np.exp(  46.5485 + -29.6495*tf.T9i + -46.88*tf.T913i + -0.536629*tf.T913
                  + 0.144715*tf.T9 + -0.197624*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.n_Mg24__He4_Ne21 = rate

@numba.njit()
def p_Mg24__He4_Na21(rate_eval, tf):
    # Mg24 + p --> He4 + Na21
    rate = 0.0

    # nacr 
    rate += np.exp(  42.3867 + -79.897*tf.T9i + -49.9621*tf.T913i + 5.90498*tf.T913
                  + -1.6598*tf.T9 + 0.117817*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Mg24__He4_Na21 = rate

@numba.njit()
def He4_Mg24__p_Al27(rate_eval, tf):
    # Mg24 + He4 --> p + Al27
    rate = 0.0

    # il10r
    rate += np.exp(  -26.2862 + -19.5422*tf.T9i + 5.18642*tf.T913i + -34.7936*tf.T913
                  + 168.225*tf.T9 + -115.825*tf.T953 + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  -6.44575 + -22.8216*tf.T9i + 18.0416*tf.T913
                  + -1.54137*tf.T9 + 0.0847506*tf.T953 + -1.5*tf.lnT9)
    # il10n
    rate += np.exp(  30.0397 + -18.5791*tf.T9i + -26.4162*tf.T913i
                  + -2.0*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Mg24__p_Al27 = rate

@numba.njit()
def He4_Mg24__C12_O16(rate_eval, tf):
    # Mg24 + He4 --> C12 + O16
    rate = 0.0

    # cf88r
    rate += np.exp(  49.5738 + -78.202*tf.T9i + -133.413*tf.T913i + 50.1572*tf.T913
                  + -3.15987*tf.T9 + 0.0178251*tf.T953 + -23.7027*tf.lnT9)

    rate_eval.He4_Mg24__C12_O16 = rate

@numba.njit()
def p_Mg25__n_Al25(rate_eval, tf):
    # Mg25 + p --> n + Al25
    rate = 0.0

    # ths8r
    rate += np.exp(  18.4104 + -58.7072*tf.T9i + 2.28536*tf.T913
                  + -0.38512*tf.T9 + 0.0288056*tf.T953)

    rate_eval.p_Mg25__n_Al25 = rate

@numba.njit()
def p_Mg25__He4_Na22(rate_eval, tf):
    # Mg25 + p --> He4 + Na22
    rate = 0.0

    # ths8r
    rate += np.exp(  46.3217 + -36.5117*tf.T9i + -50.0924*tf.T913i + 0.807739*tf.T913
                  + -0.956029*tf.T9 + 0.0793321*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Mg25__He4_Na22 = rate

@numba.njit()
def He4_Mg25__n_Si28(rate_eval, tf):
    # Mg25 + He4 --> n + Si28
    rate = 0.0

    # nacrn
    rate += np.exp(  45.7613 + -53.415*tf.T913i + -1.46489*tf.T913
                  + 1.7777*tf.T9 + -0.903499*tf.T953 + -0.666667*tf.lnT9)
    # nacr 
    rate += np.exp(  38.337 + -53.416*tf.T913i + 8.01209*tf.T913
                  + -2.64791*tf.T9 + 0.218637*tf.T953 + 0.333333*tf.lnT9)

    rate_eval.He4_Mg25__n_Si28 = rate

@numba.njit()
def He4_Mg25__p_Al28(rate_eval, tf):
    # Mg25 + He4 --> p + Al28
    rate = 0.0

    # cf88r
    rate += np.exp(  19.0822 + -13.995*tf.T9i + -23.271*tf.T913i + 5.6671*tf.T913
                  + -0.708839*tf.T9 + -0.0561047*tf.T953 + 1.33333*tf.lnT9)
    # cf88n
    rate += np.exp(  19.5932 + -13.995*tf.T9i + -23.271*tf.T913i + 1.22438*tf.T913
                  + 1.57567*tf.T9 + -0.209975*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Mg25__p_Al28 = rate

@numba.njit()
def p_Mg26__n_Al26(rate_eval, tf):
    # Mg26 + p --> n + Al26
    rate = 0.0

    # ol11n
    rate += np.exp(  13.5366 + -55.5463*tf.T9i + 7.36773*tf.T913
                  + -2.42424*tf.T9 + 0.313743*tf.T953)
    # ol11r
    rate += np.exp(  15.9629 + -55.7499*tf.T9i + 3.87515*tf.T913
                  + 0.228327*tf.T9 + -0.045872*tf.T953 + -1.5*tf.lnT9)
    # ol11r
    rate += np.exp(  14.8989 + -55.6072*tf.T9i
                  + -1.5*tf.lnT9)

    rate_eval.p_Mg26__n_Al26 = rate

@numba.njit()
def p_Mg26__He4_Na23(rate_eval, tf):
    # Mg26 + p --> He4 + Na23
    rate = 0.0

    # ths8r
    rate += np.exp(  47.1157 + -21.128*tf.T9i + -50.2042*tf.T913i + 1.76141*tf.T913
                  + -1.36813*tf.T9 + 0.123087*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Mg26__He4_Na23 = rate

@numba.njit()
def He4_Mg26__n_Si29(rate_eval, tf):
    # Mg26 + He4 --> n + Si29
    rate = 0.0

    # nacrn
    rate += np.exp(  45.8397 + -53.505*tf.T913i + 0.045598*tf.T913
                  + -0.194481*tf.T9 + -0.0748153*tf.T953 + -0.666667*tf.lnT9)
    # nacrr
    rate += np.exp(  15.7057 + -18.73*tf.T9i + -9.54056*tf.T913
                  + 4.71712*tf.T9 + -0.461053*tf.T953 + -1.5*tf.lnT9)

    rate_eval.He4_Mg26__n_Si29 = rate

@numba.njit()
def n_Al25__p_Mg25(rate_eval, tf):
    # Al25 + n --> p + Mg25
    rate = 0.0

    # ths8r
    rate += np.exp(  18.4104 + 2.28536*tf.T913
                  + -0.38512*tf.T9 + 0.0288056*tf.T953)

    rate_eval.n_Al25__p_Mg25 = rate

@numba.njit()
def n_Al25__He4_Na22(rate_eval, tf):
    # Al25 + n --> He4 + Na22
    rate = 0.0

    # ths8r
    rate += np.exp(  8.93927 + 2.92382*tf.T913
                  + 0.706669*tf.T9 + -0.0950292*tf.T953)

    rate_eval.n_Al25__He4_Na22 = rate

@numba.njit()
def p_Al25__He4_Mg22(rate_eval, tf):
    # Al25 + p --> He4 + Mg22
    rate = 0.0

    # ths8r
    rate += np.exp(  45.4217 + -42.4162*tf.T9i + -53.0847*tf.T913i + 0.634581*tf.T913
                  + -0.829647*tf.T9 + 0.0650313*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Al25__He4_Mg22 = rate

@numba.njit()
def He4_Al25__p_Si28(rate_eval, tf):
    # Al25 + He4 --> p + Si28
    rate = 0.0

    # ths8r
    rate += np.exp(  47.6167 + -56.3424*tf.T913i + 0.553763*tf.T913
                  + -0.84072*tf.T9 + 0.0634219*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Al25__p_Si28 = rate

@numba.njit()
def n_Al26__p_Mg26(rate_eval, tf):
    # Al26 + n --> p + Mg26
    rate = 0.0

    # ol11r
    rate += np.exp(  13.565 + -0.203581*tf.T9i + 3.87515*tf.T913
                  + 0.228327*tf.T9 + -0.045872*tf.T953 + -1.5*tf.lnT9)
    # ol11r
    rate += np.exp(  12.501 + -0.0608268*tf.T9i
                  + -1.5*tf.lnT9)
    # ol11n
    rate += np.exp(  11.1387 + 7.36773*tf.T913
                  + -2.42424*tf.T9 + 0.313743*tf.T953)

    rate_eval.n_Al26__p_Mg26 = rate

@numba.njit()
def n_Al26__He4_Na23(rate_eval, tf):
    # Al26 + n --> He4 + Na23
    rate = 0.0

    # ol11r
    rate += np.exp(  13.62 + -0.0661138*tf.T9i
                  + -1.5*tf.lnT9)
    # ol11n
    rate += np.exp(  11.6414 + 5.07134*tf.T913
                  + -0.557537*tf.T9 + 0.0451737*tf.T953)
    # ol11r
    rate += np.exp(  14.7127 + -0.41015*tf.T9i
                  + -1.5*tf.lnT9)

    rate_eval.n_Al26__He4_Na23 = rate

@numba.njit()
def p_Al26__n_Si26(rate_eval, tf):
    # Al26 + p --> n + Si26
    rate = 0.0

    # ths8r
    rate += np.exp(  15.9928 + -67.8633*tf.T9i + 2.01618*tf.T913
                  + -0.300371*tf.T9 + 0.0228631*tf.T953)

    rate_eval.p_Al26__n_Si26 = rate

@numba.njit()
def p_Al26__He4_Mg23(rate_eval, tf):
    # Al26 + p --> He4 + Mg23
    rate = 0.0

    # ths8r
    rate += np.exp(  46.4058 + -21.7293*tf.T9i + -53.203*tf.T913i + 0.71292*tf.T913
                  + -0.892548*tf.T9 + 0.0709813*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Al26__He4_Mg23 = rate

@numba.njit()
def He4_Al26__n_P29(rate_eval, tf):
    # Al26 + He4 --> n + P29
    rate = 0.0

    # ths8r
    rate += np.exp(  -42.1718 + -10.4914*tf.T9i + 41.9938*tf.T913
                  + -4.54859*tf.T9 + 0.243841*tf.T953)

    rate_eval.He4_Al26__n_P29 = rate

@numba.njit()
def He4_Al26__p_Si29(rate_eval, tf):
    # Al26 + He4 --> p + Si29
    rate = 0.0

    # ths8r
    rate += np.exp(  47.7092 + -56.4422*tf.T913i + 0.705353*tf.T913
                  + -0.957427*tf.T9 + 0.0756045*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Al26__p_Si29 = rate

@numba.njit()
def n_Al27__He4_Na24(rate_eval, tf):
    # Al27 + n --> He4 + Na24
    rate = 0.0

    # rath 
    rate += np.exp(  69.2087 + -39.0154*tf.T9i + 49.3636*tf.T913i + -124.715*tf.T913
                  + 5.57965*tf.T9 + -0.275543*tf.T953 + 65.8717*tf.lnT9)

    rate_eval.n_Al27__He4_Na24 = rate

@numba.njit()
def p_Al27__He4_Mg24(rate_eval, tf):
    # Al27 + p --> He4 + Mg24
    rate = 0.0

    # il10r
    rate += np.exp(  -7.02789 + -4.2425*tf.T9i + 18.0416*tf.T913
                  + -1.54137*tf.T9 + 0.0847506*tf.T953 + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  -26.8683 + -0.963012*tf.T9i + 5.18642*tf.T913i + -34.7936*tf.T913
                  + 168.225*tf.T9 + -115.825*tf.T953 + -1.5*tf.lnT9)
    # il10n
    rate += np.exp(  29.4576 + -26.4162*tf.T913i
                  + -2.0*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Al27__He4_Mg24 = rate

@numba.njit()
def p_Al27__C12_O16(rate_eval, tf):
    # Al27 + p --> C12 + O16
    rate = 0.0

    # cf88r
    rate += np.exp(  68.9829 + -59.8017*tf.T9i + -119.242*tf.T913i + 13.3667*tf.T913
                  + 0.295425*tf.T9 + -0.267288*tf.T953 + -9.91729*tf.lnT9)

    rate_eval.p_Al27__C12_O16 = rate

@numba.njit()
def He4_Al27__n_P30(rate_eval, tf):
    # Al27 + He4 --> n + P30
    rate = 0.0

    # nacr 
    rate += np.exp(  6.09163 + -30.667*tf.T9i + 5.40982*tf.T913
                  + -0.265676*tf.T9 + 1.0*tf.lnT9)
    # nacr 
    rate += np.exp(  11.3074 + -30.667*tf.T9i + 0.311974*tf.T913
                  + -2.02044*tf.T9)

    rate_eval.He4_Al27__n_P30 = rate

@numba.njit()
def He4_Al27__p_Si30(rate_eval, tf):
    # Al27 + He4 --> p + Si30
    rate = 0.0

    # ths8r
    rate += np.exp(  47.4856 + -56.5351*tf.T913i + 1.60477*tf.T913
                  + -1.40594*tf.T9 + 0.127353*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Al27__p_Si30 = rate

@numba.njit()
def p_Al28__n_Si28(rate_eval, tf):
    # Al28 + p --> n + Si28
    rate = 0.0

    # rath 
    rate += np.exp(  72.7891 + -3.04903*tf.T9i + 105.263*tf.T913i + -174.93*tf.T913
                  + 9.08555*tf.T9 + -0.495072*tf.T953 + 89.3425*tf.lnT9)

    rate_eval.p_Al28__n_Si28 = rate

@numba.njit()
def p_Al28__He4_Mg25(rate_eval, tf):
    # Al28 + p --> He4 + Mg25
    rate = 0.0

    # cf88n
    rate += np.exp(  20.6553 + -23.271*tf.T913i + 1.22438*tf.T913
                  + 1.57567*tf.T9 + -0.209975*tf.T953 + -0.666667*tf.lnT9)
    # cf88r
    rate += np.exp(  20.1444 + -23.271*tf.T913i + 5.6671*tf.T913
                  + -0.708839*tf.T9 + -0.0561047*tf.T953 + 1.33333*tf.lnT9)

    rate_eval.p_Al28__He4_Mg25 = rate

@numba.njit()
def He4_Al28__n_P31(rate_eval, tf):
    # Al28 + He4 --> n + P31
    rate = 0.0

    # rath 
    rate += np.exp(  -4.29707 + -2.4128*tf.T9i + -27.1835*tf.T913i + 29.3511*tf.T913
                  + -4.90371*tf.T9 + 0.376759*tf.T953 + 0.485028*tf.lnT9)

    rate_eval.He4_Al28__n_P31 = rate

@numba.njit()
def n_Si26__p_Al26(rate_eval, tf):
    # Si26 + n --> p + Al26
    rate = 0.0

    # ths8r
    rate += np.exp(  18.3907 + 2.01618*tf.T913
                  + -0.300371*tf.T9 + 0.0228631*tf.T953)

    rate_eval.n_Si26__p_Al26 = rate

@numba.njit()
def n_Si26__He4_Mg23(rate_eval, tf):
    # Si26 + n --> He4 + Mg23
    rate = 0.0

    # ths8r
    rate += np.exp(  16.107 + 1.08923*tf.T913
                  + -0.0248723*tf.T9 + 0.00450822*tf.T953)

    rate_eval.n_Si26__He4_Mg23 = rate

@numba.njit()
def He4_Si26__p_P29(rate_eval, tf):
    # Si26 + He4 --> p + P29
    rate = 0.0

    # ths8r
    rate += np.exp(  48.8732 + -59.3013*tf.T913i + 0.480742*tf.T913
                  + -0.834505*tf.T9 + 0.0621841*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Si26__p_P29 = rate

@numba.njit()
def n_Si28__p_Al28(rate_eval, tf):
    # Si28 + n --> p + Al28
    rate = 0.0

    # rath 
    rate += np.exp(  74.735 + -47.8424*tf.T9i + 105.263*tf.T913i + -174.93*tf.T913
                  + 9.08555*tf.T9 + -0.495072*tf.T953 + 89.3425*tf.lnT9)

    rate_eval.n_Si28__p_Al28 = rate

@numba.njit()
def n_Si28__He4_Mg25(rate_eval, tf):
    # Si28 + n --> He4 + Mg25
    rate = 0.0

    # nacrn
    rate += np.exp(  48.7694 + -30.7983*tf.T9i + -53.415*tf.T913i + -1.46489*tf.T913
                  + 1.7777*tf.T9 + -0.903499*tf.T953 + -0.666667*tf.lnT9)
    # nacr 
    rate += np.exp(  41.3451 + -30.7983*tf.T9i + -53.416*tf.T913i + 8.01209*tf.T913
                  + -2.64791*tf.T9 + 0.218637*tf.T953 + 0.333333*tf.lnT9)

    rate_eval.n_Si28__He4_Mg25 = rate

@numba.njit()
def p_Si28__He4_Al25(rate_eval, tf):
    # Si28 + p --> He4 + Al25
    rate = 0.0

    # ths8r
    rate += np.exp(  50.6248 + -89.5005*tf.T9i + -56.3424*tf.T913i + 0.553763*tf.T913
                  + -0.84072*tf.T9 + 0.0634219*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Si28__He4_Al25 = rate

@numba.njit()
def He4_Si28__n_S31(rate_eval, tf):
    # Si28 + He4 --> n + S31
    rate = 0.0

    # ths8r
    rate += np.exp(  17.5194 + -93.9332*tf.T9i + 1.12938*tf.T913
                  + -0.122539*tf.T9 + 0.00424404*tf.T953)

    rate_eval.He4_Si28__n_S31 = rate

@numba.njit()
def He4_Si28__p_P31(rate_eval, tf):
    # Si28 + He4 --> p + P31
    rate = 0.0

    # il10r
    rate += np.exp(  -13.4595 + -24.112*tf.T9i
                  + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  -11.4335 + -25.6606*tf.T9i + 21.521*tf.T913
                  + -1.90355*tf.T9 + 0.092724*tf.T953 + -1.5*tf.lnT9)
    # il10n
    rate += np.exp(  60.3424 + -22.2348*tf.T9i + -31.932*tf.T913i + -77.0334*tf.T913
                  + -43.6847*tf.T9 + -4.28955*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Si28__p_P31 = rate

@numba.njit()
def He4_Si28__C12_Ne20(rate_eval, tf):
    # Si28 + He4 --> C12 + Ne20
    rate = 0.0

    # rolfr
    rate += np.exp(  -307.762 + -186.722*tf.T9i + 514.197*tf.T913i + -200.896*tf.T913
                  + -6.42713*tf.T9 + 0.758256*tf.T953 + 236.359*tf.lnT9)

    rate_eval.He4_Si28__C12_Ne20 = rate

@numba.njit()
def He4_Si28__O16_O16(rate_eval, tf):
    # Si28 + He4 --> O16 + O16
    rate = 0.0

    # cf88r
    rate += np.exp(  97.7904 + -111.595*tf.T9i + -119.324*tf.T913i + -32.2497*tf.T913
                  + 1.46214*tf.T9 + -0.200893*tf.T953 + 13.2148*tf.lnT9)

    rate_eval.He4_Si28__O16_O16 = rate

@numba.njit()
def n_Si29__He4_Mg26(rate_eval, tf):
    # Si29 + n --> He4 + Mg26
    rate = 0.0

    # nacrn
    rate += np.exp(  46.369 + -0.394553*tf.T9i + -53.505*tf.T913i + 0.045598*tf.T913
                  + -0.194481*tf.T9 + -0.0748153*tf.T953 + -0.666667*tf.lnT9)
    # nacrr
    rate += np.exp(  16.235 + -19.1246*tf.T9i + -9.54056*tf.T913
                  + 4.71712*tf.T9 + -0.461053*tf.T953 + -1.5*tf.lnT9)

    rate_eval.n_Si29__He4_Mg26 = rate

@numba.njit()
def p_Si29__n_P29(rate_eval, tf):
    # Si29 + p --> n + P29
    rate = 0.0

    # ths8r
    rate += np.exp(  18.5351 + -66.4331*tf.T9i + 1.48018*tf.T913
                  + -0.177129*tf.T9 + 0.0127163*tf.T953)

    rate_eval.p_Si29__n_P29 = rate

@numba.njit()
def p_Si29__He4_Al26(rate_eval, tf):
    # Si29 + p --> He4 + Al26
    rate = 0.0

    # ths8r
    rate += np.exp(  50.6364 + -55.9416*tf.T9i + -56.4422*tf.T913i + 0.705353*tf.T913
                  + -0.957427*tf.T9 + 0.0756045*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Si29__He4_Al26 = rate

@numba.njit()
def He4_Si29__n_S32(rate_eval, tf):
    # Si29 + He4 --> n + S32
    rate = 0.0

    # ths8r
    rate += np.exp(  -2.87932 + -17.7056*tf.T9i + 9.48125*tf.T913
                  + 0.4472*tf.T9 + -0.119237*tf.T953)

    rate_eval.He4_Si29__n_S32 = rate

@numba.njit()
def He4_Si29__p_P32(rate_eval, tf):
    # Si29 + He4 --> p + P32
    rate = 0.0

    # rath 
    rate += np.exp(  117.4 + -32.3767*tf.T9i + 160.025*tf.T913i + -285.126*tf.T913
                  + 16.3226*tf.T9 + -0.928477*tf.T953 + 138.458*tf.lnT9)

    rate_eval.He4_Si29__p_P32 = rate

@numba.njit()
def p_Si30__n_P30(rate_eval, tf):
    # Si30 + p --> n + P30
    rate = 0.0

    # ths8r
    rate += np.exp(  19.4345 + -58.1931*tf.T9i + 1.88379*tf.T913
                  + -0.330243*tf.T9 + 0.0239836*tf.T953)

    rate_eval.p_Si30__n_P30 = rate

@numba.njit()
def p_Si30__He4_Al27(rate_eval, tf):
    # Si30 + p --> He4 + Al27
    rate = 0.0

    # ths8r
    rate += np.exp(  50.5056 + -27.5284*tf.T9i + -56.5351*tf.T913i + 1.60477*tf.T913
                  + -1.40594*tf.T9 + 0.127353*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Si30__He4_Al27 = rate

@numba.njit()
def n_P29__p_Si29(rate_eval, tf):
    # P29 + n --> p + Si29
    rate = 0.0

    # ths8r
    rate += np.exp(  18.5351 + 1.48018*tf.T913
                  + -0.177129*tf.T9 + 0.0127163*tf.T953)

    rate_eval.n_P29__p_Si29 = rate

@numba.njit()
def n_P29__He4_Al26(rate_eval, tf):
    # P29 + n --> He4 + Al26
    rate = 0.0

    # ths8r
    rate += np.exp(  -39.2446 + 41.9938*tf.T913
                  + -4.54859*tf.T9 + 0.243841*tf.T953)

    rate_eval.n_P29__He4_Al26 = rate

@numba.njit()
def p_P29__He4_Si26(rate_eval, tf):
    # P29 + p --> He4 + Si26
    rate = 0.0

    # ths8r
    rate += np.exp(  49.4025 + -57.372*tf.T9i + -59.3013*tf.T913i + 0.480742*tf.T913
                  + -0.834505*tf.T9 + 0.0621841*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_P29__He4_Si26 = rate

@numba.njit()
def He4_P29__n_Cl32(rate_eval, tf):
    # P29 + He4 --> n + Cl32
    rate = 0.0

    # ths8r
    rate += np.exp(  16.1204 + -107.565*tf.T9i + 1.0254*tf.T913
                  + -0.188653*tf.T9 + 0.0264376*tf.T953)

    rate_eval.He4_P29__n_Cl32 = rate

@numba.njit()
def He4_P29__p_S32(rate_eval, tf):
    # P29 + He4 --> p + S32
    rate = 0.0

    # ths8r
    rate += np.exp(  50.2824 + -62.3802*tf.T913i + 0.459085*tf.T913
                  + -0.870169*tf.T9 + 0.0631143*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_P29__p_S32 = rate

@numba.njit()
def n_P30__p_Si30(rate_eval, tf):
    # P30 + n --> p + Si30
    rate = 0.0

    # ths8r
    rate += np.exp(  18.3359 + 1.88379*tf.T913
                  + -0.330243*tf.T9 + 0.0239836*tf.T953)

    rate_eval.n_P30__p_Si30 = rate

@numba.njit()
def n_P30__He4_Al27(rate_eval, tf):
    # P30 + n --> He4 + Al27
    rate = 0.0

    # nacr 
    rate += np.exp(  8.01303 + 0.0036935*tf.T9i + 5.40982*tf.T913
                  + -0.265676*tf.T9 + 1.0*tf.lnT9)
    # nacr 
    rate += np.exp(  13.2288 + 0.0036935*tf.T9i + 0.311974*tf.T913
                  + -2.02044*tf.T9)

    rate_eval.n_P30__He4_Al27 = rate

@numba.njit()
def n_P31__He4_Al28(rate_eval, tf):
    # P31 + n --> He4 + Al28
    rate = 0.0

    # rath 
    rate += np.exp(  -1.81069 + -24.9603*tf.T9i + -27.1835*tf.T913i + 29.3511*tf.T913
                  + -4.90371*tf.T9 + 0.376759*tf.T953 + 0.485028*tf.lnT9)

    rate_eval.n_P31__He4_Al28 = rate

@numba.njit()
def p_P31__n_S31(rate_eval, tf):
    # P31 + p --> n + S31
    rate = 0.0

    # ths8r
    rate += np.exp(  17.7037 + -71.6993*tf.T9i + 2.15006*tf.T913
                  + -0.300369*tf.T9 + 0.0229589*tf.T953)

    rate_eval.p_P31__n_S31 = rate

@numba.njit()
def p_P31__He4_Si28(rate_eval, tf):
    # P31 + p --> He4 + Si28
    rate = 0.0

    # il10r
    rate += np.exp(  -10.893 + -3.42575*tf.T9i + 21.521*tf.T913
                  + -1.90355*tf.T9 + 0.092724*tf.T953 + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  -12.919 + -1.87716*tf.T9i
                  + -1.5*tf.lnT9)
    # il10n
    rate += np.exp(  60.8829 + -31.932*tf.T913i + -77.0334*tf.T913
                  + -43.6847*tf.T9 + -4.28955*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_P31__He4_Si28 = rate

@numba.njit()
def p_P31__C12_Ne20(rate_eval, tf):
    # P31 + p --> C12 + Ne20
    rate = 0.0

    # rolfr
    rate += np.exp(  -266.452 + -156.019*tf.T9i + 361.154*tf.T913i + -92.643*tf.T913
                  + -9.98738*tf.T9 + 0.892737*tf.T953 + 161.042*tf.lnT9)

    rate_eval.p_P31__C12_Ne20 = rate

@numba.njit()
def p_P31__O16_O16(rate_eval, tf):
    # P31 + p --> O16 + O16
    rate = 0.0

    # cf88r
    rate += np.exp(  86.3501 + -88.8797*tf.T9i + -145.844*tf.T913i + 8.72612*tf.T913
                  + -0.554035*tf.T9 + -0.137562*tf.T953 + -6.88807*tf.lnT9)

    rate_eval.p_P31__O16_O16 = rate

@numba.njit()
def He4_P31__n_Cl34(rate_eval, tf):
    # P31 + He4 --> n + Cl34
    rate = 0.0

    # ths8r
    rate += np.exp(  16.9244 + -65.5365*tf.T9i + -0.466069*tf.T913
                  + 0.167169*tf.T9 + -0.00537463*tf.T953)

    rate_eval.He4_P31__n_Cl34 = rate

@numba.njit()
def He4_P31__p_S34(rate_eval, tf):
    # P31 + He4 --> p + S34
    rate = 0.0

    # ths8r
    rate += np.exp(  49.0993 + -62.5433*tf.T913i + 1.69339*tf.T913
                  + -1.27872*tf.T9 + 0.106616*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_P31__p_S34 = rate

@numba.njit()
def p_P32__n_S32(rate_eval, tf):
    # P32 + p --> n + S32
    rate = 0.0

    # rath 
    rate += np.exp(  82.294 + -3.56696*tf.T9i + 124.605*tf.T913i + -206.844*tf.T913
                  + 11.0016*tf.T9 + -0.612898*tf.T953 + 104.905*tf.lnT9)

    rate_eval.p_P32__n_S32 = rate

@numba.njit()
def p_P32__He4_Si29(rate_eval, tf):
    # P32 + p --> He4 + Si29
    rate = 0.0

    # rath 
    rate += np.exp(  118.233 + -3.89928*tf.T9i + 160.025*tf.T913i + -285.126*tf.T913
                  + 16.3226*tf.T9 + -0.928477*tf.T953 + 138.458*tf.lnT9)

    rate_eval.p_P32__He4_Si29 = rate

@numba.njit()
def He4_P32__n_Cl35(rate_eval, tf):
    # P32 + He4 --> n + Cl35
    rate = 0.0

    # ths8r
    rate += np.exp(  -33.0327 + -10.8859*tf.T9i + 35.1367*tf.T913
                  + -3.30169*tf.T9 + 0.150699*tf.T953)

    rate_eval.He4_P32__n_Cl35 = rate

@numba.njit()
def n_S31__p_P31(rate_eval, tf):
    # S31 + n --> p + P31
    rate = 0.0

    # ths8r
    rate += np.exp(  17.7037 + 2.15006*tf.T913
                  + -0.300369*tf.T9 + 0.0229589*tf.T953)

    rate_eval.n_S31__p_P31 = rate

@numba.njit()
def n_S31__He4_Si28(rate_eval, tf):
    # S31 + n --> He4 + Si28
    rate = 0.0

    # ths8r
    rate += np.exp(  18.0599 + 1.12938*tf.T913
                  + -0.122539*tf.T9 + 0.00424404*tf.T953)

    rate_eval.n_S31__He4_Si28 = rate

@numba.njit()
def n_S31__C12_Ne20(rate_eval, tf):
    # S31 + n --> C12 + Ne20
    rate = 0.0

    # rolfr
    rate += np.exp(  -340.445 + -99.6981*tf.T9i + 638.957*tf.T913i + -289.04*tf.T913
                  + -3.53144*tf.T9 + 0.648991*tf.T953 + 297.721*tf.lnT9)

    rate_eval.n_S31__C12_Ne20 = rate

@numba.njit()
def n_S31__O16_O16(rate_eval, tf):
    # S31 + n --> O16 + O16
    rate = 0.0

    # cf88r
    rate += np.exp(  78.6364 + -17.7811*tf.T9i + -120.83*tf.T913i + -7.72334*tf.T913
                  + -2.27939*tf.T9 + 0.167655*tf.T953 + 7.62001*tf.lnT9)

    rate_eval.n_S31__O16_O16 = rate

@numba.njit()
def He4_S31__p_Cl34(rate_eval, tf):
    # S31 + He4 --> p + Cl34
    rate = 0.0

    # ths8r
    rate += np.exp(  50.2997 + -65.2934*tf.T913i + 1.2687*tf.T913
                  + -1.08925*tf.T9 + 0.0898609*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_S31__p_Cl34 = rate

@numba.njit()
def n_S32__p_P32(rate_eval, tf):
    # S32 + n --> p + P32
    rate = 0.0

    # rath 
    rate += np.exp(  83.3926 + -14.3359*tf.T9i + 124.605*tf.T913i + -206.844*tf.T913
                  + 11.0016*tf.T9 + -0.612898*tf.T953 + 104.905*tf.lnT9)

    rate_eval.n_S32__p_P32 = rate

@numba.njit()
def n_S32__He4_Si29(rate_eval, tf):
    # S32 + n --> He4 + Si29
    rate = 0.0

    # ths8r
    rate += np.exp(  -0.947538 + 9.48125*tf.T913
                  + 0.4472*tf.T9 + -0.119237*tf.T953)

    rate_eval.n_S32__He4_Si29 = rate

@numba.njit()
def p_S32__n_Cl32(rate_eval, tf):
    # S32 + p --> n + Cl32
    rate = 0.0

    # ths8r
    rate += np.exp(  19.4852 + -156.293*tf.T9i + 1.81434*tf.T913
                  + -0.244818*tf.T9 + 0.0169659*tf.T953)

    rate_eval.p_S32__n_Cl32 = rate

@numba.njit()
def p_S32__He4_P29(rate_eval, tf):
    # S32 + p --> He4 + P29
    rate = 0.0

    # ths8r
    rate += np.exp(  52.2142 + -48.7275*tf.T9i + -62.3802*tf.T913i + 0.459085*tf.T913
                  + -0.870169*tf.T9 + 0.0631143*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_S32__He4_P29 = rate

@numba.njit()
def He4_S32__n_Ar35(rate_eval, tf):
    # S32 + He4 --> n + Ar35
    rate = 0.0

    # ths8r
    rate += np.exp(  18.3732 + -99.9692*tf.T9i + 0.502191*tf.T913
                  + -0.0974026*tf.T9 + 0.0167829*tf.T953)

    rate_eval.He4_S32__n_Ar35 = rate

@numba.njit()
def He4_S32__p_Cl35(rate_eval, tf):
    # S32 + He4 --> p + Cl35
    rate = 0.0

    # il10r
    rate += np.exp(  -57.395 + -22.1894*tf.T9i + 25.5338*tf.T913
                  + 6.45824*tf.T9 + -0.950294*tf.T953 + -1.5*tf.lnT9)
    # il10n
    rate += np.exp(  32.2544 + -21.6564*tf.T9i + -30.9147*tf.T913i + -1.2345*tf.T913
                  + 22.5118*tf.T9 + -33.0589*tf.T953 + -0.666667*tf.lnT9)
    # il10r
    rate += np.exp(  2.42563 + -27.6662*tf.T9i + 5.33756*tf.T913
                  + 1.64418*tf.T9 + -0.246167*tf.T953 + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  -0.877602 + -25.5914*tf.T9i
                  + -1.5*tf.lnT9)

    rate_eval.He4_S32__p_Cl35 = rate

@numba.njit()
def p_S34__n_Cl34(rate_eval, tf):
    # S34 + p --> n + Cl34
    rate = 0.0

    # ths8r
    rate += np.exp(  18.4445 + -72.8108*tf.T9i + 0.847427*tf.T913
                  + 0.000125282*tf.T9 + -0.00573985*tf.T953)

    rate_eval.p_S34__n_Cl34 = rate

@numba.njit()
def p_S34__He4_P31(rate_eval, tf):
    # S34 + p --> He4 + P31
    rate = 0.0

    # ths8r
    rate += np.exp(  51.0402 + -7.27426*tf.T9i + -62.5433*tf.T913i + 1.69339*tf.T913
                  + -1.27872*tf.T9 + 0.106616*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_S34__He4_P31 = rate

@numba.njit()
def n_Cl32__p_S32(rate_eval, tf):
    # Cl32 + n --> p + S32
    rate = 0.0

    # ths8r
    rate += np.exp(  18.3866 + 1.81434*tf.T913
                  + -0.244818*tf.T9 + 0.0169659*tf.T953)

    rate_eval.n_Cl32__p_S32 = rate

@numba.njit()
def n_Cl32__He4_P29(rate_eval, tf):
    # Cl32 + n --> He4 + P29
    rate = 0.0

    # ths8r
    rate += np.exp(  16.9536 + 1.0254*tf.T913
                  + -0.188653*tf.T9 + 0.0264376*tf.T953)

    rate_eval.n_Cl32__He4_P29 = rate

@numba.njit()
def He4_Cl32__p_Ar35(rate_eval, tf):
    # Cl32 + He4 --> p + Ar35
    rate = 0.0

    # ths8r
    rate += np.exp(  52.3543 + -68.0682*tf.T913i + 0.0983684*tf.T913
                  + -0.783598*tf.T9 + 0.0556457*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Cl32__p_Ar35 = rate

@numba.njit()
def n_Cl34__p_S34(rate_eval, tf):
    # Cl34 + n --> p + S34
    rate = 0.0

    # ths8r
    rate += np.exp(  18.4445 + 0.847427*tf.T913
                  + 0.000125282*tf.T9 + -0.00573985*tf.T953)

    rate_eval.n_Cl34__p_S34 = rate

@numba.njit()
def n_Cl34__He4_P31(rate_eval, tf):
    # Cl34 + n --> He4 + P31
    rate = 0.0

    # ths8r
    rate += np.exp(  18.8653 + -0.466069*tf.T913
                  + 0.167169*tf.T9 + -0.00537463*tf.T953)

    rate_eval.n_Cl34__He4_P31 = rate

@numba.njit()
def p_Cl34__He4_S31(rate_eval, tf):
    # Cl34 + p --> He4 + S31
    rate = 0.0

    # ths8r
    rate += np.exp(  52.2406 + -6.16284*tf.T9i + -65.2934*tf.T913i + 1.2687*tf.T913
                  + -1.08925*tf.T9 + 0.0898609*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Cl34__He4_S31 = rate

@numba.njit()
def n_Cl35__He4_P32(rate_eval, tf):
    # Cl35 + n --> He4 + P32
    rate = 0.0

    # ths8r
    rate += np.exp(  -32.0685 + 35.1367*tf.T913
                  + -3.30169*tf.T9 + 0.150699*tf.T953)

    rate_eval.n_Cl35__He4_P32 = rate

@numba.njit()
def p_Cl35__n_Ar35(rate_eval, tf):
    # Cl35 + p --> n + Ar35
    rate = 0.0

    # ths8r
    rate += np.exp(  19.5533 + -78.3126*tf.T9i + 0.0970999*tf.T913
                  + 0.139648*tf.T9 + -0.0149442*tf.T953)

    rate_eval.p_Cl35__n_Ar35 = rate

@numba.njit()
def p_Cl35__He4_S32(rate_eval, tf):
    # Cl35 + p --> He4 + S32
    rate = 0.0

    # il10r
    rate += np.exp(  -1.01202 + -3.93495*tf.T9i
                  + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  -57.5294 + -0.532931*tf.T9i + 25.5338*tf.T913
                  + 6.45824*tf.T9 + -0.950294*tf.T953 + -1.5*tf.lnT9)
    # il10n
    rate += np.exp(  32.12 + -30.9147*tf.T913i + -1.2345*tf.T913
                  + 22.5118*tf.T9 + -33.0589*tf.T953 + -0.666667*tf.lnT9)
    # il10r
    rate += np.exp(  2.29121 + -6.00976*tf.T9i + 5.33756*tf.T913
                  + 1.64418*tf.T9 + -0.246167*tf.T953 + -1.5*tf.lnT9)

    rate_eval.p_Cl35__He4_S32 = rate

@numba.njit()
def n_Ar35__p_Cl35(rate_eval, tf):
    # Ar35 + n --> p + Cl35
    rate = 0.0

    # ths8r
    rate += np.exp(  19.5533 + 0.0970999*tf.T913
                  + 0.139648*tf.T9 + -0.0149442*tf.T953)

    rate_eval.n_Ar35__p_Cl35 = rate

@numba.njit()
def n_Ar35__He4_S32(rate_eval, tf):
    # Ar35 + n --> He4 + S32
    rate = 0.0

    # ths8r
    rate += np.exp(  18.2388 + 0.502191*tf.T913
                  + -0.0974026*tf.T9 + 0.0167829*tf.T953)

    rate_eval.n_Ar35__He4_S32 = rate

@numba.njit()
def p_Ar35__He4_Cl32(rate_eval, tf):
    # Ar35 + p --> He4 + Cl32
    rate = 0.0

    # ths8r
    rate += np.exp(  53.3185 + -56.3235*tf.T9i + -68.0682*tf.T913i + 0.0983684*tf.T913
                  + -0.783598*tf.T9 + 0.0556457*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Ar35__He4_Cl32 = rate

@numba.njit()
def p_d__n_p_p(rate_eval, tf):
    # d + p --> n + p + p
    rate = 0.0

    # cf88n
    rate += np.exp(  17.3271 + -25.82*tf.T9i + -3.72*tf.T913i + 0.946313*tf.T913
                  + 0.105406*tf.T9 + -0.0149431*tf.T953)

    rate_eval.p_d__n_p_p = rate

@numba.njit()
def t_t__n_n_He4(rate_eval, tf):
    # t + t --> n + n + He4
    rate = 0.0

    # cf88r
    rate += np.exp(  21.2361 + -4.872*tf.T913i + -1.72398*tf.T913
                  + 0.684775*tf.T9 + -0.0702582*tf.T953 + 0.333333*tf.lnT9)
    # cf88n
    rate += np.exp(  21.2361 + -4.872*tf.T913i + -0.0328579*tf.T913
                  + -1.13588*tf.T9 + 0.250064*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.t_t__n_n_He4 = rate

@numba.njit()
def t_He3__n_p_He4(rate_eval, tf):
    # He3 + t --> n + p + He4
    rate = 0.0

    # cf88n
    rate += np.exp(  22.7658 + -7.733*tf.T913i + -0.118902*tf.T913
                  + -0.267393*tf.T9 + 0.0275387*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.t_He3__n_p_He4 = rate

@numba.njit()
def He3_He3__p_p_He4(rate_eval, tf):
    # He3 + He3 --> p + p + He4
    rate = 0.0

    # nacrn
    rate += np.exp(  24.7788 + -12.277*tf.T913i + -0.103699*tf.T913
                  + -0.0649967*tf.T9 + 0.0168191*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He3_He3__p_p_He4 = rate

@numba.njit()
def d_Li7__n_He4_He4(rate_eval, tf):
    # Li7 + d --> n + He4 + He4
    rate = 0.0

    # cf88n
    rate += np.exp(  26.4 + -10.259*tf.T913i
                  + -0.666667*tf.lnT9)

    rate_eval.d_Li7__n_He4_He4 = rate

@numba.njit()
def d_Be7__p_He4_He4(rate_eval, tf):
    # Be7 + d --> p + He4 + He4
    rate = 0.0

    # cf88n
    rate += np.exp(  27.6987 + -12.428*tf.T913i
                  + -0.666667*tf.lnT9)

    rate_eval.d_Be7__p_He4_He4 = rate

@numba.njit()
def p_Be9__d_He4_He4(rate_eval, tf):
    # Be9 + p --> d + He4 + He4
    rate = 0.0

    # cf88r
    rate += np.exp(  20.5607 + -5.8*tf.T9i
                  + -0.75*tf.lnT9)
    # cf88r
    rate += np.exp(  20.1768 + -3.046*tf.T9i
                  + -1.0*tf.lnT9)
    # cf88n
    rate += np.exp(  26.0751 + -10.359*tf.T913i + 0.103955*tf.T913
                  + 4.4262*tf.T9 + -5.95664*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Be9__d_He4_He4 = rate

@numba.njit()
def n_B8__p_He4_He4(rate_eval, tf):
    # B8 + n --> p + He4 + He4
    rate = 0.0

    #  wagn
    rate += np.exp(  19.812 + -1.13869e-12*tf.T9i + 1.5368e-10*tf.T913i + -5.36043e-10*tf.T913
                  + 5.80628e-11*tf.T9 + -5.23778e-12*tf.T953 + 1.77211e-10*tf.lnT9)

    rate_eval.n_B8__p_He4_He4 = rate

@numba.njit()
def n_C11__He4_He4_He4(rate_eval, tf):
    # C11 + n --> 3 He4
    rate = 0.0

    # bb92r
    rate += np.exp(  22.7448 + -5.58*tf.T9i
                  + -1.5*tf.lnT9)
    # bb92n
    rate += np.exp(  18.9275)
    # bb92r
    rate += np.exp(  14.7197 + -0.917*tf.T9i
                  + -1.5*tf.lnT9)

    rate_eval.n_C11__He4_He4_He4 = rate

@numba.njit()
def t_Li7__n_n_He4_He4(rate_eval, tf):
    # Li7 + t --> n + n + He4 + He4
    rate = 0.0

    # mafon
    rate += np.exp(  27.5043 + -5.31692e-12*tf.T9i + -11.333*tf.T913i + -2.24192e-09*tf.T913
                  + 2.21773e-10*tf.T9 + -1.83941e-11*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.t_Li7__n_n_He4_He4 = rate

@numba.njit()
def He3_Li7__n_p_He4_He4(rate_eval, tf):
    # Li7 + He3 --> n + p + He4 + He4
    rate = 0.0

    # mafon
    rate += np.exp(  30.038 + -4.24733e-12*tf.T9i + -17.989*tf.T913i + -1.57523e-09*tf.T913
                  + 1.45934e-10*tf.T9 + -1.15341e-11*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He3_Li7__n_p_He4_He4 = rate

@numba.njit()
def t_Be7__n_p_He4_He4(rate_eval, tf):
    # Be7 + t --> n + p + He4 + He4
    rate = 0.0

    # mafon
    rate += np.exp(  28.6992 + -6.9004e-12*tf.T9i + -13.792*tf.T913i + -2.92021e-09*tf.T913
                  + 2.89378e-10*tf.T9 + -2.40287e-11*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.t_Be7__n_p_He4_He4 = rate

@numba.njit()
def He3_Be7__p_p_He4_He4(rate_eval, tf):
    # Be7 + He3 --> p + p + He4 + He4
    rate = 0.0

    # mafon
    rate += np.exp(  31.7435 + -5.45213e-12*tf.T9i + -21.793*tf.T913i + -1.98126e-09*tf.T913
                  + 1.84204e-10*tf.T9 + -1.46403e-11*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He3_Be7__p_p_He4_He4 = rate

@numba.njit()
def p_Be9__n_p_He4_He4(rate_eval, tf):
    # Be9 + p --> n + p + He4 + He4
    rate = 0.0

    # cf88r
    rate += np.exp(  20.7431 + -26.725*tf.T9i + 1.40505e-06*tf.T913i + -1.47128e-06*tf.T913
                  + 6.89313e-08*tf.T9 + -3.55179e-09*tf.T953 + -1.5*tf.lnT9)
    # cf88n
    rate += np.exp(  14.6035 + -21.4742*tf.T9i + -0.634849*tf.T913i + 4.82033*tf.T913
                  + -0.257317*tf.T9 + 0.0134206*tf.T953 + -1.08885*tf.lnT9)

    rate_eval.p_Be9__n_p_He4_He4 = rate

@numba.njit()
def n_n_He4__He6(rate_eval, tf):
    # n + n + He4 --> He6
    rate = 0.0

    # cf88r
    rate += np.exp(  -23.9322 + -9.585*tf.T9i + 0.694279*tf.T913i + -3.33326*tf.T913
                  + 0.507932*tf.T9 + -0.0427342*tf.T953 + -1.0*tf.lnT9)

    rate_eval.n_n_He4__He6 = rate

@numba.njit()
def n_p_He4__Li6(rate_eval, tf):
    # n + p + He4 --> Li6
    rate = 0.0

    # cf88r
    rate += np.exp(  -12.2851 + -19.353*tf.T9i + 1.44987*tf.T913i + -1.42759*tf.T913
                  + 0.0454035*tf.T9 + 0.00471161*tf.T953 + -1.0*tf.lnT9)

    rate_eval.n_p_He4__Li6 = rate

@numba.njit()
def n_He4_He4__Be9(rate_eval, tf):
    # n + He4 + He4 --> Be9
    rate = 0.0

    # ac12r
    rate += np.exp(  -6.81178 + -1.01953*tf.T9i + -1.56673*tf.T913i + -5.43497*tf.T913
                  + 0.673807*tf.T9 + -0.041014*tf.T953 + -1.5*tf.lnT9)
    # ac12n
    rate += np.exp(  -8.22898 + -13.3317*tf.T913i + 13.2237*tf.T913
                  + -9.06339*tf.T9 + -0.666667*tf.lnT9)

    rate_eval.n_He4_He4__Be9 = rate

@numba.njit()
def He4_He4_He4__C12(rate_eval, tf):
    # 3 He4 --> C12
    rate = 0.0

    # fy05r
    rate += np.exp(  -11.7884 + -1.02446*tf.T9i + -23.57*tf.T913i + 20.4886*tf.T913
                  + -12.9882*tf.T9 + -20.0*tf.T953 + -2.16667*tf.lnT9)
    # fy05n
    rate += np.exp(  -0.971052 + -37.06*tf.T913i + 29.3493*tf.T913
                  + -115.507*tf.T9 + -10.0*tf.T953 + -1.33333*tf.lnT9)
    # fy05r
    rate += np.exp(  -24.3505 + -4.12656*tf.T9i + -13.49*tf.T913i + 21.4259*tf.T913
                  + -1.34769*tf.T9 + 0.0879816*tf.T953 + -13.1653*tf.lnT9)

    rate_eval.He4_He4_He4__C12 = rate

@numba.njit()
def n_p_p__p_d(rate_eval, tf):
    # n + p + p --> p + d
    rate = 0.0

    # cf88n
    rate += np.exp(  -4.24034 + -3.72*tf.T913i + 0.946313*tf.T913
                  + 0.105406*tf.T9 + -0.0149431*tf.T953 + -1.5*tf.lnT9)

    rate_eval.n_p_p__p_d = rate

@numba.njit()
def n_n_He4__t_t(rate_eval, tf):
    # n + n + He4 --> t + t
    rate = 0.0

    # cf88n
    rate += np.exp(  -0.560128 + -131.502*tf.T9i + -4.872*tf.T913i + -0.0328579*tf.T913
                  + -1.13588*tf.T9 + 0.250064*tf.T953 + -2.16667*tf.lnT9)
    # cf88r
    rate += np.exp(  -0.560128 + -131.502*tf.T9i + -4.872*tf.T913i + -1.72398*tf.T913
                  + 0.684775*tf.T9 + -0.0702582*tf.T953 + -1.16667*tf.lnT9)

    rate_eval.n_n_He4__t_t = rate

@numba.njit()
def n_p_He4__t_He3(rate_eval, tf):
    # n + p + He4 --> t + He3
    rate = 0.0

    # cf88n
    rate += np.exp(  0.969572 + -140.368*tf.T9i + -7.733*tf.T913i + -0.118902*tf.T913
                  + -0.267393*tf.T9 + 0.0275387*tf.T953 + -2.16667*tf.lnT9)

    rate_eval.n_p_He4__t_He3 = rate

@numba.njit()
def p_p_He4__He3_He3(rate_eval, tf):
    # p + p + He4 --> He3 + He3
    rate = 0.0

    # nacrn
    rate += np.exp(  2.98257 + -149.222*tf.T9i + -12.277*tf.T913i + -0.103699*tf.T913
                  + -0.0649967*tf.T9 + 0.0168191*tf.T953 + -2.16667*tf.lnT9)

    rate_eval.p_p_He4__He3_He3 = rate

@numba.njit()
def n_He4_He4__d_Li7(rate_eval, tf):
    # n + He4 + He4 --> d + Li7
    rate = 0.0

    # cf88n
    rate += np.exp(  5.67199 + -175.472*tf.T9i + -10.259*tf.T913i
                  + -2.16667*tf.lnT9)

    rate_eval.n_He4_He4__d_Li7 = rate

@numba.njit()
def p_He4_He4__n_B8(rate_eval, tf):
    # p + He4 + He4 --> n + B8
    rate = 0.0

    #  wagn
    rate += np.exp(  -1.93853 + -218.783*tf.T9i + 1.5368e-10*tf.T913i + -5.36043e-10*tf.T913
                  + 5.80628e-11*tf.T9 + -5.23778e-12*tf.T953 + -1.5*tf.lnT9)

    rate_eval.p_He4_He4__n_B8 = rate

@numba.njit()
def p_He4_He4__d_Be7(rate_eval, tf):
    # p + He4 + He4 --> d + Be7
    rate = 0.0

    # cf88n
    rate += np.exp(  6.97069 + -194.561*tf.T9i + -12.428*tf.T913i
                  + -2.16667*tf.lnT9)

    rate_eval.p_He4_He4__d_Be7 = rate

@numba.njit()
def d_He4_He4__p_Be9(rate_eval, tf):
    # d + He4 + He4 --> p + Be9
    rate = 0.0

    # cf88n
    rate += np.exp(  2.83369 + -7.55453*tf.T9i + -10.359*tf.T913i + 0.103955*tf.T913
                  + 4.4262*tf.T9 + -5.95664*tf.T953 + -2.16667*tf.lnT9)
    # cf88r
    rate += np.exp(  -2.68071 + -13.3545*tf.T9i
                  + -2.25*tf.lnT9)
    # cf88r
    rate += np.exp(  -3.06461 + -10.6005*tf.T9i
                  + -2.5*tf.lnT9)

    rate_eval.d_He4_He4__p_Be9 = rate

@numba.njit()
def He4_He4_He4__n_C11(rate_eval, tf):
    # 3 He4 --> n + C11
    rate = 0.0

    # bb92r
    rate += np.exp(  0.961896 + -138.412*tf.T9i
                  + -3.0*tf.lnT9)
    # bb92n
    rate += np.exp(  -2.8554 + -132.832*tf.T9i
                  + -1.5*tf.lnT9)
    # bb92r
    rate += np.exp(  -7.0632 + -133.749*tf.T9i
                  + -3.0*tf.lnT9)

    rate_eval.He4_He4_He4__n_C11 = rate

@numba.njit()
def n_n_He4_He4__t_Li7(rate_eval, tf):
    # n + n + He4 + He4 --> t + Li7
    rate = 0.0

    # mafon
    rate += np.exp(  -17.4199 + -102.867*tf.T9i + -11.333*tf.T913i + -2.24192e-09*tf.T913
                  + 2.21773e-10*tf.T9 + -1.83941e-11*tf.T953 + -3.66667*tf.lnT9)

    rate_eval.n_n_He4_He4__t_Li7 = rate

@numba.njit()
def n_p_He4_He4__He3_Li7(rate_eval, tf):
    # n + p + He4 + He4 --> He3 + Li7
    rate = 0.0

    # mafon
    rate += np.exp(  -14.8862 + -111.725*tf.T9i + -17.989*tf.T913i + -1.57523e-09*tf.T913
                  + 1.45934e-10*tf.T9 + -1.15341e-11*tf.T953 + -3.66667*tf.lnT9)

    rate_eval.n_p_He4_He4__He3_Li7 = rate

@numba.njit()
def n_p_He4_He4__t_Be7(rate_eval, tf):
    # n + p + He4 + He4 --> t + Be7
    rate = 0.0

    # mafon
    rate += np.exp(  -16.225 + -121.949*tf.T9i + -13.792*tf.T913i + -2.92021e-09*tf.T913
                  + 2.89378e-10*tf.T9 + -2.40287e-11*tf.T953 + -3.66667*tf.lnT9)

    rate_eval.n_p_He4_He4__t_Be7 = rate

@numba.njit()
def n_p_He4_He4__p_Be9(rate_eval, tf):
    # n + p + He4 + He4 --> p + Be9
    rate = 0.0

    # cf88r
    rate += np.exp(  -25.452 + -8.47112*tf.T9i + 1.40505e-06*tf.T913i + -1.47128e-06*tf.T913
                  + 6.89313e-08*tf.T9 + -3.55179e-09*tf.T953 + -4.5*tf.lnT9)
    # cf88n
    rate += np.exp(  -31.5916 + -3.22032*tf.T9i + -0.634849*tf.T913i + 4.82033*tf.T913
                  + -0.257317*tf.T9 + 0.0134206*tf.T953 + -4.08885*tf.lnT9)

    rate_eval.n_p_He4_He4__p_Be9 = rate

@numba.njit()
def p_p_He4_He4__He3_Be7(rate_eval, tf):
    # p + p + He4 + He4 --> He3 + Be7
    rate = 0.0

    # mafon
    rate += np.exp(  -13.1807 + -130.807*tf.T9i + -21.793*tf.T913i + -1.98126e-09*tf.T913
                  + 1.84204e-10*tf.T9 + -1.46403e-11*tf.T953 + -3.66667*tf.lnT9)

    rate_eval.p_p_He4_He4__He3_Be7 = rate

def rhs(t, Y, rho, T, screen_func=None):
    return rhs_eq(t, Y, rho, T, screen_func)

@numba.njit()
def rhs_eq(t, Y, rho, T, screen_func):

    tf = Tfactors(T)
    rate_eval = RateEval()

    # reaclib rates
    n__p__weak__wc12(rate_eval, tf)
    t__He3__weak__wc12(rate_eval, tf)
    He3__t__weak__electron_capture(rate_eval, tf)
    He6__Li6__weak__wc12(rate_eval, tf)
    Be7__Li7__weak__electron_capture(rate_eval, tf)
    C14__N14__weak__wc12(rate_eval, tf)
    N13__C13__weak__wc12(rate_eval, tf)
    O14__N14__weak__wc12(rate_eval, tf)
    O15__N15__weak__wc12(rate_eval, tf)
    F17__O17__weak__wc12(rate_eval, tf)
    F18__O18__weak__wc12(rate_eval, tf)
    Ne18__F18__weak__wc12(rate_eval, tf)
    Ne19__F19__weak__wc12(rate_eval, tf)
    Na21__Ne21__weak__wc12(rate_eval, tf)
    Na24__Mg24__weak__wc12(rate_eval, tf)
    Mg22__Na22__weak__wc12(rate_eval, tf)
    Mg23__Na23__weak__wc12(rate_eval, tf)
    Al25__Mg25__weak__wc12(rate_eval, tf)
    Al26__Mg26__weak__wc12(rate_eval, tf)
    Al28__Si28__weak__wc12(rate_eval, tf)
    Si26__Al26__weak__wc12(rate_eval, tf)
    P29__Si29__weak__wc12(rate_eval, tf)
    P30__Si30__weak__wc12(rate_eval, tf)
    P32__S32__weak__wc12(rate_eval, tf)
    S31__P31__weak__wc12(rate_eval, tf)
    Cl32__S32__weak__wc12(rate_eval, tf)
    Cl34__S34__weak__wc12(rate_eval, tf)
    Ar35__Cl35__weak__wc12(rate_eval, tf)
    d__n_p(rate_eval, tf)
    t__n_d(rate_eval, tf)
    He3__p_d(rate_eval, tf)
    He4__n_He3(rate_eval, tf)
    He4__p_t(rate_eval, tf)
    He4__d_d(rate_eval, tf)
    Li6__He4_d(rate_eval, tf)
    Li7__n_Li6(rate_eval, tf)
    Li7__He4_t(rate_eval, tf)
    Be7__p_Li6(rate_eval, tf)
    Be7__He4_He3(rate_eval, tf)
    B8__p_Be7(rate_eval, tf)
    B8__He4_He4__weak__wc12(rate_eval, tf)
    C11__He4_Be7(rate_eval, tf)
    C12__n_C11(rate_eval, tf)
    C13__n_C12(rate_eval, tf)
    C14__n_C13(rate_eval, tf)
    N13__p_C12(rate_eval, tf)
    N14__n_N13(rate_eval, tf)
    N14__p_C13(rate_eval, tf)
    N15__n_N14(rate_eval, tf)
    N15__p_C14(rate_eval, tf)
    O14__p_N13(rate_eval, tf)
    O15__n_O14(rate_eval, tf)
    O15__p_N14(rate_eval, tf)
    O16__n_O15(rate_eval, tf)
    O16__p_N15(rate_eval, tf)
    O16__He4_C12(rate_eval, tf)
    O17__n_O16(rate_eval, tf)
    O18__n_O17(rate_eval, tf)
    O18__He4_C14(rate_eval, tf)
    F15__p_O14(rate_eval, tf)
    F16__n_F15(rate_eval, tf)
    F16__p_O15(rate_eval, tf)
    F17__n_F16(rate_eval, tf)
    F17__p_O16(rate_eval, tf)
    F18__n_F17(rate_eval, tf)
    F18__p_O17(rate_eval, tf)
    F18__He4_N14(rate_eval, tf)
    F19__n_F18(rate_eval, tf)
    F19__p_O18(rate_eval, tf)
    F19__He4_N15(rate_eval, tf)
    Ne18__p_F17(rate_eval, tf)
    Ne18__He4_O14(rate_eval, tf)
    Ne19__n_Ne18(rate_eval, tf)
    Ne19__p_F18(rate_eval, tf)
    Ne19__He4_O15(rate_eval, tf)
    Ne20__n_Ne19(rate_eval, tf)
    Ne20__p_F19(rate_eval, tf)
    Ne20__He4_O16(rate_eval, tf)
    Ne21__n_Ne20(rate_eval, tf)
    Ne21__He4_O17(rate_eval, tf)
    Na21__p_Ne20(rate_eval, tf)
    Na21__He4_F17(rate_eval, tf)
    Na22__n_Na21(rate_eval, tf)
    Na22__p_Ne21(rate_eval, tf)
    Na22__He4_F18(rate_eval, tf)
    Na23__n_Na22(rate_eval, tf)
    Na23__He4_F19(rate_eval, tf)
    Na24__n_Na23(rate_eval, tf)
    Mg22__p_Na21(rate_eval, tf)
    Mg22__He4_Ne18(rate_eval, tf)
    Mg23__n_Mg22(rate_eval, tf)
    Mg23__p_Na22(rate_eval, tf)
    Mg23__He4_Ne19(rate_eval, tf)
    Mg24__n_Mg23(rate_eval, tf)
    Mg24__p_Na23(rate_eval, tf)
    Mg24__He4_Ne20(rate_eval, tf)
    Mg25__n_Mg24(rate_eval, tf)
    Mg25__p_Na24(rate_eval, tf)
    Mg25__He4_Ne21(rate_eval, tf)
    Mg26__n_Mg25(rate_eval, tf)
    Al25__p_Mg24(rate_eval, tf)
    Al25__He4_Na21(rate_eval, tf)
    Al26__n_Al25(rate_eval, tf)
    Al26__p_Mg25(rate_eval, tf)
    Al26__He4_Na22(rate_eval, tf)
    Al27__n_Al26(rate_eval, tf)
    Al27__p_Mg26(rate_eval, tf)
    Al27__He4_Na23(rate_eval, tf)
    Al28__n_Al27(rate_eval, tf)
    Al28__He4_Na24(rate_eval, tf)
    Si26__p_Al25(rate_eval, tf)
    Si26__He4_Mg22(rate_eval, tf)
    Si28__p_Al27(rate_eval, tf)
    Si28__He4_Mg24(rate_eval, tf)
    Si29__n_Si28(rate_eval, tf)
    Si29__p_Al28(rate_eval, tf)
    Si29__He4_Mg25(rate_eval, tf)
    Si30__n_Si29(rate_eval, tf)
    Si30__He4_Mg26(rate_eval, tf)
    P29__p_Si28(rate_eval, tf)
    P29__He4_Al25(rate_eval, tf)
    P30__n_P29(rate_eval, tf)
    P30__p_Si29(rate_eval, tf)
    P30__He4_Al26(rate_eval, tf)
    P31__n_P30(rate_eval, tf)
    P31__p_Si30(rate_eval, tf)
    P31__He4_Al27(rate_eval, tf)
    P32__n_P31(rate_eval, tf)
    P32__He4_Al28(rate_eval, tf)
    S31__p_P30(rate_eval, tf)
    S32__n_S31(rate_eval, tf)
    S32__p_P31(rate_eval, tf)
    S32__He4_Si28(rate_eval, tf)
    S34__He4_Si30(rate_eval, tf)
    Cl32__p_S31(rate_eval, tf)
    Cl32__p_P31__weak__wc12(rate_eval, tf)
    Cl32__He4_Si28__weak__wc12(rate_eval, tf)
    Cl34__He4_P30(rate_eval, tf)
    Cl35__n_Cl34(rate_eval, tf)
    Cl35__p_S34(rate_eval, tf)
    Cl35__He4_P31(rate_eval, tf)
    Ar35__p_Cl34(rate_eval, tf)
    Ar35__He4_S31(rate_eval, tf)
    He6__n_n_He4(rate_eval, tf)
    Li6__n_p_He4(rate_eval, tf)
    Be9__n_He4_He4(rate_eval, tf)
    C12__He4_He4_He4(rate_eval, tf)
    n_p__d(rate_eval, tf)
    p_p__d__weak__bet_pos_(rate_eval, tf)
    p_p__d__weak__electron_capture(rate_eval, tf)
    n_d__t(rate_eval, tf)
    p_d__He3(rate_eval, tf)
    d_d__He4(rate_eval, tf)
    He4_d__Li6(rate_eval, tf)
    p_t__He4(rate_eval, tf)
    He4_t__Li7(rate_eval, tf)
    n_He3__He4(rate_eval, tf)
    p_He3__He4__weak__bet_pos_(rate_eval, tf)
    He4_He3__Be7(rate_eval, tf)
    n_Li6__Li7(rate_eval, tf)
    p_Li6__Be7(rate_eval, tf)
    p_Be7__B8(rate_eval, tf)
    He4_Be7__C11(rate_eval, tf)
    n_C11__C12(rate_eval, tf)
    n_C12__C13(rate_eval, tf)
    p_C12__N13(rate_eval, tf)
    He4_C12__O16(rate_eval, tf)
    n_C13__C14(rate_eval, tf)
    p_C13__N14(rate_eval, tf)
    p_C14__N15(rate_eval, tf)
    He4_C14__O18(rate_eval, tf)
    n_N13__N14(rate_eval, tf)
    p_N13__O14(rate_eval, tf)
    n_N14__N15(rate_eval, tf)
    p_N14__O15(rate_eval, tf)
    He4_N14__F18(rate_eval, tf)
    p_N15__O16(rate_eval, tf)
    He4_N15__F19(rate_eval, tf)
    n_O14__O15(rate_eval, tf)
    p_O14__F15(rate_eval, tf)
    He4_O14__Ne18(rate_eval, tf)
    n_O15__O16(rate_eval, tf)
    p_O15__F16(rate_eval, tf)
    He4_O15__Ne19(rate_eval, tf)
    n_O16__O17(rate_eval, tf)
    p_O16__F17(rate_eval, tf)
    He4_O16__Ne20(rate_eval, tf)
    n_O17__O18(rate_eval, tf)
    p_O17__F18(rate_eval, tf)
    He4_O17__Ne21(rate_eval, tf)
    p_O18__F19(rate_eval, tf)
    n_F15__F16(rate_eval, tf)
    n_F16__F17(rate_eval, tf)
    n_F17__F18(rate_eval, tf)
    p_F17__Ne18(rate_eval, tf)
    He4_F17__Na21(rate_eval, tf)
    n_F18__F19(rate_eval, tf)
    p_F18__Ne19(rate_eval, tf)
    He4_F18__Na22(rate_eval, tf)
    p_F19__Ne20(rate_eval, tf)
    He4_F19__Na23(rate_eval, tf)
    n_Ne18__Ne19(rate_eval, tf)
    He4_Ne18__Mg22(rate_eval, tf)
    n_Ne19__Ne20(rate_eval, tf)
    He4_Ne19__Mg23(rate_eval, tf)
    n_Ne20__Ne21(rate_eval, tf)
    p_Ne20__Na21(rate_eval, tf)
    He4_Ne20__Mg24(rate_eval, tf)
    p_Ne21__Na22(rate_eval, tf)
    He4_Ne21__Mg25(rate_eval, tf)
    n_Na21__Na22(rate_eval, tf)
    p_Na21__Mg22(rate_eval, tf)
    He4_Na21__Al25(rate_eval, tf)
    n_Na22__Na23(rate_eval, tf)
    p_Na22__Mg23(rate_eval, tf)
    He4_Na22__Al26(rate_eval, tf)
    n_Na23__Na24(rate_eval, tf)
    p_Na23__Mg24(rate_eval, tf)
    He4_Na23__Al27(rate_eval, tf)
    p_Na24__Mg25(rate_eval, tf)
    He4_Na24__Al28(rate_eval, tf)
    n_Mg22__Mg23(rate_eval, tf)
    He4_Mg22__Si26(rate_eval, tf)
    n_Mg23__Mg24(rate_eval, tf)
    n_Mg24__Mg25(rate_eval, tf)
    p_Mg24__Al25(rate_eval, tf)
    He4_Mg24__Si28(rate_eval, tf)
    n_Mg25__Mg26(rate_eval, tf)
    p_Mg25__Al26(rate_eval, tf)
    He4_Mg25__Si29(rate_eval, tf)
    p_Mg26__Al27(rate_eval, tf)
    He4_Mg26__Si30(rate_eval, tf)
    n_Al25__Al26(rate_eval, tf)
    p_Al25__Si26(rate_eval, tf)
    He4_Al25__P29(rate_eval, tf)
    n_Al26__Al27(rate_eval, tf)
    He4_Al26__P30(rate_eval, tf)
    n_Al27__Al28(rate_eval, tf)
    p_Al27__Si28(rate_eval, tf)
    He4_Al27__P31(rate_eval, tf)
    p_Al28__Si29(rate_eval, tf)
    He4_Al28__P32(rate_eval, tf)
    n_Si28__Si29(rate_eval, tf)
    p_Si28__P29(rate_eval, tf)
    He4_Si28__S32(rate_eval, tf)
    n_Si29__Si30(rate_eval, tf)
    p_Si29__P30(rate_eval, tf)
    p_Si30__P31(rate_eval, tf)
    He4_Si30__S34(rate_eval, tf)
    n_P29__P30(rate_eval, tf)
    n_P30__P31(rate_eval, tf)
    p_P30__S31(rate_eval, tf)
    He4_P30__Cl34(rate_eval, tf)
    n_P31__P32(rate_eval, tf)
    p_P31__S32(rate_eval, tf)
    He4_P31__Cl35(rate_eval, tf)
    n_S31__S32(rate_eval, tf)
    p_S31__Cl32(rate_eval, tf)
    He4_S31__Ar35(rate_eval, tf)
    p_S34__Cl35(rate_eval, tf)
    n_Cl34__Cl35(rate_eval, tf)
    p_Cl34__Ar35(rate_eval, tf)
    d_d__n_He3(rate_eval, tf)
    d_d__p_t(rate_eval, tf)
    p_t__n_He3(rate_eval, tf)
    p_t__d_d(rate_eval, tf)
    d_t__n_He4(rate_eval, tf)
    He4_t__n_Li6(rate_eval, tf)
    n_He3__p_t(rate_eval, tf)
    n_He3__d_d(rate_eval, tf)
    d_He3__p_He4(rate_eval, tf)
    t_He3__d_He4(rate_eval, tf)
    He4_He3__p_Li6(rate_eval, tf)
    n_He4__d_t(rate_eval, tf)
    p_He4__d_He3(rate_eval, tf)
    d_He4__t_He3(rate_eval, tf)
    He4_He4__n_Be7(rate_eval, tf)
    He4_He4__p_Li7(rate_eval, tf)
    n_Li6__He4_t(rate_eval, tf)
    p_Li6__He4_He3(rate_eval, tf)
    d_Li6__n_Be7(rate_eval, tf)
    d_Li6__p_Li7(rate_eval, tf)
    He4_Li6__p_Be9(rate_eval, tf)
    p_Li7__n_Be7(rate_eval, tf)
    p_Li7__d_Li6(rate_eval, tf)
    p_Li7__He4_He4(rate_eval, tf)
    t_Li7__n_Be9(rate_eval, tf)
    n_Be7__p_Li7(rate_eval, tf)
    n_Be7__d_Li6(rate_eval, tf)
    n_Be7__He4_He4(rate_eval, tf)
    n_Be9__t_Li7(rate_eval, tf)
    p_Be9__He4_Li6(rate_eval, tf)
    He4_Be9__n_C12(rate_eval, tf)
    He4_B8__p_C11(rate_eval, tf)
    p_C11__He4_B8(rate_eval, tf)
    He4_C11__n_O14(rate_eval, tf)
    He4_C11__p_N14(rate_eval, tf)
    n_C12__He4_Be9(rate_eval, tf)
    He4_C12__n_O15(rate_eval, tf)
    He4_C12__p_N15(rate_eval, tf)
    C12_C12__n_Mg23(rate_eval, tf)
    C12_C12__p_Na23(rate_eval, tf)
    C12_C12__He4_Ne20(rate_eval, tf)
    p_C13__n_N13(rate_eval, tf)
    d_C13__n_N14(rate_eval, tf)
    He4_C13__n_O16(rate_eval, tf)
    p_C14__n_N14(rate_eval, tf)
    d_C14__n_N15(rate_eval, tf)
    He4_C14__n_O17(rate_eval, tf)
    n_N13__p_C13(rate_eval, tf)
    He4_N13__n_F16(rate_eval, tf)
    He4_N13__p_O16(rate_eval, tf)
    n_N14__p_C14(rate_eval, tf)
    n_N14__d_C13(rate_eval, tf)
    p_N14__n_O14(rate_eval, tf)
    p_N14__He4_C11(rate_eval, tf)
    He4_N14__n_F17(rate_eval, tf)
    He4_N14__p_O17(rate_eval, tf)
    n_N15__d_C14(rate_eval, tf)
    p_N15__n_O15(rate_eval, tf)
    p_N15__He4_C12(rate_eval, tf)
    He4_N15__n_F18(rate_eval, tf)
    He4_N15__p_O18(rate_eval, tf)
    n_O14__p_N14(rate_eval, tf)
    n_O14__He4_C11(rate_eval, tf)
    He4_O14__p_F17(rate_eval, tf)
    n_O15__p_N15(rate_eval, tf)
    n_O15__He4_C12(rate_eval, tf)
    p_O15__n_F15(rate_eval, tf)
    He4_O15__n_Ne18(rate_eval, tf)
    He4_O15__p_F18(rate_eval, tf)
    n_O16__He4_C13(rate_eval, tf)
    p_O16__n_F16(rate_eval, tf)
    p_O16__He4_N13(rate_eval, tf)
    He4_O16__n_Ne19(rate_eval, tf)
    He4_O16__p_F19(rate_eval, tf)
    C12_O16__p_Al27(rate_eval, tf)
    C12_O16__He4_Mg24(rate_eval, tf)
    O16_O16__n_S31(rate_eval, tf)
    O16_O16__p_P31(rate_eval, tf)
    O16_O16__He4_Si28(rate_eval, tf)
    n_O17__He4_C14(rate_eval, tf)
    p_O17__n_F17(rate_eval, tf)
    p_O17__He4_N14(rate_eval, tf)
    He4_O17__n_Ne20(rate_eval, tf)
    p_O18__n_F18(rate_eval, tf)
    p_O18__He4_N15(rate_eval, tf)
    He4_O18__n_Ne21(rate_eval, tf)
    n_F15__p_O15(rate_eval, tf)
    He4_F15__p_Ne18(rate_eval, tf)
    n_F16__p_O16(rate_eval, tf)
    n_F16__He4_N13(rate_eval, tf)
    He4_F16__p_Ne19(rate_eval, tf)
    n_F17__p_O17(rate_eval, tf)
    n_F17__He4_N14(rate_eval, tf)
    p_F17__He4_O14(rate_eval, tf)
    He4_F17__p_Ne20(rate_eval, tf)
    n_F18__p_O18(rate_eval, tf)
    n_F18__He4_N15(rate_eval, tf)
    p_F18__n_Ne18(rate_eval, tf)
    p_F18__He4_O15(rate_eval, tf)
    He4_F18__n_Na21(rate_eval, tf)
    He4_F18__p_Ne21(rate_eval, tf)
    p_F19__n_Ne19(rate_eval, tf)
    p_F19__He4_O16(rate_eval, tf)
    He4_F19__n_Na22(rate_eval, tf)
    n_Ne18__p_F18(rate_eval, tf)
    n_Ne18__He4_O15(rate_eval, tf)
    p_Ne18__He4_F15(rate_eval, tf)
    He4_Ne18__p_Na21(rate_eval, tf)
    n_Ne19__p_F19(rate_eval, tf)
    n_Ne19__He4_O16(rate_eval, tf)
    p_Ne19__He4_F16(rate_eval, tf)
    He4_Ne19__n_Mg22(rate_eval, tf)
    He4_Ne19__p_Na22(rate_eval, tf)
    n_Ne20__He4_O17(rate_eval, tf)
    p_Ne20__He4_F17(rate_eval, tf)
    He4_Ne20__n_Mg23(rate_eval, tf)
    He4_Ne20__p_Na23(rate_eval, tf)
    He4_Ne20__C12_C12(rate_eval, tf)
    C12_Ne20__n_S31(rate_eval, tf)
    C12_Ne20__p_P31(rate_eval, tf)
    C12_Ne20__He4_Si28(rate_eval, tf)
    n_Ne21__He4_O18(rate_eval, tf)
    p_Ne21__n_Na21(rate_eval, tf)
    p_Ne21__He4_F18(rate_eval, tf)
    He4_Ne21__n_Mg24(rate_eval, tf)
    He4_Ne21__p_Na24(rate_eval, tf)
    n_Na21__p_Ne21(rate_eval, tf)
    n_Na21__He4_F18(rate_eval, tf)
    p_Na21__He4_Ne18(rate_eval, tf)
    He4_Na21__p_Mg24(rate_eval, tf)
    n_Na22__He4_F19(rate_eval, tf)
    p_Na22__n_Mg22(rate_eval, tf)
    p_Na22__He4_Ne19(rate_eval, tf)
    He4_Na22__n_Al25(rate_eval, tf)
    He4_Na22__p_Mg25(rate_eval, tf)
    p_Na23__n_Mg23(rate_eval, tf)
    p_Na23__He4_Ne20(rate_eval, tf)
    p_Na23__C12_C12(rate_eval, tf)
    He4_Na23__n_Al26(rate_eval, tf)
    He4_Na23__p_Mg26(rate_eval, tf)
    p_Na24__n_Mg24(rate_eval, tf)
    p_Na24__He4_Ne21(rate_eval, tf)
    He4_Na24__n_Al27(rate_eval, tf)
    n_Mg22__p_Na22(rate_eval, tf)
    n_Mg22__He4_Ne19(rate_eval, tf)
    He4_Mg22__p_Al25(rate_eval, tf)
    n_Mg23__p_Na23(rate_eval, tf)
    n_Mg23__He4_Ne20(rate_eval, tf)
    n_Mg23__C12_C12(rate_eval, tf)
    He4_Mg23__n_Si26(rate_eval, tf)
    He4_Mg23__p_Al26(rate_eval, tf)
    n_Mg24__p_Na24(rate_eval, tf)
    n_Mg24__He4_Ne21(rate_eval, tf)
    p_Mg24__He4_Na21(rate_eval, tf)
    He4_Mg24__p_Al27(rate_eval, tf)
    He4_Mg24__C12_O16(rate_eval, tf)
    p_Mg25__n_Al25(rate_eval, tf)
    p_Mg25__He4_Na22(rate_eval, tf)
    He4_Mg25__n_Si28(rate_eval, tf)
    He4_Mg25__p_Al28(rate_eval, tf)
    p_Mg26__n_Al26(rate_eval, tf)
    p_Mg26__He4_Na23(rate_eval, tf)
    He4_Mg26__n_Si29(rate_eval, tf)
    n_Al25__p_Mg25(rate_eval, tf)
    n_Al25__He4_Na22(rate_eval, tf)
    p_Al25__He4_Mg22(rate_eval, tf)
    He4_Al25__p_Si28(rate_eval, tf)
    n_Al26__p_Mg26(rate_eval, tf)
    n_Al26__He4_Na23(rate_eval, tf)
    p_Al26__n_Si26(rate_eval, tf)
    p_Al26__He4_Mg23(rate_eval, tf)
    He4_Al26__n_P29(rate_eval, tf)
    He4_Al26__p_Si29(rate_eval, tf)
    n_Al27__He4_Na24(rate_eval, tf)
    p_Al27__He4_Mg24(rate_eval, tf)
    p_Al27__C12_O16(rate_eval, tf)
    He4_Al27__n_P30(rate_eval, tf)
    He4_Al27__p_Si30(rate_eval, tf)
    p_Al28__n_Si28(rate_eval, tf)
    p_Al28__He4_Mg25(rate_eval, tf)
    He4_Al28__n_P31(rate_eval, tf)
    n_Si26__p_Al26(rate_eval, tf)
    n_Si26__He4_Mg23(rate_eval, tf)
    He4_Si26__p_P29(rate_eval, tf)
    n_Si28__p_Al28(rate_eval, tf)
    n_Si28__He4_Mg25(rate_eval, tf)
    p_Si28__He4_Al25(rate_eval, tf)
    He4_Si28__n_S31(rate_eval, tf)
    He4_Si28__p_P31(rate_eval, tf)
    He4_Si28__C12_Ne20(rate_eval, tf)
    He4_Si28__O16_O16(rate_eval, tf)
    n_Si29__He4_Mg26(rate_eval, tf)
    p_Si29__n_P29(rate_eval, tf)
    p_Si29__He4_Al26(rate_eval, tf)
    He4_Si29__n_S32(rate_eval, tf)
    He4_Si29__p_P32(rate_eval, tf)
    p_Si30__n_P30(rate_eval, tf)
    p_Si30__He4_Al27(rate_eval, tf)
    n_P29__p_Si29(rate_eval, tf)
    n_P29__He4_Al26(rate_eval, tf)
    p_P29__He4_Si26(rate_eval, tf)
    He4_P29__n_Cl32(rate_eval, tf)
    He4_P29__p_S32(rate_eval, tf)
    n_P30__p_Si30(rate_eval, tf)
    n_P30__He4_Al27(rate_eval, tf)
    n_P31__He4_Al28(rate_eval, tf)
    p_P31__n_S31(rate_eval, tf)
    p_P31__He4_Si28(rate_eval, tf)
    p_P31__C12_Ne20(rate_eval, tf)
    p_P31__O16_O16(rate_eval, tf)
    He4_P31__n_Cl34(rate_eval, tf)
    He4_P31__p_S34(rate_eval, tf)
    p_P32__n_S32(rate_eval, tf)
    p_P32__He4_Si29(rate_eval, tf)
    He4_P32__n_Cl35(rate_eval, tf)
    n_S31__p_P31(rate_eval, tf)
    n_S31__He4_Si28(rate_eval, tf)
    n_S31__C12_Ne20(rate_eval, tf)
    n_S31__O16_O16(rate_eval, tf)
    He4_S31__p_Cl34(rate_eval, tf)
    n_S32__p_P32(rate_eval, tf)
    n_S32__He4_Si29(rate_eval, tf)
    p_S32__n_Cl32(rate_eval, tf)
    p_S32__He4_P29(rate_eval, tf)
    He4_S32__n_Ar35(rate_eval, tf)
    He4_S32__p_Cl35(rate_eval, tf)
    p_S34__n_Cl34(rate_eval, tf)
    p_S34__He4_P31(rate_eval, tf)
    n_Cl32__p_S32(rate_eval, tf)
    n_Cl32__He4_P29(rate_eval, tf)
    He4_Cl32__p_Ar35(rate_eval, tf)
    n_Cl34__p_S34(rate_eval, tf)
    n_Cl34__He4_P31(rate_eval, tf)
    p_Cl34__He4_S31(rate_eval, tf)
    n_Cl35__He4_P32(rate_eval, tf)
    p_Cl35__n_Ar35(rate_eval, tf)
    p_Cl35__He4_S32(rate_eval, tf)
    n_Ar35__p_Cl35(rate_eval, tf)
    n_Ar35__He4_S32(rate_eval, tf)
    p_Ar35__He4_Cl32(rate_eval, tf)
    p_d__n_p_p(rate_eval, tf)
    t_t__n_n_He4(rate_eval, tf)
    t_He3__n_p_He4(rate_eval, tf)
    He3_He3__p_p_He4(rate_eval, tf)
    d_Li7__n_He4_He4(rate_eval, tf)
    d_Be7__p_He4_He4(rate_eval, tf)
    p_Be9__d_He4_He4(rate_eval, tf)
    n_B8__p_He4_He4(rate_eval, tf)
    n_C11__He4_He4_He4(rate_eval, tf)
    t_Li7__n_n_He4_He4(rate_eval, tf)
    He3_Li7__n_p_He4_He4(rate_eval, tf)
    t_Be7__n_p_He4_He4(rate_eval, tf)
    He3_Be7__p_p_He4_He4(rate_eval, tf)
    p_Be9__n_p_He4_He4(rate_eval, tf)
    n_n_He4__He6(rate_eval, tf)
    n_p_He4__Li6(rate_eval, tf)
    n_He4_He4__Be9(rate_eval, tf)
    He4_He4_He4__C12(rate_eval, tf)
    n_p_p__p_d(rate_eval, tf)
    n_n_He4__t_t(rate_eval, tf)
    n_p_He4__t_He3(rate_eval, tf)
    p_p_He4__He3_He3(rate_eval, tf)
    n_He4_He4__d_Li7(rate_eval, tf)
    p_He4_He4__n_B8(rate_eval, tf)
    p_He4_He4__d_Be7(rate_eval, tf)
    d_He4_He4__p_Be9(rate_eval, tf)
    He4_He4_He4__n_C11(rate_eval, tf)
    n_n_He4_He4__t_Li7(rate_eval, tf)
    n_p_He4_He4__He3_Li7(rate_eval, tf)
    n_p_He4_He4__t_Be7(rate_eval, tf)
    n_p_He4_He4__p_Be9(rate_eval, tf)
    p_p_He4_He4__He3_Be7(rate_eval, tf)

    if screen_func is not None:
        plasma_state = PlasmaState(T, rho, Y, Z)

        scn_fac = ScreenFactors(1, 1, 1, 1)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_p__d__weak__bet_pos_ *= scor
        rate_eval.p_p__d__weak__electron_capture *= scor
        rate_eval.n_p_p__p_d *= scor
        rate_eval.p_p_He4_He4__He3_Be7 *= scor

        scn_fac = ScreenFactors(1, 1, 1, 2)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_d__He3 *= scor
        rate_eval.p_d__n_p_p *= scor

        scn_fac = ScreenFactors(1, 2, 1, 2)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.d_d__He4 *= scor
        rate_eval.d_d__n_He3 *= scor
        rate_eval.d_d__p_t *= scor

        scn_fac = ScreenFactors(1, 2, 2, 4)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_d__Li6 *= scor
        rate_eval.d_He4__t_He3 *= scor

        scn_fac = ScreenFactors(1, 1, 1, 3)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_t__He4 *= scor
        rate_eval.p_t__n_He3 *= scor
        rate_eval.p_t__d_d *= scor

        scn_fac = ScreenFactors(1, 3, 2, 4)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_t__Li7 *= scor
        rate_eval.He4_t__n_Li6 *= scor

        scn_fac = ScreenFactors(1, 1, 2, 3)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_He3__He4__weak__bet_pos_ *= scor

        scn_fac = ScreenFactors(2, 4, 2, 3)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_He3__Be7 *= scor
        rate_eval.He4_He3__p_Li6 *= scor

        scn_fac = ScreenFactors(1, 1, 3, 6)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Li6__Be7 *= scor
        rate_eval.p_Li6__He4_He3 *= scor

        scn_fac = ScreenFactors(1, 1, 4, 7)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Be7__B8 *= scor

        scn_fac = ScreenFactors(2, 4, 4, 7)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Be7__C11 *= scor

        scn_fac = ScreenFactors(1, 1, 6, 12)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_C12__N13 *= scor

        scn_fac = ScreenFactors(2, 4, 6, 12)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_C12__O16 *= scor
        rate_eval.He4_C12__n_O15 *= scor
        rate_eval.He4_C12__p_N15 *= scor

        scn_fac = ScreenFactors(1, 1, 6, 13)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_C13__N14 *= scor
        rate_eval.p_C13__n_N13 *= scor

        scn_fac = ScreenFactors(1, 1, 6, 14)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_C14__N15 *= scor
        rate_eval.p_C14__n_N14 *= scor

        scn_fac = ScreenFactors(2, 4, 6, 14)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_C14__O18 *= scor
        rate_eval.He4_C14__n_O17 *= scor

        scn_fac = ScreenFactors(1, 1, 7, 13)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_N13__O14 *= scor

        scn_fac = ScreenFactors(1, 1, 7, 14)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_N14__O15 *= scor
        rate_eval.p_N14__n_O14 *= scor
        rate_eval.p_N14__He4_C11 *= scor

        scn_fac = ScreenFactors(2, 4, 7, 14)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_N14__F18 *= scor
        rate_eval.He4_N14__n_F17 *= scor
        rate_eval.He4_N14__p_O17 *= scor

        scn_fac = ScreenFactors(1, 1, 7, 15)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_N15__O16 *= scor
        rate_eval.p_N15__n_O15 *= scor
        rate_eval.p_N15__He4_C12 *= scor

        scn_fac = ScreenFactors(2, 4, 7, 15)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_N15__F19 *= scor
        rate_eval.He4_N15__n_F18 *= scor
        rate_eval.He4_N15__p_O18 *= scor

        scn_fac = ScreenFactors(1, 1, 8, 14)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_O14__F15 *= scor

        scn_fac = ScreenFactors(2, 4, 8, 14)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_O14__Ne18 *= scor
        rate_eval.He4_O14__p_F17 *= scor

        scn_fac = ScreenFactors(1, 1, 8, 15)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_O15__F16 *= scor
        rate_eval.p_O15__n_F15 *= scor

        scn_fac = ScreenFactors(2, 4, 8, 15)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_O15__Ne19 *= scor
        rate_eval.He4_O15__n_Ne18 *= scor
        rate_eval.He4_O15__p_F18 *= scor

        scn_fac = ScreenFactors(1, 1, 8, 16)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_O16__F17 *= scor
        rate_eval.p_O16__n_F16 *= scor
        rate_eval.p_O16__He4_N13 *= scor

        scn_fac = ScreenFactors(2, 4, 8, 16)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_O16__Ne20 *= scor
        rate_eval.He4_O16__n_Ne19 *= scor
        rate_eval.He4_O16__p_F19 *= scor

        scn_fac = ScreenFactors(1, 1, 8, 17)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_O17__F18 *= scor
        rate_eval.p_O17__n_F17 *= scor
        rate_eval.p_O17__He4_N14 *= scor

        scn_fac = ScreenFactors(2, 4, 8, 17)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_O17__Ne21 *= scor
        rate_eval.He4_O17__n_Ne20 *= scor

        scn_fac = ScreenFactors(1, 1, 8, 18)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_O18__F19 *= scor
        rate_eval.p_O18__n_F18 *= scor
        rate_eval.p_O18__He4_N15 *= scor

        scn_fac = ScreenFactors(1, 1, 9, 17)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_F17__Ne18 *= scor
        rate_eval.p_F17__He4_O14 *= scor

        scn_fac = ScreenFactors(2, 4, 9, 17)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_F17__Na21 *= scor
        rate_eval.He4_F17__p_Ne20 *= scor

        scn_fac = ScreenFactors(1, 1, 9, 18)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_F18__Ne19 *= scor
        rate_eval.p_F18__n_Ne18 *= scor
        rate_eval.p_F18__He4_O15 *= scor

        scn_fac = ScreenFactors(2, 4, 9, 18)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_F18__Na22 *= scor
        rate_eval.He4_F18__n_Na21 *= scor
        rate_eval.He4_F18__p_Ne21 *= scor

        scn_fac = ScreenFactors(1, 1, 9, 19)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_F19__Ne20 *= scor
        rate_eval.p_F19__n_Ne19 *= scor
        rate_eval.p_F19__He4_O16 *= scor

        scn_fac = ScreenFactors(2, 4, 9, 19)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_F19__Na23 *= scor
        rate_eval.He4_F19__n_Na22 *= scor

        scn_fac = ScreenFactors(2, 4, 10, 18)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Ne18__Mg22 *= scor
        rate_eval.He4_Ne18__p_Na21 *= scor

        scn_fac = ScreenFactors(2, 4, 10, 19)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Ne19__Mg23 *= scor
        rate_eval.He4_Ne19__n_Mg22 *= scor
        rate_eval.He4_Ne19__p_Na22 *= scor

        scn_fac = ScreenFactors(1, 1, 10, 20)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Ne20__Na21 *= scor
        rate_eval.p_Ne20__He4_F17 *= scor

        scn_fac = ScreenFactors(2, 4, 10, 20)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Ne20__Mg24 *= scor
        rate_eval.He4_Ne20__n_Mg23 *= scor
        rate_eval.He4_Ne20__p_Na23 *= scor
        rate_eval.He4_Ne20__C12_C12 *= scor

        scn_fac = ScreenFactors(1, 1, 10, 21)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Ne21__Na22 *= scor
        rate_eval.p_Ne21__n_Na21 *= scor
        rate_eval.p_Ne21__He4_F18 *= scor

        scn_fac = ScreenFactors(2, 4, 10, 21)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Ne21__Mg25 *= scor
        rate_eval.He4_Ne21__n_Mg24 *= scor
        rate_eval.He4_Ne21__p_Na24 *= scor

        scn_fac = ScreenFactors(1, 1, 11, 21)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Na21__Mg22 *= scor
        rate_eval.p_Na21__He4_Ne18 *= scor

        scn_fac = ScreenFactors(2, 4, 11, 21)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Na21__Al25 *= scor
        rate_eval.He4_Na21__p_Mg24 *= scor

        scn_fac = ScreenFactors(1, 1, 11, 22)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Na22__Mg23 *= scor
        rate_eval.p_Na22__n_Mg22 *= scor
        rate_eval.p_Na22__He4_Ne19 *= scor

        scn_fac = ScreenFactors(2, 4, 11, 22)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Na22__Al26 *= scor
        rate_eval.He4_Na22__n_Al25 *= scor
        rate_eval.He4_Na22__p_Mg25 *= scor

        scn_fac = ScreenFactors(1, 1, 11, 23)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Na23__Mg24 *= scor
        rate_eval.p_Na23__n_Mg23 *= scor
        rate_eval.p_Na23__He4_Ne20 *= scor
        rate_eval.p_Na23__C12_C12 *= scor

        scn_fac = ScreenFactors(2, 4, 11, 23)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Na23__Al27 *= scor
        rate_eval.He4_Na23__n_Al26 *= scor
        rate_eval.He4_Na23__p_Mg26 *= scor

        scn_fac = ScreenFactors(1, 1, 11, 24)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Na24__Mg25 *= scor
        rate_eval.p_Na24__n_Mg24 *= scor
        rate_eval.p_Na24__He4_Ne21 *= scor

        scn_fac = ScreenFactors(2, 4, 11, 24)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Na24__Al28 *= scor
        rate_eval.He4_Na24__n_Al27 *= scor

        scn_fac = ScreenFactors(2, 4, 12, 22)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Mg22__Si26 *= scor
        rate_eval.He4_Mg22__p_Al25 *= scor

        scn_fac = ScreenFactors(1, 1, 12, 24)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Mg24__Al25 *= scor
        rate_eval.p_Mg24__He4_Na21 *= scor

        scn_fac = ScreenFactors(2, 4, 12, 24)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Mg24__Si28 *= scor
        rate_eval.He4_Mg24__p_Al27 *= scor
        rate_eval.He4_Mg24__C12_O16 *= scor

        scn_fac = ScreenFactors(1, 1, 12, 25)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Mg25__Al26 *= scor
        rate_eval.p_Mg25__n_Al25 *= scor
        rate_eval.p_Mg25__He4_Na22 *= scor

        scn_fac = ScreenFactors(2, 4, 12, 25)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Mg25__Si29 *= scor
        rate_eval.He4_Mg25__n_Si28 *= scor
        rate_eval.He4_Mg25__p_Al28 *= scor

        scn_fac = ScreenFactors(1, 1, 12, 26)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Mg26__Al27 *= scor
        rate_eval.p_Mg26__n_Al26 *= scor
        rate_eval.p_Mg26__He4_Na23 *= scor

        scn_fac = ScreenFactors(2, 4, 12, 26)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Mg26__Si30 *= scor
        rate_eval.He4_Mg26__n_Si29 *= scor

        scn_fac = ScreenFactors(1, 1, 13, 25)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Al25__Si26 *= scor
        rate_eval.p_Al25__He4_Mg22 *= scor

        scn_fac = ScreenFactors(2, 4, 13, 25)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Al25__P29 *= scor
        rate_eval.He4_Al25__p_Si28 *= scor

        scn_fac = ScreenFactors(2, 4, 13, 26)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Al26__P30 *= scor
        rate_eval.He4_Al26__n_P29 *= scor
        rate_eval.He4_Al26__p_Si29 *= scor

        scn_fac = ScreenFactors(1, 1, 13, 27)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Al27__Si28 *= scor
        rate_eval.p_Al27__He4_Mg24 *= scor
        rate_eval.p_Al27__C12_O16 *= scor

        scn_fac = ScreenFactors(2, 4, 13, 27)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Al27__P31 *= scor
        rate_eval.He4_Al27__n_P30 *= scor
        rate_eval.He4_Al27__p_Si30 *= scor

        scn_fac = ScreenFactors(1, 1, 13, 28)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Al28__Si29 *= scor
        rate_eval.p_Al28__n_Si28 *= scor
        rate_eval.p_Al28__He4_Mg25 *= scor

        scn_fac = ScreenFactors(2, 4, 13, 28)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Al28__P32 *= scor
        rate_eval.He4_Al28__n_P31 *= scor

        scn_fac = ScreenFactors(1, 1, 14, 28)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Si28__P29 *= scor
        rate_eval.p_Si28__He4_Al25 *= scor

        scn_fac = ScreenFactors(2, 4, 14, 28)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Si28__S32 *= scor
        rate_eval.He4_Si28__n_S31 *= scor
        rate_eval.He4_Si28__p_P31 *= scor
        rate_eval.He4_Si28__C12_Ne20 *= scor
        rate_eval.He4_Si28__O16_O16 *= scor

        scn_fac = ScreenFactors(1, 1, 14, 29)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Si29__P30 *= scor
        rate_eval.p_Si29__n_P29 *= scor
        rate_eval.p_Si29__He4_Al26 *= scor

        scn_fac = ScreenFactors(1, 1, 14, 30)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Si30__P31 *= scor
        rate_eval.p_Si30__n_P30 *= scor
        rate_eval.p_Si30__He4_Al27 *= scor

        scn_fac = ScreenFactors(2, 4, 14, 30)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Si30__S34 *= scor

        scn_fac = ScreenFactors(1, 1, 15, 30)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_P30__S31 *= scor

        scn_fac = ScreenFactors(2, 4, 15, 30)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_P30__Cl34 *= scor

        scn_fac = ScreenFactors(1, 1, 15, 31)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_P31__S32 *= scor
        rate_eval.p_P31__n_S31 *= scor
        rate_eval.p_P31__He4_Si28 *= scor
        rate_eval.p_P31__C12_Ne20 *= scor
        rate_eval.p_P31__O16_O16 *= scor

        scn_fac = ScreenFactors(2, 4, 15, 31)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_P31__Cl35 *= scor
        rate_eval.He4_P31__n_Cl34 *= scor
        rate_eval.He4_P31__p_S34 *= scor

        scn_fac = ScreenFactors(1, 1, 16, 31)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_S31__Cl32 *= scor

        scn_fac = ScreenFactors(2, 4, 16, 31)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_S31__Ar35 *= scor
        rate_eval.He4_S31__p_Cl34 *= scor

        scn_fac = ScreenFactors(1, 1, 16, 34)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_S34__Cl35 *= scor
        rate_eval.p_S34__n_Cl34 *= scor
        rate_eval.p_S34__He4_P31 *= scor

        scn_fac = ScreenFactors(1, 1, 17, 34)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Cl34__Ar35 *= scor
        rate_eval.p_Cl34__He4_S31 *= scor

        scn_fac = ScreenFactors(1, 2, 1, 3)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.d_t__n_He4 *= scor

        scn_fac = ScreenFactors(1, 2, 2, 3)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.d_He3__p_He4 *= scor

        scn_fac = ScreenFactors(1, 3, 2, 3)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.t_He3__d_He4 *= scor
        rate_eval.t_He3__n_p_He4 *= scor

        scn_fac = ScreenFactors(1, 1, 2, 4)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_He4__d_He3 *= scor
        rate_eval.n_p_He4__Li6 *= scor
        rate_eval.n_p_He4__t_He3 *= scor

        scn_fac = ScreenFactors(2, 4, 2, 4)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_He4__n_Be7 *= scor
        rate_eval.He4_He4__p_Li7 *= scor
        rate_eval.n_He4_He4__Be9 *= scor
        rate_eval.n_He4_He4__d_Li7 *= scor
        rate_eval.n_n_He4_He4__t_Li7 *= scor

        scn_fac = ScreenFactors(1, 2, 3, 6)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.d_Li6__n_Be7 *= scor
        rate_eval.d_Li6__p_Li7 *= scor

        scn_fac = ScreenFactors(2, 4, 3, 6)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Li6__p_Be9 *= scor

        scn_fac = ScreenFactors(1, 1, 3, 7)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Li7__n_Be7 *= scor
        rate_eval.p_Li7__d_Li6 *= scor
        rate_eval.p_Li7__He4_He4 *= scor

        scn_fac = ScreenFactors(1, 3, 3, 7)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.t_Li7__n_Be9 *= scor
        rate_eval.t_Li7__n_n_He4_He4 *= scor

        scn_fac = ScreenFactors(1, 1, 4, 9)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Be9__He4_Li6 *= scor
        rate_eval.p_Be9__d_He4_He4 *= scor
        rate_eval.p_Be9__n_p_He4_He4 *= scor

        scn_fac = ScreenFactors(2, 4, 4, 9)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Be9__n_C12 *= scor

        scn_fac = ScreenFactors(2, 4, 5, 8)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_B8__p_C11 *= scor

        scn_fac = ScreenFactors(1, 1, 6, 11)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_C11__He4_B8 *= scor

        scn_fac = ScreenFactors(2, 4, 6, 11)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_C11__n_O14 *= scor
        rate_eval.He4_C11__p_N14 *= scor

        scn_fac = ScreenFactors(6, 12, 6, 12)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.C12_C12__n_Mg23 *= scor
        rate_eval.C12_C12__p_Na23 *= scor
        rate_eval.C12_C12__He4_Ne20 *= scor

        scn_fac = ScreenFactors(1, 2, 6, 13)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.d_C13__n_N14 *= scor

        scn_fac = ScreenFactors(2, 4, 6, 13)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_C13__n_O16 *= scor

        scn_fac = ScreenFactors(1, 2, 6, 14)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.d_C14__n_N15 *= scor

        scn_fac = ScreenFactors(2, 4, 7, 13)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_N13__n_F16 *= scor
        rate_eval.He4_N13__p_O16 *= scor

        scn_fac = ScreenFactors(6, 12, 8, 16)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.C12_O16__p_Al27 *= scor
        rate_eval.C12_O16__He4_Mg24 *= scor

        scn_fac = ScreenFactors(8, 16, 8, 16)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.O16_O16__n_S31 *= scor
        rate_eval.O16_O16__p_P31 *= scor
        rate_eval.O16_O16__He4_Si28 *= scor

        scn_fac = ScreenFactors(2, 4, 8, 18)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_O18__n_Ne21 *= scor

        scn_fac = ScreenFactors(2, 4, 9, 15)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_F15__p_Ne18 *= scor

        scn_fac = ScreenFactors(2, 4, 9, 16)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_F16__p_Ne19 *= scor

        scn_fac = ScreenFactors(1, 1, 10, 18)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Ne18__He4_F15 *= scor

        scn_fac = ScreenFactors(1, 1, 10, 19)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Ne19__He4_F16 *= scor

        scn_fac = ScreenFactors(6, 12, 10, 20)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.C12_Ne20__n_S31 *= scor
        rate_eval.C12_Ne20__p_P31 *= scor
        rate_eval.C12_Ne20__He4_Si28 *= scor

        scn_fac = ScreenFactors(2, 4, 12, 23)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Mg23__n_Si26 *= scor
        rate_eval.He4_Mg23__p_Al26 *= scor

        scn_fac = ScreenFactors(1, 1, 13, 26)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Al26__n_Si26 *= scor
        rate_eval.p_Al26__He4_Mg23 *= scor

        scn_fac = ScreenFactors(2, 4, 14, 26)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Si26__p_P29 *= scor

        scn_fac = ScreenFactors(2, 4, 14, 29)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Si29__n_S32 *= scor
        rate_eval.He4_Si29__p_P32 *= scor

        scn_fac = ScreenFactors(1, 1, 15, 29)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_P29__He4_Si26 *= scor

        scn_fac = ScreenFactors(2, 4, 15, 29)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_P29__n_Cl32 *= scor
        rate_eval.He4_P29__p_S32 *= scor

        scn_fac = ScreenFactors(1, 1, 15, 32)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_P32__n_S32 *= scor
        rate_eval.p_P32__He4_Si29 *= scor

        scn_fac = ScreenFactors(2, 4, 15, 32)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_P32__n_Cl35 *= scor

        scn_fac = ScreenFactors(1, 1, 16, 32)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_S32__n_Cl32 *= scor
        rate_eval.p_S32__He4_P29 *= scor

        scn_fac = ScreenFactors(2, 4, 16, 32)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_S32__n_Ar35 *= scor
        rate_eval.He4_S32__p_Cl35 *= scor

        scn_fac = ScreenFactors(2, 4, 17, 32)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Cl32__p_Ar35 *= scor

        scn_fac = ScreenFactors(1, 1, 17, 35)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Cl35__n_Ar35 *= scor
        rate_eval.p_Cl35__He4_S32 *= scor

        scn_fac = ScreenFactors(1, 1, 18, 35)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Ar35__He4_Cl32 *= scor

        scn_fac = ScreenFactors(1, 3, 1, 3)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.t_t__n_n_He4 *= scor

        scn_fac = ScreenFactors(2, 3, 2, 3)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He3_He3__p_p_He4 *= scor

        scn_fac = ScreenFactors(1, 2, 3, 7)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.d_Li7__n_He4_He4 *= scor

        scn_fac = ScreenFactors(1, 2, 4, 7)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.d_Be7__p_He4_He4 *= scor

        scn_fac = ScreenFactors(2, 3, 3, 7)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He3_Li7__n_p_He4_He4 *= scor

        scn_fac = ScreenFactors(1, 3, 4, 7)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.t_Be7__n_p_He4_He4 *= scor

        scn_fac = ScreenFactors(2, 3, 4, 7)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He3_Be7__p_p_He4_He4 *= scor

        scn_fac = ScreenFactors(2, 4, 2, 4)
        scor = screen_func(plasma_state, scn_fac)
        scn_fac2 = ScreenFactors(2, 4, 4, 8)
        scor2 = screen_func(plasma_state, scn_fac2)
        rate_eval.He4_He4_He4__C12 *= scor * scor2
        rate_eval.He4_He4_He4__n_C11 *= scor * scor2

        scn_fac = ScreenFactors(1, 1, 1, 1)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_p_He4__He3_He3 *= scor

        scn_fac = ScreenFactors(1, 1, 2, 4)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_He4_He4__n_B8 *= scor
        rate_eval.p_He4_He4__d_Be7 *= scor
        rate_eval.n_p_He4_He4__He3_Li7 *= scor
        rate_eval.n_p_He4_He4__t_Be7 *= scor
        rate_eval.n_p_He4_He4__p_Be9 *= scor

        scn_fac = ScreenFactors(1, 2, 2, 4)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.d_He4_He4__p_Be9 *= scor

    dYdt = np.zeros((nnuc), dtype=np.float64)

    dYdt[jn] = (
          -Y[jn]*rate_eval.n__p__weak__wc12  +
          ( -rho*Y[jn]*Y[jp]*rate_eval.n_p__d +Y[jd]*rate_eval.d__n_p ) +
          ( -rho*Y[jn]*Y[jd]*rate_eval.n_d__t +Y[jt]*rate_eval.t__n_d ) +
          ( -rho*Y[jn]*Y[jhe3]*rate_eval.n_He3__He4 +Y[jhe4]*rate_eval.He4__n_He3 ) +
          ( -rho*Y[jn]*Y[jli6]*rate_eval.n_Li6__Li7 +Y[jli7]*rate_eval.Li7__n_Li6 ) +
          ( -rho*Y[jn]*Y[jc11]*rate_eval.n_C11__C12 +Y[jc12]*rate_eval.C12__n_C11 ) +
          ( -rho*Y[jn]*Y[jc12]*rate_eval.n_C12__C13 +Y[jc13]*rate_eval.C13__n_C12 ) +
          ( -rho*Y[jn]*Y[jc13]*rate_eval.n_C13__C14 +Y[jc14]*rate_eval.C14__n_C13 ) +
          ( -rho*Y[jn]*Y[jn13]*rate_eval.n_N13__N14 +Y[jn14]*rate_eval.N14__n_N13 ) +
          ( -rho*Y[jn]*Y[jn14]*rate_eval.n_N14__N15 +Y[jn15]*rate_eval.N15__n_N14 ) +
          ( -rho*Y[jn]*Y[jo14]*rate_eval.n_O14__O15 +Y[jo15]*rate_eval.O15__n_O14 ) +
          ( -rho*Y[jn]*Y[jo15]*rate_eval.n_O15__O16 +Y[jo16]*rate_eval.O16__n_O15 ) +
          ( -rho*Y[jn]*Y[jo16]*rate_eval.n_O16__O17 +Y[jo17]*rate_eval.O17__n_O16 ) +
          ( -rho*Y[jn]*Y[jo17]*rate_eval.n_O17__O18 +Y[jo18]*rate_eval.O18__n_O17 ) +
          ( -rho*Y[jn]*Y[jf15]*rate_eval.n_F15__F16 +Y[jf16]*rate_eval.F16__n_F15 ) +
          ( -rho*Y[jn]*Y[jf16]*rate_eval.n_F16__F17 +Y[jf17]*rate_eval.F17__n_F16 ) +
          ( -rho*Y[jn]*Y[jf17]*rate_eval.n_F17__F18 +Y[jf18]*rate_eval.F18__n_F17 ) +
          ( -rho*Y[jn]*Y[jf18]*rate_eval.n_F18__F19 +Y[jf19]*rate_eval.F19__n_F18 ) +
          ( -rho*Y[jn]*Y[jne18]*rate_eval.n_Ne18__Ne19 +Y[jne19]*rate_eval.Ne19__n_Ne18 ) +
          ( -rho*Y[jn]*Y[jne19]*rate_eval.n_Ne19__Ne20 +Y[jne20]*rate_eval.Ne20__n_Ne19 ) +
          ( -rho*Y[jn]*Y[jne20]*rate_eval.n_Ne20__Ne21 +Y[jne21]*rate_eval.Ne21__n_Ne20 ) +
          ( -rho*Y[jn]*Y[jna21]*rate_eval.n_Na21__Na22 +Y[jna22]*rate_eval.Na22__n_Na21 ) +
          ( -rho*Y[jn]*Y[jna22]*rate_eval.n_Na22__Na23 +Y[jna23]*rate_eval.Na23__n_Na22 ) +
          ( -rho*Y[jn]*Y[jna23]*rate_eval.n_Na23__Na24 +Y[jna24]*rate_eval.Na24__n_Na23 ) +
          ( -rho*Y[jn]*Y[jmg22]*rate_eval.n_Mg22__Mg23 +Y[jmg23]*rate_eval.Mg23__n_Mg22 ) +
          ( -rho*Y[jn]*Y[jmg23]*rate_eval.n_Mg23__Mg24 +Y[jmg24]*rate_eval.Mg24__n_Mg23 ) +
          ( -rho*Y[jn]*Y[jmg24]*rate_eval.n_Mg24__Mg25 +Y[jmg25]*rate_eval.Mg25__n_Mg24 ) +
          ( -rho*Y[jn]*Y[jmg25]*rate_eval.n_Mg25__Mg26 +Y[jmg26]*rate_eval.Mg26__n_Mg25 ) +
          ( -rho*Y[jn]*Y[jal25]*rate_eval.n_Al25__Al26 +Y[jal26]*rate_eval.Al26__n_Al25 ) +
          ( -rho*Y[jn]*Y[jal26]*rate_eval.n_Al26__Al27 +Y[jal27]*rate_eval.Al27__n_Al26 ) +
          ( -rho*Y[jn]*Y[jal27]*rate_eval.n_Al27__Al28 +Y[jal28]*rate_eval.Al28__n_Al27 ) +
          ( -rho*Y[jn]*Y[jsi28]*rate_eval.n_Si28__Si29 +Y[jsi29]*rate_eval.Si29__n_Si28 ) +
          ( -rho*Y[jn]*Y[jsi29]*rate_eval.n_Si29__Si30 +Y[jsi30]*rate_eval.Si30__n_Si29 ) +
          ( -rho*Y[jn]*Y[jp29]*rate_eval.n_P29__P30 +Y[jp30]*rate_eval.P30__n_P29 ) +
          ( -rho*Y[jn]*Y[jp30]*rate_eval.n_P30__P31 +Y[jp31]*rate_eval.P31__n_P30 ) +
          ( -rho*Y[jn]*Y[jp31]*rate_eval.n_P31__P32 +Y[jp32]*rate_eval.P32__n_P31 ) +
          ( -rho*Y[jn]*Y[js31]*rate_eval.n_S31__S32 +Y[js32]*rate_eval.S32__n_S31 ) +
          ( -rho*Y[jn]*Y[jcl34]*rate_eval.n_Cl34__Cl35 +Y[jcl35]*rate_eval.Cl35__n_Cl34 ) +
          ( +5.00000000000000e-01*rho*Y[jd]**2*rate_eval.d_d__n_He3 -rho*Y[jn]*Y[jhe3]*rate_eval.n_He3__d_d ) +
          ( +rho*Y[jd]*Y[jt]*rate_eval.d_t__n_He4 -rho*Y[jn]*Y[jhe4]*rate_eval.n_He4__d_t ) +
          ( -rho*Y[jn]*Y[jhe3]*rate_eval.n_He3__p_t +rho*Y[jp]*Y[jt]*rate_eval.p_t__n_He3 ) +
          ( -rho*Y[jn]*Y[jli6]*rate_eval.n_Li6__He4_t +rho*Y[jt]*Y[jhe4]*rate_eval.He4_t__n_Li6 ) +
          ( +rho*Y[jd]*Y[jli6]*rate_eval.d_Li6__n_Be7 -rho*Y[jn]*Y[jbe7]*rate_eval.n_Be7__d_Li6 ) +
          ( +rho*Y[jt]*Y[jli7]*rate_eval.t_Li7__n_Be9 -rho*Y[jn]*Y[jbe9]*rate_eval.n_Be9__t_Li7 ) +
          ( -rho*Y[jn]*Y[jbe7]*rate_eval.n_Be7__p_Li7 +rho*Y[jp]*Y[jli7]*rate_eval.p_Li7__n_Be7 ) +
          ( -rho*Y[jn]*Y[jbe7]*rate_eval.n_Be7__He4_He4 +5.00000000000000e-01*rho*Y[jhe4]**2*rate_eval.He4_He4__n_Be7 ) +
          ( +rho*Y[jhe4]*Y[jbe9]*rate_eval.He4_Be9__n_C12 -rho*Y[jn]*Y[jc12]*rate_eval.n_C12__He4_Be9 ) +
          ( +rho*Y[jd]*Y[jc13]*rate_eval.d_C13__n_N14 -rho*Y[jn]*Y[jn14]*rate_eval.n_N14__d_C13 ) +
          ( +rho*Y[jhe4]*Y[jc13]*rate_eval.He4_C13__n_O16 -rho*Y[jn]*Y[jo16]*rate_eval.n_O16__He4_C13 ) +
          ( +rho*Y[jd]*Y[jc14]*rate_eval.d_C14__n_N15 -rho*Y[jn]*Y[jn15]*rate_eval.n_N15__d_C14 ) +
          ( -rho*Y[jn]*Y[jn13]*rate_eval.n_N13__p_C13 +rho*Y[jp]*Y[jc13]*rate_eval.p_C13__n_N13 ) +
          ( -rho*Y[jn]*Y[jn14]*rate_eval.n_N14__p_C14 +rho*Y[jp]*Y[jc14]*rate_eval.p_C14__n_N14 ) +
          ( -rho*Y[jn]*Y[jo14]*rate_eval.n_O14__p_N14 +rho*Y[jp]*Y[jn14]*rate_eval.p_N14__n_O14 ) +
          ( -rho*Y[jn]*Y[jo14]*rate_eval.n_O14__He4_C11 +rho*Y[jhe4]*Y[jc11]*rate_eval.He4_C11__n_O14 ) +
          ( -rho*Y[jn]*Y[jo15]*rate_eval.n_O15__p_N15 +rho*Y[jp]*Y[jn15]*rate_eval.p_N15__n_O15 ) +
          ( -rho*Y[jn]*Y[jo15]*rate_eval.n_O15__He4_C12 +rho*Y[jhe4]*Y[jc12]*rate_eval.He4_C12__n_O15 ) +
          ( +5.00000000000000e-01*rho*Y[jo16]**2*rate_eval.O16_O16__n_S31 -rho*Y[jn]*Y[js31]*rate_eval.n_S31__O16_O16 ) +
          ( -rho*Y[jn]*Y[jo17]*rate_eval.n_O17__He4_C14 +rho*Y[jhe4]*Y[jc14]*rate_eval.He4_C14__n_O17 ) +
          ( +rho*Y[jhe4]*Y[jo17]*rate_eval.He4_O17__n_Ne20 -rho*Y[jn]*Y[jne20]*rate_eval.n_Ne20__He4_O17 ) +
          ( -rho*Y[jn]*Y[jf15]*rate_eval.n_F15__p_O15 +rho*Y[jp]*Y[jo15]*rate_eval.p_O15__n_F15 ) +
          ( -rho*Y[jn]*Y[jf16]*rate_eval.n_F16__p_O16 +rho*Y[jp]*Y[jo16]*rate_eval.p_O16__n_F16 ) +
          ( -rho*Y[jn]*Y[jf16]*rate_eval.n_F16__He4_N13 +rho*Y[jhe4]*Y[jn13]*rate_eval.He4_N13__n_F16 ) +
          ( -rho*Y[jn]*Y[jf17]*rate_eval.n_F17__p_O17 +rho*Y[jp]*Y[jo17]*rate_eval.p_O17__n_F17 ) +
          ( -rho*Y[jn]*Y[jf17]*rate_eval.n_F17__He4_N14 +rho*Y[jhe4]*Y[jn14]*rate_eval.He4_N14__n_F17 ) +
          ( -rho*Y[jn]*Y[jf18]*rate_eval.n_F18__p_O18 +rho*Y[jp]*Y[jo18]*rate_eval.p_O18__n_F18 ) +
          ( -rho*Y[jn]*Y[jf18]*rate_eval.n_F18__He4_N15 +rho*Y[jhe4]*Y[jn15]*rate_eval.He4_N15__n_F18 ) +
          ( -rho*Y[jn]*Y[jne18]*rate_eval.n_Ne18__p_F18 +rho*Y[jp]*Y[jf18]*rate_eval.p_F18__n_Ne18 ) +
          ( -rho*Y[jn]*Y[jne18]*rate_eval.n_Ne18__He4_O15 +rho*Y[jhe4]*Y[jo15]*rate_eval.He4_O15__n_Ne18 ) +
          ( -rho*Y[jn]*Y[jne19]*rate_eval.n_Ne19__p_F19 +rho*Y[jp]*Y[jf19]*rate_eval.p_F19__n_Ne19 ) +
          ( -rho*Y[jn]*Y[jne19]*rate_eval.n_Ne19__He4_O16 +rho*Y[jhe4]*Y[jo16]*rate_eval.He4_O16__n_Ne19 ) +
          ( +rho*Y[jc12]*Y[jne20]*rate_eval.C12_Ne20__n_S31 -rho*Y[jn]*Y[js31]*rate_eval.n_S31__C12_Ne20 ) +
          ( -rho*Y[jn]*Y[jne21]*rate_eval.n_Ne21__He4_O18 +rho*Y[jhe4]*Y[jo18]*rate_eval.He4_O18__n_Ne21 ) +
          ( +rho*Y[jhe4]*Y[jne21]*rate_eval.He4_Ne21__n_Mg24 -rho*Y[jn]*Y[jmg24]*rate_eval.n_Mg24__He4_Ne21 ) +
          ( -rho*Y[jn]*Y[jna21]*rate_eval.n_Na21__p_Ne21 +rho*Y[jp]*Y[jne21]*rate_eval.p_Ne21__n_Na21 ) +
          ( -rho*Y[jn]*Y[jna21]*rate_eval.n_Na21__He4_F18 +rho*Y[jhe4]*Y[jf18]*rate_eval.He4_F18__n_Na21 ) +
          ( -rho*Y[jn]*Y[jna22]*rate_eval.n_Na22__He4_F19 +rho*Y[jhe4]*Y[jf19]*rate_eval.He4_F19__n_Na22 ) +
          ( +rho*Y[jp]*Y[jna24]*rate_eval.p_Na24__n_Mg24 -rho*Y[jn]*Y[jmg24]*rate_eval.n_Mg24__p_Na24 ) +
          ( +rho*Y[jhe4]*Y[jna24]*rate_eval.He4_Na24__n_Al27 -rho*Y[jn]*Y[jal27]*rate_eval.n_Al27__He4_Na24 ) +
          ( -rho*Y[jn]*Y[jmg22]*rate_eval.n_Mg22__p_Na22 +rho*Y[jp]*Y[jna22]*rate_eval.p_Na22__n_Mg22 ) +
          ( -rho*Y[jn]*Y[jmg22]*rate_eval.n_Mg22__He4_Ne19 +rho*Y[jhe4]*Y[jne19]*rate_eval.He4_Ne19__n_Mg22 ) +
          ( -rho*Y[jn]*Y[jmg23]*rate_eval.n_Mg23__p_Na23 +rho*Y[jp]*Y[jna23]*rate_eval.p_Na23__n_Mg23 ) +
          ( -rho*Y[jn]*Y[jmg23]*rate_eval.n_Mg23__He4_Ne20 +rho*Y[jhe4]*Y[jne20]*rate_eval.He4_Ne20__n_Mg23 ) +
          ( -rho*Y[jn]*Y[jmg23]*rate_eval.n_Mg23__C12_C12 +5.00000000000000e-01*rho*Y[jc12]**2*rate_eval.C12_C12__n_Mg23 ) +
          ( +rho*Y[jhe4]*Y[jmg25]*rate_eval.He4_Mg25__n_Si28 -rho*Y[jn]*Y[jsi28]*rate_eval.n_Si28__He4_Mg25 ) +
          ( +rho*Y[jhe4]*Y[jmg26]*rate_eval.He4_Mg26__n_Si29 -rho*Y[jn]*Y[jsi29]*rate_eval.n_Si29__He4_Mg26 ) +
          ( -rho*Y[jn]*Y[jal25]*rate_eval.n_Al25__p_Mg25 +rho*Y[jp]*Y[jmg25]*rate_eval.p_Mg25__n_Al25 ) +
          ( -rho*Y[jn]*Y[jal25]*rate_eval.n_Al25__He4_Na22 +rho*Y[jhe4]*Y[jna22]*rate_eval.He4_Na22__n_Al25 ) +
          ( -rho*Y[jn]*Y[jal26]*rate_eval.n_Al26__p_Mg26 +rho*Y[jp]*Y[jmg26]*rate_eval.p_Mg26__n_Al26 ) +
          ( -rho*Y[jn]*Y[jal26]*rate_eval.n_Al26__He4_Na23 +rho*Y[jhe4]*Y[jna23]*rate_eval.He4_Na23__n_Al26 ) +
          ( +rho*Y[jp]*Y[jal28]*rate_eval.p_Al28__n_Si28 -rho*Y[jn]*Y[jsi28]*rate_eval.n_Si28__p_Al28 ) +
          ( +rho*Y[jhe4]*Y[jal28]*rate_eval.He4_Al28__n_P31 -rho*Y[jn]*Y[jp31]*rate_eval.n_P31__He4_Al28 ) +
          ( -rho*Y[jn]*Y[jsi26]*rate_eval.n_Si26__p_Al26 +rho*Y[jp]*Y[jal26]*rate_eval.p_Al26__n_Si26 ) +
          ( -rho*Y[jn]*Y[jsi26]*rate_eval.n_Si26__He4_Mg23 +rho*Y[jhe4]*Y[jmg23]*rate_eval.He4_Mg23__n_Si26 ) +
          ( -rho*Y[jn]*Y[jp29]*rate_eval.n_P29__p_Si29 +rho*Y[jp]*Y[jsi29]*rate_eval.p_Si29__n_P29 ) +
          ( -rho*Y[jn]*Y[jp29]*rate_eval.n_P29__He4_Al26 +rho*Y[jhe4]*Y[jal26]*rate_eval.He4_Al26__n_P29 ) +
          ( -rho*Y[jn]*Y[jp30]*rate_eval.n_P30__p_Si30 +rho*Y[jp]*Y[jsi30]*rate_eval.p_Si30__n_P30 ) +
          ( -rho*Y[jn]*Y[jp30]*rate_eval.n_P30__He4_Al27 +rho*Y[jhe4]*Y[jal27]*rate_eval.He4_Al27__n_P30 ) +
          ( +rho*Y[jp]*Y[jp32]*rate_eval.p_P32__n_S32 -rho*Y[jn]*Y[js32]*rate_eval.n_S32__p_P32 ) +
          ( -rho*Y[jn]*Y[js31]*rate_eval.n_S31__p_P31 +rho*Y[jp]*Y[jp31]*rate_eval.p_P31__n_S31 ) +
          ( -rho*Y[jn]*Y[js31]*rate_eval.n_S31__He4_Si28 +rho*Y[jhe4]*Y[jsi28]*rate_eval.He4_Si28__n_S31 ) +
          ( -rho*Y[jn]*Y[js32]*rate_eval.n_S32__He4_Si29 +rho*Y[jhe4]*Y[jsi29]*rate_eval.He4_Si29__n_S32 ) +
          ( -rho*Y[jn]*Y[jcl32]*rate_eval.n_Cl32__p_S32 +rho*Y[jp]*Y[js32]*rate_eval.p_S32__n_Cl32 ) +
          ( -rho*Y[jn]*Y[jcl32]*rate_eval.n_Cl32__He4_P29 +rho*Y[jhe4]*Y[jp29]*rate_eval.He4_P29__n_Cl32 ) +
          ( -rho*Y[jn]*Y[jcl34]*rate_eval.n_Cl34__p_S34 +rho*Y[jp]*Y[js34]*rate_eval.p_S34__n_Cl34 ) +
          ( -rho*Y[jn]*Y[jcl34]*rate_eval.n_Cl34__He4_P31 +rho*Y[jhe4]*Y[jp31]*rate_eval.He4_P31__n_Cl34 ) +
          ( -rho*Y[jn]*Y[jcl35]*rate_eval.n_Cl35__He4_P32 +rho*Y[jhe4]*Y[jp32]*rate_eval.He4_P32__n_Cl35 ) +
          ( -rho*Y[jn]*Y[jar35]*rate_eval.n_Ar35__p_Cl35 +rho*Y[jp]*Y[jcl35]*rate_eval.p_Cl35__n_Ar35 ) +
          ( -rho*Y[jn]*Y[jar35]*rate_eval.n_Ar35__He4_S32 +rho*Y[jhe4]*Y[js32]*rate_eval.He4_S32__n_Ar35 ) +
          ( + 2*5.00000000000000e-01*rho*Y[jt]**2*rate_eval.t_t__n_n_He4 + -2*5.00000000000000e-01*rho**2*Y[jn]**2*Y[jhe4]*rate_eval.n_n_He4__t_t ) +
          ( +rho*Y[jt]*Y[jhe3]*rate_eval.t_He3__n_p_He4 -rho**2*Y[jn]*Y[jp]*Y[jhe4]*rate_eval.n_p_He4__t_He3 ) +
          ( +rho*Y[jd]*Y[jli7]*rate_eval.d_Li7__n_He4_He4 -5.00000000000000e-01*rho**2*Y[jn]*Y[jhe4]**2*rate_eval.n_He4_He4__d_Li7 ) +
          ( -rho*Y[jn]*Y[jb8]*rate_eval.n_B8__p_He4_He4 +5.00000000000000e-01*rho**2*Y[jp]*Y[jhe4]**2*rate_eval.p_He4_He4__n_B8 ) +
          ( -rho*Y[jn]*Y[jc11]*rate_eval.n_C11__He4_He4_He4 +1.66666666666667e-01*rho**2*Y[jhe4]**3*rate_eval.He4_He4_He4__n_C11 ) +
          ( + 2*rho*Y[jt]*Y[jli7]*rate_eval.t_Li7__n_n_He4_He4 + -2*2.50000000000000e-01*rho**3*Y[jn]**2*Y[jhe4]**2*rate_eval.n_n_He4_He4__t_Li7 ) +
          ( +rho*Y[jhe3]*Y[jli7]*rate_eval.He3_Li7__n_p_He4_He4 -5.00000000000000e-01*rho**3*Y[jn]*Y[jp]*Y[jhe4]**2*rate_eval.n_p_He4_He4__He3_Li7 ) +
          ( +rho*Y[jt]*Y[jbe7]*rate_eval.t_Be7__n_p_He4_He4 -5.00000000000000e-01*rho**3*Y[jn]*Y[jp]*Y[jhe4]**2*rate_eval.n_p_He4_He4__t_Be7 ) +
          ( + -2*5.00000000000000e-01*rho**2*Y[jn]**2*Y[jhe4]*rate_eval.n_n_He4__He6 + 2*Y[jhe6]*rate_eval.He6__n_n_He4 ) +
          ( -rho**2*Y[jn]*Y[jp]*Y[jhe4]*rate_eval.n_p_He4__Li6 +Y[jli6]*rate_eval.Li6__n_p_He4 ) +
          ( -5.00000000000000e-01*rho**2*Y[jn]*Y[jhe4]**2*rate_eval.n_He4_He4__Be9 +Y[jbe9]*rate_eval.Be9__n_He4_He4 ) +
          ( -5.00000000000000e-01*rho**2*Y[jn]*Y[jp]**2*rate_eval.n_p_p__p_d +rho*Y[jp]*Y[jd]*rate_eval.p_d__n_p_p ) +
          ( -5.00000000000000e-01*rho**3*Y[jn]*Y[jp]*Y[jhe4]**2*rate_eval.n_p_He4_He4__p_Be9 +rho*Y[jp]*Y[jbe9]*rate_eval.p_Be9__n_p_He4_He4 )
       )

    dYdt[jp] = (
          +Y[jn]*rate_eval.n__p__weak__wc12  +
          ( +Y[jf15]*rate_eval.F15__p_O14 -rho*Y[jp]*Y[jo14]*rate_eval.p_O14__F15 ) +
          ( +Y[jf16]*rate_eval.F16__p_O15 -rho*Y[jp]*Y[jo15]*rate_eval.p_O15__F16 ) +
          +Y[jcl32]*rate_eval.Cl32__p_P31__weak__wc12  +
          ( -rho*Y[jn]*Y[jp]*rate_eval.n_p__d +Y[jd]*rate_eval.d__n_p ) +
          + -2*5.00000000000000e-01*rho*Y[jp]**2*rate_eval.p_p__d__weak__bet_pos_  +
          + -2*5.00000000000000e-01*rho**2*ye(Y)*Y[jp]**2*rate_eval.p_p__d__weak__electron_capture  +
          ( -rho*Y[jp]*Y[jd]*rate_eval.p_d__He3 +Y[jhe3]*rate_eval.He3__p_d ) +
          ( -rho*Y[jp]*Y[jt]*rate_eval.p_t__He4 +Y[jhe4]*rate_eval.He4__p_t ) +
          -rho*Y[jp]*Y[jhe3]*rate_eval.p_He3__He4__weak__bet_pos_  +
          ( -rho*Y[jp]*Y[jli6]*rate_eval.p_Li6__Be7 +Y[jbe7]*rate_eval.Be7__p_Li6 ) +
          ( -rho*Y[jp]*Y[jbe7]*rate_eval.p_Be7__B8 +Y[jb8]*rate_eval.B8__p_Be7 ) +
          ( -rho*Y[jp]*Y[jc12]*rate_eval.p_C12__N13 +Y[jn13]*rate_eval.N13__p_C12 ) +
          ( -rho*Y[jp]*Y[jc13]*rate_eval.p_C13__N14 +Y[jn14]*rate_eval.N14__p_C13 ) +
          ( -rho*Y[jp]*Y[jc14]*rate_eval.p_C14__N15 +Y[jn15]*rate_eval.N15__p_C14 ) +
          ( -rho*Y[jp]*Y[jn13]*rate_eval.p_N13__O14 +Y[jo14]*rate_eval.O14__p_N13 ) +
          ( -rho*Y[jp]*Y[jn14]*rate_eval.p_N14__O15 +Y[jo15]*rate_eval.O15__p_N14 ) +
          ( -rho*Y[jp]*Y[jn15]*rate_eval.p_N15__O16 +Y[jo16]*rate_eval.O16__p_N15 ) +
          ( -rho*Y[jp]*Y[jo16]*rate_eval.p_O16__F17 +Y[jf17]*rate_eval.F17__p_O16 ) +
          ( -rho*Y[jp]*Y[jo17]*rate_eval.p_O17__F18 +Y[jf18]*rate_eval.F18__p_O17 ) +
          ( -rho*Y[jp]*Y[jo18]*rate_eval.p_O18__F19 +Y[jf19]*rate_eval.F19__p_O18 ) +
          ( -rho*Y[jp]*Y[jf17]*rate_eval.p_F17__Ne18 +Y[jne18]*rate_eval.Ne18__p_F17 ) +
          ( -rho*Y[jp]*Y[jf18]*rate_eval.p_F18__Ne19 +Y[jne19]*rate_eval.Ne19__p_F18 ) +
          ( -rho*Y[jp]*Y[jf19]*rate_eval.p_F19__Ne20 +Y[jne20]*rate_eval.Ne20__p_F19 ) +
          ( -rho*Y[jp]*Y[jne20]*rate_eval.p_Ne20__Na21 +Y[jna21]*rate_eval.Na21__p_Ne20 ) +
          ( -rho*Y[jp]*Y[jne21]*rate_eval.p_Ne21__Na22 +Y[jna22]*rate_eval.Na22__p_Ne21 ) +
          ( -rho*Y[jp]*Y[jna21]*rate_eval.p_Na21__Mg22 +Y[jmg22]*rate_eval.Mg22__p_Na21 ) +
          ( -rho*Y[jp]*Y[jna22]*rate_eval.p_Na22__Mg23 +Y[jmg23]*rate_eval.Mg23__p_Na22 ) +
          ( -rho*Y[jp]*Y[jna23]*rate_eval.p_Na23__Mg24 +Y[jmg24]*rate_eval.Mg24__p_Na23 ) +
          ( -rho*Y[jp]*Y[jna24]*rate_eval.p_Na24__Mg25 +Y[jmg25]*rate_eval.Mg25__p_Na24 ) +
          ( -rho*Y[jp]*Y[jmg24]*rate_eval.p_Mg24__Al25 +Y[jal25]*rate_eval.Al25__p_Mg24 ) +
          ( -rho*Y[jp]*Y[jmg25]*rate_eval.p_Mg25__Al26 +Y[jal26]*rate_eval.Al26__p_Mg25 ) +
          ( -rho*Y[jp]*Y[jmg26]*rate_eval.p_Mg26__Al27 +Y[jal27]*rate_eval.Al27__p_Mg26 ) +
          ( -rho*Y[jp]*Y[jal25]*rate_eval.p_Al25__Si26 +Y[jsi26]*rate_eval.Si26__p_Al25 ) +
          ( -rho*Y[jp]*Y[jal27]*rate_eval.p_Al27__Si28 +Y[jsi28]*rate_eval.Si28__p_Al27 ) +
          ( -rho*Y[jp]*Y[jal28]*rate_eval.p_Al28__Si29 +Y[jsi29]*rate_eval.Si29__p_Al28 ) +
          ( -rho*Y[jp]*Y[jsi28]*rate_eval.p_Si28__P29 +Y[jp29]*rate_eval.P29__p_Si28 ) +
          ( -rho*Y[jp]*Y[jsi29]*rate_eval.p_Si29__P30 +Y[jp30]*rate_eval.P30__p_Si29 ) +
          ( -rho*Y[jp]*Y[jsi30]*rate_eval.p_Si30__P31 +Y[jp31]*rate_eval.P31__p_Si30 ) +
          ( -rho*Y[jp]*Y[jp30]*rate_eval.p_P30__S31 +Y[js31]*rate_eval.S31__p_P30 ) +
          ( -rho*Y[jp]*Y[jp31]*rate_eval.p_P31__S32 +Y[js32]*rate_eval.S32__p_P31 ) +
          ( -rho*Y[jp]*Y[js31]*rate_eval.p_S31__Cl32 +Y[jcl32]*rate_eval.Cl32__p_S31 ) +
          ( -rho*Y[jp]*Y[js34]*rate_eval.p_S34__Cl35 +Y[jcl35]*rate_eval.Cl35__p_S34 ) +
          ( -rho*Y[jp]*Y[jcl34]*rate_eval.p_Cl34__Ar35 +Y[jar35]*rate_eval.Ar35__p_Cl34 ) +
          ( +5.00000000000000e-01*rho*Y[jd]**2*rate_eval.d_d__p_t -rho*Y[jp]*Y[jt]*rate_eval.p_t__d_d ) +
          ( +rho*Y[jn]*Y[jhe3]*rate_eval.n_He3__p_t -rho*Y[jp]*Y[jt]*rate_eval.p_t__n_He3 ) +
          ( +rho*Y[jd]*Y[jhe3]*rate_eval.d_He3__p_He4 -rho*Y[jp]*Y[jhe4]*rate_eval.p_He4__d_He3 ) +
          ( -rho*Y[jp]*Y[jli6]*rate_eval.p_Li6__He4_He3 +rho*Y[jhe3]*Y[jhe4]*rate_eval.He4_He3__p_Li6 ) +
          ( +rho*Y[jd]*Y[jli6]*rate_eval.d_Li6__p_Li7 -rho*Y[jp]*Y[jli7]*rate_eval.p_Li7__d_Li6 ) +
          ( -rho*Y[jp]*Y[jli7]*rate_eval.p_Li7__He4_He4 +5.00000000000000e-01*rho*Y[jhe4]**2*rate_eval.He4_He4__p_Li7 ) +
          ( +rho*Y[jn]*Y[jbe7]*rate_eval.n_Be7__p_Li7 -rho*Y[jp]*Y[jli7]*rate_eval.p_Li7__n_Be7 ) +
          ( -rho*Y[jp]*Y[jbe9]*rate_eval.p_Be9__He4_Li6 +rho*Y[jhe4]*Y[jli6]*rate_eval.He4_Li6__p_Be9 ) +
          ( +rho*Y[jhe4]*Y[jb8]*rate_eval.He4_B8__p_C11 -rho*Y[jp]*Y[jc11]*rate_eval.p_C11__He4_B8 ) +
          ( +rho*Y[jhe4]*Y[jc11]*rate_eval.He4_C11__p_N14 -rho*Y[jp]*Y[jn14]*rate_eval.p_N14__He4_C11 ) +
          ( +5.00000000000000e-01*rho*Y[jc12]**2*rate_eval.C12_C12__p_Na23 -rho*Y[jp]*Y[jna23]*rate_eval.p_Na23__C12_C12 ) +
          ( +rho*Y[jn]*Y[jn13]*rate_eval.n_N13__p_C13 -rho*Y[jp]*Y[jc13]*rate_eval.p_C13__n_N13 ) +
          ( +rho*Y[jhe4]*Y[jn13]*rate_eval.He4_N13__p_O16 -rho*Y[jp]*Y[jo16]*rate_eval.p_O16__He4_N13 ) +
          ( +rho*Y[jn]*Y[jn14]*rate_eval.n_N14__p_C14 -rho*Y[jp]*Y[jc14]*rate_eval.p_C14__n_N14 ) +
          ( -rho*Y[jp]*Y[jn15]*rate_eval.p_N15__He4_C12 +rho*Y[jhe4]*Y[jc12]*rate_eval.He4_C12__p_N15 ) +
          ( +rho*Y[jn]*Y[jo14]*rate_eval.n_O14__p_N14 -rho*Y[jp]*Y[jn14]*rate_eval.p_N14__n_O14 ) +
          ( +rho*Y[jhe4]*Y[jo14]*rate_eval.He4_O14__p_F17 -rho*Y[jp]*Y[jf17]*rate_eval.p_F17__He4_O14 ) +
          ( +rho*Y[jn]*Y[jo15]*rate_eval.n_O15__p_N15 -rho*Y[jp]*Y[jn15]*rate_eval.p_N15__n_O15 ) +
          ( +rho*Y[jc12]*Y[jo16]*rate_eval.C12_O16__p_Al27 -rho*Y[jp]*Y[jal27]*rate_eval.p_Al27__C12_O16 ) +
          ( +5.00000000000000e-01*rho*Y[jo16]**2*rate_eval.O16_O16__p_P31 -rho*Y[jp]*Y[jp31]*rate_eval.p_P31__O16_O16 ) +
          ( -rho*Y[jp]*Y[jo17]*rate_eval.p_O17__He4_N14 +rho*Y[jhe4]*Y[jn14]*rate_eval.He4_N14__p_O17 ) +
          ( -rho*Y[jp]*Y[jo18]*rate_eval.p_O18__He4_N15 +rho*Y[jhe4]*Y[jn15]*rate_eval.He4_N15__p_O18 ) +
          ( +rho*Y[jn]*Y[jf15]*rate_eval.n_F15__p_O15 -rho*Y[jp]*Y[jo15]*rate_eval.p_O15__n_F15 ) +
          ( +rho*Y[jhe4]*Y[jf15]*rate_eval.He4_F15__p_Ne18 -rho*Y[jp]*Y[jne18]*rate_eval.p_Ne18__He4_F15 ) +
          ( +rho*Y[jn]*Y[jf16]*rate_eval.n_F16__p_O16 -rho*Y[jp]*Y[jo16]*rate_eval.p_O16__n_F16 ) +
          ( +rho*Y[jhe4]*Y[jf16]*rate_eval.He4_F16__p_Ne19 -rho*Y[jp]*Y[jne19]*rate_eval.p_Ne19__He4_F16 ) +
          ( +rho*Y[jn]*Y[jf17]*rate_eval.n_F17__p_O17 -rho*Y[jp]*Y[jo17]*rate_eval.p_O17__n_F17 ) +
          ( +rho*Y[jhe4]*Y[jf17]*rate_eval.He4_F17__p_Ne20 -rho*Y[jp]*Y[jne20]*rate_eval.p_Ne20__He4_F17 ) +
          ( +rho*Y[jn]*Y[jf18]*rate_eval.n_F18__p_O18 -rho*Y[jp]*Y[jo18]*rate_eval.p_O18__n_F18 ) +
          ( -rho*Y[jp]*Y[jf18]*rate_eval.p_F18__He4_O15 +rho*Y[jhe4]*Y[jo15]*rate_eval.He4_O15__p_F18 ) +
          ( +rho*Y[jhe4]*Y[jf18]*rate_eval.He4_F18__p_Ne21 -rho*Y[jp]*Y[jne21]*rate_eval.p_Ne21__He4_F18 ) +
          ( -rho*Y[jp]*Y[jf19]*rate_eval.p_F19__He4_O16 +rho*Y[jhe4]*Y[jo16]*rate_eval.He4_O16__p_F19 ) +
          ( +rho*Y[jn]*Y[jne18]*rate_eval.n_Ne18__p_F18 -rho*Y[jp]*Y[jf18]*rate_eval.p_F18__n_Ne18 ) +
          ( +rho*Y[jhe4]*Y[jne18]*rate_eval.He4_Ne18__p_Na21 -rho*Y[jp]*Y[jna21]*rate_eval.p_Na21__He4_Ne18 ) +
          ( +rho*Y[jn]*Y[jne19]*rate_eval.n_Ne19__p_F19 -rho*Y[jp]*Y[jf19]*rate_eval.p_F19__n_Ne19 ) +
          ( +rho*Y[jhe4]*Y[jne19]*rate_eval.He4_Ne19__p_Na22 -rho*Y[jp]*Y[jna22]*rate_eval.p_Na22__He4_Ne19 ) +
          ( +rho*Y[jc12]*Y[jne20]*rate_eval.C12_Ne20__p_P31 -rho*Y[jp]*Y[jp31]*rate_eval.p_P31__C12_Ne20 ) +
          ( +rho*Y[jn]*Y[jna21]*rate_eval.n_Na21__p_Ne21 -rho*Y[jp]*Y[jne21]*rate_eval.p_Ne21__n_Na21 ) +
          ( +rho*Y[jhe4]*Y[jna21]*rate_eval.He4_Na21__p_Mg24 -rho*Y[jp]*Y[jmg24]*rate_eval.p_Mg24__He4_Na21 ) +
          ( +rho*Y[jhe4]*Y[jna22]*rate_eval.He4_Na22__p_Mg25 -rho*Y[jp]*Y[jmg25]*rate_eval.p_Mg25__He4_Na22 ) +
          ( -rho*Y[jp]*Y[jna23]*rate_eval.p_Na23__He4_Ne20 +rho*Y[jhe4]*Y[jne20]*rate_eval.He4_Ne20__p_Na23 ) +
          ( +rho*Y[jhe4]*Y[jna23]*rate_eval.He4_Na23__p_Mg26 -rho*Y[jp]*Y[jmg26]*rate_eval.p_Mg26__He4_Na23 ) +
          ( -rho*Y[jp]*Y[jna24]*rate_eval.p_Na24__n_Mg24 +rho*Y[jn]*Y[jmg24]*rate_eval.n_Mg24__p_Na24 ) +
          ( -rho*Y[jp]*Y[jna24]*rate_eval.p_Na24__He4_Ne21 +rho*Y[jhe4]*Y[jne21]*rate_eval.He4_Ne21__p_Na24 ) +
          ( +rho*Y[jn]*Y[jmg22]*rate_eval.n_Mg22__p_Na22 -rho*Y[jp]*Y[jna22]*rate_eval.p_Na22__n_Mg22 ) +
          ( +rho*Y[jhe4]*Y[jmg22]*rate_eval.He4_Mg22__p_Al25 -rho*Y[jp]*Y[jal25]*rate_eval.p_Al25__He4_Mg22 ) +
          ( +rho*Y[jn]*Y[jmg23]*rate_eval.n_Mg23__p_Na23 -rho*Y[jp]*Y[jna23]*rate_eval.p_Na23__n_Mg23 ) +
          ( +rho*Y[jhe4]*Y[jmg23]*rate_eval.He4_Mg23__p_Al26 -rho*Y[jp]*Y[jal26]*rate_eval.p_Al26__He4_Mg23 ) +
          ( +rho*Y[jn]*Y[jal25]*rate_eval.n_Al25__p_Mg25 -rho*Y[jp]*Y[jmg25]*rate_eval.p_Mg25__n_Al25 ) +
          ( +rho*Y[jhe4]*Y[jal25]*rate_eval.He4_Al25__p_Si28 -rho*Y[jp]*Y[jsi28]*rate_eval.p_Si28__He4_Al25 ) +
          ( +rho*Y[jn]*Y[jal26]*rate_eval.n_Al26__p_Mg26 -rho*Y[jp]*Y[jmg26]*rate_eval.p_Mg26__n_Al26 ) +
          ( +rho*Y[jhe4]*Y[jal26]*rate_eval.He4_Al26__p_Si29 -rho*Y[jp]*Y[jsi29]*rate_eval.p_Si29__He4_Al26 ) +
          ( -rho*Y[jp]*Y[jal27]*rate_eval.p_Al27__He4_Mg24 +rho*Y[jhe4]*Y[jmg24]*rate_eval.He4_Mg24__p_Al27 ) +
          ( +rho*Y[jhe4]*Y[jal27]*rate_eval.He4_Al27__p_Si30 -rho*Y[jp]*Y[jsi30]*rate_eval.p_Si30__He4_Al27 ) +
          ( -rho*Y[jp]*Y[jal28]*rate_eval.p_Al28__n_Si28 +rho*Y[jn]*Y[jsi28]*rate_eval.n_Si28__p_Al28 ) +
          ( -rho*Y[jp]*Y[jal28]*rate_eval.p_Al28__He4_Mg25 +rho*Y[jhe4]*Y[jmg25]*rate_eval.He4_Mg25__p_Al28 ) +
          ( +rho*Y[jn]*Y[jsi26]*rate_eval.n_Si26__p_Al26 -rho*Y[jp]*Y[jal26]*rate_eval.p_Al26__n_Si26 ) +
          ( +rho*Y[jhe4]*Y[jsi26]*rate_eval.He4_Si26__p_P29 -rho*Y[jp]*Y[jp29]*rate_eval.p_P29__He4_Si26 ) +
          ( +rho*Y[jn]*Y[jp29]*rate_eval.n_P29__p_Si29 -rho*Y[jp]*Y[jsi29]*rate_eval.p_Si29__n_P29 ) +
          ( +rho*Y[jhe4]*Y[jp29]*rate_eval.He4_P29__p_S32 -rho*Y[jp]*Y[js32]*rate_eval.p_S32__He4_P29 ) +
          ( +rho*Y[jn]*Y[jp30]*rate_eval.n_P30__p_Si30 -rho*Y[jp]*Y[jsi30]*rate_eval.p_Si30__n_P30 ) +
          ( -rho*Y[jp]*Y[jp31]*rate_eval.p_P31__He4_Si28 +rho*Y[jhe4]*Y[jsi28]*rate_eval.He4_Si28__p_P31 ) +
          ( +rho*Y[jhe4]*Y[jp31]*rate_eval.He4_P31__p_S34 -rho*Y[jp]*Y[js34]*rate_eval.p_S34__He4_P31 ) +
          ( -rho*Y[jp]*Y[jp32]*rate_eval.p_P32__n_S32 +rho*Y[jn]*Y[js32]*rate_eval.n_S32__p_P32 ) +
          ( -rho*Y[jp]*Y[jp32]*rate_eval.p_P32__He4_Si29 +rho*Y[jhe4]*Y[jsi29]*rate_eval.He4_Si29__p_P32 ) +
          ( +rho*Y[jn]*Y[js31]*rate_eval.n_S31__p_P31 -rho*Y[jp]*Y[jp31]*rate_eval.p_P31__n_S31 ) +
          ( +rho*Y[jhe4]*Y[js31]*rate_eval.He4_S31__p_Cl34 -rho*Y[jp]*Y[jcl34]*rate_eval.p_Cl34__He4_S31 ) +
          ( +rho*Y[jn]*Y[jcl32]*rate_eval.n_Cl32__p_S32 -rho*Y[jp]*Y[js32]*rate_eval.p_S32__n_Cl32 ) +
          ( +rho*Y[jhe4]*Y[jcl32]*rate_eval.He4_Cl32__p_Ar35 -rho*Y[jp]*Y[jar35]*rate_eval.p_Ar35__He4_Cl32 ) +
          ( +rho*Y[jn]*Y[jcl34]*rate_eval.n_Cl34__p_S34 -rho*Y[jp]*Y[js34]*rate_eval.p_S34__n_Cl34 ) +
          ( -rho*Y[jp]*Y[jcl35]*rate_eval.p_Cl35__He4_S32 +rho*Y[jhe4]*Y[js32]*rate_eval.He4_S32__p_Cl35 ) +
          ( +rho*Y[jn]*Y[jar35]*rate_eval.n_Ar35__p_Cl35 -rho*Y[jp]*Y[jcl35]*rate_eval.p_Cl35__n_Ar35 ) +
          ( +rho*Y[jt]*Y[jhe3]*rate_eval.t_He3__n_p_He4 -rho**2*Y[jn]*Y[jp]*Y[jhe4]*rate_eval.n_p_He4__t_He3 ) +
          ( + 2*5.00000000000000e-01*rho*Y[jhe3]**2*rate_eval.He3_He3__p_p_He4 + -2*5.00000000000000e-01*rho**2*Y[jp]**2*Y[jhe4]*rate_eval.p_p_He4__He3_He3 ) +
          ( +rho*Y[jd]*Y[jbe7]*rate_eval.d_Be7__p_He4_He4 -5.00000000000000e-01*rho**2*Y[jp]*Y[jhe4]**2*rate_eval.p_He4_He4__d_Be7 ) +
          ( -rho*Y[jp]*Y[jbe9]*rate_eval.p_Be9__d_He4_He4 +5.00000000000000e-01*rho**2*Y[jd]*Y[jhe4]**2*rate_eval.d_He4_He4__p_Be9 ) +
          ( +rho*Y[jn]*Y[jb8]*rate_eval.n_B8__p_He4_He4 -5.00000000000000e-01*rho**2*Y[jp]*Y[jhe4]**2*rate_eval.p_He4_He4__n_B8 ) +
          ( +rho*Y[jhe3]*Y[jli7]*rate_eval.He3_Li7__n_p_He4_He4 -5.00000000000000e-01*rho**3*Y[jn]*Y[jp]*Y[jhe4]**2*rate_eval.n_p_He4_He4__He3_Li7 ) +
          ( +rho*Y[jt]*Y[jbe7]*rate_eval.t_Be7__n_p_He4_He4 -5.00000000000000e-01*rho**3*Y[jn]*Y[jp]*Y[jhe4]**2*rate_eval.n_p_He4_He4__t_Be7 ) +
          ( + 2*rho*Y[jhe3]*Y[jbe7]*rate_eval.He3_Be7__p_p_He4_He4 + -2*2.50000000000000e-01*rho**3*Y[jp]**2*Y[jhe4]**2*rate_eval.p_p_He4_He4__He3_Be7 ) +
          ( -rho**2*Y[jn]*Y[jp]*Y[jhe4]*rate_eval.n_p_He4__Li6 +Y[jli6]*rate_eval.Li6__n_p_He4 ) +
          ( -5.00000000000000e-01*rho**2*Y[jn]*Y[jp]**2*rate_eval.n_p_p__p_d +rho*Y[jp]*Y[jd]*rate_eval.p_d__n_p_p ) +
          ( + 0*5.00000000000000e-01*rho**3*Y[jn]*Y[jp]*Y[jhe4]**2*rate_eval.n_p_He4_He4__p_Be9 + 0*rho*Y[jp]*Y[jbe9]*rate_eval.p_Be9__n_p_He4_He4 )
       )

    dYdt[jd] = (
          ( +rho*Y[jn]*Y[jp]*rate_eval.n_p__d -Y[jd]*rate_eval.d__n_p ) +
          +5.00000000000000e-01*rho*Y[jp]**2*rate_eval.p_p__d__weak__bet_pos_  +
          +5.00000000000000e-01*rho**2*ye(Y)*Y[jp]**2*rate_eval.p_p__d__weak__electron_capture  +
          ( -rho*Y[jn]*Y[jd]*rate_eval.n_d__t +Y[jt]*rate_eval.t__n_d ) +
          ( -rho*Y[jp]*Y[jd]*rate_eval.p_d__He3 +Y[jhe3]*rate_eval.He3__p_d ) +
          ( + -2*5.00000000000000e-01*rho*Y[jd]**2*rate_eval.d_d__He4 + 2*Y[jhe4]*rate_eval.He4__d_d ) +
          ( -rho*Y[jd]*Y[jhe4]*rate_eval.He4_d__Li6 +Y[jli6]*rate_eval.Li6__He4_d ) +
          ( + -2*5.00000000000000e-01*rho*Y[jd]**2*rate_eval.d_d__n_He3 + 2*rho*Y[jn]*Y[jhe3]*rate_eval.n_He3__d_d ) +
          ( + -2*5.00000000000000e-01*rho*Y[jd]**2*rate_eval.d_d__p_t + 2*rho*Y[jp]*Y[jt]*rate_eval.p_t__d_d ) +
          ( -rho*Y[jd]*Y[jt]*rate_eval.d_t__n_He4 +rho*Y[jn]*Y[jhe4]*rate_eval.n_He4__d_t ) +
          ( -rho*Y[jd]*Y[jhe3]*rate_eval.d_He3__p_He4 +rho*Y[jp]*Y[jhe4]*rate_eval.p_He4__d_He3 ) +
          ( +rho*Y[jt]*Y[jhe3]*rate_eval.t_He3__d_He4 -rho*Y[jd]*Y[jhe4]*rate_eval.d_He4__t_He3 ) +
          ( -rho*Y[jd]*Y[jli6]*rate_eval.d_Li6__n_Be7 +rho*Y[jn]*Y[jbe7]*rate_eval.n_Be7__d_Li6 ) +
          ( -rho*Y[jd]*Y[jli6]*rate_eval.d_Li6__p_Li7 +rho*Y[jp]*Y[jli7]*rate_eval.p_Li7__d_Li6 ) +
          ( -rho*Y[jd]*Y[jc13]*rate_eval.d_C13__n_N14 +rho*Y[jn]*Y[jn14]*rate_eval.n_N14__d_C13 ) +
          ( -rho*Y[jd]*Y[jc14]*rate_eval.d_C14__n_N15 +rho*Y[jn]*Y[jn15]*rate_eval.n_N15__d_C14 ) +
          ( -rho*Y[jd]*Y[jli7]*rate_eval.d_Li7__n_He4_He4 +5.00000000000000e-01*rho**2*Y[jn]*Y[jhe4]**2*rate_eval.n_He4_He4__d_Li7 ) +
          ( -rho*Y[jd]*Y[jbe7]*rate_eval.d_Be7__p_He4_He4 +5.00000000000000e-01*rho**2*Y[jp]*Y[jhe4]**2*rate_eval.p_He4_He4__d_Be7 ) +
          ( +rho*Y[jp]*Y[jbe9]*rate_eval.p_Be9__d_He4_He4 -5.00000000000000e-01*rho**2*Y[jd]*Y[jhe4]**2*rate_eval.d_He4_He4__p_Be9 ) +
          ( +5.00000000000000e-01*rho**2*Y[jn]*Y[jp]**2*rate_eval.n_p_p__p_d -rho*Y[jp]*Y[jd]*rate_eval.p_d__n_p_p )
       )

    dYdt[jt] = (
          ( -Y[jt]*rate_eval.t__He3__weak__wc12 +rho*ye(Y)*Y[jhe3]*rate_eval.He3__t__weak__electron_capture ) +
          ( +rho*Y[jn]*Y[jd]*rate_eval.n_d__t -Y[jt]*rate_eval.t__n_d ) +
          ( -rho*Y[jp]*Y[jt]*rate_eval.p_t__He4 +Y[jhe4]*rate_eval.He4__p_t ) +
          ( -rho*Y[jt]*Y[jhe4]*rate_eval.He4_t__Li7 +Y[jli7]*rate_eval.Li7__He4_t ) +
          ( +5.00000000000000e-01*rho*Y[jd]**2*rate_eval.d_d__p_t -rho*Y[jp]*Y[jt]*rate_eval.p_t__d_d ) +
          ( -rho*Y[jd]*Y[jt]*rate_eval.d_t__n_He4 +rho*Y[jn]*Y[jhe4]*rate_eval.n_He4__d_t ) +
          ( +rho*Y[jn]*Y[jhe3]*rate_eval.n_He3__p_t -rho*Y[jp]*Y[jt]*rate_eval.p_t__n_He3 ) +
          ( -rho*Y[jt]*Y[jhe3]*rate_eval.t_He3__d_He4 +rho*Y[jd]*Y[jhe4]*rate_eval.d_He4__t_He3 ) +
          ( +rho*Y[jn]*Y[jli6]*rate_eval.n_Li6__He4_t -rho*Y[jt]*Y[jhe4]*rate_eval.He4_t__n_Li6 ) +
          ( -rho*Y[jt]*Y[jli7]*rate_eval.t_Li7__n_Be9 +rho*Y[jn]*Y[jbe9]*rate_eval.n_Be9__t_Li7 ) +
          ( + -2*5.00000000000000e-01*rho*Y[jt]**2*rate_eval.t_t__n_n_He4 + 2*5.00000000000000e-01*rho**2*Y[jn]**2*Y[jhe4]*rate_eval.n_n_He4__t_t ) +
          ( -rho*Y[jt]*Y[jhe3]*rate_eval.t_He3__n_p_He4 +rho**2*Y[jn]*Y[jp]*Y[jhe4]*rate_eval.n_p_He4__t_He3 ) +
          ( -rho*Y[jt]*Y[jli7]*rate_eval.t_Li7__n_n_He4_He4 +2.50000000000000e-01*rho**3*Y[jn]**2*Y[jhe4]**2*rate_eval.n_n_He4_He4__t_Li7 ) +
          ( -rho*Y[jt]*Y[jbe7]*rate_eval.t_Be7__n_p_He4_He4 +5.00000000000000e-01*rho**3*Y[jn]*Y[jp]*Y[jhe4]**2*rate_eval.n_p_He4_He4__t_Be7 )
       )

    dYdt[jhe3] = (
          ( +Y[jt]*rate_eval.t__He3__weak__wc12 -rho*ye(Y)*Y[jhe3]*rate_eval.He3__t__weak__electron_capture ) +
          ( +rho*Y[jp]*Y[jd]*rate_eval.p_d__He3 -Y[jhe3]*rate_eval.He3__p_d ) +
          ( -rho*Y[jn]*Y[jhe3]*rate_eval.n_He3__He4 +Y[jhe4]*rate_eval.He4__n_He3 ) +
          -rho*Y[jp]*Y[jhe3]*rate_eval.p_He3__He4__weak__bet_pos_  +
          ( -rho*Y[jhe3]*Y[jhe4]*rate_eval.He4_He3__Be7 +Y[jbe7]*rate_eval.Be7__He4_He3 ) +
          ( +5.00000000000000e-01*rho*Y[jd]**2*rate_eval.d_d__n_He3 -rho*Y[jn]*Y[jhe3]*rate_eval.n_He3__d_d ) +
          ( -rho*Y[jn]*Y[jhe3]*rate_eval.n_He3__p_t +rho*Y[jp]*Y[jt]*rate_eval.p_t__n_He3 ) +
          ( -rho*Y[jd]*Y[jhe3]*rate_eval.d_He3__p_He4 +rho*Y[jp]*Y[jhe4]*rate_eval.p_He4__d_He3 ) +
          ( -rho*Y[jt]*Y[jhe3]*rate_eval.t_He3__d_He4 +rho*Y[jd]*Y[jhe4]*rate_eval.d_He4__t_He3 ) +
          ( +rho*Y[jp]*Y[jli6]*rate_eval.p_Li6__He4_He3 -rho*Y[jhe3]*Y[jhe4]*rate_eval.He4_He3__p_Li6 ) +
          ( -rho*Y[jt]*Y[jhe3]*rate_eval.t_He3__n_p_He4 +rho**2*Y[jn]*Y[jp]*Y[jhe4]*rate_eval.n_p_He4__t_He3 ) +
          ( + -2*5.00000000000000e-01*rho*Y[jhe3]**2*rate_eval.He3_He3__p_p_He4 + 2*5.00000000000000e-01*rho**2*Y[jp]**2*Y[jhe4]*rate_eval.p_p_He4__He3_He3 ) +
          ( -rho*Y[jhe3]*Y[jli7]*rate_eval.He3_Li7__n_p_He4_He4 +5.00000000000000e-01*rho**3*Y[jn]*Y[jp]*Y[jhe4]**2*rate_eval.n_p_He4_He4__He3_Li7 ) +
          ( -rho*Y[jhe3]*Y[jbe7]*rate_eval.He3_Be7__p_p_He4_He4 +2.50000000000000e-01*rho**3*Y[jp]**2*Y[jhe4]**2*rate_eval.p_p_He4_He4__He3_Be7 )
       )

    dYdt[jhe4] = (
          + 2*Y[jb8]*rate_eval.B8__He4_He4__weak__wc12  +
          +Y[jcl32]*rate_eval.Cl32__He4_Si28__weak__wc12  +
          ( +5.00000000000000e-01*rho*Y[jd]**2*rate_eval.d_d__He4 -Y[jhe4]*rate_eval.He4__d_d ) +
          ( -rho*Y[jd]*Y[jhe4]*rate_eval.He4_d__Li6 +Y[jli6]*rate_eval.Li6__He4_d ) +
          ( +rho*Y[jp]*Y[jt]*rate_eval.p_t__He4 -Y[jhe4]*rate_eval.He4__p_t ) +
          ( -rho*Y[jt]*Y[jhe4]*rate_eval.He4_t__Li7 +Y[jli7]*rate_eval.Li7__He4_t ) +
          ( +rho*Y[jn]*Y[jhe3]*rate_eval.n_He3__He4 -Y[jhe4]*rate_eval.He4__n_He3 ) +
          +rho*Y[jp]*Y[jhe3]*rate_eval.p_He3__He4__weak__bet_pos_  +
          ( -rho*Y[jhe3]*Y[jhe4]*rate_eval.He4_He3__Be7 +Y[jbe7]*rate_eval.Be7__He4_He3 ) +
          ( -rho*Y[jhe4]*Y[jbe7]*rate_eval.He4_Be7__C11 +Y[jc11]*rate_eval.C11__He4_Be7 ) +
          ( -rho*Y[jhe4]*Y[jc12]*rate_eval.He4_C12__O16 +Y[jo16]*rate_eval.O16__He4_C12 ) +
          ( -rho*Y[jhe4]*Y[jc14]*rate_eval.He4_C14__O18 +Y[jo18]*rate_eval.O18__He4_C14 ) +
          ( -rho*Y[jhe4]*Y[jn14]*rate_eval.He4_N14__F18 +Y[jf18]*rate_eval.F18__He4_N14 ) +
          ( -rho*Y[jhe4]*Y[jn15]*rate_eval.He4_N15__F19 +Y[jf19]*rate_eval.F19__He4_N15 ) +
          ( -rho*Y[jhe4]*Y[jo14]*rate_eval.He4_O14__Ne18 +Y[jne18]*rate_eval.Ne18__He4_O14 ) +
          ( -rho*Y[jhe4]*Y[jo15]*rate_eval.He4_O15__Ne19 +Y[jne19]*rate_eval.Ne19__He4_O15 ) +
          ( -rho*Y[jhe4]*Y[jo16]*rate_eval.He4_O16__Ne20 +Y[jne20]*rate_eval.Ne20__He4_O16 ) +
          ( -rho*Y[jhe4]*Y[jo17]*rate_eval.He4_O17__Ne21 +Y[jne21]*rate_eval.Ne21__He4_O17 ) +
          ( -rho*Y[jhe4]*Y[jf17]*rate_eval.He4_F17__Na21 +Y[jna21]*rate_eval.Na21__He4_F17 ) +
          ( -rho*Y[jhe4]*Y[jf18]*rate_eval.He4_F18__Na22 +Y[jna22]*rate_eval.Na22__He4_F18 ) +
          ( -rho*Y[jhe4]*Y[jf19]*rate_eval.He4_F19__Na23 +Y[jna23]*rate_eval.Na23__He4_F19 ) +
          ( -rho*Y[jhe4]*Y[jne18]*rate_eval.He4_Ne18__Mg22 +Y[jmg22]*rate_eval.Mg22__He4_Ne18 ) +
          ( -rho*Y[jhe4]*Y[jne19]*rate_eval.He4_Ne19__Mg23 +Y[jmg23]*rate_eval.Mg23__He4_Ne19 ) +
          ( -rho*Y[jhe4]*Y[jne20]*rate_eval.He4_Ne20__Mg24 +Y[jmg24]*rate_eval.Mg24__He4_Ne20 ) +
          ( -rho*Y[jhe4]*Y[jne21]*rate_eval.He4_Ne21__Mg25 +Y[jmg25]*rate_eval.Mg25__He4_Ne21 ) +
          ( -rho*Y[jhe4]*Y[jna21]*rate_eval.He4_Na21__Al25 +Y[jal25]*rate_eval.Al25__He4_Na21 ) +
          ( -rho*Y[jhe4]*Y[jna22]*rate_eval.He4_Na22__Al26 +Y[jal26]*rate_eval.Al26__He4_Na22 ) +
          ( -rho*Y[jhe4]*Y[jna23]*rate_eval.He4_Na23__Al27 +Y[jal27]*rate_eval.Al27__He4_Na23 ) +
          ( -rho*Y[jhe4]*Y[jna24]*rate_eval.He4_Na24__Al28 +Y[jal28]*rate_eval.Al28__He4_Na24 ) +
          ( -rho*Y[jhe4]*Y[jmg22]*rate_eval.He4_Mg22__Si26 +Y[jsi26]*rate_eval.Si26__He4_Mg22 ) +
          ( -rho*Y[jhe4]*Y[jmg24]*rate_eval.He4_Mg24__Si28 +Y[jsi28]*rate_eval.Si28__He4_Mg24 ) +
          ( -rho*Y[jhe4]*Y[jmg25]*rate_eval.He4_Mg25__Si29 +Y[jsi29]*rate_eval.Si29__He4_Mg25 ) +
          ( -rho*Y[jhe4]*Y[jmg26]*rate_eval.He4_Mg26__Si30 +Y[jsi30]*rate_eval.Si30__He4_Mg26 ) +
          ( -rho*Y[jhe4]*Y[jal25]*rate_eval.He4_Al25__P29 +Y[jp29]*rate_eval.P29__He4_Al25 ) +
          ( -rho*Y[jhe4]*Y[jal26]*rate_eval.He4_Al26__P30 +Y[jp30]*rate_eval.P30__He4_Al26 ) +
          ( -rho*Y[jhe4]*Y[jal27]*rate_eval.He4_Al27__P31 +Y[jp31]*rate_eval.P31__He4_Al27 ) +
          ( -rho*Y[jhe4]*Y[jal28]*rate_eval.He4_Al28__P32 +Y[jp32]*rate_eval.P32__He4_Al28 ) +
          ( -rho*Y[jhe4]*Y[jsi28]*rate_eval.He4_Si28__S32 +Y[js32]*rate_eval.S32__He4_Si28 ) +
          ( -rho*Y[jhe4]*Y[jsi30]*rate_eval.He4_Si30__S34 +Y[js34]*rate_eval.S34__He4_Si30 ) +
          ( -rho*Y[jhe4]*Y[jp30]*rate_eval.He4_P30__Cl34 +Y[jcl34]*rate_eval.Cl34__He4_P30 ) +
          ( -rho*Y[jhe4]*Y[jp31]*rate_eval.He4_P31__Cl35 +Y[jcl35]*rate_eval.Cl35__He4_P31 ) +
          ( -rho*Y[jhe4]*Y[js31]*rate_eval.He4_S31__Ar35 +Y[jar35]*rate_eval.Ar35__He4_S31 ) +
          ( +rho*Y[jd]*Y[jt]*rate_eval.d_t__n_He4 -rho*Y[jn]*Y[jhe4]*rate_eval.n_He4__d_t ) +
          ( +rho*Y[jd]*Y[jhe3]*rate_eval.d_He3__p_He4 -rho*Y[jp]*Y[jhe4]*rate_eval.p_He4__d_He3 ) +
          ( +rho*Y[jt]*Y[jhe3]*rate_eval.t_He3__d_He4 -rho*Y[jd]*Y[jhe4]*rate_eval.d_He4__t_He3 ) +
          ( +rho*Y[jn]*Y[jli6]*rate_eval.n_Li6__He4_t -rho*Y[jt]*Y[jhe4]*rate_eval.He4_t__n_Li6 ) +
          ( +rho*Y[jp]*Y[jli6]*rate_eval.p_Li6__He4_He3 -rho*Y[jhe3]*Y[jhe4]*rate_eval.He4_He3__p_Li6 ) +
          ( + 2*rho*Y[jp]*Y[jli7]*rate_eval.p_Li7__He4_He4 + -2*5.00000000000000e-01*rho*Y[jhe4]**2*rate_eval.He4_He4__p_Li7 ) +
          ( + 2*rho*Y[jn]*Y[jbe7]*rate_eval.n_Be7__He4_He4 + -2*5.00000000000000e-01*rho*Y[jhe4]**2*rate_eval.He4_He4__n_Be7 ) +
          ( +rho*Y[jp]*Y[jbe9]*rate_eval.p_Be9__He4_Li6 -rho*Y[jhe4]*Y[jli6]*rate_eval.He4_Li6__p_Be9 ) +
          ( -rho*Y[jhe4]*Y[jbe9]*rate_eval.He4_Be9__n_C12 +rho*Y[jn]*Y[jc12]*rate_eval.n_C12__He4_Be9 ) +
          ( -rho*Y[jhe4]*Y[jb8]*rate_eval.He4_B8__p_C11 +rho*Y[jp]*Y[jc11]*rate_eval.p_C11__He4_B8 ) +
          ( -rho*Y[jhe4]*Y[jc11]*rate_eval.He4_C11__p_N14 +rho*Y[jp]*Y[jn14]*rate_eval.p_N14__He4_C11 ) +
          ( +5.00000000000000e-01*rho*Y[jc12]**2*rate_eval.C12_C12__He4_Ne20 -rho*Y[jhe4]*Y[jne20]*rate_eval.He4_Ne20__C12_C12 ) +
          ( -rho*Y[jhe4]*Y[jc13]*rate_eval.He4_C13__n_O16 +rho*Y[jn]*Y[jo16]*rate_eval.n_O16__He4_C13 ) +
          ( -rho*Y[jhe4]*Y[jn13]*rate_eval.He4_N13__p_O16 +rho*Y[jp]*Y[jo16]*rate_eval.p_O16__He4_N13 ) +
          ( +rho*Y[jp]*Y[jn15]*rate_eval.p_N15__He4_C12 -rho*Y[jhe4]*Y[jc12]*rate_eval.He4_C12__p_N15 ) +
          ( +rho*Y[jn]*Y[jo14]*rate_eval.n_O14__He4_C11 -rho*Y[jhe4]*Y[jc11]*rate_eval.He4_C11__n_O14 ) +
          ( -rho*Y[jhe4]*Y[jo14]*rate_eval.He4_O14__p_F17 +rho*Y[jp]*Y[jf17]*rate_eval.p_F17__He4_O14 ) +
          ( +rho*Y[jn]*Y[jo15]*rate_eval.n_O15__He4_C12 -rho*Y[jhe4]*Y[jc12]*rate_eval.He4_C12__n_O15 ) +
          ( +rho*Y[jc12]*Y[jo16]*rate_eval.C12_O16__He4_Mg24 -rho*Y[jhe4]*Y[jmg24]*rate_eval.He4_Mg24__C12_O16 ) +
          ( +5.00000000000000e-01*rho*Y[jo16]**2*rate_eval.O16_O16__He4_Si28 -rho*Y[jhe4]*Y[jsi28]*rate_eval.He4_Si28__O16_O16 ) +
          ( +rho*Y[jn]*Y[jo17]*rate_eval.n_O17__He4_C14 -rho*Y[jhe4]*Y[jc14]*rate_eval.He4_C14__n_O17 ) +
          ( +rho*Y[jp]*Y[jo17]*rate_eval.p_O17__He4_N14 -rho*Y[jhe4]*Y[jn14]*rate_eval.He4_N14__p_O17 ) +
          ( -rho*Y[jhe4]*Y[jo17]*rate_eval.He4_O17__n_Ne20 +rho*Y[jn]*Y[jne20]*rate_eval.n_Ne20__He4_O17 ) +
          ( +rho*Y[jp]*Y[jo18]*rate_eval.p_O18__He4_N15 -rho*Y[jhe4]*Y[jn15]*rate_eval.He4_N15__p_O18 ) +
          ( -rho*Y[jhe4]*Y[jf15]*rate_eval.He4_F15__p_Ne18 +rho*Y[jp]*Y[jne18]*rate_eval.p_Ne18__He4_F15 ) +
          ( +rho*Y[jn]*Y[jf16]*rate_eval.n_F16__He4_N13 -rho*Y[jhe4]*Y[jn13]*rate_eval.He4_N13__n_F16 ) +
          ( -rho*Y[jhe4]*Y[jf16]*rate_eval.He4_F16__p_Ne19 +rho*Y[jp]*Y[jne19]*rate_eval.p_Ne19__He4_F16 ) +
          ( +rho*Y[jn]*Y[jf17]*rate_eval.n_F17__He4_N14 -rho*Y[jhe4]*Y[jn14]*rate_eval.He4_N14__n_F17 ) +
          ( -rho*Y[jhe4]*Y[jf17]*rate_eval.He4_F17__p_Ne20 +rho*Y[jp]*Y[jne20]*rate_eval.p_Ne20__He4_F17 ) +
          ( +rho*Y[jn]*Y[jf18]*rate_eval.n_F18__He4_N15 -rho*Y[jhe4]*Y[jn15]*rate_eval.He4_N15__n_F18 ) +
          ( +rho*Y[jp]*Y[jf18]*rate_eval.p_F18__He4_O15 -rho*Y[jhe4]*Y[jo15]*rate_eval.He4_O15__p_F18 ) +
          ( -rho*Y[jhe4]*Y[jf18]*rate_eval.He4_F18__p_Ne21 +rho*Y[jp]*Y[jne21]*rate_eval.p_Ne21__He4_F18 ) +
          ( +rho*Y[jp]*Y[jf19]*rate_eval.p_F19__He4_O16 -rho*Y[jhe4]*Y[jo16]*rate_eval.He4_O16__p_F19 ) +
          ( +rho*Y[jn]*Y[jne18]*rate_eval.n_Ne18__He4_O15 -rho*Y[jhe4]*Y[jo15]*rate_eval.He4_O15__n_Ne18 ) +
          ( -rho*Y[jhe4]*Y[jne18]*rate_eval.He4_Ne18__p_Na21 +rho*Y[jp]*Y[jna21]*rate_eval.p_Na21__He4_Ne18 ) +
          ( +rho*Y[jn]*Y[jne19]*rate_eval.n_Ne19__He4_O16 -rho*Y[jhe4]*Y[jo16]*rate_eval.He4_O16__n_Ne19 ) +
          ( -rho*Y[jhe4]*Y[jne19]*rate_eval.He4_Ne19__p_Na22 +rho*Y[jp]*Y[jna22]*rate_eval.p_Na22__He4_Ne19 ) +
          ( +rho*Y[jc12]*Y[jne20]*rate_eval.C12_Ne20__He4_Si28 -rho*Y[jhe4]*Y[jsi28]*rate_eval.He4_Si28__C12_Ne20 ) +
          ( +rho*Y[jn]*Y[jne21]*rate_eval.n_Ne21__He4_O18 -rho*Y[jhe4]*Y[jo18]*rate_eval.He4_O18__n_Ne21 ) +
          ( -rho*Y[jhe4]*Y[jne21]*rate_eval.He4_Ne21__n_Mg24 +rho*Y[jn]*Y[jmg24]*rate_eval.n_Mg24__He4_Ne21 ) +
          ( +rho*Y[jn]*Y[jna21]*rate_eval.n_Na21__He4_F18 -rho*Y[jhe4]*Y[jf18]*rate_eval.He4_F18__n_Na21 ) +
          ( -rho*Y[jhe4]*Y[jna21]*rate_eval.He4_Na21__p_Mg24 +rho*Y[jp]*Y[jmg24]*rate_eval.p_Mg24__He4_Na21 ) +
          ( +rho*Y[jn]*Y[jna22]*rate_eval.n_Na22__He4_F19 -rho*Y[jhe4]*Y[jf19]*rate_eval.He4_F19__n_Na22 ) +
          ( -rho*Y[jhe4]*Y[jna22]*rate_eval.He4_Na22__p_Mg25 +rho*Y[jp]*Y[jmg25]*rate_eval.p_Mg25__He4_Na22 ) +
          ( +rho*Y[jp]*Y[jna23]*rate_eval.p_Na23__He4_Ne20 -rho*Y[jhe4]*Y[jne20]*rate_eval.He4_Ne20__p_Na23 ) +
          ( -rho*Y[jhe4]*Y[jna23]*rate_eval.He4_Na23__p_Mg26 +rho*Y[jp]*Y[jmg26]*rate_eval.p_Mg26__He4_Na23 ) +
          ( +rho*Y[jp]*Y[jna24]*rate_eval.p_Na24__He4_Ne21 -rho*Y[jhe4]*Y[jne21]*rate_eval.He4_Ne21__p_Na24 ) +
          ( -rho*Y[jhe4]*Y[jna24]*rate_eval.He4_Na24__n_Al27 +rho*Y[jn]*Y[jal27]*rate_eval.n_Al27__He4_Na24 ) +
          ( +rho*Y[jn]*Y[jmg22]*rate_eval.n_Mg22__He4_Ne19 -rho*Y[jhe4]*Y[jne19]*rate_eval.He4_Ne19__n_Mg22 ) +
          ( -rho*Y[jhe4]*Y[jmg22]*rate_eval.He4_Mg22__p_Al25 +rho*Y[jp]*Y[jal25]*rate_eval.p_Al25__He4_Mg22 ) +
          ( +rho*Y[jn]*Y[jmg23]*rate_eval.n_Mg23__He4_Ne20 -rho*Y[jhe4]*Y[jne20]*rate_eval.He4_Ne20__n_Mg23 ) +
          ( -rho*Y[jhe4]*Y[jmg23]*rate_eval.He4_Mg23__p_Al26 +rho*Y[jp]*Y[jal26]*rate_eval.p_Al26__He4_Mg23 ) +
          ( -rho*Y[jhe4]*Y[jmg25]*rate_eval.He4_Mg25__n_Si28 +rho*Y[jn]*Y[jsi28]*rate_eval.n_Si28__He4_Mg25 ) +
          ( -rho*Y[jhe4]*Y[jmg26]*rate_eval.He4_Mg26__n_Si29 +rho*Y[jn]*Y[jsi29]*rate_eval.n_Si29__He4_Mg26 ) +
          ( +rho*Y[jn]*Y[jal25]*rate_eval.n_Al25__He4_Na22 -rho*Y[jhe4]*Y[jna22]*rate_eval.He4_Na22__n_Al25 ) +
          ( -rho*Y[jhe4]*Y[jal25]*rate_eval.He4_Al25__p_Si28 +rho*Y[jp]*Y[jsi28]*rate_eval.p_Si28__He4_Al25 ) +
          ( +rho*Y[jn]*Y[jal26]*rate_eval.n_Al26__He4_Na23 -rho*Y[jhe4]*Y[jna23]*rate_eval.He4_Na23__n_Al26 ) +
          ( -rho*Y[jhe4]*Y[jal26]*rate_eval.He4_Al26__p_Si29 +rho*Y[jp]*Y[jsi29]*rate_eval.p_Si29__He4_Al26 ) +
          ( +rho*Y[jp]*Y[jal27]*rate_eval.p_Al27__He4_Mg24 -rho*Y[jhe4]*Y[jmg24]*rate_eval.He4_Mg24__p_Al27 ) +
          ( -rho*Y[jhe4]*Y[jal27]*rate_eval.He4_Al27__p_Si30 +rho*Y[jp]*Y[jsi30]*rate_eval.p_Si30__He4_Al27 ) +
          ( +rho*Y[jp]*Y[jal28]*rate_eval.p_Al28__He4_Mg25 -rho*Y[jhe4]*Y[jmg25]*rate_eval.He4_Mg25__p_Al28 ) +
          ( -rho*Y[jhe4]*Y[jal28]*rate_eval.He4_Al28__n_P31 +rho*Y[jn]*Y[jp31]*rate_eval.n_P31__He4_Al28 ) +
          ( +rho*Y[jn]*Y[jsi26]*rate_eval.n_Si26__He4_Mg23 -rho*Y[jhe4]*Y[jmg23]*rate_eval.He4_Mg23__n_Si26 ) +
          ( -rho*Y[jhe4]*Y[jsi26]*rate_eval.He4_Si26__p_P29 +rho*Y[jp]*Y[jp29]*rate_eval.p_P29__He4_Si26 ) +
          ( +rho*Y[jn]*Y[jp29]*rate_eval.n_P29__He4_Al26 -rho*Y[jhe4]*Y[jal26]*rate_eval.He4_Al26__n_P29 ) +
          ( -rho*Y[jhe4]*Y[jp29]*rate_eval.He4_P29__p_S32 +rho*Y[jp]*Y[js32]*rate_eval.p_S32__He4_P29 ) +
          ( +rho*Y[jn]*Y[jp30]*rate_eval.n_P30__He4_Al27 -rho*Y[jhe4]*Y[jal27]*rate_eval.He4_Al27__n_P30 ) +
          ( +rho*Y[jp]*Y[jp31]*rate_eval.p_P31__He4_Si28 -rho*Y[jhe4]*Y[jsi28]*rate_eval.He4_Si28__p_P31 ) +
          ( -rho*Y[jhe4]*Y[jp31]*rate_eval.He4_P31__p_S34 +rho*Y[jp]*Y[js34]*rate_eval.p_S34__He4_P31 ) +
          ( +rho*Y[jp]*Y[jp32]*rate_eval.p_P32__He4_Si29 -rho*Y[jhe4]*Y[jsi29]*rate_eval.He4_Si29__p_P32 ) +
          ( +rho*Y[jn]*Y[js31]*rate_eval.n_S31__He4_Si28 -rho*Y[jhe4]*Y[jsi28]*rate_eval.He4_Si28__n_S31 ) +
          ( -rho*Y[jhe4]*Y[js31]*rate_eval.He4_S31__p_Cl34 +rho*Y[jp]*Y[jcl34]*rate_eval.p_Cl34__He4_S31 ) +
          ( +rho*Y[jn]*Y[js32]*rate_eval.n_S32__He4_Si29 -rho*Y[jhe4]*Y[jsi29]*rate_eval.He4_Si29__n_S32 ) +
          ( +rho*Y[jn]*Y[jcl32]*rate_eval.n_Cl32__He4_P29 -rho*Y[jhe4]*Y[jp29]*rate_eval.He4_P29__n_Cl32 ) +
          ( -rho*Y[jhe4]*Y[jcl32]*rate_eval.He4_Cl32__p_Ar35 +rho*Y[jp]*Y[jar35]*rate_eval.p_Ar35__He4_Cl32 ) +
          ( +rho*Y[jn]*Y[jcl34]*rate_eval.n_Cl34__He4_P31 -rho*Y[jhe4]*Y[jp31]*rate_eval.He4_P31__n_Cl34 ) +
          ( +rho*Y[jn]*Y[jcl35]*rate_eval.n_Cl35__He4_P32 -rho*Y[jhe4]*Y[jp32]*rate_eval.He4_P32__n_Cl35 ) +
          ( +rho*Y[jp]*Y[jcl35]*rate_eval.p_Cl35__He4_S32 -rho*Y[jhe4]*Y[js32]*rate_eval.He4_S32__p_Cl35 ) +
          ( +rho*Y[jn]*Y[jar35]*rate_eval.n_Ar35__He4_S32 -rho*Y[jhe4]*Y[js32]*rate_eval.He4_S32__n_Ar35 ) +
          ( +5.00000000000000e-01*rho*Y[jt]**2*rate_eval.t_t__n_n_He4 -5.00000000000000e-01*rho**2*Y[jn]**2*Y[jhe4]*rate_eval.n_n_He4__t_t ) +
          ( +rho*Y[jt]*Y[jhe3]*rate_eval.t_He3__n_p_He4 -rho**2*Y[jn]*Y[jp]*Y[jhe4]*rate_eval.n_p_He4__t_He3 ) +
          ( +5.00000000000000e-01*rho*Y[jhe3]**2*rate_eval.He3_He3__p_p_He4 -5.00000000000000e-01*rho**2*Y[jp]**2*Y[jhe4]*rate_eval.p_p_He4__He3_He3 ) +
          ( + 2*rho*Y[jd]*Y[jli7]*rate_eval.d_Li7__n_He4_He4 + -2*5.00000000000000e-01*rho**2*Y[jn]*Y[jhe4]**2*rate_eval.n_He4_He4__d_Li7 ) +
          ( + 2*rho*Y[jd]*Y[jbe7]*rate_eval.d_Be7__p_He4_He4 + -2*5.00000000000000e-01*rho**2*Y[jp]*Y[jhe4]**2*rate_eval.p_He4_He4__d_Be7 ) +
          ( + 2*rho*Y[jp]*Y[jbe9]*rate_eval.p_Be9__d_He4_He4 + -2*5.00000000000000e-01*rho**2*Y[jd]*Y[jhe4]**2*rate_eval.d_He4_He4__p_Be9 ) +
          ( + 2*rho*Y[jn]*Y[jb8]*rate_eval.n_B8__p_He4_He4 + -2*5.00000000000000e-01*rho**2*Y[jp]*Y[jhe4]**2*rate_eval.p_He4_He4__n_B8 ) +
          ( + 3*rho*Y[jn]*Y[jc11]*rate_eval.n_C11__He4_He4_He4 + -3*1.66666666666667e-01*rho**2*Y[jhe4]**3*rate_eval.He4_He4_He4__n_C11 ) +
          ( + 2*rho*Y[jt]*Y[jli7]*rate_eval.t_Li7__n_n_He4_He4 + -2*2.50000000000000e-01*rho**3*Y[jn]**2*Y[jhe4]**2*rate_eval.n_n_He4_He4__t_Li7 ) +
          ( + 2*rho*Y[jhe3]*Y[jli7]*rate_eval.He3_Li7__n_p_He4_He4 + -2*5.00000000000000e-01*rho**3*Y[jn]*Y[jp]*Y[jhe4]**2*rate_eval.n_p_He4_He4__He3_Li7 ) +
          ( + 2*rho*Y[jt]*Y[jbe7]*rate_eval.t_Be7__n_p_He4_He4 + -2*5.00000000000000e-01*rho**3*Y[jn]*Y[jp]*Y[jhe4]**2*rate_eval.n_p_He4_He4__t_Be7 ) +
          ( + 2*rho*Y[jhe3]*Y[jbe7]*rate_eval.He3_Be7__p_p_He4_He4 + -2*2.50000000000000e-01*rho**3*Y[jp]**2*Y[jhe4]**2*rate_eval.p_p_He4_He4__He3_Be7 ) +
          ( -5.00000000000000e-01*rho**2*Y[jn]**2*Y[jhe4]*rate_eval.n_n_He4__He6 +Y[jhe6]*rate_eval.He6__n_n_He4 ) +
          ( -rho**2*Y[jn]*Y[jp]*Y[jhe4]*rate_eval.n_p_He4__Li6 +Y[jli6]*rate_eval.Li6__n_p_He4 ) +
          ( + -2*5.00000000000000e-01*rho**2*Y[jn]*Y[jhe4]**2*rate_eval.n_He4_He4__Be9 + 2*Y[jbe9]*rate_eval.Be9__n_He4_He4 ) +
          ( + -3*1.66666666666667e-01*rho**2*Y[jhe4]**3*rate_eval.He4_He4_He4__C12 + 3*Y[jc12]*rate_eval.C12__He4_He4_He4 ) +
          ( + -2*5.00000000000000e-01*rho**3*Y[jn]*Y[jp]*Y[jhe4]**2*rate_eval.n_p_He4_He4__p_Be9 + 2*rho*Y[jp]*Y[jbe9]*rate_eval.p_Be9__n_p_He4_He4 )
       )

    dYdt[jhe6] = (
          -Y[jhe6]*rate_eval.He6__Li6__weak__wc12  +
          ( +5.00000000000000e-01*rho**2*Y[jn]**2*Y[jhe4]*rate_eval.n_n_He4__He6 -Y[jhe6]*rate_eval.He6__n_n_He4 )
       )

    dYdt[jli6] = (
          +Y[jhe6]*rate_eval.He6__Li6__weak__wc12  +
          ( +rho*Y[jd]*Y[jhe4]*rate_eval.He4_d__Li6 -Y[jli6]*rate_eval.Li6__He4_d ) +
          ( -rho*Y[jn]*Y[jli6]*rate_eval.n_Li6__Li7 +Y[jli7]*rate_eval.Li7__n_Li6 ) +
          ( -rho*Y[jp]*Y[jli6]*rate_eval.p_Li6__Be7 +Y[jbe7]*rate_eval.Be7__p_Li6 ) +
          ( -rho*Y[jn]*Y[jli6]*rate_eval.n_Li6__He4_t +rho*Y[jt]*Y[jhe4]*rate_eval.He4_t__n_Li6 ) +
          ( -rho*Y[jp]*Y[jli6]*rate_eval.p_Li6__He4_He3 +rho*Y[jhe3]*Y[jhe4]*rate_eval.He4_He3__p_Li6 ) +
          ( -rho*Y[jd]*Y[jli6]*rate_eval.d_Li6__n_Be7 +rho*Y[jn]*Y[jbe7]*rate_eval.n_Be7__d_Li6 ) +
          ( -rho*Y[jd]*Y[jli6]*rate_eval.d_Li6__p_Li7 +rho*Y[jp]*Y[jli7]*rate_eval.p_Li7__d_Li6 ) +
          ( +rho*Y[jp]*Y[jbe9]*rate_eval.p_Be9__He4_Li6 -rho*Y[jhe4]*Y[jli6]*rate_eval.He4_Li6__p_Be9 ) +
          ( +rho**2*Y[jn]*Y[jp]*Y[jhe4]*rate_eval.n_p_He4__Li6 -Y[jli6]*rate_eval.Li6__n_p_He4 )
       )

    dYdt[jli7] = (
          +rho*ye(Y)*Y[jbe7]*rate_eval.Be7__Li7__weak__electron_capture  +
          ( +rho*Y[jt]*Y[jhe4]*rate_eval.He4_t__Li7 -Y[jli7]*rate_eval.Li7__He4_t ) +
          ( +rho*Y[jn]*Y[jli6]*rate_eval.n_Li6__Li7 -Y[jli7]*rate_eval.Li7__n_Li6 ) +
          ( +rho*Y[jd]*Y[jli6]*rate_eval.d_Li6__p_Li7 -rho*Y[jp]*Y[jli7]*rate_eval.p_Li7__d_Li6 ) +
          ( -rho*Y[jp]*Y[jli7]*rate_eval.p_Li7__He4_He4 +5.00000000000000e-01*rho*Y[jhe4]**2*rate_eval.He4_He4__p_Li7 ) +
          ( -rho*Y[jt]*Y[jli7]*rate_eval.t_Li7__n_Be9 +rho*Y[jn]*Y[jbe9]*rate_eval.n_Be9__t_Li7 ) +
          ( +rho*Y[jn]*Y[jbe7]*rate_eval.n_Be7__p_Li7 -rho*Y[jp]*Y[jli7]*rate_eval.p_Li7__n_Be7 ) +
          ( -rho*Y[jd]*Y[jli7]*rate_eval.d_Li7__n_He4_He4 +5.00000000000000e-01*rho**2*Y[jn]*Y[jhe4]**2*rate_eval.n_He4_He4__d_Li7 ) +
          ( -rho*Y[jt]*Y[jli7]*rate_eval.t_Li7__n_n_He4_He4 +2.50000000000000e-01*rho**3*Y[jn]**2*Y[jhe4]**2*rate_eval.n_n_He4_He4__t_Li7 ) +
          ( -rho*Y[jhe3]*Y[jli7]*rate_eval.He3_Li7__n_p_He4_He4 +5.00000000000000e-01*rho**3*Y[jn]*Y[jp]*Y[jhe4]**2*rate_eval.n_p_He4_He4__He3_Li7 )
       )

    dYdt[jbe7] = (
          -rho*ye(Y)*Y[jbe7]*rate_eval.Be7__Li7__weak__electron_capture  +
          ( +rho*Y[jhe3]*Y[jhe4]*rate_eval.He4_He3__Be7 -Y[jbe7]*rate_eval.Be7__He4_He3 ) +
          ( +rho*Y[jp]*Y[jli6]*rate_eval.p_Li6__Be7 -Y[jbe7]*rate_eval.Be7__p_Li6 ) +
          ( -rho*Y[jp]*Y[jbe7]*rate_eval.p_Be7__B8 +Y[jb8]*rate_eval.B8__p_Be7 ) +
          ( -rho*Y[jhe4]*Y[jbe7]*rate_eval.He4_Be7__C11 +Y[jc11]*rate_eval.C11__He4_Be7 ) +
          ( +rho*Y[jd]*Y[jli6]*rate_eval.d_Li6__n_Be7 -rho*Y[jn]*Y[jbe7]*rate_eval.n_Be7__d_Li6 ) +
          ( -rho*Y[jn]*Y[jbe7]*rate_eval.n_Be7__p_Li7 +rho*Y[jp]*Y[jli7]*rate_eval.p_Li7__n_Be7 ) +
          ( -rho*Y[jn]*Y[jbe7]*rate_eval.n_Be7__He4_He4 +5.00000000000000e-01*rho*Y[jhe4]**2*rate_eval.He4_He4__n_Be7 ) +
          ( -rho*Y[jd]*Y[jbe7]*rate_eval.d_Be7__p_He4_He4 +5.00000000000000e-01*rho**2*Y[jp]*Y[jhe4]**2*rate_eval.p_He4_He4__d_Be7 ) +
          ( -rho*Y[jt]*Y[jbe7]*rate_eval.t_Be7__n_p_He4_He4 +5.00000000000000e-01*rho**3*Y[jn]*Y[jp]*Y[jhe4]**2*rate_eval.n_p_He4_He4__t_Be7 ) +
          ( -rho*Y[jhe3]*Y[jbe7]*rate_eval.He3_Be7__p_p_He4_He4 +2.50000000000000e-01*rho**3*Y[jp]**2*Y[jhe4]**2*rate_eval.p_p_He4_He4__He3_Be7 )
       )

    dYdt[jbe9] = (
          ( +rho*Y[jt]*Y[jli7]*rate_eval.t_Li7__n_Be9 -rho*Y[jn]*Y[jbe9]*rate_eval.n_Be9__t_Li7 ) +
          ( -rho*Y[jp]*Y[jbe9]*rate_eval.p_Be9__He4_Li6 +rho*Y[jhe4]*Y[jli6]*rate_eval.He4_Li6__p_Be9 ) +
          ( -rho*Y[jhe4]*Y[jbe9]*rate_eval.He4_Be9__n_C12 +rho*Y[jn]*Y[jc12]*rate_eval.n_C12__He4_Be9 ) +
          ( -rho*Y[jp]*Y[jbe9]*rate_eval.p_Be9__d_He4_He4 +5.00000000000000e-01*rho**2*Y[jd]*Y[jhe4]**2*rate_eval.d_He4_He4__p_Be9 ) +
          ( +5.00000000000000e-01*rho**2*Y[jn]*Y[jhe4]**2*rate_eval.n_He4_He4__Be9 -Y[jbe9]*rate_eval.Be9__n_He4_He4 ) +
          ( +5.00000000000000e-01*rho**3*Y[jn]*Y[jp]*Y[jhe4]**2*rate_eval.n_p_He4_He4__p_Be9 -rho*Y[jp]*Y[jbe9]*rate_eval.p_Be9__n_p_He4_He4 )
       )

    dYdt[jb8] = (
          -Y[jb8]*rate_eval.B8__He4_He4__weak__wc12  +
          ( +rho*Y[jp]*Y[jbe7]*rate_eval.p_Be7__B8 -Y[jb8]*rate_eval.B8__p_Be7 ) +
          ( -rho*Y[jhe4]*Y[jb8]*rate_eval.He4_B8__p_C11 +rho*Y[jp]*Y[jc11]*rate_eval.p_C11__He4_B8 ) +
          ( -rho*Y[jn]*Y[jb8]*rate_eval.n_B8__p_He4_He4 +5.00000000000000e-01*rho**2*Y[jp]*Y[jhe4]**2*rate_eval.p_He4_He4__n_B8 )
       )

    dYdt[jc11] = (
          ( +rho*Y[jhe4]*Y[jbe7]*rate_eval.He4_Be7__C11 -Y[jc11]*rate_eval.C11__He4_Be7 ) +
          ( -rho*Y[jn]*Y[jc11]*rate_eval.n_C11__C12 +Y[jc12]*rate_eval.C12__n_C11 ) +
          ( +rho*Y[jhe4]*Y[jb8]*rate_eval.He4_B8__p_C11 -rho*Y[jp]*Y[jc11]*rate_eval.p_C11__He4_B8 ) +
          ( -rho*Y[jhe4]*Y[jc11]*rate_eval.He4_C11__p_N14 +rho*Y[jp]*Y[jn14]*rate_eval.p_N14__He4_C11 ) +
          ( +rho*Y[jn]*Y[jo14]*rate_eval.n_O14__He4_C11 -rho*Y[jhe4]*Y[jc11]*rate_eval.He4_C11__n_O14 ) +
          ( -rho*Y[jn]*Y[jc11]*rate_eval.n_C11__He4_He4_He4 +1.66666666666667e-01*rho**2*Y[jhe4]**3*rate_eval.He4_He4_He4__n_C11 )
       )

    dYdt[jc12] = (
          ( +rho*Y[jn]*Y[jc11]*rate_eval.n_C11__C12 -Y[jc12]*rate_eval.C12__n_C11 ) +
          ( -rho*Y[jn]*Y[jc12]*rate_eval.n_C12__C13 +Y[jc13]*rate_eval.C13__n_C12 ) +
          ( -rho*Y[jp]*Y[jc12]*rate_eval.p_C12__N13 +Y[jn13]*rate_eval.N13__p_C12 ) +
          ( -rho*Y[jhe4]*Y[jc12]*rate_eval.He4_C12__O16 +Y[jo16]*rate_eval.O16__He4_C12 ) +
          ( +rho*Y[jhe4]*Y[jbe9]*rate_eval.He4_Be9__n_C12 -rho*Y[jn]*Y[jc12]*rate_eval.n_C12__He4_Be9 ) +
          ( + -2*5.00000000000000e-01*rho*Y[jc12]**2*rate_eval.C12_C12__p_Na23 + 2*rho*Y[jp]*Y[jna23]*rate_eval.p_Na23__C12_C12 ) +
          ( + -2*5.00000000000000e-01*rho*Y[jc12]**2*rate_eval.C12_C12__He4_Ne20 + 2*rho*Y[jhe4]*Y[jne20]*rate_eval.He4_Ne20__C12_C12 ) +
          ( +rho*Y[jp]*Y[jn15]*rate_eval.p_N15__He4_C12 -rho*Y[jhe4]*Y[jc12]*rate_eval.He4_C12__p_N15 ) +
          ( +rho*Y[jn]*Y[jo15]*rate_eval.n_O15__He4_C12 -rho*Y[jhe4]*Y[jc12]*rate_eval.He4_C12__n_O15 ) +
          ( -rho*Y[jc12]*Y[jo16]*rate_eval.C12_O16__p_Al27 +rho*Y[jp]*Y[jal27]*rate_eval.p_Al27__C12_O16 ) +
          ( -rho*Y[jc12]*Y[jo16]*rate_eval.C12_O16__He4_Mg24 +rho*Y[jhe4]*Y[jmg24]*rate_eval.He4_Mg24__C12_O16 ) +
          ( -rho*Y[jc12]*Y[jne20]*rate_eval.C12_Ne20__n_S31 +rho*Y[jn]*Y[js31]*rate_eval.n_S31__C12_Ne20 ) +
          ( -rho*Y[jc12]*Y[jne20]*rate_eval.C12_Ne20__p_P31 +rho*Y[jp]*Y[jp31]*rate_eval.p_P31__C12_Ne20 ) +
          ( -rho*Y[jc12]*Y[jne20]*rate_eval.C12_Ne20__He4_Si28 +rho*Y[jhe4]*Y[jsi28]*rate_eval.He4_Si28__C12_Ne20 ) +
          ( + 2*rho*Y[jn]*Y[jmg23]*rate_eval.n_Mg23__C12_C12 + -2*5.00000000000000e-01*rho*Y[jc12]**2*rate_eval.C12_C12__n_Mg23 ) +
          ( +1.66666666666667e-01*rho**2*Y[jhe4]**3*rate_eval.He4_He4_He4__C12 -Y[jc12]*rate_eval.C12__He4_He4_He4 )
       )

    dYdt[jc13] = (
          +Y[jn13]*rate_eval.N13__C13__weak__wc12  +
          ( +rho*Y[jn]*Y[jc12]*rate_eval.n_C12__C13 -Y[jc13]*rate_eval.C13__n_C12 ) +
          ( -rho*Y[jn]*Y[jc13]*rate_eval.n_C13__C14 +Y[jc14]*rate_eval.C14__n_C13 ) +
          ( -rho*Y[jp]*Y[jc13]*rate_eval.p_C13__N14 +Y[jn14]*rate_eval.N14__p_C13 ) +
          ( -rho*Y[jd]*Y[jc13]*rate_eval.d_C13__n_N14 +rho*Y[jn]*Y[jn14]*rate_eval.n_N14__d_C13 ) +
          ( -rho*Y[jhe4]*Y[jc13]*rate_eval.He4_C13__n_O16 +rho*Y[jn]*Y[jo16]*rate_eval.n_O16__He4_C13 ) +
          ( +rho*Y[jn]*Y[jn13]*rate_eval.n_N13__p_C13 -rho*Y[jp]*Y[jc13]*rate_eval.p_C13__n_N13 )
       )

    dYdt[jc14] = (
          -Y[jc14]*rate_eval.C14__N14__weak__wc12  +
          ( +rho*Y[jn]*Y[jc13]*rate_eval.n_C13__C14 -Y[jc14]*rate_eval.C14__n_C13 ) +
          ( -rho*Y[jp]*Y[jc14]*rate_eval.p_C14__N15 +Y[jn15]*rate_eval.N15__p_C14 ) +
          ( -rho*Y[jhe4]*Y[jc14]*rate_eval.He4_C14__O18 +Y[jo18]*rate_eval.O18__He4_C14 ) +
          ( -rho*Y[jd]*Y[jc14]*rate_eval.d_C14__n_N15 +rho*Y[jn]*Y[jn15]*rate_eval.n_N15__d_C14 ) +
          ( +rho*Y[jn]*Y[jn14]*rate_eval.n_N14__p_C14 -rho*Y[jp]*Y[jc14]*rate_eval.p_C14__n_N14 ) +
          ( +rho*Y[jn]*Y[jo17]*rate_eval.n_O17__He4_C14 -rho*Y[jhe4]*Y[jc14]*rate_eval.He4_C14__n_O17 )
       )

    dYdt[jn13] = (
          -Y[jn13]*rate_eval.N13__C13__weak__wc12  +
          ( +rho*Y[jp]*Y[jc12]*rate_eval.p_C12__N13 -Y[jn13]*rate_eval.N13__p_C12 ) +
          ( -rho*Y[jn]*Y[jn13]*rate_eval.n_N13__N14 +Y[jn14]*rate_eval.N14__n_N13 ) +
          ( -rho*Y[jp]*Y[jn13]*rate_eval.p_N13__O14 +Y[jo14]*rate_eval.O14__p_N13 ) +
          ( -rho*Y[jn]*Y[jn13]*rate_eval.n_N13__p_C13 +rho*Y[jp]*Y[jc13]*rate_eval.p_C13__n_N13 ) +
          ( -rho*Y[jhe4]*Y[jn13]*rate_eval.He4_N13__p_O16 +rho*Y[jp]*Y[jo16]*rate_eval.p_O16__He4_N13 ) +
          ( +rho*Y[jn]*Y[jf16]*rate_eval.n_F16__He4_N13 -rho*Y[jhe4]*Y[jn13]*rate_eval.He4_N13__n_F16 )
       )

    dYdt[jn14] = (
          +Y[jc14]*rate_eval.C14__N14__weak__wc12  +
          +Y[jo14]*rate_eval.O14__N14__weak__wc12  +
          ( +rho*Y[jp]*Y[jc13]*rate_eval.p_C13__N14 -Y[jn14]*rate_eval.N14__p_C13 ) +
          ( +rho*Y[jn]*Y[jn13]*rate_eval.n_N13__N14 -Y[jn14]*rate_eval.N14__n_N13 ) +
          ( -rho*Y[jn]*Y[jn14]*rate_eval.n_N14__N15 +Y[jn15]*rate_eval.N15__n_N14 ) +
          ( -rho*Y[jp]*Y[jn14]*rate_eval.p_N14__O15 +Y[jo15]*rate_eval.O15__p_N14 ) +
          ( -rho*Y[jhe4]*Y[jn14]*rate_eval.He4_N14__F18 +Y[jf18]*rate_eval.F18__He4_N14 ) +
          ( +rho*Y[jhe4]*Y[jc11]*rate_eval.He4_C11__p_N14 -rho*Y[jp]*Y[jn14]*rate_eval.p_N14__He4_C11 ) +
          ( +rho*Y[jd]*Y[jc13]*rate_eval.d_C13__n_N14 -rho*Y[jn]*Y[jn14]*rate_eval.n_N14__d_C13 ) +
          ( -rho*Y[jn]*Y[jn14]*rate_eval.n_N14__p_C14 +rho*Y[jp]*Y[jc14]*rate_eval.p_C14__n_N14 ) +
          ( +rho*Y[jn]*Y[jo14]*rate_eval.n_O14__p_N14 -rho*Y[jp]*Y[jn14]*rate_eval.p_N14__n_O14 ) +
          ( +rho*Y[jp]*Y[jo17]*rate_eval.p_O17__He4_N14 -rho*Y[jhe4]*Y[jn14]*rate_eval.He4_N14__p_O17 ) +
          ( +rho*Y[jn]*Y[jf17]*rate_eval.n_F17__He4_N14 -rho*Y[jhe4]*Y[jn14]*rate_eval.He4_N14__n_F17 )
       )

    dYdt[jn15] = (
          +Y[jo15]*rate_eval.O15__N15__weak__wc12  +
          ( +rho*Y[jp]*Y[jc14]*rate_eval.p_C14__N15 -Y[jn15]*rate_eval.N15__p_C14 ) +
          ( +rho*Y[jn]*Y[jn14]*rate_eval.n_N14__N15 -Y[jn15]*rate_eval.N15__n_N14 ) +
          ( -rho*Y[jp]*Y[jn15]*rate_eval.p_N15__O16 +Y[jo16]*rate_eval.O16__p_N15 ) +
          ( -rho*Y[jhe4]*Y[jn15]*rate_eval.He4_N15__F19 +Y[jf19]*rate_eval.F19__He4_N15 ) +
          ( +rho*Y[jd]*Y[jc14]*rate_eval.d_C14__n_N15 -rho*Y[jn]*Y[jn15]*rate_eval.n_N15__d_C14 ) +
          ( -rho*Y[jp]*Y[jn15]*rate_eval.p_N15__He4_C12 +rho*Y[jhe4]*Y[jc12]*rate_eval.He4_C12__p_N15 ) +
          ( +rho*Y[jn]*Y[jo15]*rate_eval.n_O15__p_N15 -rho*Y[jp]*Y[jn15]*rate_eval.p_N15__n_O15 ) +
          ( +rho*Y[jp]*Y[jo18]*rate_eval.p_O18__He4_N15 -rho*Y[jhe4]*Y[jn15]*rate_eval.He4_N15__p_O18 ) +
          ( +rho*Y[jn]*Y[jf18]*rate_eval.n_F18__He4_N15 -rho*Y[jhe4]*Y[jn15]*rate_eval.He4_N15__n_F18 )
       )

    dYdt[jo14] = (
          -Y[jo14]*rate_eval.O14__N14__weak__wc12  +
          ( +Y[jf15]*rate_eval.F15__p_O14 -rho*Y[jp]*Y[jo14]*rate_eval.p_O14__F15 ) +
          ( +rho*Y[jp]*Y[jn13]*rate_eval.p_N13__O14 -Y[jo14]*rate_eval.O14__p_N13 ) +
          ( -rho*Y[jn]*Y[jo14]*rate_eval.n_O14__O15 +Y[jo15]*rate_eval.O15__n_O14 ) +
          ( -rho*Y[jhe4]*Y[jo14]*rate_eval.He4_O14__Ne18 +Y[jne18]*rate_eval.Ne18__He4_O14 ) +
          ( -rho*Y[jn]*Y[jo14]*rate_eval.n_O14__p_N14 +rho*Y[jp]*Y[jn14]*rate_eval.p_N14__n_O14 ) +
          ( -rho*Y[jn]*Y[jo14]*rate_eval.n_O14__He4_C11 +rho*Y[jhe4]*Y[jc11]*rate_eval.He4_C11__n_O14 ) +
          ( -rho*Y[jhe4]*Y[jo14]*rate_eval.He4_O14__p_F17 +rho*Y[jp]*Y[jf17]*rate_eval.p_F17__He4_O14 )
       )

    dYdt[jo15] = (
          -Y[jo15]*rate_eval.O15__N15__weak__wc12  +
          ( +Y[jf16]*rate_eval.F16__p_O15 -rho*Y[jp]*Y[jo15]*rate_eval.p_O15__F16 ) +
          ( +rho*Y[jp]*Y[jn14]*rate_eval.p_N14__O15 -Y[jo15]*rate_eval.O15__p_N14 ) +
          ( +rho*Y[jn]*Y[jo14]*rate_eval.n_O14__O15 -Y[jo15]*rate_eval.O15__n_O14 ) +
          ( -rho*Y[jn]*Y[jo15]*rate_eval.n_O15__O16 +Y[jo16]*rate_eval.O16__n_O15 ) +
          ( -rho*Y[jhe4]*Y[jo15]*rate_eval.He4_O15__Ne19 +Y[jne19]*rate_eval.Ne19__He4_O15 ) +
          ( -rho*Y[jn]*Y[jo15]*rate_eval.n_O15__p_N15 +rho*Y[jp]*Y[jn15]*rate_eval.p_N15__n_O15 ) +
          ( -rho*Y[jn]*Y[jo15]*rate_eval.n_O15__He4_C12 +rho*Y[jhe4]*Y[jc12]*rate_eval.He4_C12__n_O15 ) +
          ( +rho*Y[jn]*Y[jf15]*rate_eval.n_F15__p_O15 -rho*Y[jp]*Y[jo15]*rate_eval.p_O15__n_F15 ) +
          ( +rho*Y[jp]*Y[jf18]*rate_eval.p_F18__He4_O15 -rho*Y[jhe4]*Y[jo15]*rate_eval.He4_O15__p_F18 ) +
          ( +rho*Y[jn]*Y[jne18]*rate_eval.n_Ne18__He4_O15 -rho*Y[jhe4]*Y[jo15]*rate_eval.He4_O15__n_Ne18 )
       )

    dYdt[jo16] = (
          ( +rho*Y[jhe4]*Y[jc12]*rate_eval.He4_C12__O16 -Y[jo16]*rate_eval.O16__He4_C12 ) +
          ( +rho*Y[jp]*Y[jn15]*rate_eval.p_N15__O16 -Y[jo16]*rate_eval.O16__p_N15 ) +
          ( +rho*Y[jn]*Y[jo15]*rate_eval.n_O15__O16 -Y[jo16]*rate_eval.O16__n_O15 ) +
          ( -rho*Y[jn]*Y[jo16]*rate_eval.n_O16__O17 +Y[jo17]*rate_eval.O17__n_O16 ) +
          ( -rho*Y[jp]*Y[jo16]*rate_eval.p_O16__F17 +Y[jf17]*rate_eval.F17__p_O16 ) +
          ( -rho*Y[jhe4]*Y[jo16]*rate_eval.He4_O16__Ne20 +Y[jne20]*rate_eval.Ne20__He4_O16 ) +
          ( +rho*Y[jhe4]*Y[jc13]*rate_eval.He4_C13__n_O16 -rho*Y[jn]*Y[jo16]*rate_eval.n_O16__He4_C13 ) +
          ( +rho*Y[jhe4]*Y[jn13]*rate_eval.He4_N13__p_O16 -rho*Y[jp]*Y[jo16]*rate_eval.p_O16__He4_N13 ) +
          ( -rho*Y[jc12]*Y[jo16]*rate_eval.C12_O16__p_Al27 +rho*Y[jp]*Y[jal27]*rate_eval.p_Al27__C12_O16 ) +
          ( -rho*Y[jc12]*Y[jo16]*rate_eval.C12_O16__He4_Mg24 +rho*Y[jhe4]*Y[jmg24]*rate_eval.He4_Mg24__C12_O16 ) +
          ( + -2*5.00000000000000e-01*rho*Y[jo16]**2*rate_eval.O16_O16__n_S31 + 2*rho*Y[jn]*Y[js31]*rate_eval.n_S31__O16_O16 ) +
          ( + -2*5.00000000000000e-01*rho*Y[jo16]**2*rate_eval.O16_O16__p_P31 + 2*rho*Y[jp]*Y[jp31]*rate_eval.p_P31__O16_O16 ) +
          ( + -2*5.00000000000000e-01*rho*Y[jo16]**2*rate_eval.O16_O16__He4_Si28 + 2*rho*Y[jhe4]*Y[jsi28]*rate_eval.He4_Si28__O16_O16 ) +
          ( +rho*Y[jn]*Y[jf16]*rate_eval.n_F16__p_O16 -rho*Y[jp]*Y[jo16]*rate_eval.p_O16__n_F16 ) +
          ( +rho*Y[jp]*Y[jf19]*rate_eval.p_F19__He4_O16 -rho*Y[jhe4]*Y[jo16]*rate_eval.He4_O16__p_F19 ) +
          ( +rho*Y[jn]*Y[jne19]*rate_eval.n_Ne19__He4_O16 -rho*Y[jhe4]*Y[jo16]*rate_eval.He4_O16__n_Ne19 )
       )

    dYdt[jo17] = (
          +Y[jf17]*rate_eval.F17__O17__weak__wc12  +
          ( +rho*Y[jn]*Y[jo16]*rate_eval.n_O16__O17 -Y[jo17]*rate_eval.O17__n_O16 ) +
          ( -rho*Y[jn]*Y[jo17]*rate_eval.n_O17__O18 +Y[jo18]*rate_eval.O18__n_O17 ) +
          ( -rho*Y[jp]*Y[jo17]*rate_eval.p_O17__F18 +Y[jf18]*rate_eval.F18__p_O17 ) +
          ( -rho*Y[jhe4]*Y[jo17]*rate_eval.He4_O17__Ne21 +Y[jne21]*rate_eval.Ne21__He4_O17 ) +
          ( -rho*Y[jn]*Y[jo17]*rate_eval.n_O17__He4_C14 +rho*Y[jhe4]*Y[jc14]*rate_eval.He4_C14__n_O17 ) +
          ( -rho*Y[jp]*Y[jo17]*rate_eval.p_O17__He4_N14 +rho*Y[jhe4]*Y[jn14]*rate_eval.He4_N14__p_O17 ) +
          ( -rho*Y[jhe4]*Y[jo17]*rate_eval.He4_O17__n_Ne20 +rho*Y[jn]*Y[jne20]*rate_eval.n_Ne20__He4_O17 ) +
          ( +rho*Y[jn]*Y[jf17]*rate_eval.n_F17__p_O17 -rho*Y[jp]*Y[jo17]*rate_eval.p_O17__n_F17 )
       )

    dYdt[jo18] = (
          +Y[jf18]*rate_eval.F18__O18__weak__wc12  +
          ( +rho*Y[jhe4]*Y[jc14]*rate_eval.He4_C14__O18 -Y[jo18]*rate_eval.O18__He4_C14 ) +
          ( +rho*Y[jn]*Y[jo17]*rate_eval.n_O17__O18 -Y[jo18]*rate_eval.O18__n_O17 ) +
          ( -rho*Y[jp]*Y[jo18]*rate_eval.p_O18__F19 +Y[jf19]*rate_eval.F19__p_O18 ) +
          ( -rho*Y[jp]*Y[jo18]*rate_eval.p_O18__He4_N15 +rho*Y[jhe4]*Y[jn15]*rate_eval.He4_N15__p_O18 ) +
          ( +rho*Y[jn]*Y[jf18]*rate_eval.n_F18__p_O18 -rho*Y[jp]*Y[jo18]*rate_eval.p_O18__n_F18 ) +
          ( +rho*Y[jn]*Y[jne21]*rate_eval.n_Ne21__He4_O18 -rho*Y[jhe4]*Y[jo18]*rate_eval.He4_O18__n_Ne21 )
       )

    dYdt[jf15] = (
          ( -Y[jf15]*rate_eval.F15__p_O14 +rho*Y[jp]*Y[jo14]*rate_eval.p_O14__F15 ) +
          ( -rho*Y[jn]*Y[jf15]*rate_eval.n_F15__F16 +Y[jf16]*rate_eval.F16__n_F15 ) +
          ( -rho*Y[jn]*Y[jf15]*rate_eval.n_F15__p_O15 +rho*Y[jp]*Y[jo15]*rate_eval.p_O15__n_F15 ) +
          ( -rho*Y[jhe4]*Y[jf15]*rate_eval.He4_F15__p_Ne18 +rho*Y[jp]*Y[jne18]*rate_eval.p_Ne18__He4_F15 )
       )

    dYdt[jf16] = (
          ( -Y[jf16]*rate_eval.F16__p_O15 +rho*Y[jp]*Y[jo15]*rate_eval.p_O15__F16 ) +
          ( +rho*Y[jn]*Y[jf15]*rate_eval.n_F15__F16 -Y[jf16]*rate_eval.F16__n_F15 ) +
          ( -rho*Y[jn]*Y[jf16]*rate_eval.n_F16__F17 +Y[jf17]*rate_eval.F17__n_F16 ) +
          ( -rho*Y[jn]*Y[jf16]*rate_eval.n_F16__p_O16 +rho*Y[jp]*Y[jo16]*rate_eval.p_O16__n_F16 ) +
          ( -rho*Y[jn]*Y[jf16]*rate_eval.n_F16__He4_N13 +rho*Y[jhe4]*Y[jn13]*rate_eval.He4_N13__n_F16 ) +
          ( -rho*Y[jhe4]*Y[jf16]*rate_eval.He4_F16__p_Ne19 +rho*Y[jp]*Y[jne19]*rate_eval.p_Ne19__He4_F16 )
       )

    dYdt[jf17] = (
          -Y[jf17]*rate_eval.F17__O17__weak__wc12  +
          ( +rho*Y[jp]*Y[jo16]*rate_eval.p_O16__F17 -Y[jf17]*rate_eval.F17__p_O16 ) +
          ( +rho*Y[jn]*Y[jf16]*rate_eval.n_F16__F17 -Y[jf17]*rate_eval.F17__n_F16 ) +
          ( -rho*Y[jn]*Y[jf17]*rate_eval.n_F17__F18 +Y[jf18]*rate_eval.F18__n_F17 ) +
          ( -rho*Y[jp]*Y[jf17]*rate_eval.p_F17__Ne18 +Y[jne18]*rate_eval.Ne18__p_F17 ) +
          ( -rho*Y[jhe4]*Y[jf17]*rate_eval.He4_F17__Na21 +Y[jna21]*rate_eval.Na21__He4_F17 ) +
          ( +rho*Y[jhe4]*Y[jo14]*rate_eval.He4_O14__p_F17 -rho*Y[jp]*Y[jf17]*rate_eval.p_F17__He4_O14 ) +
          ( -rho*Y[jn]*Y[jf17]*rate_eval.n_F17__p_O17 +rho*Y[jp]*Y[jo17]*rate_eval.p_O17__n_F17 ) +
          ( -rho*Y[jn]*Y[jf17]*rate_eval.n_F17__He4_N14 +rho*Y[jhe4]*Y[jn14]*rate_eval.He4_N14__n_F17 ) +
          ( -rho*Y[jhe4]*Y[jf17]*rate_eval.He4_F17__p_Ne20 +rho*Y[jp]*Y[jne20]*rate_eval.p_Ne20__He4_F17 )
       )

    dYdt[jf18] = (
          -Y[jf18]*rate_eval.F18__O18__weak__wc12  +
          +Y[jne18]*rate_eval.Ne18__F18__weak__wc12  +
          ( +rho*Y[jhe4]*Y[jn14]*rate_eval.He4_N14__F18 -Y[jf18]*rate_eval.F18__He4_N14 ) +
          ( +rho*Y[jp]*Y[jo17]*rate_eval.p_O17__F18 -Y[jf18]*rate_eval.F18__p_O17 ) +
          ( +rho*Y[jn]*Y[jf17]*rate_eval.n_F17__F18 -Y[jf18]*rate_eval.F18__n_F17 ) +
          ( -rho*Y[jn]*Y[jf18]*rate_eval.n_F18__F19 +Y[jf19]*rate_eval.F19__n_F18 ) +
          ( -rho*Y[jp]*Y[jf18]*rate_eval.p_F18__Ne19 +Y[jne19]*rate_eval.Ne19__p_F18 ) +
          ( -rho*Y[jhe4]*Y[jf18]*rate_eval.He4_F18__Na22 +Y[jna22]*rate_eval.Na22__He4_F18 ) +
          ( -rho*Y[jn]*Y[jf18]*rate_eval.n_F18__p_O18 +rho*Y[jp]*Y[jo18]*rate_eval.p_O18__n_F18 ) +
          ( -rho*Y[jn]*Y[jf18]*rate_eval.n_F18__He4_N15 +rho*Y[jhe4]*Y[jn15]*rate_eval.He4_N15__n_F18 ) +
          ( -rho*Y[jp]*Y[jf18]*rate_eval.p_F18__He4_O15 +rho*Y[jhe4]*Y[jo15]*rate_eval.He4_O15__p_F18 ) +
          ( -rho*Y[jhe4]*Y[jf18]*rate_eval.He4_F18__p_Ne21 +rho*Y[jp]*Y[jne21]*rate_eval.p_Ne21__He4_F18 ) +
          ( +rho*Y[jn]*Y[jne18]*rate_eval.n_Ne18__p_F18 -rho*Y[jp]*Y[jf18]*rate_eval.p_F18__n_Ne18 ) +
          ( +rho*Y[jn]*Y[jna21]*rate_eval.n_Na21__He4_F18 -rho*Y[jhe4]*Y[jf18]*rate_eval.He4_F18__n_Na21 )
       )

    dYdt[jf19] = (
          +Y[jne19]*rate_eval.Ne19__F19__weak__wc12  +
          ( +rho*Y[jhe4]*Y[jn15]*rate_eval.He4_N15__F19 -Y[jf19]*rate_eval.F19__He4_N15 ) +
          ( +rho*Y[jp]*Y[jo18]*rate_eval.p_O18__F19 -Y[jf19]*rate_eval.F19__p_O18 ) +
          ( +rho*Y[jn]*Y[jf18]*rate_eval.n_F18__F19 -Y[jf19]*rate_eval.F19__n_F18 ) +
          ( -rho*Y[jp]*Y[jf19]*rate_eval.p_F19__Ne20 +Y[jne20]*rate_eval.Ne20__p_F19 ) +
          ( -rho*Y[jhe4]*Y[jf19]*rate_eval.He4_F19__Na23 +Y[jna23]*rate_eval.Na23__He4_F19 ) +
          ( -rho*Y[jp]*Y[jf19]*rate_eval.p_F19__He4_O16 +rho*Y[jhe4]*Y[jo16]*rate_eval.He4_O16__p_F19 ) +
          ( +rho*Y[jn]*Y[jne19]*rate_eval.n_Ne19__p_F19 -rho*Y[jp]*Y[jf19]*rate_eval.p_F19__n_Ne19 ) +
          ( +rho*Y[jn]*Y[jna22]*rate_eval.n_Na22__He4_F19 -rho*Y[jhe4]*Y[jf19]*rate_eval.He4_F19__n_Na22 )
       )

    dYdt[jne18] = (
          -Y[jne18]*rate_eval.Ne18__F18__weak__wc12  +
          ( +rho*Y[jhe4]*Y[jo14]*rate_eval.He4_O14__Ne18 -Y[jne18]*rate_eval.Ne18__He4_O14 ) +
          ( +rho*Y[jp]*Y[jf17]*rate_eval.p_F17__Ne18 -Y[jne18]*rate_eval.Ne18__p_F17 ) +
          ( -rho*Y[jn]*Y[jne18]*rate_eval.n_Ne18__Ne19 +Y[jne19]*rate_eval.Ne19__n_Ne18 ) +
          ( -rho*Y[jhe4]*Y[jne18]*rate_eval.He4_Ne18__Mg22 +Y[jmg22]*rate_eval.Mg22__He4_Ne18 ) +
          ( +rho*Y[jhe4]*Y[jf15]*rate_eval.He4_F15__p_Ne18 -rho*Y[jp]*Y[jne18]*rate_eval.p_Ne18__He4_F15 ) +
          ( -rho*Y[jn]*Y[jne18]*rate_eval.n_Ne18__p_F18 +rho*Y[jp]*Y[jf18]*rate_eval.p_F18__n_Ne18 ) +
          ( -rho*Y[jn]*Y[jne18]*rate_eval.n_Ne18__He4_O15 +rho*Y[jhe4]*Y[jo15]*rate_eval.He4_O15__n_Ne18 ) +
          ( -rho*Y[jhe4]*Y[jne18]*rate_eval.He4_Ne18__p_Na21 +rho*Y[jp]*Y[jna21]*rate_eval.p_Na21__He4_Ne18 )
       )

    dYdt[jne19] = (
          -Y[jne19]*rate_eval.Ne19__F19__weak__wc12  +
          ( +rho*Y[jhe4]*Y[jo15]*rate_eval.He4_O15__Ne19 -Y[jne19]*rate_eval.Ne19__He4_O15 ) +
          ( +rho*Y[jp]*Y[jf18]*rate_eval.p_F18__Ne19 -Y[jne19]*rate_eval.Ne19__p_F18 ) +
          ( +rho*Y[jn]*Y[jne18]*rate_eval.n_Ne18__Ne19 -Y[jne19]*rate_eval.Ne19__n_Ne18 ) +
          ( -rho*Y[jn]*Y[jne19]*rate_eval.n_Ne19__Ne20 +Y[jne20]*rate_eval.Ne20__n_Ne19 ) +
          ( -rho*Y[jhe4]*Y[jne19]*rate_eval.He4_Ne19__Mg23 +Y[jmg23]*rate_eval.Mg23__He4_Ne19 ) +
          ( +rho*Y[jhe4]*Y[jf16]*rate_eval.He4_F16__p_Ne19 -rho*Y[jp]*Y[jne19]*rate_eval.p_Ne19__He4_F16 ) +
          ( -rho*Y[jn]*Y[jne19]*rate_eval.n_Ne19__p_F19 +rho*Y[jp]*Y[jf19]*rate_eval.p_F19__n_Ne19 ) +
          ( -rho*Y[jn]*Y[jne19]*rate_eval.n_Ne19__He4_O16 +rho*Y[jhe4]*Y[jo16]*rate_eval.He4_O16__n_Ne19 ) +
          ( -rho*Y[jhe4]*Y[jne19]*rate_eval.He4_Ne19__p_Na22 +rho*Y[jp]*Y[jna22]*rate_eval.p_Na22__He4_Ne19 ) +
          ( +rho*Y[jn]*Y[jmg22]*rate_eval.n_Mg22__He4_Ne19 -rho*Y[jhe4]*Y[jne19]*rate_eval.He4_Ne19__n_Mg22 )
       )

    dYdt[jne20] = (
          ( +rho*Y[jhe4]*Y[jo16]*rate_eval.He4_O16__Ne20 -Y[jne20]*rate_eval.Ne20__He4_O16 ) +
          ( +rho*Y[jp]*Y[jf19]*rate_eval.p_F19__Ne20 -Y[jne20]*rate_eval.Ne20__p_F19 ) +
          ( +rho*Y[jn]*Y[jne19]*rate_eval.n_Ne19__Ne20 -Y[jne20]*rate_eval.Ne20__n_Ne19 ) +
          ( -rho*Y[jn]*Y[jne20]*rate_eval.n_Ne20__Ne21 +Y[jne21]*rate_eval.Ne21__n_Ne20 ) +
          ( -rho*Y[jp]*Y[jne20]*rate_eval.p_Ne20__Na21 +Y[jna21]*rate_eval.Na21__p_Ne20 ) +
          ( -rho*Y[jhe4]*Y[jne20]*rate_eval.He4_Ne20__Mg24 +Y[jmg24]*rate_eval.Mg24__He4_Ne20 ) +
          ( +5.00000000000000e-01*rho*Y[jc12]**2*rate_eval.C12_C12__He4_Ne20 -rho*Y[jhe4]*Y[jne20]*rate_eval.He4_Ne20__C12_C12 ) +
          ( +rho*Y[jhe4]*Y[jo17]*rate_eval.He4_O17__n_Ne20 -rho*Y[jn]*Y[jne20]*rate_eval.n_Ne20__He4_O17 ) +
          ( +rho*Y[jhe4]*Y[jf17]*rate_eval.He4_F17__p_Ne20 -rho*Y[jp]*Y[jne20]*rate_eval.p_Ne20__He4_F17 ) +
          ( -rho*Y[jc12]*Y[jne20]*rate_eval.C12_Ne20__n_S31 +rho*Y[jn]*Y[js31]*rate_eval.n_S31__C12_Ne20 ) +
          ( -rho*Y[jc12]*Y[jne20]*rate_eval.C12_Ne20__p_P31 +rho*Y[jp]*Y[jp31]*rate_eval.p_P31__C12_Ne20 ) +
          ( -rho*Y[jc12]*Y[jne20]*rate_eval.C12_Ne20__He4_Si28 +rho*Y[jhe4]*Y[jsi28]*rate_eval.He4_Si28__C12_Ne20 ) +
          ( +rho*Y[jp]*Y[jna23]*rate_eval.p_Na23__He4_Ne20 -rho*Y[jhe4]*Y[jne20]*rate_eval.He4_Ne20__p_Na23 ) +
          ( +rho*Y[jn]*Y[jmg23]*rate_eval.n_Mg23__He4_Ne20 -rho*Y[jhe4]*Y[jne20]*rate_eval.He4_Ne20__n_Mg23 )
       )

    dYdt[jne21] = (
          +Y[jna21]*rate_eval.Na21__Ne21__weak__wc12  +
          ( +rho*Y[jhe4]*Y[jo17]*rate_eval.He4_O17__Ne21 -Y[jne21]*rate_eval.Ne21__He4_O17 ) +
          ( +rho*Y[jn]*Y[jne20]*rate_eval.n_Ne20__Ne21 -Y[jne21]*rate_eval.Ne21__n_Ne20 ) +
          ( -rho*Y[jp]*Y[jne21]*rate_eval.p_Ne21__Na22 +Y[jna22]*rate_eval.Na22__p_Ne21 ) +
          ( -rho*Y[jhe4]*Y[jne21]*rate_eval.He4_Ne21__Mg25 +Y[jmg25]*rate_eval.Mg25__He4_Ne21 ) +
          ( +rho*Y[jhe4]*Y[jf18]*rate_eval.He4_F18__p_Ne21 -rho*Y[jp]*Y[jne21]*rate_eval.p_Ne21__He4_F18 ) +
          ( -rho*Y[jn]*Y[jne21]*rate_eval.n_Ne21__He4_O18 +rho*Y[jhe4]*Y[jo18]*rate_eval.He4_O18__n_Ne21 ) +
          ( -rho*Y[jhe4]*Y[jne21]*rate_eval.He4_Ne21__n_Mg24 +rho*Y[jn]*Y[jmg24]*rate_eval.n_Mg24__He4_Ne21 ) +
          ( +rho*Y[jn]*Y[jna21]*rate_eval.n_Na21__p_Ne21 -rho*Y[jp]*Y[jne21]*rate_eval.p_Ne21__n_Na21 ) +
          ( +rho*Y[jp]*Y[jna24]*rate_eval.p_Na24__He4_Ne21 -rho*Y[jhe4]*Y[jne21]*rate_eval.He4_Ne21__p_Na24 )
       )

    dYdt[jna21] = (
          -Y[jna21]*rate_eval.Na21__Ne21__weak__wc12  +
          ( +rho*Y[jhe4]*Y[jf17]*rate_eval.He4_F17__Na21 -Y[jna21]*rate_eval.Na21__He4_F17 ) +
          ( +rho*Y[jp]*Y[jne20]*rate_eval.p_Ne20__Na21 -Y[jna21]*rate_eval.Na21__p_Ne20 ) +
          ( -rho*Y[jn]*Y[jna21]*rate_eval.n_Na21__Na22 +Y[jna22]*rate_eval.Na22__n_Na21 ) +
          ( -rho*Y[jp]*Y[jna21]*rate_eval.p_Na21__Mg22 +Y[jmg22]*rate_eval.Mg22__p_Na21 ) +
          ( -rho*Y[jhe4]*Y[jna21]*rate_eval.He4_Na21__Al25 +Y[jal25]*rate_eval.Al25__He4_Na21 ) +
          ( +rho*Y[jhe4]*Y[jne18]*rate_eval.He4_Ne18__p_Na21 -rho*Y[jp]*Y[jna21]*rate_eval.p_Na21__He4_Ne18 ) +
          ( -rho*Y[jn]*Y[jna21]*rate_eval.n_Na21__p_Ne21 +rho*Y[jp]*Y[jne21]*rate_eval.p_Ne21__n_Na21 ) +
          ( -rho*Y[jn]*Y[jna21]*rate_eval.n_Na21__He4_F18 +rho*Y[jhe4]*Y[jf18]*rate_eval.He4_F18__n_Na21 ) +
          ( -rho*Y[jhe4]*Y[jna21]*rate_eval.He4_Na21__p_Mg24 +rho*Y[jp]*Y[jmg24]*rate_eval.p_Mg24__He4_Na21 )
       )

    dYdt[jna22] = (
          +Y[jmg22]*rate_eval.Mg22__Na22__weak__wc12  +
          ( +rho*Y[jhe4]*Y[jf18]*rate_eval.He4_F18__Na22 -Y[jna22]*rate_eval.Na22__He4_F18 ) +
          ( +rho*Y[jp]*Y[jne21]*rate_eval.p_Ne21__Na22 -Y[jna22]*rate_eval.Na22__p_Ne21 ) +
          ( +rho*Y[jn]*Y[jna21]*rate_eval.n_Na21__Na22 -Y[jna22]*rate_eval.Na22__n_Na21 ) +
          ( -rho*Y[jn]*Y[jna22]*rate_eval.n_Na22__Na23 +Y[jna23]*rate_eval.Na23__n_Na22 ) +
          ( -rho*Y[jp]*Y[jna22]*rate_eval.p_Na22__Mg23 +Y[jmg23]*rate_eval.Mg23__p_Na22 ) +
          ( -rho*Y[jhe4]*Y[jna22]*rate_eval.He4_Na22__Al26 +Y[jal26]*rate_eval.Al26__He4_Na22 ) +
          ( +rho*Y[jhe4]*Y[jne19]*rate_eval.He4_Ne19__p_Na22 -rho*Y[jp]*Y[jna22]*rate_eval.p_Na22__He4_Ne19 ) +
          ( -rho*Y[jn]*Y[jna22]*rate_eval.n_Na22__He4_F19 +rho*Y[jhe4]*Y[jf19]*rate_eval.He4_F19__n_Na22 ) +
          ( -rho*Y[jhe4]*Y[jna22]*rate_eval.He4_Na22__p_Mg25 +rho*Y[jp]*Y[jmg25]*rate_eval.p_Mg25__He4_Na22 ) +
          ( +rho*Y[jn]*Y[jmg22]*rate_eval.n_Mg22__p_Na22 -rho*Y[jp]*Y[jna22]*rate_eval.p_Na22__n_Mg22 ) +
          ( +rho*Y[jn]*Y[jal25]*rate_eval.n_Al25__He4_Na22 -rho*Y[jhe4]*Y[jna22]*rate_eval.He4_Na22__n_Al25 )
       )

    dYdt[jna23] = (
          +Y[jmg23]*rate_eval.Mg23__Na23__weak__wc12  +
          ( +rho*Y[jhe4]*Y[jf19]*rate_eval.He4_F19__Na23 -Y[jna23]*rate_eval.Na23__He4_F19 ) +
          ( +rho*Y[jn]*Y[jna22]*rate_eval.n_Na22__Na23 -Y[jna23]*rate_eval.Na23__n_Na22 ) +
          ( -rho*Y[jn]*Y[jna23]*rate_eval.n_Na23__Na24 +Y[jna24]*rate_eval.Na24__n_Na23 ) +
          ( -rho*Y[jp]*Y[jna23]*rate_eval.p_Na23__Mg24 +Y[jmg24]*rate_eval.Mg24__p_Na23 ) +
          ( -rho*Y[jhe4]*Y[jna23]*rate_eval.He4_Na23__Al27 +Y[jal27]*rate_eval.Al27__He4_Na23 ) +
          ( +5.00000000000000e-01*rho*Y[jc12]**2*rate_eval.C12_C12__p_Na23 -rho*Y[jp]*Y[jna23]*rate_eval.p_Na23__C12_C12 ) +
          ( -rho*Y[jp]*Y[jna23]*rate_eval.p_Na23__He4_Ne20 +rho*Y[jhe4]*Y[jne20]*rate_eval.He4_Ne20__p_Na23 ) +
          ( -rho*Y[jhe4]*Y[jna23]*rate_eval.He4_Na23__p_Mg26 +rho*Y[jp]*Y[jmg26]*rate_eval.p_Mg26__He4_Na23 ) +
          ( +rho*Y[jn]*Y[jmg23]*rate_eval.n_Mg23__p_Na23 -rho*Y[jp]*Y[jna23]*rate_eval.p_Na23__n_Mg23 ) +
          ( +rho*Y[jn]*Y[jal26]*rate_eval.n_Al26__He4_Na23 -rho*Y[jhe4]*Y[jna23]*rate_eval.He4_Na23__n_Al26 )
       )

    dYdt[jna24] = (
          -Y[jna24]*rate_eval.Na24__Mg24__weak__wc12  +
          ( +rho*Y[jn]*Y[jna23]*rate_eval.n_Na23__Na24 -Y[jna24]*rate_eval.Na24__n_Na23 ) +
          ( -rho*Y[jp]*Y[jna24]*rate_eval.p_Na24__Mg25 +Y[jmg25]*rate_eval.Mg25__p_Na24 ) +
          ( -rho*Y[jhe4]*Y[jna24]*rate_eval.He4_Na24__Al28 +Y[jal28]*rate_eval.Al28__He4_Na24 ) +
          ( -rho*Y[jp]*Y[jna24]*rate_eval.p_Na24__n_Mg24 +rho*Y[jn]*Y[jmg24]*rate_eval.n_Mg24__p_Na24 ) +
          ( -rho*Y[jp]*Y[jna24]*rate_eval.p_Na24__He4_Ne21 +rho*Y[jhe4]*Y[jne21]*rate_eval.He4_Ne21__p_Na24 ) +
          ( -rho*Y[jhe4]*Y[jna24]*rate_eval.He4_Na24__n_Al27 +rho*Y[jn]*Y[jal27]*rate_eval.n_Al27__He4_Na24 )
       )

    dYdt[jmg22] = (
          -Y[jmg22]*rate_eval.Mg22__Na22__weak__wc12  +
          ( +rho*Y[jhe4]*Y[jne18]*rate_eval.He4_Ne18__Mg22 -Y[jmg22]*rate_eval.Mg22__He4_Ne18 ) +
          ( +rho*Y[jp]*Y[jna21]*rate_eval.p_Na21__Mg22 -Y[jmg22]*rate_eval.Mg22__p_Na21 ) +
          ( -rho*Y[jn]*Y[jmg22]*rate_eval.n_Mg22__Mg23 +Y[jmg23]*rate_eval.Mg23__n_Mg22 ) +
          ( -rho*Y[jhe4]*Y[jmg22]*rate_eval.He4_Mg22__Si26 +Y[jsi26]*rate_eval.Si26__He4_Mg22 ) +
          ( -rho*Y[jn]*Y[jmg22]*rate_eval.n_Mg22__p_Na22 +rho*Y[jp]*Y[jna22]*rate_eval.p_Na22__n_Mg22 ) +
          ( -rho*Y[jn]*Y[jmg22]*rate_eval.n_Mg22__He4_Ne19 +rho*Y[jhe4]*Y[jne19]*rate_eval.He4_Ne19__n_Mg22 ) +
          ( -rho*Y[jhe4]*Y[jmg22]*rate_eval.He4_Mg22__p_Al25 +rho*Y[jp]*Y[jal25]*rate_eval.p_Al25__He4_Mg22 )
       )

    dYdt[jmg23] = (
          -Y[jmg23]*rate_eval.Mg23__Na23__weak__wc12  +
          ( +rho*Y[jhe4]*Y[jne19]*rate_eval.He4_Ne19__Mg23 -Y[jmg23]*rate_eval.Mg23__He4_Ne19 ) +
          ( +rho*Y[jp]*Y[jna22]*rate_eval.p_Na22__Mg23 -Y[jmg23]*rate_eval.Mg23__p_Na22 ) +
          ( +rho*Y[jn]*Y[jmg22]*rate_eval.n_Mg22__Mg23 -Y[jmg23]*rate_eval.Mg23__n_Mg22 ) +
          ( -rho*Y[jn]*Y[jmg23]*rate_eval.n_Mg23__Mg24 +Y[jmg24]*rate_eval.Mg24__n_Mg23 ) +
          ( -rho*Y[jn]*Y[jmg23]*rate_eval.n_Mg23__p_Na23 +rho*Y[jp]*Y[jna23]*rate_eval.p_Na23__n_Mg23 ) +
          ( -rho*Y[jn]*Y[jmg23]*rate_eval.n_Mg23__He4_Ne20 +rho*Y[jhe4]*Y[jne20]*rate_eval.He4_Ne20__n_Mg23 ) +
          ( -rho*Y[jn]*Y[jmg23]*rate_eval.n_Mg23__C12_C12 +5.00000000000000e-01*rho*Y[jc12]**2*rate_eval.C12_C12__n_Mg23 ) +
          ( -rho*Y[jhe4]*Y[jmg23]*rate_eval.He4_Mg23__p_Al26 +rho*Y[jp]*Y[jal26]*rate_eval.p_Al26__He4_Mg23 ) +
          ( +rho*Y[jn]*Y[jsi26]*rate_eval.n_Si26__He4_Mg23 -rho*Y[jhe4]*Y[jmg23]*rate_eval.He4_Mg23__n_Si26 )
       )

    dYdt[jmg24] = (
          +Y[jna24]*rate_eval.Na24__Mg24__weak__wc12  +
          ( +rho*Y[jhe4]*Y[jne20]*rate_eval.He4_Ne20__Mg24 -Y[jmg24]*rate_eval.Mg24__He4_Ne20 ) +
          ( +rho*Y[jp]*Y[jna23]*rate_eval.p_Na23__Mg24 -Y[jmg24]*rate_eval.Mg24__p_Na23 ) +
          ( +rho*Y[jn]*Y[jmg23]*rate_eval.n_Mg23__Mg24 -Y[jmg24]*rate_eval.Mg24__n_Mg23 ) +
          ( -rho*Y[jn]*Y[jmg24]*rate_eval.n_Mg24__Mg25 +Y[jmg25]*rate_eval.Mg25__n_Mg24 ) +
          ( -rho*Y[jp]*Y[jmg24]*rate_eval.p_Mg24__Al25 +Y[jal25]*rate_eval.Al25__p_Mg24 ) +
          ( -rho*Y[jhe4]*Y[jmg24]*rate_eval.He4_Mg24__Si28 +Y[jsi28]*rate_eval.Si28__He4_Mg24 ) +
          ( +rho*Y[jc12]*Y[jo16]*rate_eval.C12_O16__He4_Mg24 -rho*Y[jhe4]*Y[jmg24]*rate_eval.He4_Mg24__C12_O16 ) +
          ( +rho*Y[jhe4]*Y[jne21]*rate_eval.He4_Ne21__n_Mg24 -rho*Y[jn]*Y[jmg24]*rate_eval.n_Mg24__He4_Ne21 ) +
          ( +rho*Y[jhe4]*Y[jna21]*rate_eval.He4_Na21__p_Mg24 -rho*Y[jp]*Y[jmg24]*rate_eval.p_Mg24__He4_Na21 ) +
          ( +rho*Y[jp]*Y[jna24]*rate_eval.p_Na24__n_Mg24 -rho*Y[jn]*Y[jmg24]*rate_eval.n_Mg24__p_Na24 ) +
          ( +rho*Y[jp]*Y[jal27]*rate_eval.p_Al27__He4_Mg24 -rho*Y[jhe4]*Y[jmg24]*rate_eval.He4_Mg24__p_Al27 )
       )

    dYdt[jmg25] = (
          +Y[jal25]*rate_eval.Al25__Mg25__weak__wc12  +
          ( +rho*Y[jhe4]*Y[jne21]*rate_eval.He4_Ne21__Mg25 -Y[jmg25]*rate_eval.Mg25__He4_Ne21 ) +
          ( +rho*Y[jp]*Y[jna24]*rate_eval.p_Na24__Mg25 -Y[jmg25]*rate_eval.Mg25__p_Na24 ) +
          ( +rho*Y[jn]*Y[jmg24]*rate_eval.n_Mg24__Mg25 -Y[jmg25]*rate_eval.Mg25__n_Mg24 ) +
          ( -rho*Y[jn]*Y[jmg25]*rate_eval.n_Mg25__Mg26 +Y[jmg26]*rate_eval.Mg26__n_Mg25 ) +
          ( -rho*Y[jp]*Y[jmg25]*rate_eval.p_Mg25__Al26 +Y[jal26]*rate_eval.Al26__p_Mg25 ) +
          ( -rho*Y[jhe4]*Y[jmg25]*rate_eval.He4_Mg25__Si29 +Y[jsi29]*rate_eval.Si29__He4_Mg25 ) +
          ( +rho*Y[jhe4]*Y[jna22]*rate_eval.He4_Na22__p_Mg25 -rho*Y[jp]*Y[jmg25]*rate_eval.p_Mg25__He4_Na22 ) +
          ( -rho*Y[jhe4]*Y[jmg25]*rate_eval.He4_Mg25__n_Si28 +rho*Y[jn]*Y[jsi28]*rate_eval.n_Si28__He4_Mg25 ) +
          ( +rho*Y[jn]*Y[jal25]*rate_eval.n_Al25__p_Mg25 -rho*Y[jp]*Y[jmg25]*rate_eval.p_Mg25__n_Al25 ) +
          ( +rho*Y[jp]*Y[jal28]*rate_eval.p_Al28__He4_Mg25 -rho*Y[jhe4]*Y[jmg25]*rate_eval.He4_Mg25__p_Al28 )
       )

    dYdt[jmg26] = (
          +Y[jal26]*rate_eval.Al26__Mg26__weak__wc12  +
          ( +rho*Y[jn]*Y[jmg25]*rate_eval.n_Mg25__Mg26 -Y[jmg26]*rate_eval.Mg26__n_Mg25 ) +
          ( -rho*Y[jp]*Y[jmg26]*rate_eval.p_Mg26__Al27 +Y[jal27]*rate_eval.Al27__p_Mg26 ) +
          ( -rho*Y[jhe4]*Y[jmg26]*rate_eval.He4_Mg26__Si30 +Y[jsi30]*rate_eval.Si30__He4_Mg26 ) +
          ( +rho*Y[jhe4]*Y[jna23]*rate_eval.He4_Na23__p_Mg26 -rho*Y[jp]*Y[jmg26]*rate_eval.p_Mg26__He4_Na23 ) +
          ( -rho*Y[jhe4]*Y[jmg26]*rate_eval.He4_Mg26__n_Si29 +rho*Y[jn]*Y[jsi29]*rate_eval.n_Si29__He4_Mg26 ) +
          ( +rho*Y[jn]*Y[jal26]*rate_eval.n_Al26__p_Mg26 -rho*Y[jp]*Y[jmg26]*rate_eval.p_Mg26__n_Al26 )
       )

    dYdt[jal25] = (
          -Y[jal25]*rate_eval.Al25__Mg25__weak__wc12  +
          ( +rho*Y[jhe4]*Y[jna21]*rate_eval.He4_Na21__Al25 -Y[jal25]*rate_eval.Al25__He4_Na21 ) +
          ( +rho*Y[jp]*Y[jmg24]*rate_eval.p_Mg24__Al25 -Y[jal25]*rate_eval.Al25__p_Mg24 ) +
          ( -rho*Y[jn]*Y[jal25]*rate_eval.n_Al25__Al26 +Y[jal26]*rate_eval.Al26__n_Al25 ) +
          ( -rho*Y[jp]*Y[jal25]*rate_eval.p_Al25__Si26 +Y[jsi26]*rate_eval.Si26__p_Al25 ) +
          ( -rho*Y[jhe4]*Y[jal25]*rate_eval.He4_Al25__P29 +Y[jp29]*rate_eval.P29__He4_Al25 ) +
          ( +rho*Y[jhe4]*Y[jmg22]*rate_eval.He4_Mg22__p_Al25 -rho*Y[jp]*Y[jal25]*rate_eval.p_Al25__He4_Mg22 ) +
          ( -rho*Y[jn]*Y[jal25]*rate_eval.n_Al25__p_Mg25 +rho*Y[jp]*Y[jmg25]*rate_eval.p_Mg25__n_Al25 ) +
          ( -rho*Y[jn]*Y[jal25]*rate_eval.n_Al25__He4_Na22 +rho*Y[jhe4]*Y[jna22]*rate_eval.He4_Na22__n_Al25 ) +
          ( -rho*Y[jhe4]*Y[jal25]*rate_eval.He4_Al25__p_Si28 +rho*Y[jp]*Y[jsi28]*rate_eval.p_Si28__He4_Al25 )
       )

    dYdt[jal26] = (
          -Y[jal26]*rate_eval.Al26__Mg26__weak__wc12  +
          +Y[jsi26]*rate_eval.Si26__Al26__weak__wc12  +
          ( +rho*Y[jhe4]*Y[jna22]*rate_eval.He4_Na22__Al26 -Y[jal26]*rate_eval.Al26__He4_Na22 ) +
          ( +rho*Y[jp]*Y[jmg25]*rate_eval.p_Mg25__Al26 -Y[jal26]*rate_eval.Al26__p_Mg25 ) +
          ( +rho*Y[jn]*Y[jal25]*rate_eval.n_Al25__Al26 -Y[jal26]*rate_eval.Al26__n_Al25 ) +
          ( -rho*Y[jn]*Y[jal26]*rate_eval.n_Al26__Al27 +Y[jal27]*rate_eval.Al27__n_Al26 ) +
          ( -rho*Y[jhe4]*Y[jal26]*rate_eval.He4_Al26__P30 +Y[jp30]*rate_eval.P30__He4_Al26 ) +
          ( +rho*Y[jhe4]*Y[jmg23]*rate_eval.He4_Mg23__p_Al26 -rho*Y[jp]*Y[jal26]*rate_eval.p_Al26__He4_Mg23 ) +
          ( -rho*Y[jn]*Y[jal26]*rate_eval.n_Al26__p_Mg26 +rho*Y[jp]*Y[jmg26]*rate_eval.p_Mg26__n_Al26 ) +
          ( -rho*Y[jn]*Y[jal26]*rate_eval.n_Al26__He4_Na23 +rho*Y[jhe4]*Y[jna23]*rate_eval.He4_Na23__n_Al26 ) +
          ( -rho*Y[jhe4]*Y[jal26]*rate_eval.He4_Al26__p_Si29 +rho*Y[jp]*Y[jsi29]*rate_eval.p_Si29__He4_Al26 ) +
          ( +rho*Y[jn]*Y[jsi26]*rate_eval.n_Si26__p_Al26 -rho*Y[jp]*Y[jal26]*rate_eval.p_Al26__n_Si26 ) +
          ( +rho*Y[jn]*Y[jp29]*rate_eval.n_P29__He4_Al26 -rho*Y[jhe4]*Y[jal26]*rate_eval.He4_Al26__n_P29 )
       )

    dYdt[jal27] = (
          ( +rho*Y[jhe4]*Y[jna23]*rate_eval.He4_Na23__Al27 -Y[jal27]*rate_eval.Al27__He4_Na23 ) +
          ( +rho*Y[jp]*Y[jmg26]*rate_eval.p_Mg26__Al27 -Y[jal27]*rate_eval.Al27__p_Mg26 ) +
          ( +rho*Y[jn]*Y[jal26]*rate_eval.n_Al26__Al27 -Y[jal27]*rate_eval.Al27__n_Al26 ) +
          ( -rho*Y[jn]*Y[jal27]*rate_eval.n_Al27__Al28 +Y[jal28]*rate_eval.Al28__n_Al27 ) +
          ( -rho*Y[jp]*Y[jal27]*rate_eval.p_Al27__Si28 +Y[jsi28]*rate_eval.Si28__p_Al27 ) +
          ( -rho*Y[jhe4]*Y[jal27]*rate_eval.He4_Al27__P31 +Y[jp31]*rate_eval.P31__He4_Al27 ) +
          ( +rho*Y[jc12]*Y[jo16]*rate_eval.C12_O16__p_Al27 -rho*Y[jp]*Y[jal27]*rate_eval.p_Al27__C12_O16 ) +
          ( +rho*Y[jhe4]*Y[jna24]*rate_eval.He4_Na24__n_Al27 -rho*Y[jn]*Y[jal27]*rate_eval.n_Al27__He4_Na24 ) +
          ( -rho*Y[jp]*Y[jal27]*rate_eval.p_Al27__He4_Mg24 +rho*Y[jhe4]*Y[jmg24]*rate_eval.He4_Mg24__p_Al27 ) +
          ( -rho*Y[jhe4]*Y[jal27]*rate_eval.He4_Al27__p_Si30 +rho*Y[jp]*Y[jsi30]*rate_eval.p_Si30__He4_Al27 ) +
          ( +rho*Y[jn]*Y[jp30]*rate_eval.n_P30__He4_Al27 -rho*Y[jhe4]*Y[jal27]*rate_eval.He4_Al27__n_P30 )
       )

    dYdt[jal28] = (
          -Y[jal28]*rate_eval.Al28__Si28__weak__wc12  +
          ( +rho*Y[jhe4]*Y[jna24]*rate_eval.He4_Na24__Al28 -Y[jal28]*rate_eval.Al28__He4_Na24 ) +
          ( +rho*Y[jn]*Y[jal27]*rate_eval.n_Al27__Al28 -Y[jal28]*rate_eval.Al28__n_Al27 ) +
          ( -rho*Y[jp]*Y[jal28]*rate_eval.p_Al28__Si29 +Y[jsi29]*rate_eval.Si29__p_Al28 ) +
          ( -rho*Y[jhe4]*Y[jal28]*rate_eval.He4_Al28__P32 +Y[jp32]*rate_eval.P32__He4_Al28 ) +
          ( -rho*Y[jp]*Y[jal28]*rate_eval.p_Al28__n_Si28 +rho*Y[jn]*Y[jsi28]*rate_eval.n_Si28__p_Al28 ) +
          ( -rho*Y[jp]*Y[jal28]*rate_eval.p_Al28__He4_Mg25 +rho*Y[jhe4]*Y[jmg25]*rate_eval.He4_Mg25__p_Al28 ) +
          ( -rho*Y[jhe4]*Y[jal28]*rate_eval.He4_Al28__n_P31 +rho*Y[jn]*Y[jp31]*rate_eval.n_P31__He4_Al28 )
       )

    dYdt[jsi26] = (
          -Y[jsi26]*rate_eval.Si26__Al26__weak__wc12  +
          ( +rho*Y[jhe4]*Y[jmg22]*rate_eval.He4_Mg22__Si26 -Y[jsi26]*rate_eval.Si26__He4_Mg22 ) +
          ( +rho*Y[jp]*Y[jal25]*rate_eval.p_Al25__Si26 -Y[jsi26]*rate_eval.Si26__p_Al25 ) +
          ( -rho*Y[jn]*Y[jsi26]*rate_eval.n_Si26__p_Al26 +rho*Y[jp]*Y[jal26]*rate_eval.p_Al26__n_Si26 ) +
          ( -rho*Y[jn]*Y[jsi26]*rate_eval.n_Si26__He4_Mg23 +rho*Y[jhe4]*Y[jmg23]*rate_eval.He4_Mg23__n_Si26 ) +
          ( -rho*Y[jhe4]*Y[jsi26]*rate_eval.He4_Si26__p_P29 +rho*Y[jp]*Y[jp29]*rate_eval.p_P29__He4_Si26 )
       )

    dYdt[jsi28] = (
          +Y[jal28]*rate_eval.Al28__Si28__weak__wc12  +
          +Y[jcl32]*rate_eval.Cl32__He4_Si28__weak__wc12  +
          ( +rho*Y[jhe4]*Y[jmg24]*rate_eval.He4_Mg24__Si28 -Y[jsi28]*rate_eval.Si28__He4_Mg24 ) +
          ( +rho*Y[jp]*Y[jal27]*rate_eval.p_Al27__Si28 -Y[jsi28]*rate_eval.Si28__p_Al27 ) +
          ( -rho*Y[jn]*Y[jsi28]*rate_eval.n_Si28__Si29 +Y[jsi29]*rate_eval.Si29__n_Si28 ) +
          ( -rho*Y[jp]*Y[jsi28]*rate_eval.p_Si28__P29 +Y[jp29]*rate_eval.P29__p_Si28 ) +
          ( -rho*Y[jhe4]*Y[jsi28]*rate_eval.He4_Si28__S32 +Y[js32]*rate_eval.S32__He4_Si28 ) +
          ( +5.00000000000000e-01*rho*Y[jo16]**2*rate_eval.O16_O16__He4_Si28 -rho*Y[jhe4]*Y[jsi28]*rate_eval.He4_Si28__O16_O16 ) +
          ( +rho*Y[jc12]*Y[jne20]*rate_eval.C12_Ne20__He4_Si28 -rho*Y[jhe4]*Y[jsi28]*rate_eval.He4_Si28__C12_Ne20 ) +
          ( +rho*Y[jhe4]*Y[jmg25]*rate_eval.He4_Mg25__n_Si28 -rho*Y[jn]*Y[jsi28]*rate_eval.n_Si28__He4_Mg25 ) +
          ( +rho*Y[jhe4]*Y[jal25]*rate_eval.He4_Al25__p_Si28 -rho*Y[jp]*Y[jsi28]*rate_eval.p_Si28__He4_Al25 ) +
          ( +rho*Y[jp]*Y[jal28]*rate_eval.p_Al28__n_Si28 -rho*Y[jn]*Y[jsi28]*rate_eval.n_Si28__p_Al28 ) +
          ( +rho*Y[jp]*Y[jp31]*rate_eval.p_P31__He4_Si28 -rho*Y[jhe4]*Y[jsi28]*rate_eval.He4_Si28__p_P31 ) +
          ( +rho*Y[jn]*Y[js31]*rate_eval.n_S31__He4_Si28 -rho*Y[jhe4]*Y[jsi28]*rate_eval.He4_Si28__n_S31 )
       )

    dYdt[jsi29] = (
          +Y[jp29]*rate_eval.P29__Si29__weak__wc12  +
          ( +rho*Y[jhe4]*Y[jmg25]*rate_eval.He4_Mg25__Si29 -Y[jsi29]*rate_eval.Si29__He4_Mg25 ) +
          ( +rho*Y[jp]*Y[jal28]*rate_eval.p_Al28__Si29 -Y[jsi29]*rate_eval.Si29__p_Al28 ) +
          ( +rho*Y[jn]*Y[jsi28]*rate_eval.n_Si28__Si29 -Y[jsi29]*rate_eval.Si29__n_Si28 ) +
          ( -rho*Y[jn]*Y[jsi29]*rate_eval.n_Si29__Si30 +Y[jsi30]*rate_eval.Si30__n_Si29 ) +
          ( -rho*Y[jp]*Y[jsi29]*rate_eval.p_Si29__P30 +Y[jp30]*rate_eval.P30__p_Si29 ) +
          ( +rho*Y[jhe4]*Y[jmg26]*rate_eval.He4_Mg26__n_Si29 -rho*Y[jn]*Y[jsi29]*rate_eval.n_Si29__He4_Mg26 ) +
          ( +rho*Y[jhe4]*Y[jal26]*rate_eval.He4_Al26__p_Si29 -rho*Y[jp]*Y[jsi29]*rate_eval.p_Si29__He4_Al26 ) +
          ( +rho*Y[jn]*Y[jp29]*rate_eval.n_P29__p_Si29 -rho*Y[jp]*Y[jsi29]*rate_eval.p_Si29__n_P29 ) +
          ( +rho*Y[jp]*Y[jp32]*rate_eval.p_P32__He4_Si29 -rho*Y[jhe4]*Y[jsi29]*rate_eval.He4_Si29__p_P32 ) +
          ( +rho*Y[jn]*Y[js32]*rate_eval.n_S32__He4_Si29 -rho*Y[jhe4]*Y[jsi29]*rate_eval.He4_Si29__n_S32 )
       )

    dYdt[jsi30] = (
          +Y[jp30]*rate_eval.P30__Si30__weak__wc12  +
          ( +rho*Y[jhe4]*Y[jmg26]*rate_eval.He4_Mg26__Si30 -Y[jsi30]*rate_eval.Si30__He4_Mg26 ) +
          ( +rho*Y[jn]*Y[jsi29]*rate_eval.n_Si29__Si30 -Y[jsi30]*rate_eval.Si30__n_Si29 ) +
          ( -rho*Y[jp]*Y[jsi30]*rate_eval.p_Si30__P31 +Y[jp31]*rate_eval.P31__p_Si30 ) +
          ( -rho*Y[jhe4]*Y[jsi30]*rate_eval.He4_Si30__S34 +Y[js34]*rate_eval.S34__He4_Si30 ) +
          ( +rho*Y[jhe4]*Y[jal27]*rate_eval.He4_Al27__p_Si30 -rho*Y[jp]*Y[jsi30]*rate_eval.p_Si30__He4_Al27 ) +
          ( +rho*Y[jn]*Y[jp30]*rate_eval.n_P30__p_Si30 -rho*Y[jp]*Y[jsi30]*rate_eval.p_Si30__n_P30 )
       )

    dYdt[jp29] = (
          -Y[jp29]*rate_eval.P29__Si29__weak__wc12  +
          ( +rho*Y[jhe4]*Y[jal25]*rate_eval.He4_Al25__P29 -Y[jp29]*rate_eval.P29__He4_Al25 ) +
          ( +rho*Y[jp]*Y[jsi28]*rate_eval.p_Si28__P29 -Y[jp29]*rate_eval.P29__p_Si28 ) +
          ( -rho*Y[jn]*Y[jp29]*rate_eval.n_P29__P30 +Y[jp30]*rate_eval.P30__n_P29 ) +
          ( +rho*Y[jhe4]*Y[jsi26]*rate_eval.He4_Si26__p_P29 -rho*Y[jp]*Y[jp29]*rate_eval.p_P29__He4_Si26 ) +
          ( -rho*Y[jn]*Y[jp29]*rate_eval.n_P29__p_Si29 +rho*Y[jp]*Y[jsi29]*rate_eval.p_Si29__n_P29 ) +
          ( -rho*Y[jn]*Y[jp29]*rate_eval.n_P29__He4_Al26 +rho*Y[jhe4]*Y[jal26]*rate_eval.He4_Al26__n_P29 ) +
          ( -rho*Y[jhe4]*Y[jp29]*rate_eval.He4_P29__p_S32 +rho*Y[jp]*Y[js32]*rate_eval.p_S32__He4_P29 ) +
          ( +rho*Y[jn]*Y[jcl32]*rate_eval.n_Cl32__He4_P29 -rho*Y[jhe4]*Y[jp29]*rate_eval.He4_P29__n_Cl32 )
       )

    dYdt[jp30] = (
          -Y[jp30]*rate_eval.P30__Si30__weak__wc12  +
          ( +rho*Y[jhe4]*Y[jal26]*rate_eval.He4_Al26__P30 -Y[jp30]*rate_eval.P30__He4_Al26 ) +
          ( +rho*Y[jp]*Y[jsi29]*rate_eval.p_Si29__P30 -Y[jp30]*rate_eval.P30__p_Si29 ) +
          ( +rho*Y[jn]*Y[jp29]*rate_eval.n_P29__P30 -Y[jp30]*rate_eval.P30__n_P29 ) +
          ( -rho*Y[jn]*Y[jp30]*rate_eval.n_P30__P31 +Y[jp31]*rate_eval.P31__n_P30 ) +
          ( -rho*Y[jp]*Y[jp30]*rate_eval.p_P30__S31 +Y[js31]*rate_eval.S31__p_P30 ) +
          ( -rho*Y[jhe4]*Y[jp30]*rate_eval.He4_P30__Cl34 +Y[jcl34]*rate_eval.Cl34__He4_P30 ) +
          ( -rho*Y[jn]*Y[jp30]*rate_eval.n_P30__p_Si30 +rho*Y[jp]*Y[jsi30]*rate_eval.p_Si30__n_P30 ) +
          ( -rho*Y[jn]*Y[jp30]*rate_eval.n_P30__He4_Al27 +rho*Y[jhe4]*Y[jal27]*rate_eval.He4_Al27__n_P30 )
       )

    dYdt[jp31] = (
          +Y[js31]*rate_eval.S31__P31__weak__wc12  +
          +Y[jcl32]*rate_eval.Cl32__p_P31__weak__wc12  +
          ( +rho*Y[jhe4]*Y[jal27]*rate_eval.He4_Al27__P31 -Y[jp31]*rate_eval.P31__He4_Al27 ) +
          ( +rho*Y[jp]*Y[jsi30]*rate_eval.p_Si30__P31 -Y[jp31]*rate_eval.P31__p_Si30 ) +
          ( +rho*Y[jn]*Y[jp30]*rate_eval.n_P30__P31 -Y[jp31]*rate_eval.P31__n_P30 ) +
          ( -rho*Y[jn]*Y[jp31]*rate_eval.n_P31__P32 +Y[jp32]*rate_eval.P32__n_P31 ) +
          ( -rho*Y[jp]*Y[jp31]*rate_eval.p_P31__S32 +Y[js32]*rate_eval.S32__p_P31 ) +
          ( -rho*Y[jhe4]*Y[jp31]*rate_eval.He4_P31__Cl35 +Y[jcl35]*rate_eval.Cl35__He4_P31 ) +
          ( +5.00000000000000e-01*rho*Y[jo16]**2*rate_eval.O16_O16__p_P31 -rho*Y[jp]*Y[jp31]*rate_eval.p_P31__O16_O16 ) +
          ( +rho*Y[jc12]*Y[jne20]*rate_eval.C12_Ne20__p_P31 -rho*Y[jp]*Y[jp31]*rate_eval.p_P31__C12_Ne20 ) +
          ( +rho*Y[jhe4]*Y[jal28]*rate_eval.He4_Al28__n_P31 -rho*Y[jn]*Y[jp31]*rate_eval.n_P31__He4_Al28 ) +
          ( -rho*Y[jp]*Y[jp31]*rate_eval.p_P31__He4_Si28 +rho*Y[jhe4]*Y[jsi28]*rate_eval.He4_Si28__p_P31 ) +
          ( -rho*Y[jhe4]*Y[jp31]*rate_eval.He4_P31__p_S34 +rho*Y[jp]*Y[js34]*rate_eval.p_S34__He4_P31 ) +
          ( +rho*Y[jn]*Y[js31]*rate_eval.n_S31__p_P31 -rho*Y[jp]*Y[jp31]*rate_eval.p_P31__n_S31 ) +
          ( +rho*Y[jn]*Y[jcl34]*rate_eval.n_Cl34__He4_P31 -rho*Y[jhe4]*Y[jp31]*rate_eval.He4_P31__n_Cl34 )
       )

    dYdt[jp32] = (
          -Y[jp32]*rate_eval.P32__S32__weak__wc12  +
          ( +rho*Y[jhe4]*Y[jal28]*rate_eval.He4_Al28__P32 -Y[jp32]*rate_eval.P32__He4_Al28 ) +
          ( +rho*Y[jn]*Y[jp31]*rate_eval.n_P31__P32 -Y[jp32]*rate_eval.P32__n_P31 ) +
          ( -rho*Y[jp]*Y[jp32]*rate_eval.p_P32__n_S32 +rho*Y[jn]*Y[js32]*rate_eval.n_S32__p_P32 ) +
          ( -rho*Y[jp]*Y[jp32]*rate_eval.p_P32__He4_Si29 +rho*Y[jhe4]*Y[jsi29]*rate_eval.He4_Si29__p_P32 ) +
          ( +rho*Y[jn]*Y[jcl35]*rate_eval.n_Cl35__He4_P32 -rho*Y[jhe4]*Y[jp32]*rate_eval.He4_P32__n_Cl35 )
       )

    dYdt[js31] = (
          -Y[js31]*rate_eval.S31__P31__weak__wc12  +
          ( +rho*Y[jp]*Y[jp30]*rate_eval.p_P30__S31 -Y[js31]*rate_eval.S31__p_P30 ) +
          ( -rho*Y[jn]*Y[js31]*rate_eval.n_S31__S32 +Y[js32]*rate_eval.S32__n_S31 ) +
          ( -rho*Y[jp]*Y[js31]*rate_eval.p_S31__Cl32 +Y[jcl32]*rate_eval.Cl32__p_S31 ) +
          ( -rho*Y[jhe4]*Y[js31]*rate_eval.He4_S31__Ar35 +Y[jar35]*rate_eval.Ar35__He4_S31 ) +
          ( +5.00000000000000e-01*rho*Y[jo16]**2*rate_eval.O16_O16__n_S31 -rho*Y[jn]*Y[js31]*rate_eval.n_S31__O16_O16 ) +
          ( +rho*Y[jc12]*Y[jne20]*rate_eval.C12_Ne20__n_S31 -rho*Y[jn]*Y[js31]*rate_eval.n_S31__C12_Ne20 ) +
          ( -rho*Y[jn]*Y[js31]*rate_eval.n_S31__p_P31 +rho*Y[jp]*Y[jp31]*rate_eval.p_P31__n_S31 ) +
          ( -rho*Y[jn]*Y[js31]*rate_eval.n_S31__He4_Si28 +rho*Y[jhe4]*Y[jsi28]*rate_eval.He4_Si28__n_S31 ) +
          ( -rho*Y[jhe4]*Y[js31]*rate_eval.He4_S31__p_Cl34 +rho*Y[jp]*Y[jcl34]*rate_eval.p_Cl34__He4_S31 )
       )

    dYdt[js32] = (
          +Y[jp32]*rate_eval.P32__S32__weak__wc12  +
          +Y[jcl32]*rate_eval.Cl32__S32__weak__wc12  +
          ( +rho*Y[jhe4]*Y[jsi28]*rate_eval.He4_Si28__S32 -Y[js32]*rate_eval.S32__He4_Si28 ) +
          ( +rho*Y[jp]*Y[jp31]*rate_eval.p_P31__S32 -Y[js32]*rate_eval.S32__p_P31 ) +
          ( +rho*Y[jn]*Y[js31]*rate_eval.n_S31__S32 -Y[js32]*rate_eval.S32__n_S31 ) +
          ( +rho*Y[jhe4]*Y[jp29]*rate_eval.He4_P29__p_S32 -rho*Y[jp]*Y[js32]*rate_eval.p_S32__He4_P29 ) +
          ( +rho*Y[jp]*Y[jp32]*rate_eval.p_P32__n_S32 -rho*Y[jn]*Y[js32]*rate_eval.n_S32__p_P32 ) +
          ( -rho*Y[jn]*Y[js32]*rate_eval.n_S32__He4_Si29 +rho*Y[jhe4]*Y[jsi29]*rate_eval.He4_Si29__n_S32 ) +
          ( +rho*Y[jn]*Y[jcl32]*rate_eval.n_Cl32__p_S32 -rho*Y[jp]*Y[js32]*rate_eval.p_S32__n_Cl32 ) +
          ( +rho*Y[jp]*Y[jcl35]*rate_eval.p_Cl35__He4_S32 -rho*Y[jhe4]*Y[js32]*rate_eval.He4_S32__p_Cl35 ) +
          ( +rho*Y[jn]*Y[jar35]*rate_eval.n_Ar35__He4_S32 -rho*Y[jhe4]*Y[js32]*rate_eval.He4_S32__n_Ar35 )
       )

    dYdt[js34] = (
          +Y[jcl34]*rate_eval.Cl34__S34__weak__wc12  +
          ( +rho*Y[jhe4]*Y[jsi30]*rate_eval.He4_Si30__S34 -Y[js34]*rate_eval.S34__He4_Si30 ) +
          ( -rho*Y[jp]*Y[js34]*rate_eval.p_S34__Cl35 +Y[jcl35]*rate_eval.Cl35__p_S34 ) +
          ( +rho*Y[jhe4]*Y[jp31]*rate_eval.He4_P31__p_S34 -rho*Y[jp]*Y[js34]*rate_eval.p_S34__He4_P31 ) +
          ( +rho*Y[jn]*Y[jcl34]*rate_eval.n_Cl34__p_S34 -rho*Y[jp]*Y[js34]*rate_eval.p_S34__n_Cl34 )
       )

    dYdt[jcl32] = (
          -Y[jcl32]*rate_eval.Cl32__S32__weak__wc12  +
          -Y[jcl32]*rate_eval.Cl32__p_P31__weak__wc12  +
          -Y[jcl32]*rate_eval.Cl32__He4_Si28__weak__wc12  +
          ( +rho*Y[jp]*Y[js31]*rate_eval.p_S31__Cl32 -Y[jcl32]*rate_eval.Cl32__p_S31 ) +
          ( -rho*Y[jn]*Y[jcl32]*rate_eval.n_Cl32__p_S32 +rho*Y[jp]*Y[js32]*rate_eval.p_S32__n_Cl32 ) +
          ( -rho*Y[jn]*Y[jcl32]*rate_eval.n_Cl32__He4_P29 +rho*Y[jhe4]*Y[jp29]*rate_eval.He4_P29__n_Cl32 ) +
          ( -rho*Y[jhe4]*Y[jcl32]*rate_eval.He4_Cl32__p_Ar35 +rho*Y[jp]*Y[jar35]*rate_eval.p_Ar35__He4_Cl32 )
       )

    dYdt[jcl34] = (
          -Y[jcl34]*rate_eval.Cl34__S34__weak__wc12  +
          ( +rho*Y[jhe4]*Y[jp30]*rate_eval.He4_P30__Cl34 -Y[jcl34]*rate_eval.Cl34__He4_P30 ) +
          ( -rho*Y[jn]*Y[jcl34]*rate_eval.n_Cl34__Cl35 +Y[jcl35]*rate_eval.Cl35__n_Cl34 ) +
          ( -rho*Y[jp]*Y[jcl34]*rate_eval.p_Cl34__Ar35 +Y[jar35]*rate_eval.Ar35__p_Cl34 ) +
          ( +rho*Y[jhe4]*Y[js31]*rate_eval.He4_S31__p_Cl34 -rho*Y[jp]*Y[jcl34]*rate_eval.p_Cl34__He4_S31 ) +
          ( -rho*Y[jn]*Y[jcl34]*rate_eval.n_Cl34__p_S34 +rho*Y[jp]*Y[js34]*rate_eval.p_S34__n_Cl34 ) +
          ( -rho*Y[jn]*Y[jcl34]*rate_eval.n_Cl34__He4_P31 +rho*Y[jhe4]*Y[jp31]*rate_eval.He4_P31__n_Cl34 )
       )

    dYdt[jcl35] = (
          +Y[jar35]*rate_eval.Ar35__Cl35__weak__wc12  +
          ( +rho*Y[jhe4]*Y[jp31]*rate_eval.He4_P31__Cl35 -Y[jcl35]*rate_eval.Cl35__He4_P31 ) +
          ( +rho*Y[jp]*Y[js34]*rate_eval.p_S34__Cl35 -Y[jcl35]*rate_eval.Cl35__p_S34 ) +
          ( +rho*Y[jn]*Y[jcl34]*rate_eval.n_Cl34__Cl35 -Y[jcl35]*rate_eval.Cl35__n_Cl34 ) +
          ( -rho*Y[jn]*Y[jcl35]*rate_eval.n_Cl35__He4_P32 +rho*Y[jhe4]*Y[jp32]*rate_eval.He4_P32__n_Cl35 ) +
          ( -rho*Y[jp]*Y[jcl35]*rate_eval.p_Cl35__He4_S32 +rho*Y[jhe4]*Y[js32]*rate_eval.He4_S32__p_Cl35 ) +
          ( +rho*Y[jn]*Y[jar35]*rate_eval.n_Ar35__p_Cl35 -rho*Y[jp]*Y[jcl35]*rate_eval.p_Cl35__n_Ar35 )
       )

    dYdt[jar35] = (
          -Y[jar35]*rate_eval.Ar35__Cl35__weak__wc12  +
          ( +rho*Y[jhe4]*Y[js31]*rate_eval.He4_S31__Ar35 -Y[jar35]*rate_eval.Ar35__He4_S31 ) +
          ( +rho*Y[jp]*Y[jcl34]*rate_eval.p_Cl34__Ar35 -Y[jar35]*rate_eval.Ar35__p_Cl34 ) +
          ( +rho*Y[jhe4]*Y[jcl32]*rate_eval.He4_Cl32__p_Ar35 -rho*Y[jp]*Y[jar35]*rate_eval.p_Ar35__He4_Cl32 ) +
          ( -rho*Y[jn]*Y[jar35]*rate_eval.n_Ar35__p_Cl35 +rho*Y[jp]*Y[jcl35]*rate_eval.p_Cl35__n_Ar35 ) +
          ( -rho*Y[jn]*Y[jar35]*rate_eval.n_Ar35__He4_S32 +rho*Y[jhe4]*Y[js32]*rate_eval.He4_S32__n_Ar35 )
       )

    return dYdt

def jacobian(t, Y, rho, T, screen_func=None):
    return jacobian_eq(t, Y, rho, T, screen_func)

@numba.njit()
def jacobian_eq(t, Y, rho, T, screen_func):

    tf = Tfactors(T)
    rate_eval = RateEval()

    # reaclib rates
    n__p__weak__wc12(rate_eval, tf)
    t__He3__weak__wc12(rate_eval, tf)
    He3__t__weak__electron_capture(rate_eval, tf)
    He6__Li6__weak__wc12(rate_eval, tf)
    Be7__Li7__weak__electron_capture(rate_eval, tf)
    C14__N14__weak__wc12(rate_eval, tf)
    N13__C13__weak__wc12(rate_eval, tf)
    O14__N14__weak__wc12(rate_eval, tf)
    O15__N15__weak__wc12(rate_eval, tf)
    F17__O17__weak__wc12(rate_eval, tf)
    F18__O18__weak__wc12(rate_eval, tf)
    Ne18__F18__weak__wc12(rate_eval, tf)
    Ne19__F19__weak__wc12(rate_eval, tf)
    Na21__Ne21__weak__wc12(rate_eval, tf)
    Na24__Mg24__weak__wc12(rate_eval, tf)
    Mg22__Na22__weak__wc12(rate_eval, tf)
    Mg23__Na23__weak__wc12(rate_eval, tf)
    Al25__Mg25__weak__wc12(rate_eval, tf)
    Al26__Mg26__weak__wc12(rate_eval, tf)
    Al28__Si28__weak__wc12(rate_eval, tf)
    Si26__Al26__weak__wc12(rate_eval, tf)
    P29__Si29__weak__wc12(rate_eval, tf)
    P30__Si30__weak__wc12(rate_eval, tf)
    P32__S32__weak__wc12(rate_eval, tf)
    S31__P31__weak__wc12(rate_eval, tf)
    Cl32__S32__weak__wc12(rate_eval, tf)
    Cl34__S34__weak__wc12(rate_eval, tf)
    Ar35__Cl35__weak__wc12(rate_eval, tf)
    d__n_p(rate_eval, tf)
    t__n_d(rate_eval, tf)
    He3__p_d(rate_eval, tf)
    He4__n_He3(rate_eval, tf)
    He4__p_t(rate_eval, tf)
    He4__d_d(rate_eval, tf)
    Li6__He4_d(rate_eval, tf)
    Li7__n_Li6(rate_eval, tf)
    Li7__He4_t(rate_eval, tf)
    Be7__p_Li6(rate_eval, tf)
    Be7__He4_He3(rate_eval, tf)
    B8__p_Be7(rate_eval, tf)
    B8__He4_He4__weak__wc12(rate_eval, tf)
    C11__He4_Be7(rate_eval, tf)
    C12__n_C11(rate_eval, tf)
    C13__n_C12(rate_eval, tf)
    C14__n_C13(rate_eval, tf)
    N13__p_C12(rate_eval, tf)
    N14__n_N13(rate_eval, tf)
    N14__p_C13(rate_eval, tf)
    N15__n_N14(rate_eval, tf)
    N15__p_C14(rate_eval, tf)
    O14__p_N13(rate_eval, tf)
    O15__n_O14(rate_eval, tf)
    O15__p_N14(rate_eval, tf)
    O16__n_O15(rate_eval, tf)
    O16__p_N15(rate_eval, tf)
    O16__He4_C12(rate_eval, tf)
    O17__n_O16(rate_eval, tf)
    O18__n_O17(rate_eval, tf)
    O18__He4_C14(rate_eval, tf)
    F15__p_O14(rate_eval, tf)
    F16__n_F15(rate_eval, tf)
    F16__p_O15(rate_eval, tf)
    F17__n_F16(rate_eval, tf)
    F17__p_O16(rate_eval, tf)
    F18__n_F17(rate_eval, tf)
    F18__p_O17(rate_eval, tf)
    F18__He4_N14(rate_eval, tf)
    F19__n_F18(rate_eval, tf)
    F19__p_O18(rate_eval, tf)
    F19__He4_N15(rate_eval, tf)
    Ne18__p_F17(rate_eval, tf)
    Ne18__He4_O14(rate_eval, tf)
    Ne19__n_Ne18(rate_eval, tf)
    Ne19__p_F18(rate_eval, tf)
    Ne19__He4_O15(rate_eval, tf)
    Ne20__n_Ne19(rate_eval, tf)
    Ne20__p_F19(rate_eval, tf)
    Ne20__He4_O16(rate_eval, tf)
    Ne21__n_Ne20(rate_eval, tf)
    Ne21__He4_O17(rate_eval, tf)
    Na21__p_Ne20(rate_eval, tf)
    Na21__He4_F17(rate_eval, tf)
    Na22__n_Na21(rate_eval, tf)
    Na22__p_Ne21(rate_eval, tf)
    Na22__He4_F18(rate_eval, tf)
    Na23__n_Na22(rate_eval, tf)
    Na23__He4_F19(rate_eval, tf)
    Na24__n_Na23(rate_eval, tf)
    Mg22__p_Na21(rate_eval, tf)
    Mg22__He4_Ne18(rate_eval, tf)
    Mg23__n_Mg22(rate_eval, tf)
    Mg23__p_Na22(rate_eval, tf)
    Mg23__He4_Ne19(rate_eval, tf)
    Mg24__n_Mg23(rate_eval, tf)
    Mg24__p_Na23(rate_eval, tf)
    Mg24__He4_Ne20(rate_eval, tf)
    Mg25__n_Mg24(rate_eval, tf)
    Mg25__p_Na24(rate_eval, tf)
    Mg25__He4_Ne21(rate_eval, tf)
    Mg26__n_Mg25(rate_eval, tf)
    Al25__p_Mg24(rate_eval, tf)
    Al25__He4_Na21(rate_eval, tf)
    Al26__n_Al25(rate_eval, tf)
    Al26__p_Mg25(rate_eval, tf)
    Al26__He4_Na22(rate_eval, tf)
    Al27__n_Al26(rate_eval, tf)
    Al27__p_Mg26(rate_eval, tf)
    Al27__He4_Na23(rate_eval, tf)
    Al28__n_Al27(rate_eval, tf)
    Al28__He4_Na24(rate_eval, tf)
    Si26__p_Al25(rate_eval, tf)
    Si26__He4_Mg22(rate_eval, tf)
    Si28__p_Al27(rate_eval, tf)
    Si28__He4_Mg24(rate_eval, tf)
    Si29__n_Si28(rate_eval, tf)
    Si29__p_Al28(rate_eval, tf)
    Si29__He4_Mg25(rate_eval, tf)
    Si30__n_Si29(rate_eval, tf)
    Si30__He4_Mg26(rate_eval, tf)
    P29__p_Si28(rate_eval, tf)
    P29__He4_Al25(rate_eval, tf)
    P30__n_P29(rate_eval, tf)
    P30__p_Si29(rate_eval, tf)
    P30__He4_Al26(rate_eval, tf)
    P31__n_P30(rate_eval, tf)
    P31__p_Si30(rate_eval, tf)
    P31__He4_Al27(rate_eval, tf)
    P32__n_P31(rate_eval, tf)
    P32__He4_Al28(rate_eval, tf)
    S31__p_P30(rate_eval, tf)
    S32__n_S31(rate_eval, tf)
    S32__p_P31(rate_eval, tf)
    S32__He4_Si28(rate_eval, tf)
    S34__He4_Si30(rate_eval, tf)
    Cl32__p_S31(rate_eval, tf)
    Cl32__p_P31__weak__wc12(rate_eval, tf)
    Cl32__He4_Si28__weak__wc12(rate_eval, tf)
    Cl34__He4_P30(rate_eval, tf)
    Cl35__n_Cl34(rate_eval, tf)
    Cl35__p_S34(rate_eval, tf)
    Cl35__He4_P31(rate_eval, tf)
    Ar35__p_Cl34(rate_eval, tf)
    Ar35__He4_S31(rate_eval, tf)
    He6__n_n_He4(rate_eval, tf)
    Li6__n_p_He4(rate_eval, tf)
    Be9__n_He4_He4(rate_eval, tf)
    C12__He4_He4_He4(rate_eval, tf)
    n_p__d(rate_eval, tf)
    p_p__d__weak__bet_pos_(rate_eval, tf)
    p_p__d__weak__electron_capture(rate_eval, tf)
    n_d__t(rate_eval, tf)
    p_d__He3(rate_eval, tf)
    d_d__He4(rate_eval, tf)
    He4_d__Li6(rate_eval, tf)
    p_t__He4(rate_eval, tf)
    He4_t__Li7(rate_eval, tf)
    n_He3__He4(rate_eval, tf)
    p_He3__He4__weak__bet_pos_(rate_eval, tf)
    He4_He3__Be7(rate_eval, tf)
    n_Li6__Li7(rate_eval, tf)
    p_Li6__Be7(rate_eval, tf)
    p_Be7__B8(rate_eval, tf)
    He4_Be7__C11(rate_eval, tf)
    n_C11__C12(rate_eval, tf)
    n_C12__C13(rate_eval, tf)
    p_C12__N13(rate_eval, tf)
    He4_C12__O16(rate_eval, tf)
    n_C13__C14(rate_eval, tf)
    p_C13__N14(rate_eval, tf)
    p_C14__N15(rate_eval, tf)
    He4_C14__O18(rate_eval, tf)
    n_N13__N14(rate_eval, tf)
    p_N13__O14(rate_eval, tf)
    n_N14__N15(rate_eval, tf)
    p_N14__O15(rate_eval, tf)
    He4_N14__F18(rate_eval, tf)
    p_N15__O16(rate_eval, tf)
    He4_N15__F19(rate_eval, tf)
    n_O14__O15(rate_eval, tf)
    p_O14__F15(rate_eval, tf)
    He4_O14__Ne18(rate_eval, tf)
    n_O15__O16(rate_eval, tf)
    p_O15__F16(rate_eval, tf)
    He4_O15__Ne19(rate_eval, tf)
    n_O16__O17(rate_eval, tf)
    p_O16__F17(rate_eval, tf)
    He4_O16__Ne20(rate_eval, tf)
    n_O17__O18(rate_eval, tf)
    p_O17__F18(rate_eval, tf)
    He4_O17__Ne21(rate_eval, tf)
    p_O18__F19(rate_eval, tf)
    n_F15__F16(rate_eval, tf)
    n_F16__F17(rate_eval, tf)
    n_F17__F18(rate_eval, tf)
    p_F17__Ne18(rate_eval, tf)
    He4_F17__Na21(rate_eval, tf)
    n_F18__F19(rate_eval, tf)
    p_F18__Ne19(rate_eval, tf)
    He4_F18__Na22(rate_eval, tf)
    p_F19__Ne20(rate_eval, tf)
    He4_F19__Na23(rate_eval, tf)
    n_Ne18__Ne19(rate_eval, tf)
    He4_Ne18__Mg22(rate_eval, tf)
    n_Ne19__Ne20(rate_eval, tf)
    He4_Ne19__Mg23(rate_eval, tf)
    n_Ne20__Ne21(rate_eval, tf)
    p_Ne20__Na21(rate_eval, tf)
    He4_Ne20__Mg24(rate_eval, tf)
    p_Ne21__Na22(rate_eval, tf)
    He4_Ne21__Mg25(rate_eval, tf)
    n_Na21__Na22(rate_eval, tf)
    p_Na21__Mg22(rate_eval, tf)
    He4_Na21__Al25(rate_eval, tf)
    n_Na22__Na23(rate_eval, tf)
    p_Na22__Mg23(rate_eval, tf)
    He4_Na22__Al26(rate_eval, tf)
    n_Na23__Na24(rate_eval, tf)
    p_Na23__Mg24(rate_eval, tf)
    He4_Na23__Al27(rate_eval, tf)
    p_Na24__Mg25(rate_eval, tf)
    He4_Na24__Al28(rate_eval, tf)
    n_Mg22__Mg23(rate_eval, tf)
    He4_Mg22__Si26(rate_eval, tf)
    n_Mg23__Mg24(rate_eval, tf)
    n_Mg24__Mg25(rate_eval, tf)
    p_Mg24__Al25(rate_eval, tf)
    He4_Mg24__Si28(rate_eval, tf)
    n_Mg25__Mg26(rate_eval, tf)
    p_Mg25__Al26(rate_eval, tf)
    He4_Mg25__Si29(rate_eval, tf)
    p_Mg26__Al27(rate_eval, tf)
    He4_Mg26__Si30(rate_eval, tf)
    n_Al25__Al26(rate_eval, tf)
    p_Al25__Si26(rate_eval, tf)
    He4_Al25__P29(rate_eval, tf)
    n_Al26__Al27(rate_eval, tf)
    He4_Al26__P30(rate_eval, tf)
    n_Al27__Al28(rate_eval, tf)
    p_Al27__Si28(rate_eval, tf)
    He4_Al27__P31(rate_eval, tf)
    p_Al28__Si29(rate_eval, tf)
    He4_Al28__P32(rate_eval, tf)
    n_Si28__Si29(rate_eval, tf)
    p_Si28__P29(rate_eval, tf)
    He4_Si28__S32(rate_eval, tf)
    n_Si29__Si30(rate_eval, tf)
    p_Si29__P30(rate_eval, tf)
    p_Si30__P31(rate_eval, tf)
    He4_Si30__S34(rate_eval, tf)
    n_P29__P30(rate_eval, tf)
    n_P30__P31(rate_eval, tf)
    p_P30__S31(rate_eval, tf)
    He4_P30__Cl34(rate_eval, tf)
    n_P31__P32(rate_eval, tf)
    p_P31__S32(rate_eval, tf)
    He4_P31__Cl35(rate_eval, tf)
    n_S31__S32(rate_eval, tf)
    p_S31__Cl32(rate_eval, tf)
    He4_S31__Ar35(rate_eval, tf)
    p_S34__Cl35(rate_eval, tf)
    n_Cl34__Cl35(rate_eval, tf)
    p_Cl34__Ar35(rate_eval, tf)
    d_d__n_He3(rate_eval, tf)
    d_d__p_t(rate_eval, tf)
    p_t__n_He3(rate_eval, tf)
    p_t__d_d(rate_eval, tf)
    d_t__n_He4(rate_eval, tf)
    He4_t__n_Li6(rate_eval, tf)
    n_He3__p_t(rate_eval, tf)
    n_He3__d_d(rate_eval, tf)
    d_He3__p_He4(rate_eval, tf)
    t_He3__d_He4(rate_eval, tf)
    He4_He3__p_Li6(rate_eval, tf)
    n_He4__d_t(rate_eval, tf)
    p_He4__d_He3(rate_eval, tf)
    d_He4__t_He3(rate_eval, tf)
    He4_He4__n_Be7(rate_eval, tf)
    He4_He4__p_Li7(rate_eval, tf)
    n_Li6__He4_t(rate_eval, tf)
    p_Li6__He4_He3(rate_eval, tf)
    d_Li6__n_Be7(rate_eval, tf)
    d_Li6__p_Li7(rate_eval, tf)
    He4_Li6__p_Be9(rate_eval, tf)
    p_Li7__n_Be7(rate_eval, tf)
    p_Li7__d_Li6(rate_eval, tf)
    p_Li7__He4_He4(rate_eval, tf)
    t_Li7__n_Be9(rate_eval, tf)
    n_Be7__p_Li7(rate_eval, tf)
    n_Be7__d_Li6(rate_eval, tf)
    n_Be7__He4_He4(rate_eval, tf)
    n_Be9__t_Li7(rate_eval, tf)
    p_Be9__He4_Li6(rate_eval, tf)
    He4_Be9__n_C12(rate_eval, tf)
    He4_B8__p_C11(rate_eval, tf)
    p_C11__He4_B8(rate_eval, tf)
    He4_C11__n_O14(rate_eval, tf)
    He4_C11__p_N14(rate_eval, tf)
    n_C12__He4_Be9(rate_eval, tf)
    He4_C12__n_O15(rate_eval, tf)
    He4_C12__p_N15(rate_eval, tf)
    C12_C12__n_Mg23(rate_eval, tf)
    C12_C12__p_Na23(rate_eval, tf)
    C12_C12__He4_Ne20(rate_eval, tf)
    p_C13__n_N13(rate_eval, tf)
    d_C13__n_N14(rate_eval, tf)
    He4_C13__n_O16(rate_eval, tf)
    p_C14__n_N14(rate_eval, tf)
    d_C14__n_N15(rate_eval, tf)
    He4_C14__n_O17(rate_eval, tf)
    n_N13__p_C13(rate_eval, tf)
    He4_N13__n_F16(rate_eval, tf)
    He4_N13__p_O16(rate_eval, tf)
    n_N14__p_C14(rate_eval, tf)
    n_N14__d_C13(rate_eval, tf)
    p_N14__n_O14(rate_eval, tf)
    p_N14__He4_C11(rate_eval, tf)
    He4_N14__n_F17(rate_eval, tf)
    He4_N14__p_O17(rate_eval, tf)
    n_N15__d_C14(rate_eval, tf)
    p_N15__n_O15(rate_eval, tf)
    p_N15__He4_C12(rate_eval, tf)
    He4_N15__n_F18(rate_eval, tf)
    He4_N15__p_O18(rate_eval, tf)
    n_O14__p_N14(rate_eval, tf)
    n_O14__He4_C11(rate_eval, tf)
    He4_O14__p_F17(rate_eval, tf)
    n_O15__p_N15(rate_eval, tf)
    n_O15__He4_C12(rate_eval, tf)
    p_O15__n_F15(rate_eval, tf)
    He4_O15__n_Ne18(rate_eval, tf)
    He4_O15__p_F18(rate_eval, tf)
    n_O16__He4_C13(rate_eval, tf)
    p_O16__n_F16(rate_eval, tf)
    p_O16__He4_N13(rate_eval, tf)
    He4_O16__n_Ne19(rate_eval, tf)
    He4_O16__p_F19(rate_eval, tf)
    C12_O16__p_Al27(rate_eval, tf)
    C12_O16__He4_Mg24(rate_eval, tf)
    O16_O16__n_S31(rate_eval, tf)
    O16_O16__p_P31(rate_eval, tf)
    O16_O16__He4_Si28(rate_eval, tf)
    n_O17__He4_C14(rate_eval, tf)
    p_O17__n_F17(rate_eval, tf)
    p_O17__He4_N14(rate_eval, tf)
    He4_O17__n_Ne20(rate_eval, tf)
    p_O18__n_F18(rate_eval, tf)
    p_O18__He4_N15(rate_eval, tf)
    He4_O18__n_Ne21(rate_eval, tf)
    n_F15__p_O15(rate_eval, tf)
    He4_F15__p_Ne18(rate_eval, tf)
    n_F16__p_O16(rate_eval, tf)
    n_F16__He4_N13(rate_eval, tf)
    He4_F16__p_Ne19(rate_eval, tf)
    n_F17__p_O17(rate_eval, tf)
    n_F17__He4_N14(rate_eval, tf)
    p_F17__He4_O14(rate_eval, tf)
    He4_F17__p_Ne20(rate_eval, tf)
    n_F18__p_O18(rate_eval, tf)
    n_F18__He4_N15(rate_eval, tf)
    p_F18__n_Ne18(rate_eval, tf)
    p_F18__He4_O15(rate_eval, tf)
    He4_F18__n_Na21(rate_eval, tf)
    He4_F18__p_Ne21(rate_eval, tf)
    p_F19__n_Ne19(rate_eval, tf)
    p_F19__He4_O16(rate_eval, tf)
    He4_F19__n_Na22(rate_eval, tf)
    n_Ne18__p_F18(rate_eval, tf)
    n_Ne18__He4_O15(rate_eval, tf)
    p_Ne18__He4_F15(rate_eval, tf)
    He4_Ne18__p_Na21(rate_eval, tf)
    n_Ne19__p_F19(rate_eval, tf)
    n_Ne19__He4_O16(rate_eval, tf)
    p_Ne19__He4_F16(rate_eval, tf)
    He4_Ne19__n_Mg22(rate_eval, tf)
    He4_Ne19__p_Na22(rate_eval, tf)
    n_Ne20__He4_O17(rate_eval, tf)
    p_Ne20__He4_F17(rate_eval, tf)
    He4_Ne20__n_Mg23(rate_eval, tf)
    He4_Ne20__p_Na23(rate_eval, tf)
    He4_Ne20__C12_C12(rate_eval, tf)
    C12_Ne20__n_S31(rate_eval, tf)
    C12_Ne20__p_P31(rate_eval, tf)
    C12_Ne20__He4_Si28(rate_eval, tf)
    n_Ne21__He4_O18(rate_eval, tf)
    p_Ne21__n_Na21(rate_eval, tf)
    p_Ne21__He4_F18(rate_eval, tf)
    He4_Ne21__n_Mg24(rate_eval, tf)
    He4_Ne21__p_Na24(rate_eval, tf)
    n_Na21__p_Ne21(rate_eval, tf)
    n_Na21__He4_F18(rate_eval, tf)
    p_Na21__He4_Ne18(rate_eval, tf)
    He4_Na21__p_Mg24(rate_eval, tf)
    n_Na22__He4_F19(rate_eval, tf)
    p_Na22__n_Mg22(rate_eval, tf)
    p_Na22__He4_Ne19(rate_eval, tf)
    He4_Na22__n_Al25(rate_eval, tf)
    He4_Na22__p_Mg25(rate_eval, tf)
    p_Na23__n_Mg23(rate_eval, tf)
    p_Na23__He4_Ne20(rate_eval, tf)
    p_Na23__C12_C12(rate_eval, tf)
    He4_Na23__n_Al26(rate_eval, tf)
    He4_Na23__p_Mg26(rate_eval, tf)
    p_Na24__n_Mg24(rate_eval, tf)
    p_Na24__He4_Ne21(rate_eval, tf)
    He4_Na24__n_Al27(rate_eval, tf)
    n_Mg22__p_Na22(rate_eval, tf)
    n_Mg22__He4_Ne19(rate_eval, tf)
    He4_Mg22__p_Al25(rate_eval, tf)
    n_Mg23__p_Na23(rate_eval, tf)
    n_Mg23__He4_Ne20(rate_eval, tf)
    n_Mg23__C12_C12(rate_eval, tf)
    He4_Mg23__n_Si26(rate_eval, tf)
    He4_Mg23__p_Al26(rate_eval, tf)
    n_Mg24__p_Na24(rate_eval, tf)
    n_Mg24__He4_Ne21(rate_eval, tf)
    p_Mg24__He4_Na21(rate_eval, tf)
    He4_Mg24__p_Al27(rate_eval, tf)
    He4_Mg24__C12_O16(rate_eval, tf)
    p_Mg25__n_Al25(rate_eval, tf)
    p_Mg25__He4_Na22(rate_eval, tf)
    He4_Mg25__n_Si28(rate_eval, tf)
    He4_Mg25__p_Al28(rate_eval, tf)
    p_Mg26__n_Al26(rate_eval, tf)
    p_Mg26__He4_Na23(rate_eval, tf)
    He4_Mg26__n_Si29(rate_eval, tf)
    n_Al25__p_Mg25(rate_eval, tf)
    n_Al25__He4_Na22(rate_eval, tf)
    p_Al25__He4_Mg22(rate_eval, tf)
    He4_Al25__p_Si28(rate_eval, tf)
    n_Al26__p_Mg26(rate_eval, tf)
    n_Al26__He4_Na23(rate_eval, tf)
    p_Al26__n_Si26(rate_eval, tf)
    p_Al26__He4_Mg23(rate_eval, tf)
    He4_Al26__n_P29(rate_eval, tf)
    He4_Al26__p_Si29(rate_eval, tf)
    n_Al27__He4_Na24(rate_eval, tf)
    p_Al27__He4_Mg24(rate_eval, tf)
    p_Al27__C12_O16(rate_eval, tf)
    He4_Al27__n_P30(rate_eval, tf)
    He4_Al27__p_Si30(rate_eval, tf)
    p_Al28__n_Si28(rate_eval, tf)
    p_Al28__He4_Mg25(rate_eval, tf)
    He4_Al28__n_P31(rate_eval, tf)
    n_Si26__p_Al26(rate_eval, tf)
    n_Si26__He4_Mg23(rate_eval, tf)
    He4_Si26__p_P29(rate_eval, tf)
    n_Si28__p_Al28(rate_eval, tf)
    n_Si28__He4_Mg25(rate_eval, tf)
    p_Si28__He4_Al25(rate_eval, tf)
    He4_Si28__n_S31(rate_eval, tf)
    He4_Si28__p_P31(rate_eval, tf)
    He4_Si28__C12_Ne20(rate_eval, tf)
    He4_Si28__O16_O16(rate_eval, tf)
    n_Si29__He4_Mg26(rate_eval, tf)
    p_Si29__n_P29(rate_eval, tf)
    p_Si29__He4_Al26(rate_eval, tf)
    He4_Si29__n_S32(rate_eval, tf)
    He4_Si29__p_P32(rate_eval, tf)
    p_Si30__n_P30(rate_eval, tf)
    p_Si30__He4_Al27(rate_eval, tf)
    n_P29__p_Si29(rate_eval, tf)
    n_P29__He4_Al26(rate_eval, tf)
    p_P29__He4_Si26(rate_eval, tf)
    He4_P29__n_Cl32(rate_eval, tf)
    He4_P29__p_S32(rate_eval, tf)
    n_P30__p_Si30(rate_eval, tf)
    n_P30__He4_Al27(rate_eval, tf)
    n_P31__He4_Al28(rate_eval, tf)
    p_P31__n_S31(rate_eval, tf)
    p_P31__He4_Si28(rate_eval, tf)
    p_P31__C12_Ne20(rate_eval, tf)
    p_P31__O16_O16(rate_eval, tf)
    He4_P31__n_Cl34(rate_eval, tf)
    He4_P31__p_S34(rate_eval, tf)
    p_P32__n_S32(rate_eval, tf)
    p_P32__He4_Si29(rate_eval, tf)
    He4_P32__n_Cl35(rate_eval, tf)
    n_S31__p_P31(rate_eval, tf)
    n_S31__He4_Si28(rate_eval, tf)
    n_S31__C12_Ne20(rate_eval, tf)
    n_S31__O16_O16(rate_eval, tf)
    He4_S31__p_Cl34(rate_eval, tf)
    n_S32__p_P32(rate_eval, tf)
    n_S32__He4_Si29(rate_eval, tf)
    p_S32__n_Cl32(rate_eval, tf)
    p_S32__He4_P29(rate_eval, tf)
    He4_S32__n_Ar35(rate_eval, tf)
    He4_S32__p_Cl35(rate_eval, tf)
    p_S34__n_Cl34(rate_eval, tf)
    p_S34__He4_P31(rate_eval, tf)
    n_Cl32__p_S32(rate_eval, tf)
    n_Cl32__He4_P29(rate_eval, tf)
    He4_Cl32__p_Ar35(rate_eval, tf)
    n_Cl34__p_S34(rate_eval, tf)
    n_Cl34__He4_P31(rate_eval, tf)
    p_Cl34__He4_S31(rate_eval, tf)
    n_Cl35__He4_P32(rate_eval, tf)
    p_Cl35__n_Ar35(rate_eval, tf)
    p_Cl35__He4_S32(rate_eval, tf)
    n_Ar35__p_Cl35(rate_eval, tf)
    n_Ar35__He4_S32(rate_eval, tf)
    p_Ar35__He4_Cl32(rate_eval, tf)
    p_d__n_p_p(rate_eval, tf)
    t_t__n_n_He4(rate_eval, tf)
    t_He3__n_p_He4(rate_eval, tf)
    He3_He3__p_p_He4(rate_eval, tf)
    d_Li7__n_He4_He4(rate_eval, tf)
    d_Be7__p_He4_He4(rate_eval, tf)
    p_Be9__d_He4_He4(rate_eval, tf)
    n_B8__p_He4_He4(rate_eval, tf)
    n_C11__He4_He4_He4(rate_eval, tf)
    t_Li7__n_n_He4_He4(rate_eval, tf)
    He3_Li7__n_p_He4_He4(rate_eval, tf)
    t_Be7__n_p_He4_He4(rate_eval, tf)
    He3_Be7__p_p_He4_He4(rate_eval, tf)
    p_Be9__n_p_He4_He4(rate_eval, tf)
    n_n_He4__He6(rate_eval, tf)
    n_p_He4__Li6(rate_eval, tf)
    n_He4_He4__Be9(rate_eval, tf)
    He4_He4_He4__C12(rate_eval, tf)
    n_p_p__p_d(rate_eval, tf)
    n_n_He4__t_t(rate_eval, tf)
    n_p_He4__t_He3(rate_eval, tf)
    p_p_He4__He3_He3(rate_eval, tf)
    n_He4_He4__d_Li7(rate_eval, tf)
    p_He4_He4__n_B8(rate_eval, tf)
    p_He4_He4__d_Be7(rate_eval, tf)
    d_He4_He4__p_Be9(rate_eval, tf)
    He4_He4_He4__n_C11(rate_eval, tf)
    n_n_He4_He4__t_Li7(rate_eval, tf)
    n_p_He4_He4__He3_Li7(rate_eval, tf)
    n_p_He4_He4__t_Be7(rate_eval, tf)
    n_p_He4_He4__p_Be9(rate_eval, tf)
    p_p_He4_He4__He3_Be7(rate_eval, tf)

    if screen_func is not None:
        plasma_state = PlasmaState(T, rho, Y, Z)

        scn_fac = ScreenFactors(1, 1, 1, 1)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_p__d__weak__bet_pos_ *= scor
        rate_eval.p_p__d__weak__electron_capture *= scor
        rate_eval.n_p_p__p_d *= scor
        rate_eval.p_p_He4_He4__He3_Be7 *= scor

        scn_fac = ScreenFactors(1, 1, 1, 2)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_d__He3 *= scor
        rate_eval.p_d__n_p_p *= scor

        scn_fac = ScreenFactors(1, 2, 1, 2)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.d_d__He4 *= scor
        rate_eval.d_d__n_He3 *= scor
        rate_eval.d_d__p_t *= scor

        scn_fac = ScreenFactors(1, 2, 2, 4)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_d__Li6 *= scor
        rate_eval.d_He4__t_He3 *= scor

        scn_fac = ScreenFactors(1, 1, 1, 3)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_t__He4 *= scor
        rate_eval.p_t__n_He3 *= scor
        rate_eval.p_t__d_d *= scor

        scn_fac = ScreenFactors(1, 3, 2, 4)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_t__Li7 *= scor
        rate_eval.He4_t__n_Li6 *= scor

        scn_fac = ScreenFactors(1, 1, 2, 3)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_He3__He4__weak__bet_pos_ *= scor

        scn_fac = ScreenFactors(2, 4, 2, 3)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_He3__Be7 *= scor
        rate_eval.He4_He3__p_Li6 *= scor

        scn_fac = ScreenFactors(1, 1, 3, 6)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Li6__Be7 *= scor
        rate_eval.p_Li6__He4_He3 *= scor

        scn_fac = ScreenFactors(1, 1, 4, 7)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Be7__B8 *= scor

        scn_fac = ScreenFactors(2, 4, 4, 7)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Be7__C11 *= scor

        scn_fac = ScreenFactors(1, 1, 6, 12)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_C12__N13 *= scor

        scn_fac = ScreenFactors(2, 4, 6, 12)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_C12__O16 *= scor
        rate_eval.He4_C12__n_O15 *= scor
        rate_eval.He4_C12__p_N15 *= scor

        scn_fac = ScreenFactors(1, 1, 6, 13)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_C13__N14 *= scor
        rate_eval.p_C13__n_N13 *= scor

        scn_fac = ScreenFactors(1, 1, 6, 14)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_C14__N15 *= scor
        rate_eval.p_C14__n_N14 *= scor

        scn_fac = ScreenFactors(2, 4, 6, 14)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_C14__O18 *= scor
        rate_eval.He4_C14__n_O17 *= scor

        scn_fac = ScreenFactors(1, 1, 7, 13)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_N13__O14 *= scor

        scn_fac = ScreenFactors(1, 1, 7, 14)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_N14__O15 *= scor
        rate_eval.p_N14__n_O14 *= scor
        rate_eval.p_N14__He4_C11 *= scor

        scn_fac = ScreenFactors(2, 4, 7, 14)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_N14__F18 *= scor
        rate_eval.He4_N14__n_F17 *= scor
        rate_eval.He4_N14__p_O17 *= scor

        scn_fac = ScreenFactors(1, 1, 7, 15)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_N15__O16 *= scor
        rate_eval.p_N15__n_O15 *= scor
        rate_eval.p_N15__He4_C12 *= scor

        scn_fac = ScreenFactors(2, 4, 7, 15)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_N15__F19 *= scor
        rate_eval.He4_N15__n_F18 *= scor
        rate_eval.He4_N15__p_O18 *= scor

        scn_fac = ScreenFactors(1, 1, 8, 14)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_O14__F15 *= scor

        scn_fac = ScreenFactors(2, 4, 8, 14)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_O14__Ne18 *= scor
        rate_eval.He4_O14__p_F17 *= scor

        scn_fac = ScreenFactors(1, 1, 8, 15)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_O15__F16 *= scor
        rate_eval.p_O15__n_F15 *= scor

        scn_fac = ScreenFactors(2, 4, 8, 15)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_O15__Ne19 *= scor
        rate_eval.He4_O15__n_Ne18 *= scor
        rate_eval.He4_O15__p_F18 *= scor

        scn_fac = ScreenFactors(1, 1, 8, 16)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_O16__F17 *= scor
        rate_eval.p_O16__n_F16 *= scor
        rate_eval.p_O16__He4_N13 *= scor

        scn_fac = ScreenFactors(2, 4, 8, 16)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_O16__Ne20 *= scor
        rate_eval.He4_O16__n_Ne19 *= scor
        rate_eval.He4_O16__p_F19 *= scor

        scn_fac = ScreenFactors(1, 1, 8, 17)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_O17__F18 *= scor
        rate_eval.p_O17__n_F17 *= scor
        rate_eval.p_O17__He4_N14 *= scor

        scn_fac = ScreenFactors(2, 4, 8, 17)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_O17__Ne21 *= scor
        rate_eval.He4_O17__n_Ne20 *= scor

        scn_fac = ScreenFactors(1, 1, 8, 18)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_O18__F19 *= scor
        rate_eval.p_O18__n_F18 *= scor
        rate_eval.p_O18__He4_N15 *= scor

        scn_fac = ScreenFactors(1, 1, 9, 17)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_F17__Ne18 *= scor
        rate_eval.p_F17__He4_O14 *= scor

        scn_fac = ScreenFactors(2, 4, 9, 17)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_F17__Na21 *= scor
        rate_eval.He4_F17__p_Ne20 *= scor

        scn_fac = ScreenFactors(1, 1, 9, 18)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_F18__Ne19 *= scor
        rate_eval.p_F18__n_Ne18 *= scor
        rate_eval.p_F18__He4_O15 *= scor

        scn_fac = ScreenFactors(2, 4, 9, 18)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_F18__Na22 *= scor
        rate_eval.He4_F18__n_Na21 *= scor
        rate_eval.He4_F18__p_Ne21 *= scor

        scn_fac = ScreenFactors(1, 1, 9, 19)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_F19__Ne20 *= scor
        rate_eval.p_F19__n_Ne19 *= scor
        rate_eval.p_F19__He4_O16 *= scor

        scn_fac = ScreenFactors(2, 4, 9, 19)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_F19__Na23 *= scor
        rate_eval.He4_F19__n_Na22 *= scor

        scn_fac = ScreenFactors(2, 4, 10, 18)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Ne18__Mg22 *= scor
        rate_eval.He4_Ne18__p_Na21 *= scor

        scn_fac = ScreenFactors(2, 4, 10, 19)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Ne19__Mg23 *= scor
        rate_eval.He4_Ne19__n_Mg22 *= scor
        rate_eval.He4_Ne19__p_Na22 *= scor

        scn_fac = ScreenFactors(1, 1, 10, 20)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Ne20__Na21 *= scor
        rate_eval.p_Ne20__He4_F17 *= scor

        scn_fac = ScreenFactors(2, 4, 10, 20)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Ne20__Mg24 *= scor
        rate_eval.He4_Ne20__n_Mg23 *= scor
        rate_eval.He4_Ne20__p_Na23 *= scor
        rate_eval.He4_Ne20__C12_C12 *= scor

        scn_fac = ScreenFactors(1, 1, 10, 21)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Ne21__Na22 *= scor
        rate_eval.p_Ne21__n_Na21 *= scor
        rate_eval.p_Ne21__He4_F18 *= scor

        scn_fac = ScreenFactors(2, 4, 10, 21)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Ne21__Mg25 *= scor
        rate_eval.He4_Ne21__n_Mg24 *= scor
        rate_eval.He4_Ne21__p_Na24 *= scor

        scn_fac = ScreenFactors(1, 1, 11, 21)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Na21__Mg22 *= scor
        rate_eval.p_Na21__He4_Ne18 *= scor

        scn_fac = ScreenFactors(2, 4, 11, 21)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Na21__Al25 *= scor
        rate_eval.He4_Na21__p_Mg24 *= scor

        scn_fac = ScreenFactors(1, 1, 11, 22)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Na22__Mg23 *= scor
        rate_eval.p_Na22__n_Mg22 *= scor
        rate_eval.p_Na22__He4_Ne19 *= scor

        scn_fac = ScreenFactors(2, 4, 11, 22)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Na22__Al26 *= scor
        rate_eval.He4_Na22__n_Al25 *= scor
        rate_eval.He4_Na22__p_Mg25 *= scor

        scn_fac = ScreenFactors(1, 1, 11, 23)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Na23__Mg24 *= scor
        rate_eval.p_Na23__n_Mg23 *= scor
        rate_eval.p_Na23__He4_Ne20 *= scor
        rate_eval.p_Na23__C12_C12 *= scor

        scn_fac = ScreenFactors(2, 4, 11, 23)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Na23__Al27 *= scor
        rate_eval.He4_Na23__n_Al26 *= scor
        rate_eval.He4_Na23__p_Mg26 *= scor

        scn_fac = ScreenFactors(1, 1, 11, 24)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Na24__Mg25 *= scor
        rate_eval.p_Na24__n_Mg24 *= scor
        rate_eval.p_Na24__He4_Ne21 *= scor

        scn_fac = ScreenFactors(2, 4, 11, 24)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Na24__Al28 *= scor
        rate_eval.He4_Na24__n_Al27 *= scor

        scn_fac = ScreenFactors(2, 4, 12, 22)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Mg22__Si26 *= scor
        rate_eval.He4_Mg22__p_Al25 *= scor

        scn_fac = ScreenFactors(1, 1, 12, 24)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Mg24__Al25 *= scor
        rate_eval.p_Mg24__He4_Na21 *= scor

        scn_fac = ScreenFactors(2, 4, 12, 24)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Mg24__Si28 *= scor
        rate_eval.He4_Mg24__p_Al27 *= scor
        rate_eval.He4_Mg24__C12_O16 *= scor

        scn_fac = ScreenFactors(1, 1, 12, 25)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Mg25__Al26 *= scor
        rate_eval.p_Mg25__n_Al25 *= scor
        rate_eval.p_Mg25__He4_Na22 *= scor

        scn_fac = ScreenFactors(2, 4, 12, 25)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Mg25__Si29 *= scor
        rate_eval.He4_Mg25__n_Si28 *= scor
        rate_eval.He4_Mg25__p_Al28 *= scor

        scn_fac = ScreenFactors(1, 1, 12, 26)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Mg26__Al27 *= scor
        rate_eval.p_Mg26__n_Al26 *= scor
        rate_eval.p_Mg26__He4_Na23 *= scor

        scn_fac = ScreenFactors(2, 4, 12, 26)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Mg26__Si30 *= scor
        rate_eval.He4_Mg26__n_Si29 *= scor

        scn_fac = ScreenFactors(1, 1, 13, 25)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Al25__Si26 *= scor
        rate_eval.p_Al25__He4_Mg22 *= scor

        scn_fac = ScreenFactors(2, 4, 13, 25)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Al25__P29 *= scor
        rate_eval.He4_Al25__p_Si28 *= scor

        scn_fac = ScreenFactors(2, 4, 13, 26)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Al26__P30 *= scor
        rate_eval.He4_Al26__n_P29 *= scor
        rate_eval.He4_Al26__p_Si29 *= scor

        scn_fac = ScreenFactors(1, 1, 13, 27)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Al27__Si28 *= scor
        rate_eval.p_Al27__He4_Mg24 *= scor
        rate_eval.p_Al27__C12_O16 *= scor

        scn_fac = ScreenFactors(2, 4, 13, 27)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Al27__P31 *= scor
        rate_eval.He4_Al27__n_P30 *= scor
        rate_eval.He4_Al27__p_Si30 *= scor

        scn_fac = ScreenFactors(1, 1, 13, 28)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Al28__Si29 *= scor
        rate_eval.p_Al28__n_Si28 *= scor
        rate_eval.p_Al28__He4_Mg25 *= scor

        scn_fac = ScreenFactors(2, 4, 13, 28)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Al28__P32 *= scor
        rate_eval.He4_Al28__n_P31 *= scor

        scn_fac = ScreenFactors(1, 1, 14, 28)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Si28__P29 *= scor
        rate_eval.p_Si28__He4_Al25 *= scor

        scn_fac = ScreenFactors(2, 4, 14, 28)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Si28__S32 *= scor
        rate_eval.He4_Si28__n_S31 *= scor
        rate_eval.He4_Si28__p_P31 *= scor
        rate_eval.He4_Si28__C12_Ne20 *= scor
        rate_eval.He4_Si28__O16_O16 *= scor

        scn_fac = ScreenFactors(1, 1, 14, 29)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Si29__P30 *= scor
        rate_eval.p_Si29__n_P29 *= scor
        rate_eval.p_Si29__He4_Al26 *= scor

        scn_fac = ScreenFactors(1, 1, 14, 30)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Si30__P31 *= scor
        rate_eval.p_Si30__n_P30 *= scor
        rate_eval.p_Si30__He4_Al27 *= scor

        scn_fac = ScreenFactors(2, 4, 14, 30)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Si30__S34 *= scor

        scn_fac = ScreenFactors(1, 1, 15, 30)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_P30__S31 *= scor

        scn_fac = ScreenFactors(2, 4, 15, 30)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_P30__Cl34 *= scor

        scn_fac = ScreenFactors(1, 1, 15, 31)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_P31__S32 *= scor
        rate_eval.p_P31__n_S31 *= scor
        rate_eval.p_P31__He4_Si28 *= scor
        rate_eval.p_P31__C12_Ne20 *= scor
        rate_eval.p_P31__O16_O16 *= scor

        scn_fac = ScreenFactors(2, 4, 15, 31)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_P31__Cl35 *= scor
        rate_eval.He4_P31__n_Cl34 *= scor
        rate_eval.He4_P31__p_S34 *= scor

        scn_fac = ScreenFactors(1, 1, 16, 31)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_S31__Cl32 *= scor

        scn_fac = ScreenFactors(2, 4, 16, 31)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_S31__Ar35 *= scor
        rate_eval.He4_S31__p_Cl34 *= scor

        scn_fac = ScreenFactors(1, 1, 16, 34)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_S34__Cl35 *= scor
        rate_eval.p_S34__n_Cl34 *= scor
        rate_eval.p_S34__He4_P31 *= scor

        scn_fac = ScreenFactors(1, 1, 17, 34)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Cl34__Ar35 *= scor
        rate_eval.p_Cl34__He4_S31 *= scor

        scn_fac = ScreenFactors(1, 2, 1, 3)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.d_t__n_He4 *= scor

        scn_fac = ScreenFactors(1, 2, 2, 3)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.d_He3__p_He4 *= scor

        scn_fac = ScreenFactors(1, 3, 2, 3)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.t_He3__d_He4 *= scor
        rate_eval.t_He3__n_p_He4 *= scor

        scn_fac = ScreenFactors(1, 1, 2, 4)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_He4__d_He3 *= scor
        rate_eval.n_p_He4__Li6 *= scor
        rate_eval.n_p_He4__t_He3 *= scor

        scn_fac = ScreenFactors(2, 4, 2, 4)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_He4__n_Be7 *= scor
        rate_eval.He4_He4__p_Li7 *= scor
        rate_eval.n_He4_He4__Be9 *= scor
        rate_eval.n_He4_He4__d_Li7 *= scor
        rate_eval.n_n_He4_He4__t_Li7 *= scor

        scn_fac = ScreenFactors(1, 2, 3, 6)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.d_Li6__n_Be7 *= scor
        rate_eval.d_Li6__p_Li7 *= scor

        scn_fac = ScreenFactors(2, 4, 3, 6)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Li6__p_Be9 *= scor

        scn_fac = ScreenFactors(1, 1, 3, 7)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Li7__n_Be7 *= scor
        rate_eval.p_Li7__d_Li6 *= scor
        rate_eval.p_Li7__He4_He4 *= scor

        scn_fac = ScreenFactors(1, 3, 3, 7)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.t_Li7__n_Be9 *= scor
        rate_eval.t_Li7__n_n_He4_He4 *= scor

        scn_fac = ScreenFactors(1, 1, 4, 9)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Be9__He4_Li6 *= scor
        rate_eval.p_Be9__d_He4_He4 *= scor
        rate_eval.p_Be9__n_p_He4_He4 *= scor

        scn_fac = ScreenFactors(2, 4, 4, 9)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Be9__n_C12 *= scor

        scn_fac = ScreenFactors(2, 4, 5, 8)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_B8__p_C11 *= scor

        scn_fac = ScreenFactors(1, 1, 6, 11)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_C11__He4_B8 *= scor

        scn_fac = ScreenFactors(2, 4, 6, 11)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_C11__n_O14 *= scor
        rate_eval.He4_C11__p_N14 *= scor

        scn_fac = ScreenFactors(6, 12, 6, 12)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.C12_C12__n_Mg23 *= scor
        rate_eval.C12_C12__p_Na23 *= scor
        rate_eval.C12_C12__He4_Ne20 *= scor

        scn_fac = ScreenFactors(1, 2, 6, 13)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.d_C13__n_N14 *= scor

        scn_fac = ScreenFactors(2, 4, 6, 13)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_C13__n_O16 *= scor

        scn_fac = ScreenFactors(1, 2, 6, 14)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.d_C14__n_N15 *= scor

        scn_fac = ScreenFactors(2, 4, 7, 13)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_N13__n_F16 *= scor
        rate_eval.He4_N13__p_O16 *= scor

        scn_fac = ScreenFactors(6, 12, 8, 16)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.C12_O16__p_Al27 *= scor
        rate_eval.C12_O16__He4_Mg24 *= scor

        scn_fac = ScreenFactors(8, 16, 8, 16)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.O16_O16__n_S31 *= scor
        rate_eval.O16_O16__p_P31 *= scor
        rate_eval.O16_O16__He4_Si28 *= scor

        scn_fac = ScreenFactors(2, 4, 8, 18)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_O18__n_Ne21 *= scor

        scn_fac = ScreenFactors(2, 4, 9, 15)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_F15__p_Ne18 *= scor

        scn_fac = ScreenFactors(2, 4, 9, 16)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_F16__p_Ne19 *= scor

        scn_fac = ScreenFactors(1, 1, 10, 18)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Ne18__He4_F15 *= scor

        scn_fac = ScreenFactors(1, 1, 10, 19)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Ne19__He4_F16 *= scor

        scn_fac = ScreenFactors(6, 12, 10, 20)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.C12_Ne20__n_S31 *= scor
        rate_eval.C12_Ne20__p_P31 *= scor
        rate_eval.C12_Ne20__He4_Si28 *= scor

        scn_fac = ScreenFactors(2, 4, 12, 23)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Mg23__n_Si26 *= scor
        rate_eval.He4_Mg23__p_Al26 *= scor

        scn_fac = ScreenFactors(1, 1, 13, 26)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Al26__n_Si26 *= scor
        rate_eval.p_Al26__He4_Mg23 *= scor

        scn_fac = ScreenFactors(2, 4, 14, 26)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Si26__p_P29 *= scor

        scn_fac = ScreenFactors(2, 4, 14, 29)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Si29__n_S32 *= scor
        rate_eval.He4_Si29__p_P32 *= scor

        scn_fac = ScreenFactors(1, 1, 15, 29)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_P29__He4_Si26 *= scor

        scn_fac = ScreenFactors(2, 4, 15, 29)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_P29__n_Cl32 *= scor
        rate_eval.He4_P29__p_S32 *= scor

        scn_fac = ScreenFactors(1, 1, 15, 32)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_P32__n_S32 *= scor
        rate_eval.p_P32__He4_Si29 *= scor

        scn_fac = ScreenFactors(2, 4, 15, 32)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_P32__n_Cl35 *= scor

        scn_fac = ScreenFactors(1, 1, 16, 32)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_S32__n_Cl32 *= scor
        rate_eval.p_S32__He4_P29 *= scor

        scn_fac = ScreenFactors(2, 4, 16, 32)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_S32__n_Ar35 *= scor
        rate_eval.He4_S32__p_Cl35 *= scor

        scn_fac = ScreenFactors(2, 4, 17, 32)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Cl32__p_Ar35 *= scor

        scn_fac = ScreenFactors(1, 1, 17, 35)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Cl35__n_Ar35 *= scor
        rate_eval.p_Cl35__He4_S32 *= scor

        scn_fac = ScreenFactors(1, 1, 18, 35)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Ar35__He4_Cl32 *= scor

        scn_fac = ScreenFactors(1, 3, 1, 3)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.t_t__n_n_He4 *= scor

        scn_fac = ScreenFactors(2, 3, 2, 3)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He3_He3__p_p_He4 *= scor

        scn_fac = ScreenFactors(1, 2, 3, 7)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.d_Li7__n_He4_He4 *= scor

        scn_fac = ScreenFactors(1, 2, 4, 7)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.d_Be7__p_He4_He4 *= scor

        scn_fac = ScreenFactors(2, 3, 3, 7)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He3_Li7__n_p_He4_He4 *= scor

        scn_fac = ScreenFactors(1, 3, 4, 7)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.t_Be7__n_p_He4_He4 *= scor

        scn_fac = ScreenFactors(2, 3, 4, 7)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He3_Be7__p_p_He4_He4 *= scor

        scn_fac = ScreenFactors(2, 4, 2, 4)
        scor = screen_func(plasma_state, scn_fac)
        scn_fac2 = ScreenFactors(2, 4, 4, 8)
        scor2 = screen_func(plasma_state, scn_fac2)
        rate_eval.He4_He4_He4__C12 *= scor * scor2
        rate_eval.He4_He4_He4__n_C11 *= scor * scor2

        scn_fac = ScreenFactors(1, 1, 1, 1)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_p_He4__He3_He3 *= scor

        scn_fac = ScreenFactors(1, 1, 2, 4)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_He4_He4__n_B8 *= scor
        rate_eval.p_He4_He4__d_Be7 *= scor
        rate_eval.n_p_He4_He4__He3_Li7 *= scor
        rate_eval.n_p_He4_He4__t_Be7 *= scor
        rate_eval.n_p_He4_He4__p_Be9 *= scor

        scn_fac = ScreenFactors(1, 2, 2, 4)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.d_He4_He4__p_Be9 *= scor

    jac = np.zeros((nnuc, nnuc), dtype=np.float64)

    jac[jn, jn] = (
       -rate_eval.n__p__weak__wc12
       -rho*Y[jp]*rate_eval.n_p__d
       -rho*Y[jd]*rate_eval.n_d__t
       -rho*Y[jhe3]*rate_eval.n_He3__He4
       -rho*Y[jli6]*rate_eval.n_Li6__Li7
       -rho*Y[jc11]*rate_eval.n_C11__C12
       -rho*Y[jc12]*rate_eval.n_C12__C13
       -rho*Y[jc13]*rate_eval.n_C13__C14
       -rho*Y[jn13]*rate_eval.n_N13__N14
       -rho*Y[jn14]*rate_eval.n_N14__N15
       -rho*Y[jo14]*rate_eval.n_O14__O15
       -rho*Y[jo15]*rate_eval.n_O15__O16
       -rho*Y[jo16]*rate_eval.n_O16__O17
       -rho*Y[jo17]*rate_eval.n_O17__O18
       -rho*Y[jf15]*rate_eval.n_F15__F16
       -rho*Y[jf16]*rate_eval.n_F16__F17
       -rho*Y[jf17]*rate_eval.n_F17__F18
       -rho*Y[jf18]*rate_eval.n_F18__F19
       -rho*Y[jne18]*rate_eval.n_Ne18__Ne19
       -rho*Y[jne19]*rate_eval.n_Ne19__Ne20
       -rho*Y[jne20]*rate_eval.n_Ne20__Ne21
       -rho*Y[jna21]*rate_eval.n_Na21__Na22
       -rho*Y[jna22]*rate_eval.n_Na22__Na23
       -rho*Y[jna23]*rate_eval.n_Na23__Na24
       -rho*Y[jmg22]*rate_eval.n_Mg22__Mg23
       -rho*Y[jmg23]*rate_eval.n_Mg23__Mg24
       -rho*Y[jmg24]*rate_eval.n_Mg24__Mg25
       -rho*Y[jmg25]*rate_eval.n_Mg25__Mg26
       -rho*Y[jal25]*rate_eval.n_Al25__Al26
       -rho*Y[jal26]*rate_eval.n_Al26__Al27
       -rho*Y[jal27]*rate_eval.n_Al27__Al28
       -rho*Y[jsi28]*rate_eval.n_Si28__Si29
       -rho*Y[jsi29]*rate_eval.n_Si29__Si30
       -rho*Y[jp29]*rate_eval.n_P29__P30
       -rho*Y[jp30]*rate_eval.n_P30__P31
       -rho*Y[jp31]*rate_eval.n_P31__P32
       -rho*Y[js31]*rate_eval.n_S31__S32
       -rho*Y[jcl34]*rate_eval.n_Cl34__Cl35
       -rho*Y[jhe3]*rate_eval.n_He3__p_t
       -rho*Y[jhe3]*rate_eval.n_He3__d_d
       -rho*Y[jhe4]*rate_eval.n_He4__d_t
       -rho*Y[jli6]*rate_eval.n_Li6__He4_t
       -rho*Y[jbe7]*rate_eval.n_Be7__p_Li7
       -rho*Y[jbe7]*rate_eval.n_Be7__d_Li6
       -rho*Y[jbe7]*rate_eval.n_Be7__He4_He4
       -rho*Y[jbe9]*rate_eval.n_Be9__t_Li7
       -rho*Y[jc12]*rate_eval.n_C12__He4_Be9
       -rho*Y[jn13]*rate_eval.n_N13__p_C13
       -rho*Y[jn14]*rate_eval.n_N14__p_C14
       -rho*Y[jn14]*rate_eval.n_N14__d_C13
       -rho*Y[jn15]*rate_eval.n_N15__d_C14
       -rho*Y[jo14]*rate_eval.n_O14__p_N14
       -rho*Y[jo14]*rate_eval.n_O14__He4_C11
       -rho*Y[jo15]*rate_eval.n_O15__p_N15
       -rho*Y[jo15]*rate_eval.n_O15__He4_C12
       -rho*Y[jo16]*rate_eval.n_O16__He4_C13
       -rho*Y[jo17]*rate_eval.n_O17__He4_C14
       -rho*Y[jf15]*rate_eval.n_F15__p_O15
       -rho*Y[jf16]*rate_eval.n_F16__p_O16
       -rho*Y[jf16]*rate_eval.n_F16__He4_N13
       -rho*Y[jf17]*rate_eval.n_F17__p_O17
       -rho*Y[jf17]*rate_eval.n_F17__He4_N14
       -rho*Y[jf18]*rate_eval.n_F18__p_O18
       -rho*Y[jf18]*rate_eval.n_F18__He4_N15
       -rho*Y[jne18]*rate_eval.n_Ne18__p_F18
       -rho*Y[jne18]*rate_eval.n_Ne18__He4_O15
       -rho*Y[jne19]*rate_eval.n_Ne19__p_F19
       -rho*Y[jne19]*rate_eval.n_Ne19__He4_O16
       -rho*Y[jne20]*rate_eval.n_Ne20__He4_O17
       -rho*Y[jne21]*rate_eval.n_Ne21__He4_O18
       -rho*Y[jna21]*rate_eval.n_Na21__p_Ne21
       -rho*Y[jna21]*rate_eval.n_Na21__He4_F18
       -rho*Y[jna22]*rate_eval.n_Na22__He4_F19
       -rho*Y[jmg22]*rate_eval.n_Mg22__p_Na22
       -rho*Y[jmg22]*rate_eval.n_Mg22__He4_Ne19
       -rho*Y[jmg23]*rate_eval.n_Mg23__p_Na23
       -rho*Y[jmg23]*rate_eval.n_Mg23__He4_Ne20
       -rho*Y[jmg23]*rate_eval.n_Mg23__C12_C12
       -rho*Y[jmg24]*rate_eval.n_Mg24__p_Na24
       -rho*Y[jmg24]*rate_eval.n_Mg24__He4_Ne21
       -rho*Y[jal25]*rate_eval.n_Al25__p_Mg25
       -rho*Y[jal25]*rate_eval.n_Al25__He4_Na22
       -rho*Y[jal26]*rate_eval.n_Al26__p_Mg26
       -rho*Y[jal26]*rate_eval.n_Al26__He4_Na23
       -rho*Y[jal27]*rate_eval.n_Al27__He4_Na24
       -rho*Y[jsi26]*rate_eval.n_Si26__p_Al26
       -rho*Y[jsi26]*rate_eval.n_Si26__He4_Mg23
       -rho*Y[jsi28]*rate_eval.n_Si28__p_Al28
       -rho*Y[jsi28]*rate_eval.n_Si28__He4_Mg25
       -rho*Y[jsi29]*rate_eval.n_Si29__He4_Mg26
       -rho*Y[jp29]*rate_eval.n_P29__p_Si29
       -rho*Y[jp29]*rate_eval.n_P29__He4_Al26
       -rho*Y[jp30]*rate_eval.n_P30__p_Si30
       -rho*Y[jp30]*rate_eval.n_P30__He4_Al27
       -rho*Y[jp31]*rate_eval.n_P31__He4_Al28
       -rho*Y[js31]*rate_eval.n_S31__p_P31
       -rho*Y[js31]*rate_eval.n_S31__He4_Si28
       -rho*Y[js31]*rate_eval.n_S31__C12_Ne20
       -rho*Y[js31]*rate_eval.n_S31__O16_O16
       -rho*Y[js32]*rate_eval.n_S32__p_P32
       -rho*Y[js32]*rate_eval.n_S32__He4_Si29
       -rho*Y[jcl32]*rate_eval.n_Cl32__p_S32
       -rho*Y[jcl32]*rate_eval.n_Cl32__He4_P29
       -rho*Y[jcl34]*rate_eval.n_Cl34__p_S34
       -rho*Y[jcl34]*rate_eval.n_Cl34__He4_P31
       -rho*Y[jcl35]*rate_eval.n_Cl35__He4_P32
       -rho*Y[jar35]*rate_eval.n_Ar35__p_Cl35
       -rho*Y[jar35]*rate_eval.n_Ar35__He4_S32
       -rho*Y[jb8]*rate_eval.n_B8__p_He4_He4
       -rho*Y[jc11]*rate_eval.n_C11__He4_He4_He4
       -2*5.00000000000000e-01*rho**2*2*Y[jn]*Y[jhe4]*rate_eval.n_n_He4__He6
       -rho**2*Y[jp]*Y[jhe4]*rate_eval.n_p_He4__Li6
       -5.00000000000000e-01*rho**2*Y[jhe4]**2*rate_eval.n_He4_He4__Be9
       -5.00000000000000e-01*rho**2*Y[jp]**2*rate_eval.n_p_p__p_d
       -2*5.00000000000000e-01*rho**2*2*Y[jn]*Y[jhe4]*rate_eval.n_n_He4__t_t
       -rho**2*Y[jp]*Y[jhe4]*rate_eval.n_p_He4__t_He3
       -5.00000000000000e-01*rho**2*Y[jhe4]**2*rate_eval.n_He4_He4__d_Li7
       -2*2.50000000000000e-01*rho**3*2*Y[jn]*Y[jhe4]**2*rate_eval.n_n_He4_He4__t_Li7
       -5.00000000000000e-01*rho**3*Y[jp]*Y[jhe4]**2*rate_eval.n_p_He4_He4__He3_Li7
       -5.00000000000000e-01*rho**3*Y[jp]*Y[jhe4]**2*rate_eval.n_p_He4_He4__t_Be7
       -5.00000000000000e-01*rho**3*Y[jp]*Y[jhe4]**2*rate_eval.n_p_He4_He4__p_Be9
       )

    jac[jn, jp] = (
       -rho*Y[jn]*rate_eval.n_p__d
       -rho**2*Y[jn]*Y[jhe4]*rate_eval.n_p_He4__Li6
       -5.00000000000000e-01*rho**2*Y[jn]*2*Y[jp]*rate_eval.n_p_p__p_d
       -rho**2*Y[jn]*Y[jhe4]*rate_eval.n_p_He4__t_He3
       -5.00000000000000e-01*rho**3*Y[jn]*Y[jhe4]**2*rate_eval.n_p_He4_He4__He3_Li7
       -5.00000000000000e-01*rho**3*Y[jn]*Y[jhe4]**2*rate_eval.n_p_He4_He4__t_Be7
       -5.00000000000000e-01*rho**3*Y[jn]*Y[jhe4]**2*rate_eval.n_p_He4_He4__p_Be9
       +rho*Y[jt]*rate_eval.p_t__n_He3
       +rho*Y[jli7]*rate_eval.p_Li7__n_Be7
       +rho*Y[jc13]*rate_eval.p_C13__n_N13
       +rho*Y[jc14]*rate_eval.p_C14__n_N14
       +rho*Y[jn14]*rate_eval.p_N14__n_O14
       +rho*Y[jn15]*rate_eval.p_N15__n_O15
       +rho*Y[jo15]*rate_eval.p_O15__n_F15
       +rho*Y[jo16]*rate_eval.p_O16__n_F16
       +rho*Y[jo17]*rate_eval.p_O17__n_F17
       +rho*Y[jo18]*rate_eval.p_O18__n_F18
       +rho*Y[jf18]*rate_eval.p_F18__n_Ne18
       +rho*Y[jf19]*rate_eval.p_F19__n_Ne19
       +rho*Y[jne21]*rate_eval.p_Ne21__n_Na21
       +rho*Y[jna22]*rate_eval.p_Na22__n_Mg22
       +rho*Y[jna23]*rate_eval.p_Na23__n_Mg23
       +rho*Y[jna24]*rate_eval.p_Na24__n_Mg24
       +rho*Y[jmg25]*rate_eval.p_Mg25__n_Al25
       +rho*Y[jmg26]*rate_eval.p_Mg26__n_Al26
       +rho*Y[jal26]*rate_eval.p_Al26__n_Si26
       +rho*Y[jal28]*rate_eval.p_Al28__n_Si28
       +rho*Y[jsi29]*rate_eval.p_Si29__n_P29
       +rho*Y[jsi30]*rate_eval.p_Si30__n_P30
       +rho*Y[jp31]*rate_eval.p_P31__n_S31
       +rho*Y[jp32]*rate_eval.p_P32__n_S32
       +rho*Y[js32]*rate_eval.p_S32__n_Cl32
       +rho*Y[js34]*rate_eval.p_S34__n_Cl34
       +rho*Y[jcl35]*rate_eval.p_Cl35__n_Ar35
       +rho*Y[jd]*rate_eval.p_d__n_p_p
       +rho*Y[jbe9]*rate_eval.p_Be9__n_p_He4_He4
       +5.00000000000000e-01*rho**2*Y[jhe4]**2*rate_eval.p_He4_He4__n_B8
       )

    jac[jn, jd] = (
       -rho*Y[jn]*rate_eval.n_d__t
       +rate_eval.d__n_p
       +5.00000000000000e-01*rho*2*Y[jd]*rate_eval.d_d__n_He3
       +rho*Y[jt]*rate_eval.d_t__n_He4
       +rho*Y[jli6]*rate_eval.d_Li6__n_Be7
       +rho*Y[jc13]*rate_eval.d_C13__n_N14
       +rho*Y[jc14]*rate_eval.d_C14__n_N15
       +rho*Y[jp]*rate_eval.p_d__n_p_p
       +rho*Y[jli7]*rate_eval.d_Li7__n_He4_He4
       )

    jac[jn, jt] = (
       +rate_eval.t__n_d
       +rho*Y[jp]*rate_eval.p_t__n_He3
       +rho*Y[jd]*rate_eval.d_t__n_He4
       +rho*Y[jhe4]*rate_eval.He4_t__n_Li6
       +rho*Y[jli7]*rate_eval.t_Li7__n_Be9
       +2*5.00000000000000e-01*rho*2*Y[jt]*rate_eval.t_t__n_n_He4
       +rho*Y[jhe3]*rate_eval.t_He3__n_p_He4
       +2*rho*Y[jli7]*rate_eval.t_Li7__n_n_He4_He4
       +rho*Y[jbe7]*rate_eval.t_Be7__n_p_He4_He4
       )

    jac[jn, jhe3] = (
       -rho*Y[jn]*rate_eval.n_He3__He4
       -rho*Y[jn]*rate_eval.n_He3__p_t
       -rho*Y[jn]*rate_eval.n_He3__d_d
       +rho*Y[jt]*rate_eval.t_He3__n_p_He4
       +rho*Y[jli7]*rate_eval.He3_Li7__n_p_He4_He4
       )

    jac[jn, jhe4] = (
       -rho*Y[jn]*rate_eval.n_He4__d_t
       -2*5.00000000000000e-01*rho**2*Y[jn]**2*rate_eval.n_n_He4__He6
       -rho**2*Y[jn]*Y[jp]*rate_eval.n_p_He4__Li6
       -5.00000000000000e-01*rho**2*Y[jn]*2*Y[jhe4]*rate_eval.n_He4_He4__Be9
       -2*5.00000000000000e-01*rho**2*Y[jn]**2*rate_eval.n_n_He4__t_t
       -rho**2*Y[jn]*Y[jp]*rate_eval.n_p_He4__t_He3
       -5.00000000000000e-01*rho**2*Y[jn]*2*Y[jhe4]*rate_eval.n_He4_He4__d_Li7
       -2*2.50000000000000e-01*rho**3*Y[jn]**2*2*Y[jhe4]*rate_eval.n_n_He4_He4__t_Li7
       -5.00000000000000e-01*rho**3*Y[jn]*Y[jp]*2*Y[jhe4]*rate_eval.n_p_He4_He4__He3_Li7
       -5.00000000000000e-01*rho**3*Y[jn]*Y[jp]*2*Y[jhe4]*rate_eval.n_p_He4_He4__t_Be7
       -5.00000000000000e-01*rho**3*Y[jn]*Y[jp]*2*Y[jhe4]*rate_eval.n_p_He4_He4__p_Be9
       +rate_eval.He4__n_He3
       +rho*Y[jt]*rate_eval.He4_t__n_Li6
       +5.00000000000000e-01*rho*2*Y[jhe4]*rate_eval.He4_He4__n_Be7
       +rho*Y[jbe9]*rate_eval.He4_Be9__n_C12
       +rho*Y[jc11]*rate_eval.He4_C11__n_O14
       +rho*Y[jc12]*rate_eval.He4_C12__n_O15
       +rho*Y[jc13]*rate_eval.He4_C13__n_O16
       +rho*Y[jc14]*rate_eval.He4_C14__n_O17
       +rho*Y[jn13]*rate_eval.He4_N13__n_F16
       +rho*Y[jn14]*rate_eval.He4_N14__n_F17
       +rho*Y[jn15]*rate_eval.He4_N15__n_F18
       +rho*Y[jo15]*rate_eval.He4_O15__n_Ne18
       +rho*Y[jo16]*rate_eval.He4_O16__n_Ne19
       +rho*Y[jo17]*rate_eval.He4_O17__n_Ne20
       +rho*Y[jo18]*rate_eval.He4_O18__n_Ne21
       +rho*Y[jf18]*rate_eval.He4_F18__n_Na21
       +rho*Y[jf19]*rate_eval.He4_F19__n_Na22
       +rho*Y[jne19]*rate_eval.He4_Ne19__n_Mg22
       +rho*Y[jne20]*rate_eval.He4_Ne20__n_Mg23
       +rho*Y[jne21]*rate_eval.He4_Ne21__n_Mg24
       +rho*Y[jna22]*rate_eval.He4_Na22__n_Al25
       +rho*Y[jna23]*rate_eval.He4_Na23__n_Al26
       +rho*Y[jna24]*rate_eval.He4_Na24__n_Al27
       +rho*Y[jmg23]*rate_eval.He4_Mg23__n_Si26
       +rho*Y[jmg25]*rate_eval.He4_Mg25__n_Si28
       +rho*Y[jmg26]*rate_eval.He4_Mg26__n_Si29
       +rho*Y[jal26]*rate_eval.He4_Al26__n_P29
       +rho*Y[jal27]*rate_eval.He4_Al27__n_P30
       +rho*Y[jal28]*rate_eval.He4_Al28__n_P31
       +rho*Y[jsi28]*rate_eval.He4_Si28__n_S31
       +rho*Y[jsi29]*rate_eval.He4_Si29__n_S32
       +rho*Y[jp29]*rate_eval.He4_P29__n_Cl32
       +rho*Y[jp31]*rate_eval.He4_P31__n_Cl34
       +rho*Y[jp32]*rate_eval.He4_P32__n_Cl35
       +rho*Y[js32]*rate_eval.He4_S32__n_Ar35
       +5.00000000000000e-01*rho**2*Y[jp]*2*Y[jhe4]*rate_eval.p_He4_He4__n_B8
       +1.66666666666667e-01*rho**2*3*Y[jhe4]**2*rate_eval.He4_He4_He4__n_C11
       )

    jac[jn, jhe6] = (
       +2*rate_eval.He6__n_n_He4
       )

    jac[jn, jli6] = (
       -rho*Y[jn]*rate_eval.n_Li6__Li7
       -rho*Y[jn]*rate_eval.n_Li6__He4_t
       +rate_eval.Li6__n_p_He4
       +rho*Y[jd]*rate_eval.d_Li6__n_Be7
       )

    jac[jn, jli7] = (
       +rate_eval.Li7__n_Li6
       +rho*Y[jp]*rate_eval.p_Li7__n_Be7
       +rho*Y[jt]*rate_eval.t_Li7__n_Be9
       +rho*Y[jd]*rate_eval.d_Li7__n_He4_He4
       +2*rho*Y[jt]*rate_eval.t_Li7__n_n_He4_He4
       +rho*Y[jhe3]*rate_eval.He3_Li7__n_p_He4_He4
       )

    jac[jn, jbe7] = (
       -rho*Y[jn]*rate_eval.n_Be7__p_Li7
       -rho*Y[jn]*rate_eval.n_Be7__d_Li6
       -rho*Y[jn]*rate_eval.n_Be7__He4_He4
       +rho*Y[jt]*rate_eval.t_Be7__n_p_He4_He4
       )

    jac[jn, jbe9] = (
       -rho*Y[jn]*rate_eval.n_Be9__t_Li7
       +rate_eval.Be9__n_He4_He4
       +rho*Y[jhe4]*rate_eval.He4_Be9__n_C12
       +rho*Y[jp]*rate_eval.p_Be9__n_p_He4_He4
       )

    jac[jn, jb8] = (
       -rho*Y[jn]*rate_eval.n_B8__p_He4_He4
       )

    jac[jn, jc11] = (
       -rho*Y[jn]*rate_eval.n_C11__C12
       -rho*Y[jn]*rate_eval.n_C11__He4_He4_He4
       +rho*Y[jhe4]*rate_eval.He4_C11__n_O14
       )

    jac[jn, jc12] = (
       -rho*Y[jn]*rate_eval.n_C12__C13
       -rho*Y[jn]*rate_eval.n_C12__He4_Be9
       +rate_eval.C12__n_C11
       +rho*Y[jhe4]*rate_eval.He4_C12__n_O15
       +5.00000000000000e-01*rho*2*Y[jc12]*rate_eval.C12_C12__n_Mg23
       +rho*Y[jne20]*rate_eval.C12_Ne20__n_S31
       )

    jac[jn, jc13] = (
       -rho*Y[jn]*rate_eval.n_C13__C14
       +rate_eval.C13__n_C12
       +rho*Y[jp]*rate_eval.p_C13__n_N13
       +rho*Y[jd]*rate_eval.d_C13__n_N14
       +rho*Y[jhe4]*rate_eval.He4_C13__n_O16
       )

    jac[jn, jc14] = (
       +rate_eval.C14__n_C13
       +rho*Y[jp]*rate_eval.p_C14__n_N14
       +rho*Y[jd]*rate_eval.d_C14__n_N15
       +rho*Y[jhe4]*rate_eval.He4_C14__n_O17
       )

    jac[jn, jn13] = (
       -rho*Y[jn]*rate_eval.n_N13__N14
       -rho*Y[jn]*rate_eval.n_N13__p_C13
       +rho*Y[jhe4]*rate_eval.He4_N13__n_F16
       )

    jac[jn, jn14] = (
       -rho*Y[jn]*rate_eval.n_N14__N15
       -rho*Y[jn]*rate_eval.n_N14__p_C14
       -rho*Y[jn]*rate_eval.n_N14__d_C13
       +rate_eval.N14__n_N13
       +rho*Y[jp]*rate_eval.p_N14__n_O14
       +rho*Y[jhe4]*rate_eval.He4_N14__n_F17
       )

    jac[jn, jn15] = (
       -rho*Y[jn]*rate_eval.n_N15__d_C14
       +rate_eval.N15__n_N14
       +rho*Y[jp]*rate_eval.p_N15__n_O15
       +rho*Y[jhe4]*rate_eval.He4_N15__n_F18
       )

    jac[jn, jo14] = (
       -rho*Y[jn]*rate_eval.n_O14__O15
       -rho*Y[jn]*rate_eval.n_O14__p_N14
       -rho*Y[jn]*rate_eval.n_O14__He4_C11
       )

    jac[jn, jo15] = (
       -rho*Y[jn]*rate_eval.n_O15__O16
       -rho*Y[jn]*rate_eval.n_O15__p_N15
       -rho*Y[jn]*rate_eval.n_O15__He4_C12
       +rate_eval.O15__n_O14
       +rho*Y[jp]*rate_eval.p_O15__n_F15
       +rho*Y[jhe4]*rate_eval.He4_O15__n_Ne18
       )

    jac[jn, jo16] = (
       -rho*Y[jn]*rate_eval.n_O16__O17
       -rho*Y[jn]*rate_eval.n_O16__He4_C13
       +rate_eval.O16__n_O15
       +rho*Y[jp]*rate_eval.p_O16__n_F16
       +rho*Y[jhe4]*rate_eval.He4_O16__n_Ne19
       +5.00000000000000e-01*rho*2*Y[jo16]*rate_eval.O16_O16__n_S31
       )

    jac[jn, jo17] = (
       -rho*Y[jn]*rate_eval.n_O17__O18
       -rho*Y[jn]*rate_eval.n_O17__He4_C14
       +rate_eval.O17__n_O16
       +rho*Y[jp]*rate_eval.p_O17__n_F17
       +rho*Y[jhe4]*rate_eval.He4_O17__n_Ne20
       )

    jac[jn, jo18] = (
       +rate_eval.O18__n_O17
       +rho*Y[jp]*rate_eval.p_O18__n_F18
       +rho*Y[jhe4]*rate_eval.He4_O18__n_Ne21
       )

    jac[jn, jf15] = (
       -rho*Y[jn]*rate_eval.n_F15__F16
       -rho*Y[jn]*rate_eval.n_F15__p_O15
       )

    jac[jn, jf16] = (
       -rho*Y[jn]*rate_eval.n_F16__F17
       -rho*Y[jn]*rate_eval.n_F16__p_O16
       -rho*Y[jn]*rate_eval.n_F16__He4_N13
       +rate_eval.F16__n_F15
       )

    jac[jn, jf17] = (
       -rho*Y[jn]*rate_eval.n_F17__F18
       -rho*Y[jn]*rate_eval.n_F17__p_O17
       -rho*Y[jn]*rate_eval.n_F17__He4_N14
       +rate_eval.F17__n_F16
       )

    jac[jn, jf18] = (
       -rho*Y[jn]*rate_eval.n_F18__F19
       -rho*Y[jn]*rate_eval.n_F18__p_O18
       -rho*Y[jn]*rate_eval.n_F18__He4_N15
       +rate_eval.F18__n_F17
       +rho*Y[jp]*rate_eval.p_F18__n_Ne18
       +rho*Y[jhe4]*rate_eval.He4_F18__n_Na21
       )

    jac[jn, jf19] = (
       +rate_eval.F19__n_F18
       +rho*Y[jp]*rate_eval.p_F19__n_Ne19
       +rho*Y[jhe4]*rate_eval.He4_F19__n_Na22
       )

    jac[jn, jne18] = (
       -rho*Y[jn]*rate_eval.n_Ne18__Ne19
       -rho*Y[jn]*rate_eval.n_Ne18__p_F18
       -rho*Y[jn]*rate_eval.n_Ne18__He4_O15
       )

    jac[jn, jne19] = (
       -rho*Y[jn]*rate_eval.n_Ne19__Ne20
       -rho*Y[jn]*rate_eval.n_Ne19__p_F19
       -rho*Y[jn]*rate_eval.n_Ne19__He4_O16
       +rate_eval.Ne19__n_Ne18
       +rho*Y[jhe4]*rate_eval.He4_Ne19__n_Mg22
       )

    jac[jn, jne20] = (
       -rho*Y[jn]*rate_eval.n_Ne20__Ne21
       -rho*Y[jn]*rate_eval.n_Ne20__He4_O17
       +rate_eval.Ne20__n_Ne19
       +rho*Y[jhe4]*rate_eval.He4_Ne20__n_Mg23
       +rho*Y[jc12]*rate_eval.C12_Ne20__n_S31
       )

    jac[jn, jne21] = (
       -rho*Y[jn]*rate_eval.n_Ne21__He4_O18
       +rate_eval.Ne21__n_Ne20
       +rho*Y[jp]*rate_eval.p_Ne21__n_Na21
       +rho*Y[jhe4]*rate_eval.He4_Ne21__n_Mg24
       )

    jac[jn, jna21] = (
       -rho*Y[jn]*rate_eval.n_Na21__Na22
       -rho*Y[jn]*rate_eval.n_Na21__p_Ne21
       -rho*Y[jn]*rate_eval.n_Na21__He4_F18
       )

    jac[jn, jna22] = (
       -rho*Y[jn]*rate_eval.n_Na22__Na23
       -rho*Y[jn]*rate_eval.n_Na22__He4_F19
       +rate_eval.Na22__n_Na21
       +rho*Y[jp]*rate_eval.p_Na22__n_Mg22
       +rho*Y[jhe4]*rate_eval.He4_Na22__n_Al25
       )

    jac[jn, jna23] = (
       -rho*Y[jn]*rate_eval.n_Na23__Na24
       +rate_eval.Na23__n_Na22
       +rho*Y[jp]*rate_eval.p_Na23__n_Mg23
       +rho*Y[jhe4]*rate_eval.He4_Na23__n_Al26
       )

    jac[jn, jna24] = (
       +rate_eval.Na24__n_Na23
       +rho*Y[jp]*rate_eval.p_Na24__n_Mg24
       +rho*Y[jhe4]*rate_eval.He4_Na24__n_Al27
       )

    jac[jn, jmg22] = (
       -rho*Y[jn]*rate_eval.n_Mg22__Mg23
       -rho*Y[jn]*rate_eval.n_Mg22__p_Na22
       -rho*Y[jn]*rate_eval.n_Mg22__He4_Ne19
       )

    jac[jn, jmg23] = (
       -rho*Y[jn]*rate_eval.n_Mg23__Mg24
       -rho*Y[jn]*rate_eval.n_Mg23__p_Na23
       -rho*Y[jn]*rate_eval.n_Mg23__He4_Ne20
       -rho*Y[jn]*rate_eval.n_Mg23__C12_C12
       +rate_eval.Mg23__n_Mg22
       +rho*Y[jhe4]*rate_eval.He4_Mg23__n_Si26
       )

    jac[jn, jmg24] = (
       -rho*Y[jn]*rate_eval.n_Mg24__Mg25
       -rho*Y[jn]*rate_eval.n_Mg24__p_Na24
       -rho*Y[jn]*rate_eval.n_Mg24__He4_Ne21
       +rate_eval.Mg24__n_Mg23
       )

    jac[jn, jmg25] = (
       -rho*Y[jn]*rate_eval.n_Mg25__Mg26
       +rate_eval.Mg25__n_Mg24
       +rho*Y[jp]*rate_eval.p_Mg25__n_Al25
       +rho*Y[jhe4]*rate_eval.He4_Mg25__n_Si28
       )

    jac[jn, jmg26] = (
       +rate_eval.Mg26__n_Mg25
       +rho*Y[jp]*rate_eval.p_Mg26__n_Al26
       +rho*Y[jhe4]*rate_eval.He4_Mg26__n_Si29
       )

    jac[jn, jal25] = (
       -rho*Y[jn]*rate_eval.n_Al25__Al26
       -rho*Y[jn]*rate_eval.n_Al25__p_Mg25
       -rho*Y[jn]*rate_eval.n_Al25__He4_Na22
       )

    jac[jn, jal26] = (
       -rho*Y[jn]*rate_eval.n_Al26__Al27
       -rho*Y[jn]*rate_eval.n_Al26__p_Mg26
       -rho*Y[jn]*rate_eval.n_Al26__He4_Na23
       +rate_eval.Al26__n_Al25
       +rho*Y[jp]*rate_eval.p_Al26__n_Si26
       +rho*Y[jhe4]*rate_eval.He4_Al26__n_P29
       )

    jac[jn, jal27] = (
       -rho*Y[jn]*rate_eval.n_Al27__Al28
       -rho*Y[jn]*rate_eval.n_Al27__He4_Na24
       +rate_eval.Al27__n_Al26
       +rho*Y[jhe4]*rate_eval.He4_Al27__n_P30
       )

    jac[jn, jal28] = (
       +rate_eval.Al28__n_Al27
       +rho*Y[jp]*rate_eval.p_Al28__n_Si28
       +rho*Y[jhe4]*rate_eval.He4_Al28__n_P31
       )

    jac[jn, jsi26] = (
       -rho*Y[jn]*rate_eval.n_Si26__p_Al26
       -rho*Y[jn]*rate_eval.n_Si26__He4_Mg23
       )

    jac[jn, jsi28] = (
       -rho*Y[jn]*rate_eval.n_Si28__Si29
       -rho*Y[jn]*rate_eval.n_Si28__p_Al28
       -rho*Y[jn]*rate_eval.n_Si28__He4_Mg25
       +rho*Y[jhe4]*rate_eval.He4_Si28__n_S31
       )

    jac[jn, jsi29] = (
       -rho*Y[jn]*rate_eval.n_Si29__Si30
       -rho*Y[jn]*rate_eval.n_Si29__He4_Mg26
       +rate_eval.Si29__n_Si28
       +rho*Y[jp]*rate_eval.p_Si29__n_P29
       +rho*Y[jhe4]*rate_eval.He4_Si29__n_S32
       )

    jac[jn, jsi30] = (
       +rate_eval.Si30__n_Si29
       +rho*Y[jp]*rate_eval.p_Si30__n_P30
       )

    jac[jn, jp29] = (
       -rho*Y[jn]*rate_eval.n_P29__P30
       -rho*Y[jn]*rate_eval.n_P29__p_Si29
       -rho*Y[jn]*rate_eval.n_P29__He4_Al26
       +rho*Y[jhe4]*rate_eval.He4_P29__n_Cl32
       )

    jac[jn, jp30] = (
       -rho*Y[jn]*rate_eval.n_P30__P31
       -rho*Y[jn]*rate_eval.n_P30__p_Si30
       -rho*Y[jn]*rate_eval.n_P30__He4_Al27
       +rate_eval.P30__n_P29
       )

    jac[jn, jp31] = (
       -rho*Y[jn]*rate_eval.n_P31__P32
       -rho*Y[jn]*rate_eval.n_P31__He4_Al28
       +rate_eval.P31__n_P30
       +rho*Y[jp]*rate_eval.p_P31__n_S31
       +rho*Y[jhe4]*rate_eval.He4_P31__n_Cl34
       )

    jac[jn, jp32] = (
       +rate_eval.P32__n_P31
       +rho*Y[jp]*rate_eval.p_P32__n_S32
       +rho*Y[jhe4]*rate_eval.He4_P32__n_Cl35
       )

    jac[jn, js31] = (
       -rho*Y[jn]*rate_eval.n_S31__S32
       -rho*Y[jn]*rate_eval.n_S31__p_P31
       -rho*Y[jn]*rate_eval.n_S31__He4_Si28
       -rho*Y[jn]*rate_eval.n_S31__C12_Ne20
       -rho*Y[jn]*rate_eval.n_S31__O16_O16
       )

    jac[jn, js32] = (
       -rho*Y[jn]*rate_eval.n_S32__p_P32
       -rho*Y[jn]*rate_eval.n_S32__He4_Si29
       +rate_eval.S32__n_S31
       +rho*Y[jp]*rate_eval.p_S32__n_Cl32
       +rho*Y[jhe4]*rate_eval.He4_S32__n_Ar35
       )

    jac[jn, js34] = (
       +rho*Y[jp]*rate_eval.p_S34__n_Cl34
       )

    jac[jn, jcl32] = (
       -rho*Y[jn]*rate_eval.n_Cl32__p_S32
       -rho*Y[jn]*rate_eval.n_Cl32__He4_P29
       )

    jac[jn, jcl34] = (
       -rho*Y[jn]*rate_eval.n_Cl34__Cl35
       -rho*Y[jn]*rate_eval.n_Cl34__p_S34
       -rho*Y[jn]*rate_eval.n_Cl34__He4_P31
       )

    jac[jn, jcl35] = (
       -rho*Y[jn]*rate_eval.n_Cl35__He4_P32
       +rate_eval.Cl35__n_Cl34
       +rho*Y[jp]*rate_eval.p_Cl35__n_Ar35
       )

    jac[jn, jar35] = (
       -rho*Y[jn]*rate_eval.n_Ar35__p_Cl35
       -rho*Y[jn]*rate_eval.n_Ar35__He4_S32
       )

    jac[jp, jn] = (
       -rho*Y[jp]*rate_eval.n_p__d
       -rho**2*Y[jp]*Y[jhe4]*rate_eval.n_p_He4__Li6
       -2*5.00000000000000e-01*rho**2*Y[jp]**2*rate_eval.n_p_p__p_d
       -rho**2*Y[jp]*Y[jhe4]*rate_eval.n_p_He4__t_He3
       -5.00000000000000e-01*rho**3*Y[jp]*Y[jhe4]**2*rate_eval.n_p_He4_He4__He3_Li7
       -5.00000000000000e-01*rho**3*Y[jp]*Y[jhe4]**2*rate_eval.n_p_He4_He4__t_Be7
       -5.00000000000000e-01*rho**3*Y[jp]*Y[jhe4]**2*rate_eval.n_p_He4_He4__p_Be9
       +rate_eval.n__p__weak__wc12
       +rho*Y[jhe3]*rate_eval.n_He3__p_t
       +rho*Y[jbe7]*rate_eval.n_Be7__p_Li7
       +rho*Y[jn13]*rate_eval.n_N13__p_C13
       +rho*Y[jn14]*rate_eval.n_N14__p_C14
       +rho*Y[jo14]*rate_eval.n_O14__p_N14
       +rho*Y[jo15]*rate_eval.n_O15__p_N15
       +rho*Y[jf15]*rate_eval.n_F15__p_O15
       +rho*Y[jf16]*rate_eval.n_F16__p_O16
       +rho*Y[jf17]*rate_eval.n_F17__p_O17
       +rho*Y[jf18]*rate_eval.n_F18__p_O18
       +rho*Y[jne18]*rate_eval.n_Ne18__p_F18
       +rho*Y[jne19]*rate_eval.n_Ne19__p_F19
       +rho*Y[jna21]*rate_eval.n_Na21__p_Ne21
       +rho*Y[jmg22]*rate_eval.n_Mg22__p_Na22
       +rho*Y[jmg23]*rate_eval.n_Mg23__p_Na23
       +rho*Y[jmg24]*rate_eval.n_Mg24__p_Na24
       +rho*Y[jal25]*rate_eval.n_Al25__p_Mg25
       +rho*Y[jal26]*rate_eval.n_Al26__p_Mg26
       +rho*Y[jsi26]*rate_eval.n_Si26__p_Al26
       +rho*Y[jsi28]*rate_eval.n_Si28__p_Al28
       +rho*Y[jp29]*rate_eval.n_P29__p_Si29
       +rho*Y[jp30]*rate_eval.n_P30__p_Si30
       +rho*Y[js31]*rate_eval.n_S31__p_P31
       +rho*Y[js32]*rate_eval.n_S32__p_P32
       +rho*Y[jcl32]*rate_eval.n_Cl32__p_S32
       +rho*Y[jcl34]*rate_eval.n_Cl34__p_S34
       +rho*Y[jar35]*rate_eval.n_Ar35__p_Cl35
       +rho*Y[jb8]*rate_eval.n_B8__p_He4_He4
       +5.00000000000000e-01*rho**2*Y[jp]**2*rate_eval.n_p_p__p_d
       +5.00000000000000e-01*rho**3*Y[jp]*Y[jhe4]**2*rate_eval.n_p_He4_He4__p_Be9
       )

    jac[jp, jp] = (
       -rho*Y[jn]*rate_eval.n_p__d
       -2*5.00000000000000e-01*rho*2*Y[jp]*rate_eval.p_p__d__weak__bet_pos_
       -2*5.00000000000000e-01*rho**2*ye(Y)*2*Y[jp]*rate_eval.p_p__d__weak__electron_capture
       -rho*Y[jd]*rate_eval.p_d__He3
       -rho*Y[jt]*rate_eval.p_t__He4
       -rho*Y[jhe3]*rate_eval.p_He3__He4__weak__bet_pos_
       -rho*Y[jli6]*rate_eval.p_Li6__Be7
       -rho*Y[jbe7]*rate_eval.p_Be7__B8
       -rho*Y[jc12]*rate_eval.p_C12__N13
       -rho*Y[jc13]*rate_eval.p_C13__N14
       -rho*Y[jc14]*rate_eval.p_C14__N15
       -rho*Y[jn13]*rate_eval.p_N13__O14
       -rho*Y[jn14]*rate_eval.p_N14__O15
       -rho*Y[jn15]*rate_eval.p_N15__O16
       -rho*Y[jo14]*rate_eval.p_O14__F15
       -rho*Y[jo15]*rate_eval.p_O15__F16
       -rho*Y[jo16]*rate_eval.p_O16__F17
       -rho*Y[jo17]*rate_eval.p_O17__F18
       -rho*Y[jo18]*rate_eval.p_O18__F19
       -rho*Y[jf17]*rate_eval.p_F17__Ne18
       -rho*Y[jf18]*rate_eval.p_F18__Ne19
       -rho*Y[jf19]*rate_eval.p_F19__Ne20
       -rho*Y[jne20]*rate_eval.p_Ne20__Na21
       -rho*Y[jne21]*rate_eval.p_Ne21__Na22
       -rho*Y[jna21]*rate_eval.p_Na21__Mg22
       -rho*Y[jna22]*rate_eval.p_Na22__Mg23
       -rho*Y[jna23]*rate_eval.p_Na23__Mg24
       -rho*Y[jna24]*rate_eval.p_Na24__Mg25
       -rho*Y[jmg24]*rate_eval.p_Mg24__Al25
       -rho*Y[jmg25]*rate_eval.p_Mg25__Al26
       -rho*Y[jmg26]*rate_eval.p_Mg26__Al27
       -rho*Y[jal25]*rate_eval.p_Al25__Si26
       -rho*Y[jal27]*rate_eval.p_Al27__Si28
       -rho*Y[jal28]*rate_eval.p_Al28__Si29
       -rho*Y[jsi28]*rate_eval.p_Si28__P29
       -rho*Y[jsi29]*rate_eval.p_Si29__P30
       -rho*Y[jsi30]*rate_eval.p_Si30__P31
       -rho*Y[jp30]*rate_eval.p_P30__S31
       -rho*Y[jp31]*rate_eval.p_P31__S32
       -rho*Y[js31]*rate_eval.p_S31__Cl32
       -rho*Y[js34]*rate_eval.p_S34__Cl35
       -rho*Y[jcl34]*rate_eval.p_Cl34__Ar35
       -rho*Y[jt]*rate_eval.p_t__n_He3
       -rho*Y[jt]*rate_eval.p_t__d_d
       -rho*Y[jhe4]*rate_eval.p_He4__d_He3
       -rho*Y[jli6]*rate_eval.p_Li6__He4_He3
       -rho*Y[jli7]*rate_eval.p_Li7__n_Be7
       -rho*Y[jli7]*rate_eval.p_Li7__d_Li6
       -rho*Y[jli7]*rate_eval.p_Li7__He4_He4
       -rho*Y[jbe9]*rate_eval.p_Be9__He4_Li6
       -rho*Y[jc11]*rate_eval.p_C11__He4_B8
       -rho*Y[jc13]*rate_eval.p_C13__n_N13
       -rho*Y[jc14]*rate_eval.p_C14__n_N14
       -rho*Y[jn14]*rate_eval.p_N14__n_O14
       -rho*Y[jn14]*rate_eval.p_N14__He4_C11
       -rho*Y[jn15]*rate_eval.p_N15__n_O15
       -rho*Y[jn15]*rate_eval.p_N15__He4_C12
       -rho*Y[jo15]*rate_eval.p_O15__n_F15
       -rho*Y[jo16]*rate_eval.p_O16__n_F16
       -rho*Y[jo16]*rate_eval.p_O16__He4_N13
       -rho*Y[jo17]*rate_eval.p_O17__n_F17
       -rho*Y[jo17]*rate_eval.p_O17__He4_N14
       -rho*Y[jo18]*rate_eval.p_O18__n_F18
       -rho*Y[jo18]*rate_eval.p_O18__He4_N15
       -rho*Y[jf17]*rate_eval.p_F17__He4_O14
       -rho*Y[jf18]*rate_eval.p_F18__n_Ne18
       -rho*Y[jf18]*rate_eval.p_F18__He4_O15
       -rho*Y[jf19]*rate_eval.p_F19__n_Ne19
       -rho*Y[jf19]*rate_eval.p_F19__He4_O16
       -rho*Y[jne18]*rate_eval.p_Ne18__He4_F15
       -rho*Y[jne19]*rate_eval.p_Ne19__He4_F16
       -rho*Y[jne20]*rate_eval.p_Ne20__He4_F17
       -rho*Y[jne21]*rate_eval.p_Ne21__n_Na21
       -rho*Y[jne21]*rate_eval.p_Ne21__He4_F18
       -rho*Y[jna21]*rate_eval.p_Na21__He4_Ne18
       -rho*Y[jna22]*rate_eval.p_Na22__n_Mg22
       -rho*Y[jna22]*rate_eval.p_Na22__He4_Ne19
       -rho*Y[jna23]*rate_eval.p_Na23__n_Mg23
       -rho*Y[jna23]*rate_eval.p_Na23__He4_Ne20
       -rho*Y[jna23]*rate_eval.p_Na23__C12_C12
       -rho*Y[jna24]*rate_eval.p_Na24__n_Mg24
       -rho*Y[jna24]*rate_eval.p_Na24__He4_Ne21
       -rho*Y[jmg24]*rate_eval.p_Mg24__He4_Na21
       -rho*Y[jmg25]*rate_eval.p_Mg25__n_Al25
       -rho*Y[jmg25]*rate_eval.p_Mg25__He4_Na22
       -rho*Y[jmg26]*rate_eval.p_Mg26__n_Al26
       -rho*Y[jmg26]*rate_eval.p_Mg26__He4_Na23
       -rho*Y[jal25]*rate_eval.p_Al25__He4_Mg22
       -rho*Y[jal26]*rate_eval.p_Al26__n_Si26
       -rho*Y[jal26]*rate_eval.p_Al26__He4_Mg23
       -rho*Y[jal27]*rate_eval.p_Al27__He4_Mg24
       -rho*Y[jal27]*rate_eval.p_Al27__C12_O16
       -rho*Y[jal28]*rate_eval.p_Al28__n_Si28
       -rho*Y[jal28]*rate_eval.p_Al28__He4_Mg25
       -rho*Y[jsi28]*rate_eval.p_Si28__He4_Al25
       -rho*Y[jsi29]*rate_eval.p_Si29__n_P29
       -rho*Y[jsi29]*rate_eval.p_Si29__He4_Al26
       -rho*Y[jsi30]*rate_eval.p_Si30__n_P30
       -rho*Y[jsi30]*rate_eval.p_Si30__He4_Al27
       -rho*Y[jp29]*rate_eval.p_P29__He4_Si26
       -rho*Y[jp31]*rate_eval.p_P31__n_S31
       -rho*Y[jp31]*rate_eval.p_P31__He4_Si28
       -rho*Y[jp31]*rate_eval.p_P31__C12_Ne20
       -rho*Y[jp31]*rate_eval.p_P31__O16_O16
       -rho*Y[jp32]*rate_eval.p_P32__n_S32
       -rho*Y[jp32]*rate_eval.p_P32__He4_Si29
       -rho*Y[js32]*rate_eval.p_S32__n_Cl32
       -rho*Y[js32]*rate_eval.p_S32__He4_P29
       -rho*Y[js34]*rate_eval.p_S34__n_Cl34
       -rho*Y[js34]*rate_eval.p_S34__He4_P31
       -rho*Y[jcl34]*rate_eval.p_Cl34__He4_S31
       -rho*Y[jcl35]*rate_eval.p_Cl35__n_Ar35
       -rho*Y[jcl35]*rate_eval.p_Cl35__He4_S32
       -rho*Y[jar35]*rate_eval.p_Ar35__He4_Cl32
       -rho*Y[jd]*rate_eval.p_d__n_p_p
       -rho*Y[jbe9]*rate_eval.p_Be9__d_He4_He4
       -rho*Y[jbe9]*rate_eval.p_Be9__n_p_He4_He4
       -rho**2*Y[jn]*Y[jhe4]*rate_eval.n_p_He4__Li6
       -2*5.00000000000000e-01*rho**2*Y[jn]*2*Y[jp]*rate_eval.n_p_p__p_d
       -rho**2*Y[jn]*Y[jhe4]*rate_eval.n_p_He4__t_He3
       -2*5.00000000000000e-01*rho**2*2*Y[jp]*Y[jhe4]*rate_eval.p_p_He4__He3_He3
       -5.00000000000000e-01*rho**2*Y[jhe4]**2*rate_eval.p_He4_He4__n_B8
       -5.00000000000000e-01*rho**2*Y[jhe4]**2*rate_eval.p_He4_He4__d_Be7
       -5.00000000000000e-01*rho**3*Y[jn]*Y[jhe4]**2*rate_eval.n_p_He4_He4__He3_Li7
       -5.00000000000000e-01*rho**3*Y[jn]*Y[jhe4]**2*rate_eval.n_p_He4_He4__t_Be7
       -5.00000000000000e-01*rho**3*Y[jn]*Y[jhe4]**2*rate_eval.n_p_He4_He4__p_Be9
       -2*2.50000000000000e-01*rho**3*2*Y[jp]*Y[jhe4]**2*rate_eval.p_p_He4_He4__He3_Be7
       +2*rho*Y[jd]*rate_eval.p_d__n_p_p
       +rho*Y[jbe9]*rate_eval.p_Be9__n_p_He4_He4
       +5.00000000000000e-01*rho**2*Y[jn]*2*Y[jp]*rate_eval.n_p_p__p_d
       +5.00000000000000e-01*rho**3*Y[jn]*Y[jhe4]**2*rate_eval.n_p_He4_He4__p_Be9
       )

    jac[jp, jd] = (
       -rho*Y[jp]*rate_eval.p_d__He3
       -rho*Y[jp]*rate_eval.p_d__n_p_p
       +rate_eval.d__n_p
       +5.00000000000000e-01*rho*2*Y[jd]*rate_eval.d_d__p_t
       +rho*Y[jhe3]*rate_eval.d_He3__p_He4
       +rho*Y[jli6]*rate_eval.d_Li6__p_Li7
       +2*rho*Y[jp]*rate_eval.p_d__n_p_p
       +rho*Y[jbe7]*rate_eval.d_Be7__p_He4_He4
       +5.00000000000000e-01*rho**2*Y[jhe4]**2*rate_eval.d_He4_He4__p_Be9
       )

    jac[jp, jt] = (
       -rho*Y[jp]*rate_eval.p_t__He4
       -rho*Y[jp]*rate_eval.p_t__n_He3
       -rho*Y[jp]*rate_eval.p_t__d_d
       +rho*Y[jhe3]*rate_eval.t_He3__n_p_He4
       +rho*Y[jbe7]*rate_eval.t_Be7__n_p_He4_He4
       )

    jac[jp, jhe3] = (
       -rho*Y[jp]*rate_eval.p_He3__He4__weak__bet_pos_
       +rate_eval.He3__p_d
       +rho*Y[jn]*rate_eval.n_He3__p_t
       +rho*Y[jd]*rate_eval.d_He3__p_He4
       +rho*Y[jhe4]*rate_eval.He4_He3__p_Li6
       +rho*Y[jt]*rate_eval.t_He3__n_p_He4
       +2*5.00000000000000e-01*rho*2*Y[jhe3]*rate_eval.He3_He3__p_p_He4
       +rho*Y[jli7]*rate_eval.He3_Li7__n_p_He4_He4
       +2*rho*Y[jbe7]*rate_eval.He3_Be7__p_p_He4_He4
       )

    jac[jp, jhe4] = (
       -rho*Y[jp]*rate_eval.p_He4__d_He3
       -rho**2*Y[jn]*Y[jp]*rate_eval.n_p_He4__Li6
       -rho**2*Y[jn]*Y[jp]*rate_eval.n_p_He4__t_He3
       -2*5.00000000000000e-01*rho**2*Y[jp]**2*rate_eval.p_p_He4__He3_He3
       -5.00000000000000e-01*rho**2*Y[jp]*2*Y[jhe4]*rate_eval.p_He4_He4__n_B8
       -5.00000000000000e-01*rho**2*Y[jp]*2*Y[jhe4]*rate_eval.p_He4_He4__d_Be7
       -5.00000000000000e-01*rho**3*Y[jn]*Y[jp]*2*Y[jhe4]*rate_eval.n_p_He4_He4__He3_Li7
       -5.00000000000000e-01*rho**3*Y[jn]*Y[jp]*2*Y[jhe4]*rate_eval.n_p_He4_He4__t_Be7
       -5.00000000000000e-01*rho**3*Y[jn]*Y[jp]*2*Y[jhe4]*rate_eval.n_p_He4_He4__p_Be9
       -2*2.50000000000000e-01*rho**3*Y[jp]**2*2*Y[jhe4]*rate_eval.p_p_He4_He4__He3_Be7
       +rate_eval.He4__p_t
       +rho*Y[jhe3]*rate_eval.He4_He3__p_Li6
       +5.00000000000000e-01*rho*2*Y[jhe4]*rate_eval.He4_He4__p_Li7
       +rho*Y[jli6]*rate_eval.He4_Li6__p_Be9
       +rho*Y[jb8]*rate_eval.He4_B8__p_C11
       +rho*Y[jc11]*rate_eval.He4_C11__p_N14
       +rho*Y[jc12]*rate_eval.He4_C12__p_N15
       +rho*Y[jn13]*rate_eval.He4_N13__p_O16
       +rho*Y[jn14]*rate_eval.He4_N14__p_O17
       +rho*Y[jn15]*rate_eval.He4_N15__p_O18
       +rho*Y[jo14]*rate_eval.He4_O14__p_F17
       +rho*Y[jo15]*rate_eval.He4_O15__p_F18
       +rho*Y[jo16]*rate_eval.He4_O16__p_F19
       +rho*Y[jf15]*rate_eval.He4_F15__p_Ne18
       +rho*Y[jf16]*rate_eval.He4_F16__p_Ne19
       +rho*Y[jf17]*rate_eval.He4_F17__p_Ne20
       +rho*Y[jf18]*rate_eval.He4_F18__p_Ne21
       +rho*Y[jne18]*rate_eval.He4_Ne18__p_Na21
       +rho*Y[jne19]*rate_eval.He4_Ne19__p_Na22
       +rho*Y[jne20]*rate_eval.He4_Ne20__p_Na23
       +rho*Y[jne21]*rate_eval.He4_Ne21__p_Na24
       +rho*Y[jna21]*rate_eval.He4_Na21__p_Mg24
       +rho*Y[jna22]*rate_eval.He4_Na22__p_Mg25
       +rho*Y[jna23]*rate_eval.He4_Na23__p_Mg26
       +rho*Y[jmg22]*rate_eval.He4_Mg22__p_Al25
       +rho*Y[jmg23]*rate_eval.He4_Mg23__p_Al26
       +rho*Y[jmg24]*rate_eval.He4_Mg24__p_Al27
       +rho*Y[jmg25]*rate_eval.He4_Mg25__p_Al28
       +rho*Y[jal25]*rate_eval.He4_Al25__p_Si28
       +rho*Y[jal26]*rate_eval.He4_Al26__p_Si29
       +rho*Y[jal27]*rate_eval.He4_Al27__p_Si30
       +rho*Y[jsi26]*rate_eval.He4_Si26__p_P29
       +rho*Y[jsi28]*rate_eval.He4_Si28__p_P31
       +rho*Y[jsi29]*rate_eval.He4_Si29__p_P32
       +rho*Y[jp29]*rate_eval.He4_P29__p_S32
       +rho*Y[jp31]*rate_eval.He4_P31__p_S34
       +rho*Y[js31]*rate_eval.He4_S31__p_Cl34
       +rho*Y[js32]*rate_eval.He4_S32__p_Cl35
       +rho*Y[jcl32]*rate_eval.He4_Cl32__p_Ar35
       +5.00000000000000e-01*rho**2*Y[jd]*2*Y[jhe4]*rate_eval.d_He4_He4__p_Be9
       +5.00000000000000e-01*rho**3*Y[jn]*Y[jp]*2*Y[jhe4]*rate_eval.n_p_He4_He4__p_Be9
       )

    jac[jp, jli6] = (
       -rho*Y[jp]*rate_eval.p_Li6__Be7
       -rho*Y[jp]*rate_eval.p_Li6__He4_He3
       +rate_eval.Li6__n_p_He4
       +rho*Y[jd]*rate_eval.d_Li6__p_Li7
       +rho*Y[jhe4]*rate_eval.He4_Li6__p_Be9
       )

    jac[jp, jli7] = (
       -rho*Y[jp]*rate_eval.p_Li7__n_Be7
       -rho*Y[jp]*rate_eval.p_Li7__d_Li6
       -rho*Y[jp]*rate_eval.p_Li7__He4_He4
       +rho*Y[jhe3]*rate_eval.He3_Li7__n_p_He4_He4
       )

    jac[jp, jbe7] = (
       -rho*Y[jp]*rate_eval.p_Be7__B8
       +rate_eval.Be7__p_Li6
       +rho*Y[jn]*rate_eval.n_Be7__p_Li7
       +rho*Y[jd]*rate_eval.d_Be7__p_He4_He4
       +rho*Y[jt]*rate_eval.t_Be7__n_p_He4_He4
       +2*rho*Y[jhe3]*rate_eval.He3_Be7__p_p_He4_He4
       )

    jac[jp, jbe9] = (
       -rho*Y[jp]*rate_eval.p_Be9__He4_Li6
       -rho*Y[jp]*rate_eval.p_Be9__d_He4_He4
       -rho*Y[jp]*rate_eval.p_Be9__n_p_He4_He4
       +rho*Y[jp]*rate_eval.p_Be9__n_p_He4_He4
       )

    jac[jp, jb8] = (
       +rate_eval.B8__p_Be7
       +rho*Y[jhe4]*rate_eval.He4_B8__p_C11
       +rho*Y[jn]*rate_eval.n_B8__p_He4_He4
       )

    jac[jp, jc11] = (
       -rho*Y[jp]*rate_eval.p_C11__He4_B8
       +rho*Y[jhe4]*rate_eval.He4_C11__p_N14
       )

    jac[jp, jc12] = (
       -rho*Y[jp]*rate_eval.p_C12__N13
       +rho*Y[jhe4]*rate_eval.He4_C12__p_N15
       +5.00000000000000e-01*rho*2*Y[jc12]*rate_eval.C12_C12__p_Na23
       +rho*Y[jo16]*rate_eval.C12_O16__p_Al27
       +rho*Y[jne20]*rate_eval.C12_Ne20__p_P31
       )

    jac[jp, jc13] = (
       -rho*Y[jp]*rate_eval.p_C13__N14
       -rho*Y[jp]*rate_eval.p_C13__n_N13
       )

    jac[jp, jc14] = (
       -rho*Y[jp]*rate_eval.p_C14__N15
       -rho*Y[jp]*rate_eval.p_C14__n_N14
       )

    jac[jp, jn13] = (
       -rho*Y[jp]*rate_eval.p_N13__O14
       +rate_eval.N13__p_C12
       +rho*Y[jn]*rate_eval.n_N13__p_C13
       +rho*Y[jhe4]*rate_eval.He4_N13__p_O16
       )

    jac[jp, jn14] = (
       -rho*Y[jp]*rate_eval.p_N14__O15
       -rho*Y[jp]*rate_eval.p_N14__n_O14
       -rho*Y[jp]*rate_eval.p_N14__He4_C11
       +rate_eval.N14__p_C13
       +rho*Y[jn]*rate_eval.n_N14__p_C14
       +rho*Y[jhe4]*rate_eval.He4_N14__p_O17
       )

    jac[jp, jn15] = (
       -rho*Y[jp]*rate_eval.p_N15__O16
       -rho*Y[jp]*rate_eval.p_N15__n_O15
       -rho*Y[jp]*rate_eval.p_N15__He4_C12
       +rate_eval.N15__p_C14
       +rho*Y[jhe4]*rate_eval.He4_N15__p_O18
       )

    jac[jp, jo14] = (
       -rho*Y[jp]*rate_eval.p_O14__F15
       +rate_eval.O14__p_N13
       +rho*Y[jn]*rate_eval.n_O14__p_N14
       +rho*Y[jhe4]*rate_eval.He4_O14__p_F17
       )

    jac[jp, jo15] = (
       -rho*Y[jp]*rate_eval.p_O15__F16
       -rho*Y[jp]*rate_eval.p_O15__n_F15
       +rate_eval.O15__p_N14
       +rho*Y[jn]*rate_eval.n_O15__p_N15
       +rho*Y[jhe4]*rate_eval.He4_O15__p_F18
       )

    jac[jp, jo16] = (
       -rho*Y[jp]*rate_eval.p_O16__F17
       -rho*Y[jp]*rate_eval.p_O16__n_F16
       -rho*Y[jp]*rate_eval.p_O16__He4_N13
       +rate_eval.O16__p_N15
       +rho*Y[jhe4]*rate_eval.He4_O16__p_F19
       +rho*Y[jc12]*rate_eval.C12_O16__p_Al27
       +5.00000000000000e-01*rho*2*Y[jo16]*rate_eval.O16_O16__p_P31
       )

    jac[jp, jo17] = (
       -rho*Y[jp]*rate_eval.p_O17__F18
       -rho*Y[jp]*rate_eval.p_O17__n_F17
       -rho*Y[jp]*rate_eval.p_O17__He4_N14
       )

    jac[jp, jo18] = (
       -rho*Y[jp]*rate_eval.p_O18__F19
       -rho*Y[jp]*rate_eval.p_O18__n_F18
       -rho*Y[jp]*rate_eval.p_O18__He4_N15
       )

    jac[jp, jf15] = (
       +rate_eval.F15__p_O14
       +rho*Y[jn]*rate_eval.n_F15__p_O15
       +rho*Y[jhe4]*rate_eval.He4_F15__p_Ne18
       )

    jac[jp, jf16] = (
       +rate_eval.F16__p_O15
       +rho*Y[jn]*rate_eval.n_F16__p_O16
       +rho*Y[jhe4]*rate_eval.He4_F16__p_Ne19
       )

    jac[jp, jf17] = (
       -rho*Y[jp]*rate_eval.p_F17__Ne18
       -rho*Y[jp]*rate_eval.p_F17__He4_O14
       +rate_eval.F17__p_O16
       +rho*Y[jn]*rate_eval.n_F17__p_O17
       +rho*Y[jhe4]*rate_eval.He4_F17__p_Ne20
       )

    jac[jp, jf18] = (
       -rho*Y[jp]*rate_eval.p_F18__Ne19
       -rho*Y[jp]*rate_eval.p_F18__n_Ne18
       -rho*Y[jp]*rate_eval.p_F18__He4_O15
       +rate_eval.F18__p_O17
       +rho*Y[jn]*rate_eval.n_F18__p_O18
       +rho*Y[jhe4]*rate_eval.He4_F18__p_Ne21
       )

    jac[jp, jf19] = (
       -rho*Y[jp]*rate_eval.p_F19__Ne20
       -rho*Y[jp]*rate_eval.p_F19__n_Ne19
       -rho*Y[jp]*rate_eval.p_F19__He4_O16
       +rate_eval.F19__p_O18
       )

    jac[jp, jne18] = (
       -rho*Y[jp]*rate_eval.p_Ne18__He4_F15
       +rate_eval.Ne18__p_F17
       +rho*Y[jn]*rate_eval.n_Ne18__p_F18
       +rho*Y[jhe4]*rate_eval.He4_Ne18__p_Na21
       )

    jac[jp, jne19] = (
       -rho*Y[jp]*rate_eval.p_Ne19__He4_F16
       +rate_eval.Ne19__p_F18
       +rho*Y[jn]*rate_eval.n_Ne19__p_F19
       +rho*Y[jhe4]*rate_eval.He4_Ne19__p_Na22
       )

    jac[jp, jne20] = (
       -rho*Y[jp]*rate_eval.p_Ne20__Na21
       -rho*Y[jp]*rate_eval.p_Ne20__He4_F17
       +rate_eval.Ne20__p_F19
       +rho*Y[jhe4]*rate_eval.He4_Ne20__p_Na23
       +rho*Y[jc12]*rate_eval.C12_Ne20__p_P31
       )

    jac[jp, jne21] = (
       -rho*Y[jp]*rate_eval.p_Ne21__Na22
       -rho*Y[jp]*rate_eval.p_Ne21__n_Na21
       -rho*Y[jp]*rate_eval.p_Ne21__He4_F18
       +rho*Y[jhe4]*rate_eval.He4_Ne21__p_Na24
       )

    jac[jp, jna21] = (
       -rho*Y[jp]*rate_eval.p_Na21__Mg22
       -rho*Y[jp]*rate_eval.p_Na21__He4_Ne18
       +rate_eval.Na21__p_Ne20
       +rho*Y[jn]*rate_eval.n_Na21__p_Ne21
       +rho*Y[jhe4]*rate_eval.He4_Na21__p_Mg24
       )

    jac[jp, jna22] = (
       -rho*Y[jp]*rate_eval.p_Na22__Mg23
       -rho*Y[jp]*rate_eval.p_Na22__n_Mg22
       -rho*Y[jp]*rate_eval.p_Na22__He4_Ne19
       +rate_eval.Na22__p_Ne21
       +rho*Y[jhe4]*rate_eval.He4_Na22__p_Mg25
       )

    jac[jp, jna23] = (
       -rho*Y[jp]*rate_eval.p_Na23__Mg24
       -rho*Y[jp]*rate_eval.p_Na23__n_Mg23
       -rho*Y[jp]*rate_eval.p_Na23__He4_Ne20
       -rho*Y[jp]*rate_eval.p_Na23__C12_C12
       +rho*Y[jhe4]*rate_eval.He4_Na23__p_Mg26
       )

    jac[jp, jna24] = (
       -rho*Y[jp]*rate_eval.p_Na24__Mg25
       -rho*Y[jp]*rate_eval.p_Na24__n_Mg24
       -rho*Y[jp]*rate_eval.p_Na24__He4_Ne21
       )

    jac[jp, jmg22] = (
       +rate_eval.Mg22__p_Na21
       +rho*Y[jn]*rate_eval.n_Mg22__p_Na22
       +rho*Y[jhe4]*rate_eval.He4_Mg22__p_Al25
       )

    jac[jp, jmg23] = (
       +rate_eval.Mg23__p_Na22
       +rho*Y[jn]*rate_eval.n_Mg23__p_Na23
       +rho*Y[jhe4]*rate_eval.He4_Mg23__p_Al26
       )

    jac[jp, jmg24] = (
       -rho*Y[jp]*rate_eval.p_Mg24__Al25
       -rho*Y[jp]*rate_eval.p_Mg24__He4_Na21
       +rate_eval.Mg24__p_Na23
       +rho*Y[jn]*rate_eval.n_Mg24__p_Na24
       +rho*Y[jhe4]*rate_eval.He4_Mg24__p_Al27
       )

    jac[jp, jmg25] = (
       -rho*Y[jp]*rate_eval.p_Mg25__Al26
       -rho*Y[jp]*rate_eval.p_Mg25__n_Al25
       -rho*Y[jp]*rate_eval.p_Mg25__He4_Na22
       +rate_eval.Mg25__p_Na24
       +rho*Y[jhe4]*rate_eval.He4_Mg25__p_Al28
       )

    jac[jp, jmg26] = (
       -rho*Y[jp]*rate_eval.p_Mg26__Al27
       -rho*Y[jp]*rate_eval.p_Mg26__n_Al26
       -rho*Y[jp]*rate_eval.p_Mg26__He4_Na23
       )

    jac[jp, jal25] = (
       -rho*Y[jp]*rate_eval.p_Al25__Si26
       -rho*Y[jp]*rate_eval.p_Al25__He4_Mg22
       +rate_eval.Al25__p_Mg24
       +rho*Y[jn]*rate_eval.n_Al25__p_Mg25
       +rho*Y[jhe4]*rate_eval.He4_Al25__p_Si28
       )

    jac[jp, jal26] = (
       -rho*Y[jp]*rate_eval.p_Al26__n_Si26
       -rho*Y[jp]*rate_eval.p_Al26__He4_Mg23
       +rate_eval.Al26__p_Mg25
       +rho*Y[jn]*rate_eval.n_Al26__p_Mg26
       +rho*Y[jhe4]*rate_eval.He4_Al26__p_Si29
       )

    jac[jp, jal27] = (
       -rho*Y[jp]*rate_eval.p_Al27__Si28
       -rho*Y[jp]*rate_eval.p_Al27__He4_Mg24
       -rho*Y[jp]*rate_eval.p_Al27__C12_O16
       +rate_eval.Al27__p_Mg26
       +rho*Y[jhe4]*rate_eval.He4_Al27__p_Si30
       )

    jac[jp, jal28] = (
       -rho*Y[jp]*rate_eval.p_Al28__Si29
       -rho*Y[jp]*rate_eval.p_Al28__n_Si28
       -rho*Y[jp]*rate_eval.p_Al28__He4_Mg25
       )

    jac[jp, jsi26] = (
       +rate_eval.Si26__p_Al25
       +rho*Y[jn]*rate_eval.n_Si26__p_Al26
       +rho*Y[jhe4]*rate_eval.He4_Si26__p_P29
       )

    jac[jp, jsi28] = (
       -rho*Y[jp]*rate_eval.p_Si28__P29
       -rho*Y[jp]*rate_eval.p_Si28__He4_Al25
       +rate_eval.Si28__p_Al27
       +rho*Y[jn]*rate_eval.n_Si28__p_Al28
       +rho*Y[jhe4]*rate_eval.He4_Si28__p_P31
       )

    jac[jp, jsi29] = (
       -rho*Y[jp]*rate_eval.p_Si29__P30
       -rho*Y[jp]*rate_eval.p_Si29__n_P29
       -rho*Y[jp]*rate_eval.p_Si29__He4_Al26
       +rate_eval.Si29__p_Al28
       +rho*Y[jhe4]*rate_eval.He4_Si29__p_P32
       )

    jac[jp, jsi30] = (
       -rho*Y[jp]*rate_eval.p_Si30__P31
       -rho*Y[jp]*rate_eval.p_Si30__n_P30
       -rho*Y[jp]*rate_eval.p_Si30__He4_Al27
       )

    jac[jp, jp29] = (
       -rho*Y[jp]*rate_eval.p_P29__He4_Si26
       +rate_eval.P29__p_Si28
       +rho*Y[jn]*rate_eval.n_P29__p_Si29
       +rho*Y[jhe4]*rate_eval.He4_P29__p_S32
       )

    jac[jp, jp30] = (
       -rho*Y[jp]*rate_eval.p_P30__S31
       +rate_eval.P30__p_Si29
       +rho*Y[jn]*rate_eval.n_P30__p_Si30
       )

    jac[jp, jp31] = (
       -rho*Y[jp]*rate_eval.p_P31__S32
       -rho*Y[jp]*rate_eval.p_P31__n_S31
       -rho*Y[jp]*rate_eval.p_P31__He4_Si28
       -rho*Y[jp]*rate_eval.p_P31__C12_Ne20
       -rho*Y[jp]*rate_eval.p_P31__O16_O16
       +rate_eval.P31__p_Si30
       +rho*Y[jhe4]*rate_eval.He4_P31__p_S34
       )

    jac[jp, jp32] = (
       -rho*Y[jp]*rate_eval.p_P32__n_S32
       -rho*Y[jp]*rate_eval.p_P32__He4_Si29
       )

    jac[jp, js31] = (
       -rho*Y[jp]*rate_eval.p_S31__Cl32
       +rate_eval.S31__p_P30
       +rho*Y[jn]*rate_eval.n_S31__p_P31
       +rho*Y[jhe4]*rate_eval.He4_S31__p_Cl34
       )

    jac[jp, js32] = (
       -rho*Y[jp]*rate_eval.p_S32__n_Cl32
       -rho*Y[jp]*rate_eval.p_S32__He4_P29
       +rate_eval.S32__p_P31
       +rho*Y[jn]*rate_eval.n_S32__p_P32
       +rho*Y[jhe4]*rate_eval.He4_S32__p_Cl35
       )

    jac[jp, js34] = (
       -rho*Y[jp]*rate_eval.p_S34__Cl35
       -rho*Y[jp]*rate_eval.p_S34__n_Cl34
       -rho*Y[jp]*rate_eval.p_S34__He4_P31
       )

    jac[jp, jcl32] = (
       +rate_eval.Cl32__p_S31
       +rate_eval.Cl32__p_P31__weak__wc12
       +rho*Y[jn]*rate_eval.n_Cl32__p_S32
       +rho*Y[jhe4]*rate_eval.He4_Cl32__p_Ar35
       )

    jac[jp, jcl34] = (
       -rho*Y[jp]*rate_eval.p_Cl34__Ar35
       -rho*Y[jp]*rate_eval.p_Cl34__He4_S31
       +rho*Y[jn]*rate_eval.n_Cl34__p_S34
       )

    jac[jp, jcl35] = (
       -rho*Y[jp]*rate_eval.p_Cl35__n_Ar35
       -rho*Y[jp]*rate_eval.p_Cl35__He4_S32
       +rate_eval.Cl35__p_S34
       )

    jac[jp, jar35] = (
       -rho*Y[jp]*rate_eval.p_Ar35__He4_Cl32
       +rate_eval.Ar35__p_Cl34
       +rho*Y[jn]*rate_eval.n_Ar35__p_Cl35
       )

    jac[jd, jn] = (
       -rho*Y[jd]*rate_eval.n_d__t
       +rho*Y[jp]*rate_eval.n_p__d
       +2*rho*Y[jhe3]*rate_eval.n_He3__d_d
       +rho*Y[jhe4]*rate_eval.n_He4__d_t
       +rho*Y[jbe7]*rate_eval.n_Be7__d_Li6
       +rho*Y[jn14]*rate_eval.n_N14__d_C13
       +rho*Y[jn15]*rate_eval.n_N15__d_C14
       +5.00000000000000e-01*rho**2*Y[jp]**2*rate_eval.n_p_p__p_d
       +5.00000000000000e-01*rho**2*Y[jhe4]**2*rate_eval.n_He4_He4__d_Li7
       )

    jac[jd, jp] = (
       -rho*Y[jd]*rate_eval.p_d__He3
       -rho*Y[jd]*rate_eval.p_d__n_p_p
       +rho*Y[jn]*rate_eval.n_p__d
       +5.00000000000000e-01*rho*2*Y[jp]*rate_eval.p_p__d__weak__bet_pos_
       +5.00000000000000e-01*rho**2*ye(Y)*2*Y[jp]*rate_eval.p_p__d__weak__electron_capture
       +2*rho*Y[jt]*rate_eval.p_t__d_d
       +rho*Y[jhe4]*rate_eval.p_He4__d_He3
       +rho*Y[jli7]*rate_eval.p_Li7__d_Li6
       +rho*Y[jbe9]*rate_eval.p_Be9__d_He4_He4
       +5.00000000000000e-01*rho**2*Y[jn]*2*Y[jp]*rate_eval.n_p_p__p_d
       +5.00000000000000e-01*rho**2*Y[jhe4]**2*rate_eval.p_He4_He4__d_Be7
       )

    jac[jd, jd] = (
       -rate_eval.d__n_p
       -rho*Y[jn]*rate_eval.n_d__t
       -rho*Y[jp]*rate_eval.p_d__He3
       -2*5.00000000000000e-01*rho*2*Y[jd]*rate_eval.d_d__He4
       -rho*Y[jhe4]*rate_eval.He4_d__Li6
       -2*5.00000000000000e-01*rho*2*Y[jd]*rate_eval.d_d__n_He3
       -2*5.00000000000000e-01*rho*2*Y[jd]*rate_eval.d_d__p_t
       -rho*Y[jt]*rate_eval.d_t__n_He4
       -rho*Y[jhe3]*rate_eval.d_He3__p_He4
       -rho*Y[jhe4]*rate_eval.d_He4__t_He3
       -rho*Y[jli6]*rate_eval.d_Li6__n_Be7
       -rho*Y[jli6]*rate_eval.d_Li6__p_Li7
       -rho*Y[jc13]*rate_eval.d_C13__n_N14
       -rho*Y[jc14]*rate_eval.d_C14__n_N15
       -rho*Y[jp]*rate_eval.p_d__n_p_p
       -rho*Y[jli7]*rate_eval.d_Li7__n_He4_He4
       -rho*Y[jbe7]*rate_eval.d_Be7__p_He4_He4
       -5.00000000000000e-01*rho**2*Y[jhe4]**2*rate_eval.d_He4_He4__p_Be9
       )

    jac[jd, jt] = (
       -rho*Y[jd]*rate_eval.d_t__n_He4
       +rate_eval.t__n_d
       +2*rho*Y[jp]*rate_eval.p_t__d_d
       +rho*Y[jhe3]*rate_eval.t_He3__d_He4
       )

    jac[jd, jhe3] = (
       -rho*Y[jd]*rate_eval.d_He3__p_He4
       +rate_eval.He3__p_d
       +2*rho*Y[jn]*rate_eval.n_He3__d_d
       +rho*Y[jt]*rate_eval.t_He3__d_He4
       )

    jac[jd, jhe4] = (
       -rho*Y[jd]*rate_eval.He4_d__Li6
       -rho*Y[jd]*rate_eval.d_He4__t_He3
       -5.00000000000000e-01*rho**2*Y[jd]*2*Y[jhe4]*rate_eval.d_He4_He4__p_Be9
       +2*rate_eval.He4__d_d
       +rho*Y[jn]*rate_eval.n_He4__d_t
       +rho*Y[jp]*rate_eval.p_He4__d_He3
       +5.00000000000000e-01*rho**2*Y[jn]*2*Y[jhe4]*rate_eval.n_He4_He4__d_Li7
       +5.00000000000000e-01*rho**2*Y[jp]*2*Y[jhe4]*rate_eval.p_He4_He4__d_Be7
       )

    jac[jd, jli6] = (
       -rho*Y[jd]*rate_eval.d_Li6__n_Be7
       -rho*Y[jd]*rate_eval.d_Li6__p_Li7
       +rate_eval.Li6__He4_d
       )

    jac[jd, jli7] = (
       -rho*Y[jd]*rate_eval.d_Li7__n_He4_He4
       +rho*Y[jp]*rate_eval.p_Li7__d_Li6
       )

    jac[jd, jbe7] = (
       -rho*Y[jd]*rate_eval.d_Be7__p_He4_He4
       +rho*Y[jn]*rate_eval.n_Be7__d_Li6
       )

    jac[jd, jbe9] = (
       +rho*Y[jp]*rate_eval.p_Be9__d_He4_He4
       )

    jac[jd, jc13] = (
       -rho*Y[jd]*rate_eval.d_C13__n_N14
       )

    jac[jd, jc14] = (
       -rho*Y[jd]*rate_eval.d_C14__n_N15
       )

    jac[jd, jn14] = (
       +rho*Y[jn]*rate_eval.n_N14__d_C13
       )

    jac[jd, jn15] = (
       +rho*Y[jn]*rate_eval.n_N15__d_C14
       )

    jac[jt, jn] = (
       +rho*Y[jd]*rate_eval.n_d__t
       +rho*Y[jhe3]*rate_eval.n_He3__p_t
       +rho*Y[jhe4]*rate_eval.n_He4__d_t
       +rho*Y[jli6]*rate_eval.n_Li6__He4_t
       +rho*Y[jbe9]*rate_eval.n_Be9__t_Li7
       +2*5.00000000000000e-01*rho**2*2*Y[jn]*Y[jhe4]*rate_eval.n_n_He4__t_t
       +rho**2*Y[jp]*Y[jhe4]*rate_eval.n_p_He4__t_He3
       +2.50000000000000e-01*rho**3*2*Y[jn]*Y[jhe4]**2*rate_eval.n_n_He4_He4__t_Li7
       +5.00000000000000e-01*rho**3*Y[jp]*Y[jhe4]**2*rate_eval.n_p_He4_He4__t_Be7
       )

    jac[jt, jp] = (
       -rho*Y[jt]*rate_eval.p_t__He4
       -rho*Y[jt]*rate_eval.p_t__n_He3
       -rho*Y[jt]*rate_eval.p_t__d_d
       +rho**2*Y[jn]*Y[jhe4]*rate_eval.n_p_He4__t_He3
       +5.00000000000000e-01*rho**3*Y[jn]*Y[jhe4]**2*rate_eval.n_p_He4_He4__t_Be7
       )

    jac[jt, jd] = (
       -rho*Y[jt]*rate_eval.d_t__n_He4
       +rho*Y[jn]*rate_eval.n_d__t
       +5.00000000000000e-01*rho*2*Y[jd]*rate_eval.d_d__p_t
       +rho*Y[jhe4]*rate_eval.d_He4__t_He3
       )

    jac[jt, jt] = (
       -rate_eval.t__He3__weak__wc12
       -rate_eval.t__n_d
       -rho*Y[jp]*rate_eval.p_t__He4
       -rho*Y[jhe4]*rate_eval.He4_t__Li7
       -rho*Y[jp]*rate_eval.p_t__n_He3
       -rho*Y[jp]*rate_eval.p_t__d_d
       -rho*Y[jd]*rate_eval.d_t__n_He4
       -rho*Y[jhe4]*rate_eval.He4_t__n_Li6
       -rho*Y[jhe3]*rate_eval.t_He3__d_He4
       -rho*Y[jli7]*rate_eval.t_Li7__n_Be9
       -2*5.00000000000000e-01*rho*2*Y[jt]*rate_eval.t_t__n_n_He4
       -rho*Y[jhe3]*rate_eval.t_He3__n_p_He4
       -rho*Y[jli7]*rate_eval.t_Li7__n_n_He4_He4
       -rho*Y[jbe7]*rate_eval.t_Be7__n_p_He4_He4
       )

    jac[jt, jhe3] = (
       -rho*Y[jt]*rate_eval.t_He3__d_He4
       -rho*Y[jt]*rate_eval.t_He3__n_p_He4
       +rho*ye(Y)*rate_eval.He3__t__weak__electron_capture
       +rho*Y[jn]*rate_eval.n_He3__p_t
       )

    jac[jt, jhe4] = (
       -rho*Y[jt]*rate_eval.He4_t__Li7
       -rho*Y[jt]*rate_eval.He4_t__n_Li6
       +rate_eval.He4__p_t
       +rho*Y[jn]*rate_eval.n_He4__d_t
       +rho*Y[jd]*rate_eval.d_He4__t_He3
       +2*5.00000000000000e-01*rho**2*Y[jn]**2*rate_eval.n_n_He4__t_t
       +rho**2*Y[jn]*Y[jp]*rate_eval.n_p_He4__t_He3
       +2.50000000000000e-01*rho**3*Y[jn]**2*2*Y[jhe4]*rate_eval.n_n_He4_He4__t_Li7
       +5.00000000000000e-01*rho**3*Y[jn]*Y[jp]*2*Y[jhe4]*rate_eval.n_p_He4_He4__t_Be7
       )

    jac[jt, jli6] = (
       +rho*Y[jn]*rate_eval.n_Li6__He4_t
       )

    jac[jt, jli7] = (
       -rho*Y[jt]*rate_eval.t_Li7__n_Be9
       -rho*Y[jt]*rate_eval.t_Li7__n_n_He4_He4
       +rate_eval.Li7__He4_t
       )

    jac[jt, jbe7] = (
       -rho*Y[jt]*rate_eval.t_Be7__n_p_He4_He4
       )

    jac[jt, jbe9] = (
       +rho*Y[jn]*rate_eval.n_Be9__t_Li7
       )

    jac[jhe3, jn] = (
       -rho*Y[jhe3]*rate_eval.n_He3__He4
       -rho*Y[jhe3]*rate_eval.n_He3__p_t
       -rho*Y[jhe3]*rate_eval.n_He3__d_d
       +rho**2*Y[jp]*Y[jhe4]*rate_eval.n_p_He4__t_He3
       +5.00000000000000e-01*rho**3*Y[jp]*Y[jhe4]**2*rate_eval.n_p_He4_He4__He3_Li7
       )

    jac[jhe3, jp] = (
       -rho*Y[jhe3]*rate_eval.p_He3__He4__weak__bet_pos_
       +rho*Y[jd]*rate_eval.p_d__He3
       +rho*Y[jt]*rate_eval.p_t__n_He3
       +rho*Y[jhe4]*rate_eval.p_He4__d_He3
       +rho*Y[jli6]*rate_eval.p_Li6__He4_He3
       +rho**2*Y[jn]*Y[jhe4]*rate_eval.n_p_He4__t_He3
       +2*5.00000000000000e-01*rho**2*2*Y[jp]*Y[jhe4]*rate_eval.p_p_He4__He3_He3
       +5.00000000000000e-01*rho**3*Y[jn]*Y[jhe4]**2*rate_eval.n_p_He4_He4__He3_Li7
       +2.50000000000000e-01*rho**3*2*Y[jp]*Y[jhe4]**2*rate_eval.p_p_He4_He4__He3_Be7
       )

    jac[jhe3, jd] = (
       -rho*Y[jhe3]*rate_eval.d_He3__p_He4
       +rho*Y[jp]*rate_eval.p_d__He3
       +5.00000000000000e-01*rho*2*Y[jd]*rate_eval.d_d__n_He3
       +rho*Y[jhe4]*rate_eval.d_He4__t_He3
       )

    jac[jhe3, jt] = (
       -rho*Y[jhe3]*rate_eval.t_He3__d_He4
       -rho*Y[jhe3]*rate_eval.t_He3__n_p_He4
       +rate_eval.t__He3__weak__wc12
       +rho*Y[jp]*rate_eval.p_t__n_He3
       )

    jac[jhe3, jhe3] = (
       -rho*ye(Y)*rate_eval.He3__t__weak__electron_capture
       -rate_eval.He3__p_d
       -rho*Y[jn]*rate_eval.n_He3__He4
       -rho*Y[jp]*rate_eval.p_He3__He4__weak__bet_pos_
       -rho*Y[jhe4]*rate_eval.He4_He3__Be7
       -rho*Y[jn]*rate_eval.n_He3__p_t
       -rho*Y[jn]*rate_eval.n_He3__d_d
       -rho*Y[jd]*rate_eval.d_He3__p_He4
       -rho*Y[jt]*rate_eval.t_He3__d_He4
       -rho*Y[jhe4]*rate_eval.He4_He3__p_Li6
       -rho*Y[jt]*rate_eval.t_He3__n_p_He4
       -2*5.00000000000000e-01*rho*2*Y[jhe3]*rate_eval.He3_He3__p_p_He4
       -rho*Y[jli7]*rate_eval.He3_Li7__n_p_He4_He4
       -rho*Y[jbe7]*rate_eval.He3_Be7__p_p_He4_He4
       )

    jac[jhe3, jhe4] = (
       -rho*Y[jhe3]*rate_eval.He4_He3__Be7
       -rho*Y[jhe3]*rate_eval.He4_He3__p_Li6
       +rate_eval.He4__n_He3
       +rho*Y[jp]*rate_eval.p_He4__d_He3
       +rho*Y[jd]*rate_eval.d_He4__t_He3
       +rho**2*Y[jn]*Y[jp]*rate_eval.n_p_He4__t_He3
       +2*5.00000000000000e-01*rho**2*Y[jp]**2*rate_eval.p_p_He4__He3_He3
       +5.00000000000000e-01*rho**3*Y[jn]*Y[jp]*2*Y[jhe4]*rate_eval.n_p_He4_He4__He3_Li7
       +2.50000000000000e-01*rho**3*Y[jp]**2*2*Y[jhe4]*rate_eval.p_p_He4_He4__He3_Be7
       )

    jac[jhe3, jli6] = (
       +rho*Y[jp]*rate_eval.p_Li6__He4_He3
       )

    jac[jhe3, jli7] = (
       -rho*Y[jhe3]*rate_eval.He3_Li7__n_p_He4_He4
       )

    jac[jhe3, jbe7] = (
       -rho*Y[jhe3]*rate_eval.He3_Be7__p_p_He4_He4
       +rate_eval.Be7__He4_He3
       )

    jac[jhe4, jn] = (
       -rho*Y[jhe4]*rate_eval.n_He4__d_t
       -5.00000000000000e-01*rho**2*2*Y[jn]*Y[jhe4]*rate_eval.n_n_He4__He6
       -rho**2*Y[jp]*Y[jhe4]*rate_eval.n_p_He4__Li6
       -2*5.00000000000000e-01*rho**2*Y[jhe4]**2*rate_eval.n_He4_He4__Be9
       -5.00000000000000e-01*rho**2*2*Y[jn]*Y[jhe4]*rate_eval.n_n_He4__t_t
       -rho**2*Y[jp]*Y[jhe4]*rate_eval.n_p_He4__t_He3
       -2*5.00000000000000e-01*rho**2*Y[jhe4]**2*rate_eval.n_He4_He4__d_Li7
       -2*2.50000000000000e-01*rho**3*2*Y[jn]*Y[jhe4]**2*rate_eval.n_n_He4_He4__t_Li7
       -2*5.00000000000000e-01*rho**3*Y[jp]*Y[jhe4]**2*rate_eval.n_p_He4_He4__He3_Li7
       -2*5.00000000000000e-01*rho**3*Y[jp]*Y[jhe4]**2*rate_eval.n_p_He4_He4__t_Be7
       -2*5.00000000000000e-01*rho**3*Y[jp]*Y[jhe4]**2*rate_eval.n_p_He4_He4__p_Be9
       +rho*Y[jhe3]*rate_eval.n_He3__He4
       +rho*Y[jli6]*rate_eval.n_Li6__He4_t
       +2*rho*Y[jbe7]*rate_eval.n_Be7__He4_He4
       +rho*Y[jc12]*rate_eval.n_C12__He4_Be9
       +rho*Y[jo14]*rate_eval.n_O14__He4_C11
       +rho*Y[jo15]*rate_eval.n_O15__He4_C12
       +rho*Y[jo16]*rate_eval.n_O16__He4_C13
       +rho*Y[jo17]*rate_eval.n_O17__He4_C14
       +rho*Y[jf16]*rate_eval.n_F16__He4_N13
       +rho*Y[jf17]*rate_eval.n_F17__He4_N14
       +rho*Y[jf18]*rate_eval.n_F18__He4_N15
       +rho*Y[jne18]*rate_eval.n_Ne18__He4_O15
       +rho*Y[jne19]*rate_eval.n_Ne19__He4_O16
       +rho*Y[jne20]*rate_eval.n_Ne20__He4_O17
       +rho*Y[jne21]*rate_eval.n_Ne21__He4_O18
       +rho*Y[jna21]*rate_eval.n_Na21__He4_F18
       +rho*Y[jna22]*rate_eval.n_Na22__He4_F19
       +rho*Y[jmg22]*rate_eval.n_Mg22__He4_Ne19
       +rho*Y[jmg23]*rate_eval.n_Mg23__He4_Ne20
       +rho*Y[jmg24]*rate_eval.n_Mg24__He4_Ne21
       +rho*Y[jal25]*rate_eval.n_Al25__He4_Na22
       +rho*Y[jal26]*rate_eval.n_Al26__He4_Na23
       +rho*Y[jal27]*rate_eval.n_Al27__He4_Na24
       +rho*Y[jsi26]*rate_eval.n_Si26__He4_Mg23
       +rho*Y[jsi28]*rate_eval.n_Si28__He4_Mg25
       +rho*Y[jsi29]*rate_eval.n_Si29__He4_Mg26
       +rho*Y[jp29]*rate_eval.n_P29__He4_Al26
       +rho*Y[jp30]*rate_eval.n_P30__He4_Al27
       +rho*Y[jp31]*rate_eval.n_P31__He4_Al28
       +rho*Y[js31]*rate_eval.n_S31__He4_Si28
       +rho*Y[js32]*rate_eval.n_S32__He4_Si29
       +rho*Y[jcl32]*rate_eval.n_Cl32__He4_P29
       +rho*Y[jcl34]*rate_eval.n_Cl34__He4_P31
       +rho*Y[jcl35]*rate_eval.n_Cl35__He4_P32
       +rho*Y[jar35]*rate_eval.n_Ar35__He4_S32
       +2*rho*Y[jb8]*rate_eval.n_B8__p_He4_He4
       +3*rho*Y[jc11]*rate_eval.n_C11__He4_He4_He4
       )

    jac[jhe4, jp] = (
       -rho*Y[jhe4]*rate_eval.p_He4__d_He3
       -rho**2*Y[jn]*Y[jhe4]*rate_eval.n_p_He4__Li6
       -rho**2*Y[jn]*Y[jhe4]*rate_eval.n_p_He4__t_He3
       -5.00000000000000e-01*rho**2*2*Y[jp]*Y[jhe4]*rate_eval.p_p_He4__He3_He3
       -2*5.00000000000000e-01*rho**2*Y[jhe4]**2*rate_eval.p_He4_He4__n_B8
       -2*5.00000000000000e-01*rho**2*Y[jhe4]**2*rate_eval.p_He4_He4__d_Be7
       -2*5.00000000000000e-01*rho**3*Y[jn]*Y[jhe4]**2*rate_eval.n_p_He4_He4__He3_Li7
       -2*5.00000000000000e-01*rho**3*Y[jn]*Y[jhe4]**2*rate_eval.n_p_He4_He4__t_Be7
       -2*5.00000000000000e-01*rho**3*Y[jn]*Y[jhe4]**2*rate_eval.n_p_He4_He4__p_Be9
       -2*2.50000000000000e-01*rho**3*2*Y[jp]*Y[jhe4]**2*rate_eval.p_p_He4_He4__He3_Be7
       +rho*Y[jt]*rate_eval.p_t__He4
       +rho*Y[jhe3]*rate_eval.p_He3__He4__weak__bet_pos_
       +rho*Y[jli6]*rate_eval.p_Li6__He4_He3
       +2*rho*Y[jli7]*rate_eval.p_Li7__He4_He4
       +rho*Y[jbe9]*rate_eval.p_Be9__He4_Li6
       +rho*Y[jc11]*rate_eval.p_C11__He4_B8
       +rho*Y[jn14]*rate_eval.p_N14__He4_C11
       +rho*Y[jn15]*rate_eval.p_N15__He4_C12
       +rho*Y[jo16]*rate_eval.p_O16__He4_N13
       +rho*Y[jo17]*rate_eval.p_O17__He4_N14
       +rho*Y[jo18]*rate_eval.p_O18__He4_N15
       +rho*Y[jf17]*rate_eval.p_F17__He4_O14
       +rho*Y[jf18]*rate_eval.p_F18__He4_O15
       +rho*Y[jf19]*rate_eval.p_F19__He4_O16
       +rho*Y[jne18]*rate_eval.p_Ne18__He4_F15
       +rho*Y[jne19]*rate_eval.p_Ne19__He4_F16
       +rho*Y[jne20]*rate_eval.p_Ne20__He4_F17
       +rho*Y[jne21]*rate_eval.p_Ne21__He4_F18
       +rho*Y[jna21]*rate_eval.p_Na21__He4_Ne18
       +rho*Y[jna22]*rate_eval.p_Na22__He4_Ne19
       +rho*Y[jna23]*rate_eval.p_Na23__He4_Ne20
       +rho*Y[jna24]*rate_eval.p_Na24__He4_Ne21
       +rho*Y[jmg24]*rate_eval.p_Mg24__He4_Na21
       +rho*Y[jmg25]*rate_eval.p_Mg25__He4_Na22
       +rho*Y[jmg26]*rate_eval.p_Mg26__He4_Na23
       +rho*Y[jal25]*rate_eval.p_Al25__He4_Mg22
       +rho*Y[jal26]*rate_eval.p_Al26__He4_Mg23
       +rho*Y[jal27]*rate_eval.p_Al27__He4_Mg24
       +rho*Y[jal28]*rate_eval.p_Al28__He4_Mg25
       +rho*Y[jsi28]*rate_eval.p_Si28__He4_Al25
       +rho*Y[jsi29]*rate_eval.p_Si29__He4_Al26
       +rho*Y[jsi30]*rate_eval.p_Si30__He4_Al27
       +rho*Y[jp29]*rate_eval.p_P29__He4_Si26
       +rho*Y[jp31]*rate_eval.p_P31__He4_Si28
       +rho*Y[jp32]*rate_eval.p_P32__He4_Si29
       +rho*Y[js32]*rate_eval.p_S32__He4_P29
       +rho*Y[js34]*rate_eval.p_S34__He4_P31
       +rho*Y[jcl34]*rate_eval.p_Cl34__He4_S31
       +rho*Y[jcl35]*rate_eval.p_Cl35__He4_S32
       +rho*Y[jar35]*rate_eval.p_Ar35__He4_Cl32
       +2*rho*Y[jbe9]*rate_eval.p_Be9__d_He4_He4
       +2*rho*Y[jbe9]*rate_eval.p_Be9__n_p_He4_He4
       )

    jac[jhe4, jd] = (
       -rho*Y[jhe4]*rate_eval.He4_d__Li6
       -rho*Y[jhe4]*rate_eval.d_He4__t_He3
       -2*5.00000000000000e-01*rho**2*Y[jhe4]**2*rate_eval.d_He4_He4__p_Be9
       +5.00000000000000e-01*rho*2*Y[jd]*rate_eval.d_d__He4
       +rho*Y[jt]*rate_eval.d_t__n_He4
       +rho*Y[jhe3]*rate_eval.d_He3__p_He4
       +2*rho*Y[jli7]*rate_eval.d_Li7__n_He4_He4
       +2*rho*Y[jbe7]*rate_eval.d_Be7__p_He4_He4
       )

    jac[jhe4, jt] = (
       -rho*Y[jhe4]*rate_eval.He4_t__Li7
       -rho*Y[jhe4]*rate_eval.He4_t__n_Li6
       +rho*Y[jp]*rate_eval.p_t__He4
       +rho*Y[jd]*rate_eval.d_t__n_He4
       +rho*Y[jhe3]*rate_eval.t_He3__d_He4
       +5.00000000000000e-01*rho*2*Y[jt]*rate_eval.t_t__n_n_He4
       +rho*Y[jhe3]*rate_eval.t_He3__n_p_He4
       +2*rho*Y[jli7]*rate_eval.t_Li7__n_n_He4_He4
       +2*rho*Y[jbe7]*rate_eval.t_Be7__n_p_He4_He4
       )

    jac[jhe4, jhe3] = (
       -rho*Y[jhe4]*rate_eval.He4_He3__Be7
       -rho*Y[jhe4]*rate_eval.He4_He3__p_Li6
       +rho*Y[jn]*rate_eval.n_He3__He4
       +rho*Y[jp]*rate_eval.p_He3__He4__weak__bet_pos_
       +rho*Y[jd]*rate_eval.d_He3__p_He4
       +rho*Y[jt]*rate_eval.t_He3__d_He4
       +rho*Y[jt]*rate_eval.t_He3__n_p_He4
       +5.00000000000000e-01*rho*2*Y[jhe3]*rate_eval.He3_He3__p_p_He4
       +2*rho*Y[jli7]*rate_eval.He3_Li7__n_p_He4_He4
       +2*rho*Y[jbe7]*rate_eval.He3_Be7__p_p_He4_He4
       )

    jac[jhe4, jhe4] = (
       -rate_eval.He4__n_He3
       -rate_eval.He4__p_t
       -rate_eval.He4__d_d
       -rho*Y[jd]*rate_eval.He4_d__Li6
       -rho*Y[jt]*rate_eval.He4_t__Li7
       -rho*Y[jhe3]*rate_eval.He4_He3__Be7
       -rho*Y[jbe7]*rate_eval.He4_Be7__C11
       -rho*Y[jc12]*rate_eval.He4_C12__O16
       -rho*Y[jc14]*rate_eval.He4_C14__O18
       -rho*Y[jn14]*rate_eval.He4_N14__F18
       -rho*Y[jn15]*rate_eval.He4_N15__F19
       -rho*Y[jo14]*rate_eval.He4_O14__Ne18
       -rho*Y[jo15]*rate_eval.He4_O15__Ne19
       -rho*Y[jo16]*rate_eval.He4_O16__Ne20
       -rho*Y[jo17]*rate_eval.He4_O17__Ne21
       -rho*Y[jf17]*rate_eval.He4_F17__Na21
       -rho*Y[jf18]*rate_eval.He4_F18__Na22
       -rho*Y[jf19]*rate_eval.He4_F19__Na23
       -rho*Y[jne18]*rate_eval.He4_Ne18__Mg22
       -rho*Y[jne19]*rate_eval.He4_Ne19__Mg23
       -rho*Y[jne20]*rate_eval.He4_Ne20__Mg24
       -rho*Y[jne21]*rate_eval.He4_Ne21__Mg25
       -rho*Y[jna21]*rate_eval.He4_Na21__Al25
       -rho*Y[jna22]*rate_eval.He4_Na22__Al26
       -rho*Y[jna23]*rate_eval.He4_Na23__Al27
       -rho*Y[jna24]*rate_eval.He4_Na24__Al28
       -rho*Y[jmg22]*rate_eval.He4_Mg22__Si26
       -rho*Y[jmg24]*rate_eval.He4_Mg24__Si28
       -rho*Y[jmg25]*rate_eval.He4_Mg25__Si29
       -rho*Y[jmg26]*rate_eval.He4_Mg26__Si30
       -rho*Y[jal25]*rate_eval.He4_Al25__P29
       -rho*Y[jal26]*rate_eval.He4_Al26__P30
       -rho*Y[jal27]*rate_eval.He4_Al27__P31
       -rho*Y[jal28]*rate_eval.He4_Al28__P32
       -rho*Y[jsi28]*rate_eval.He4_Si28__S32
       -rho*Y[jsi30]*rate_eval.He4_Si30__S34
       -rho*Y[jp30]*rate_eval.He4_P30__Cl34
       -rho*Y[jp31]*rate_eval.He4_P31__Cl35
       -rho*Y[js31]*rate_eval.He4_S31__Ar35
       -rho*Y[jt]*rate_eval.He4_t__n_Li6
       -rho*Y[jhe3]*rate_eval.He4_He3__p_Li6
       -rho*Y[jn]*rate_eval.n_He4__d_t
       -rho*Y[jp]*rate_eval.p_He4__d_He3
       -rho*Y[jd]*rate_eval.d_He4__t_He3
       -2*5.00000000000000e-01*rho*2*Y[jhe4]*rate_eval.He4_He4__n_Be7
       -2*5.00000000000000e-01*rho*2*Y[jhe4]*rate_eval.He4_He4__p_Li7
       -rho*Y[jli6]*rate_eval.He4_Li6__p_Be9
       -rho*Y[jbe9]*rate_eval.He4_Be9__n_C12
       -rho*Y[jb8]*rate_eval.He4_B8__p_C11
       -rho*Y[jc11]*rate_eval.He4_C11__n_O14
       -rho*Y[jc11]*rate_eval.He4_C11__p_N14
       -rho*Y[jc12]*rate_eval.He4_C12__n_O15
       -rho*Y[jc12]*rate_eval.He4_C12__p_N15
       -rho*Y[jc13]*rate_eval.He4_C13__n_O16
       -rho*Y[jc14]*rate_eval.He4_C14__n_O17
       -rho*Y[jn13]*rate_eval.He4_N13__n_F16
       -rho*Y[jn13]*rate_eval.He4_N13__p_O16
       -rho*Y[jn14]*rate_eval.He4_N14__n_F17
       -rho*Y[jn14]*rate_eval.He4_N14__p_O17
       -rho*Y[jn15]*rate_eval.He4_N15__n_F18
       -rho*Y[jn15]*rate_eval.He4_N15__p_O18
       -rho*Y[jo14]*rate_eval.He4_O14__p_F17
       -rho*Y[jo15]*rate_eval.He4_O15__n_Ne18
       -rho*Y[jo15]*rate_eval.He4_O15__p_F18
       -rho*Y[jo16]*rate_eval.He4_O16__n_Ne19
       -rho*Y[jo16]*rate_eval.He4_O16__p_F19
       -rho*Y[jo17]*rate_eval.He4_O17__n_Ne20
       -rho*Y[jo18]*rate_eval.He4_O18__n_Ne21
       -rho*Y[jf15]*rate_eval.He4_F15__p_Ne18
       -rho*Y[jf16]*rate_eval.He4_F16__p_Ne19
       -rho*Y[jf17]*rate_eval.He4_F17__p_Ne20
       -rho*Y[jf18]*rate_eval.He4_F18__n_Na21
       -rho*Y[jf18]*rate_eval.He4_F18__p_Ne21
       -rho*Y[jf19]*rate_eval.He4_F19__n_Na22
       -rho*Y[jne18]*rate_eval.He4_Ne18__p_Na21
       -rho*Y[jne19]*rate_eval.He4_Ne19__n_Mg22
       -rho*Y[jne19]*rate_eval.He4_Ne19__p_Na22
       -rho*Y[jne20]*rate_eval.He4_Ne20__n_Mg23
       -rho*Y[jne20]*rate_eval.He4_Ne20__p_Na23
       -rho*Y[jne20]*rate_eval.He4_Ne20__C12_C12
       -rho*Y[jne21]*rate_eval.He4_Ne21__n_Mg24
       -rho*Y[jne21]*rate_eval.He4_Ne21__p_Na24
       -rho*Y[jna21]*rate_eval.He4_Na21__p_Mg24
       -rho*Y[jna22]*rate_eval.He4_Na22__n_Al25
       -rho*Y[jna22]*rate_eval.He4_Na22__p_Mg25
       -rho*Y[jna23]*rate_eval.He4_Na23__n_Al26
       -rho*Y[jna23]*rate_eval.He4_Na23__p_Mg26
       -rho*Y[jna24]*rate_eval.He4_Na24__n_Al27
       -rho*Y[jmg22]*rate_eval.He4_Mg22__p_Al25
       -rho*Y[jmg23]*rate_eval.He4_Mg23__n_Si26
       -rho*Y[jmg23]*rate_eval.He4_Mg23__p_Al26
       -rho*Y[jmg24]*rate_eval.He4_Mg24__p_Al27
       -rho*Y[jmg24]*rate_eval.He4_Mg24__C12_O16
       -rho*Y[jmg25]*rate_eval.He4_Mg25__n_Si28
       -rho*Y[jmg25]*rate_eval.He4_Mg25__p_Al28
       -rho*Y[jmg26]*rate_eval.He4_Mg26__n_Si29
       -rho*Y[jal25]*rate_eval.He4_Al25__p_Si28
       -rho*Y[jal26]*rate_eval.He4_Al26__n_P29
       -rho*Y[jal26]*rate_eval.He4_Al26__p_Si29
       -rho*Y[jal27]*rate_eval.He4_Al27__n_P30
       -rho*Y[jal27]*rate_eval.He4_Al27__p_Si30
       -rho*Y[jal28]*rate_eval.He4_Al28__n_P31
       -rho*Y[jsi26]*rate_eval.He4_Si26__p_P29
       -rho*Y[jsi28]*rate_eval.He4_Si28__n_S31
       -rho*Y[jsi28]*rate_eval.He4_Si28__p_P31
       -rho*Y[jsi28]*rate_eval.He4_Si28__C12_Ne20
       -rho*Y[jsi28]*rate_eval.He4_Si28__O16_O16
       -rho*Y[jsi29]*rate_eval.He4_Si29__n_S32
       -rho*Y[jsi29]*rate_eval.He4_Si29__p_P32
       -rho*Y[jp29]*rate_eval.He4_P29__n_Cl32
       -rho*Y[jp29]*rate_eval.He4_P29__p_S32
       -rho*Y[jp31]*rate_eval.He4_P31__n_Cl34
       -rho*Y[jp31]*rate_eval.He4_P31__p_S34
       -rho*Y[jp32]*rate_eval.He4_P32__n_Cl35
       -rho*Y[js31]*rate_eval.He4_S31__p_Cl34
       -rho*Y[js32]*rate_eval.He4_S32__n_Ar35
       -rho*Y[js32]*rate_eval.He4_S32__p_Cl35
       -rho*Y[jcl32]*rate_eval.He4_Cl32__p_Ar35
       -5.00000000000000e-01*rho**2*Y[jn]**2*rate_eval.n_n_He4__He6
       -rho**2*Y[jn]*Y[jp]*rate_eval.n_p_He4__Li6
       -2*5.00000000000000e-01*rho**2*Y[jn]*2*Y[jhe4]*rate_eval.n_He4_He4__Be9
       -3*1.66666666666667e-01*rho**2*3*Y[jhe4]**2*rate_eval.He4_He4_He4__C12
       -5.00000000000000e-01*rho**2*Y[jn]**2*rate_eval.n_n_He4__t_t
       -rho**2*Y[jn]*Y[jp]*rate_eval.n_p_He4__t_He3
       -5.00000000000000e-01*rho**2*Y[jp]**2*rate_eval.p_p_He4__He3_He3
       -2*5.00000000000000e-01*rho**2*Y[jn]*2*Y[jhe4]*rate_eval.n_He4_He4__d_Li7
       -2*5.00000000000000e-01*rho**2*Y[jp]*2*Y[jhe4]*rate_eval.p_He4_He4__n_B8
       -2*5.00000000000000e-01*rho**2*Y[jp]*2*Y[jhe4]*rate_eval.p_He4_He4__d_Be7
       -2*5.00000000000000e-01*rho**2*Y[jd]*2*Y[jhe4]*rate_eval.d_He4_He4__p_Be9
       -3*1.66666666666667e-01*rho**2*3*Y[jhe4]**2*rate_eval.He4_He4_He4__n_C11
       -2*2.50000000000000e-01*rho**3*Y[jn]**2*2*Y[jhe4]*rate_eval.n_n_He4_He4__t_Li7
       -2*5.00000000000000e-01*rho**3*Y[jn]*Y[jp]*2*Y[jhe4]*rate_eval.n_p_He4_He4__He3_Li7
       -2*5.00000000000000e-01*rho**3*Y[jn]*Y[jp]*2*Y[jhe4]*rate_eval.n_p_He4_He4__t_Be7
       -2*5.00000000000000e-01*rho**3*Y[jn]*Y[jp]*2*Y[jhe4]*rate_eval.n_p_He4_He4__p_Be9
       -2*2.50000000000000e-01*rho**3*Y[jp]**2*2*Y[jhe4]*rate_eval.p_p_He4_He4__He3_Be7
       )

    jac[jhe4, jhe6] = (
       +rate_eval.He6__n_n_He4
       )

    jac[jhe4, jli6] = (
       -rho*Y[jhe4]*rate_eval.He4_Li6__p_Be9
       +rate_eval.Li6__He4_d
       +rate_eval.Li6__n_p_He4
       +rho*Y[jn]*rate_eval.n_Li6__He4_t
       +rho*Y[jp]*rate_eval.p_Li6__He4_He3
       )

    jac[jhe4, jli7] = (
       +rate_eval.Li7__He4_t
       +2*rho*Y[jp]*rate_eval.p_Li7__He4_He4
       +2*rho*Y[jd]*rate_eval.d_Li7__n_He4_He4
       +2*rho*Y[jt]*rate_eval.t_Li7__n_n_He4_He4
       +2*rho*Y[jhe3]*rate_eval.He3_Li7__n_p_He4_He4
       )

    jac[jhe4, jbe7] = (
       -rho*Y[jhe4]*rate_eval.He4_Be7__C11
       +rate_eval.Be7__He4_He3
       +2*rho*Y[jn]*rate_eval.n_Be7__He4_He4
       +2*rho*Y[jd]*rate_eval.d_Be7__p_He4_He4
       +2*rho*Y[jt]*rate_eval.t_Be7__n_p_He4_He4
       +2*rho*Y[jhe3]*rate_eval.He3_Be7__p_p_He4_He4
       )

    jac[jhe4, jbe9] = (
       -rho*Y[jhe4]*rate_eval.He4_Be9__n_C12
       +2*rate_eval.Be9__n_He4_He4
       +rho*Y[jp]*rate_eval.p_Be9__He4_Li6
       +2*rho*Y[jp]*rate_eval.p_Be9__d_He4_He4
       +2*rho*Y[jp]*rate_eval.p_Be9__n_p_He4_He4
       )

    jac[jhe4, jb8] = (
       -rho*Y[jhe4]*rate_eval.He4_B8__p_C11
       +2*rate_eval.B8__He4_He4__weak__wc12
       +2*rho*Y[jn]*rate_eval.n_B8__p_He4_He4
       )

    jac[jhe4, jc11] = (
       -rho*Y[jhe4]*rate_eval.He4_C11__n_O14
       -rho*Y[jhe4]*rate_eval.He4_C11__p_N14
       +rate_eval.C11__He4_Be7
       +rho*Y[jp]*rate_eval.p_C11__He4_B8
       +3*rho*Y[jn]*rate_eval.n_C11__He4_He4_He4
       )

    jac[jhe4, jc12] = (
       -rho*Y[jhe4]*rate_eval.He4_C12__O16
       -rho*Y[jhe4]*rate_eval.He4_C12__n_O15
       -rho*Y[jhe4]*rate_eval.He4_C12__p_N15
       +3*rate_eval.C12__He4_He4_He4
       +rho*Y[jn]*rate_eval.n_C12__He4_Be9
       +5.00000000000000e-01*rho*2*Y[jc12]*rate_eval.C12_C12__He4_Ne20
       +rho*Y[jo16]*rate_eval.C12_O16__He4_Mg24
       +rho*Y[jne20]*rate_eval.C12_Ne20__He4_Si28
       )

    jac[jhe4, jc13] = (
       -rho*Y[jhe4]*rate_eval.He4_C13__n_O16
       )

    jac[jhe4, jc14] = (
       -rho*Y[jhe4]*rate_eval.He4_C14__O18
       -rho*Y[jhe4]*rate_eval.He4_C14__n_O17
       )

    jac[jhe4, jn13] = (
       -rho*Y[jhe4]*rate_eval.He4_N13__n_F16
       -rho*Y[jhe4]*rate_eval.He4_N13__p_O16
       )

    jac[jhe4, jn14] = (
       -rho*Y[jhe4]*rate_eval.He4_N14__F18
       -rho*Y[jhe4]*rate_eval.He4_N14__n_F17
       -rho*Y[jhe4]*rate_eval.He4_N14__p_O17
       +rho*Y[jp]*rate_eval.p_N14__He4_C11
       )

    jac[jhe4, jn15] = (
       -rho*Y[jhe4]*rate_eval.He4_N15__F19
       -rho*Y[jhe4]*rate_eval.He4_N15__n_F18
       -rho*Y[jhe4]*rate_eval.He4_N15__p_O18
       +rho*Y[jp]*rate_eval.p_N15__He4_C12
       )

    jac[jhe4, jo14] = (
       -rho*Y[jhe4]*rate_eval.He4_O14__Ne18
       -rho*Y[jhe4]*rate_eval.He4_O14__p_F17
       +rho*Y[jn]*rate_eval.n_O14__He4_C11
       )

    jac[jhe4, jo15] = (
       -rho*Y[jhe4]*rate_eval.He4_O15__Ne19
       -rho*Y[jhe4]*rate_eval.He4_O15__n_Ne18
       -rho*Y[jhe4]*rate_eval.He4_O15__p_F18
       +rho*Y[jn]*rate_eval.n_O15__He4_C12
       )

    jac[jhe4, jo16] = (
       -rho*Y[jhe4]*rate_eval.He4_O16__Ne20
       -rho*Y[jhe4]*rate_eval.He4_O16__n_Ne19
       -rho*Y[jhe4]*rate_eval.He4_O16__p_F19
       +rate_eval.O16__He4_C12
       +rho*Y[jn]*rate_eval.n_O16__He4_C13
       +rho*Y[jp]*rate_eval.p_O16__He4_N13
       +rho*Y[jc12]*rate_eval.C12_O16__He4_Mg24
       +5.00000000000000e-01*rho*2*Y[jo16]*rate_eval.O16_O16__He4_Si28
       )

    jac[jhe4, jo17] = (
       -rho*Y[jhe4]*rate_eval.He4_O17__Ne21
       -rho*Y[jhe4]*rate_eval.He4_O17__n_Ne20
       +rho*Y[jn]*rate_eval.n_O17__He4_C14
       +rho*Y[jp]*rate_eval.p_O17__He4_N14
       )

    jac[jhe4, jo18] = (
       -rho*Y[jhe4]*rate_eval.He4_O18__n_Ne21
       +rate_eval.O18__He4_C14
       +rho*Y[jp]*rate_eval.p_O18__He4_N15
       )

    jac[jhe4, jf15] = (
       -rho*Y[jhe4]*rate_eval.He4_F15__p_Ne18
       )

    jac[jhe4, jf16] = (
       -rho*Y[jhe4]*rate_eval.He4_F16__p_Ne19
       +rho*Y[jn]*rate_eval.n_F16__He4_N13
       )

    jac[jhe4, jf17] = (
       -rho*Y[jhe4]*rate_eval.He4_F17__Na21
       -rho*Y[jhe4]*rate_eval.He4_F17__p_Ne20
       +rho*Y[jn]*rate_eval.n_F17__He4_N14
       +rho*Y[jp]*rate_eval.p_F17__He4_O14
       )

    jac[jhe4, jf18] = (
       -rho*Y[jhe4]*rate_eval.He4_F18__Na22
       -rho*Y[jhe4]*rate_eval.He4_F18__n_Na21
       -rho*Y[jhe4]*rate_eval.He4_F18__p_Ne21
       +rate_eval.F18__He4_N14
       +rho*Y[jn]*rate_eval.n_F18__He4_N15
       +rho*Y[jp]*rate_eval.p_F18__He4_O15
       )

    jac[jhe4, jf19] = (
       -rho*Y[jhe4]*rate_eval.He4_F19__Na23
       -rho*Y[jhe4]*rate_eval.He4_F19__n_Na22
       +rate_eval.F19__He4_N15
       +rho*Y[jp]*rate_eval.p_F19__He4_O16
       )

    jac[jhe4, jne18] = (
       -rho*Y[jhe4]*rate_eval.He4_Ne18__Mg22
       -rho*Y[jhe4]*rate_eval.He4_Ne18__p_Na21
       +rate_eval.Ne18__He4_O14
       +rho*Y[jn]*rate_eval.n_Ne18__He4_O15
       +rho*Y[jp]*rate_eval.p_Ne18__He4_F15
       )

    jac[jhe4, jne19] = (
       -rho*Y[jhe4]*rate_eval.He4_Ne19__Mg23
       -rho*Y[jhe4]*rate_eval.He4_Ne19__n_Mg22
       -rho*Y[jhe4]*rate_eval.He4_Ne19__p_Na22
       +rate_eval.Ne19__He4_O15
       +rho*Y[jn]*rate_eval.n_Ne19__He4_O16
       +rho*Y[jp]*rate_eval.p_Ne19__He4_F16
       )

    jac[jhe4, jne20] = (
       -rho*Y[jhe4]*rate_eval.He4_Ne20__Mg24
       -rho*Y[jhe4]*rate_eval.He4_Ne20__n_Mg23
       -rho*Y[jhe4]*rate_eval.He4_Ne20__p_Na23
       -rho*Y[jhe4]*rate_eval.He4_Ne20__C12_C12
       +rate_eval.Ne20__He4_O16
       +rho*Y[jn]*rate_eval.n_Ne20__He4_O17
       +rho*Y[jp]*rate_eval.p_Ne20__He4_F17
       +rho*Y[jc12]*rate_eval.C12_Ne20__He4_Si28
       )

    jac[jhe4, jne21] = (
       -rho*Y[jhe4]*rate_eval.He4_Ne21__Mg25
       -rho*Y[jhe4]*rate_eval.He4_Ne21__n_Mg24
       -rho*Y[jhe4]*rate_eval.He4_Ne21__p_Na24
       +rate_eval.Ne21__He4_O17
       +rho*Y[jn]*rate_eval.n_Ne21__He4_O18
       +rho*Y[jp]*rate_eval.p_Ne21__He4_F18
       )

    jac[jhe4, jna21] = (
       -rho*Y[jhe4]*rate_eval.He4_Na21__Al25
       -rho*Y[jhe4]*rate_eval.He4_Na21__p_Mg24
       +rate_eval.Na21__He4_F17
       +rho*Y[jn]*rate_eval.n_Na21__He4_F18
       +rho*Y[jp]*rate_eval.p_Na21__He4_Ne18
       )

    jac[jhe4, jna22] = (
       -rho*Y[jhe4]*rate_eval.He4_Na22__Al26
       -rho*Y[jhe4]*rate_eval.He4_Na22__n_Al25
       -rho*Y[jhe4]*rate_eval.He4_Na22__p_Mg25
       +rate_eval.Na22__He4_F18
       +rho*Y[jn]*rate_eval.n_Na22__He4_F19
       +rho*Y[jp]*rate_eval.p_Na22__He4_Ne19
       )

    jac[jhe4, jna23] = (
       -rho*Y[jhe4]*rate_eval.He4_Na23__Al27
       -rho*Y[jhe4]*rate_eval.He4_Na23__n_Al26
       -rho*Y[jhe4]*rate_eval.He4_Na23__p_Mg26
       +rate_eval.Na23__He4_F19
       +rho*Y[jp]*rate_eval.p_Na23__He4_Ne20
       )

    jac[jhe4, jna24] = (
       -rho*Y[jhe4]*rate_eval.He4_Na24__Al28
       -rho*Y[jhe4]*rate_eval.He4_Na24__n_Al27
       +rho*Y[jp]*rate_eval.p_Na24__He4_Ne21
       )

    jac[jhe4, jmg22] = (
       -rho*Y[jhe4]*rate_eval.He4_Mg22__Si26
       -rho*Y[jhe4]*rate_eval.He4_Mg22__p_Al25
       +rate_eval.Mg22__He4_Ne18
       +rho*Y[jn]*rate_eval.n_Mg22__He4_Ne19
       )

    jac[jhe4, jmg23] = (
       -rho*Y[jhe4]*rate_eval.He4_Mg23__n_Si26
       -rho*Y[jhe4]*rate_eval.He4_Mg23__p_Al26
       +rate_eval.Mg23__He4_Ne19
       +rho*Y[jn]*rate_eval.n_Mg23__He4_Ne20
       )

    jac[jhe4, jmg24] = (
       -rho*Y[jhe4]*rate_eval.He4_Mg24__Si28
       -rho*Y[jhe4]*rate_eval.He4_Mg24__p_Al27
       -rho*Y[jhe4]*rate_eval.He4_Mg24__C12_O16
       +rate_eval.Mg24__He4_Ne20
       +rho*Y[jn]*rate_eval.n_Mg24__He4_Ne21
       +rho*Y[jp]*rate_eval.p_Mg24__He4_Na21
       )

    jac[jhe4, jmg25] = (
       -rho*Y[jhe4]*rate_eval.He4_Mg25__Si29
       -rho*Y[jhe4]*rate_eval.He4_Mg25__n_Si28
       -rho*Y[jhe4]*rate_eval.He4_Mg25__p_Al28
       +rate_eval.Mg25__He4_Ne21
       +rho*Y[jp]*rate_eval.p_Mg25__He4_Na22
       )

    jac[jhe4, jmg26] = (
       -rho*Y[jhe4]*rate_eval.He4_Mg26__Si30
       -rho*Y[jhe4]*rate_eval.He4_Mg26__n_Si29
       +rho*Y[jp]*rate_eval.p_Mg26__He4_Na23
       )

    jac[jhe4, jal25] = (
       -rho*Y[jhe4]*rate_eval.He4_Al25__P29
       -rho*Y[jhe4]*rate_eval.He4_Al25__p_Si28
       +rate_eval.Al25__He4_Na21
       +rho*Y[jn]*rate_eval.n_Al25__He4_Na22
       +rho*Y[jp]*rate_eval.p_Al25__He4_Mg22
       )

    jac[jhe4, jal26] = (
       -rho*Y[jhe4]*rate_eval.He4_Al26__P30
       -rho*Y[jhe4]*rate_eval.He4_Al26__n_P29
       -rho*Y[jhe4]*rate_eval.He4_Al26__p_Si29
       +rate_eval.Al26__He4_Na22
       +rho*Y[jn]*rate_eval.n_Al26__He4_Na23
       +rho*Y[jp]*rate_eval.p_Al26__He4_Mg23
       )

    jac[jhe4, jal27] = (
       -rho*Y[jhe4]*rate_eval.He4_Al27__P31
       -rho*Y[jhe4]*rate_eval.He4_Al27__n_P30
       -rho*Y[jhe4]*rate_eval.He4_Al27__p_Si30
       +rate_eval.Al27__He4_Na23
       +rho*Y[jn]*rate_eval.n_Al27__He4_Na24
       +rho*Y[jp]*rate_eval.p_Al27__He4_Mg24
       )

    jac[jhe4, jal28] = (
       -rho*Y[jhe4]*rate_eval.He4_Al28__P32
       -rho*Y[jhe4]*rate_eval.He4_Al28__n_P31
       +rate_eval.Al28__He4_Na24
       +rho*Y[jp]*rate_eval.p_Al28__He4_Mg25
       )

    jac[jhe4, jsi26] = (
       -rho*Y[jhe4]*rate_eval.He4_Si26__p_P29
       +rate_eval.Si26__He4_Mg22
       +rho*Y[jn]*rate_eval.n_Si26__He4_Mg23
       )

    jac[jhe4, jsi28] = (
       -rho*Y[jhe4]*rate_eval.He4_Si28__S32
       -rho*Y[jhe4]*rate_eval.He4_Si28__n_S31
       -rho*Y[jhe4]*rate_eval.He4_Si28__p_P31
       -rho*Y[jhe4]*rate_eval.He4_Si28__C12_Ne20
       -rho*Y[jhe4]*rate_eval.He4_Si28__O16_O16
       +rate_eval.Si28__He4_Mg24
       +rho*Y[jn]*rate_eval.n_Si28__He4_Mg25
       +rho*Y[jp]*rate_eval.p_Si28__He4_Al25
       )

    jac[jhe4, jsi29] = (
       -rho*Y[jhe4]*rate_eval.He4_Si29__n_S32
       -rho*Y[jhe4]*rate_eval.He4_Si29__p_P32
       +rate_eval.Si29__He4_Mg25
       +rho*Y[jn]*rate_eval.n_Si29__He4_Mg26
       +rho*Y[jp]*rate_eval.p_Si29__He4_Al26
       )

    jac[jhe4, jsi30] = (
       -rho*Y[jhe4]*rate_eval.He4_Si30__S34
       +rate_eval.Si30__He4_Mg26
       +rho*Y[jp]*rate_eval.p_Si30__He4_Al27
       )

    jac[jhe4, jp29] = (
       -rho*Y[jhe4]*rate_eval.He4_P29__n_Cl32
       -rho*Y[jhe4]*rate_eval.He4_P29__p_S32
       +rate_eval.P29__He4_Al25
       +rho*Y[jn]*rate_eval.n_P29__He4_Al26
       +rho*Y[jp]*rate_eval.p_P29__He4_Si26
       )

    jac[jhe4, jp30] = (
       -rho*Y[jhe4]*rate_eval.He4_P30__Cl34
       +rate_eval.P30__He4_Al26
       +rho*Y[jn]*rate_eval.n_P30__He4_Al27
       )

    jac[jhe4, jp31] = (
       -rho*Y[jhe4]*rate_eval.He4_P31__Cl35
       -rho*Y[jhe4]*rate_eval.He4_P31__n_Cl34
       -rho*Y[jhe4]*rate_eval.He4_P31__p_S34
       +rate_eval.P31__He4_Al27
       +rho*Y[jn]*rate_eval.n_P31__He4_Al28
       +rho*Y[jp]*rate_eval.p_P31__He4_Si28
       )

    jac[jhe4, jp32] = (
       -rho*Y[jhe4]*rate_eval.He4_P32__n_Cl35
       +rate_eval.P32__He4_Al28
       +rho*Y[jp]*rate_eval.p_P32__He4_Si29
       )

    jac[jhe4, js31] = (
       -rho*Y[jhe4]*rate_eval.He4_S31__Ar35
       -rho*Y[jhe4]*rate_eval.He4_S31__p_Cl34
       +rho*Y[jn]*rate_eval.n_S31__He4_Si28
       )

    jac[jhe4, js32] = (
       -rho*Y[jhe4]*rate_eval.He4_S32__n_Ar35
       -rho*Y[jhe4]*rate_eval.He4_S32__p_Cl35
       +rate_eval.S32__He4_Si28
       +rho*Y[jn]*rate_eval.n_S32__He4_Si29
       +rho*Y[jp]*rate_eval.p_S32__He4_P29
       )

    jac[jhe4, js34] = (
       +rate_eval.S34__He4_Si30
       +rho*Y[jp]*rate_eval.p_S34__He4_P31
       )

    jac[jhe4, jcl32] = (
       -rho*Y[jhe4]*rate_eval.He4_Cl32__p_Ar35
       +rate_eval.Cl32__He4_Si28__weak__wc12
       +rho*Y[jn]*rate_eval.n_Cl32__He4_P29
       )

    jac[jhe4, jcl34] = (
       +rate_eval.Cl34__He4_P30
       +rho*Y[jn]*rate_eval.n_Cl34__He4_P31
       +rho*Y[jp]*rate_eval.p_Cl34__He4_S31
       )

    jac[jhe4, jcl35] = (
       +rate_eval.Cl35__He4_P31
       +rho*Y[jn]*rate_eval.n_Cl35__He4_P32
       +rho*Y[jp]*rate_eval.p_Cl35__He4_S32
       )

    jac[jhe4, jar35] = (
       +rate_eval.Ar35__He4_S31
       +rho*Y[jn]*rate_eval.n_Ar35__He4_S32
       +rho*Y[jp]*rate_eval.p_Ar35__He4_Cl32
       )

    jac[jhe6, jn] = (
       +5.00000000000000e-01*rho**2*2*Y[jn]*Y[jhe4]*rate_eval.n_n_He4__He6
       )

    jac[jhe6, jhe4] = (
       +5.00000000000000e-01*rho**2*Y[jn]**2*rate_eval.n_n_He4__He6
       )

    jac[jhe6, jhe6] = (
       -rate_eval.He6__Li6__weak__wc12
       -rate_eval.He6__n_n_He4
       )

    jac[jli6, jn] = (
       -rho*Y[jli6]*rate_eval.n_Li6__Li7
       -rho*Y[jli6]*rate_eval.n_Li6__He4_t
       +rho*Y[jbe7]*rate_eval.n_Be7__d_Li6
       +rho**2*Y[jp]*Y[jhe4]*rate_eval.n_p_He4__Li6
       )

    jac[jli6, jp] = (
       -rho*Y[jli6]*rate_eval.p_Li6__Be7
       -rho*Y[jli6]*rate_eval.p_Li6__He4_He3
       +rho*Y[jli7]*rate_eval.p_Li7__d_Li6
       +rho*Y[jbe9]*rate_eval.p_Be9__He4_Li6
       +rho**2*Y[jn]*Y[jhe4]*rate_eval.n_p_He4__Li6
       )

    jac[jli6, jd] = (
       -rho*Y[jli6]*rate_eval.d_Li6__n_Be7
       -rho*Y[jli6]*rate_eval.d_Li6__p_Li7
       +rho*Y[jhe4]*rate_eval.He4_d__Li6
       )

    jac[jli6, jt] = (
       +rho*Y[jhe4]*rate_eval.He4_t__n_Li6
       )

    jac[jli6, jhe3] = (
       +rho*Y[jhe4]*rate_eval.He4_He3__p_Li6
       )

    jac[jli6, jhe4] = (
       -rho*Y[jli6]*rate_eval.He4_Li6__p_Be9
       +rho*Y[jd]*rate_eval.He4_d__Li6
       +rho*Y[jt]*rate_eval.He4_t__n_Li6
       +rho*Y[jhe3]*rate_eval.He4_He3__p_Li6
       +rho**2*Y[jn]*Y[jp]*rate_eval.n_p_He4__Li6
       )

    jac[jli6, jhe6] = (
       +rate_eval.He6__Li6__weak__wc12
       )

    jac[jli6, jli6] = (
       -rate_eval.Li6__He4_d
       -rate_eval.Li6__n_p_He4
       -rho*Y[jn]*rate_eval.n_Li6__Li7
       -rho*Y[jp]*rate_eval.p_Li6__Be7
       -rho*Y[jn]*rate_eval.n_Li6__He4_t
       -rho*Y[jp]*rate_eval.p_Li6__He4_He3
       -rho*Y[jd]*rate_eval.d_Li6__n_Be7
       -rho*Y[jd]*rate_eval.d_Li6__p_Li7
       -rho*Y[jhe4]*rate_eval.He4_Li6__p_Be9
       )

    jac[jli6, jli7] = (
       +rate_eval.Li7__n_Li6
       +rho*Y[jp]*rate_eval.p_Li7__d_Li6
       )

    jac[jli6, jbe7] = (
       +rate_eval.Be7__p_Li6
       +rho*Y[jn]*rate_eval.n_Be7__d_Li6
       )

    jac[jli6, jbe9] = (
       +rho*Y[jp]*rate_eval.p_Be9__He4_Li6
       )

    jac[jli7, jn] = (
       +rho*Y[jli6]*rate_eval.n_Li6__Li7
       +rho*Y[jbe7]*rate_eval.n_Be7__p_Li7
       +rho*Y[jbe9]*rate_eval.n_Be9__t_Li7
       +5.00000000000000e-01*rho**2*Y[jhe4]**2*rate_eval.n_He4_He4__d_Li7
       +2.50000000000000e-01*rho**3*2*Y[jn]*Y[jhe4]**2*rate_eval.n_n_He4_He4__t_Li7
       +5.00000000000000e-01*rho**3*Y[jp]*Y[jhe4]**2*rate_eval.n_p_He4_He4__He3_Li7
       )

    jac[jli7, jp] = (
       -rho*Y[jli7]*rate_eval.p_Li7__n_Be7
       -rho*Y[jli7]*rate_eval.p_Li7__d_Li6
       -rho*Y[jli7]*rate_eval.p_Li7__He4_He4
       +5.00000000000000e-01*rho**3*Y[jn]*Y[jhe4]**2*rate_eval.n_p_He4_He4__He3_Li7
       )

    jac[jli7, jd] = (
       -rho*Y[jli7]*rate_eval.d_Li7__n_He4_He4
       +rho*Y[jli6]*rate_eval.d_Li6__p_Li7
       )

    jac[jli7, jt] = (
       -rho*Y[jli7]*rate_eval.t_Li7__n_Be9
       -rho*Y[jli7]*rate_eval.t_Li7__n_n_He4_He4
       +rho*Y[jhe4]*rate_eval.He4_t__Li7
       )

    jac[jli7, jhe3] = (
       -rho*Y[jli7]*rate_eval.He3_Li7__n_p_He4_He4
       )

    jac[jli7, jhe4] = (
       +rho*Y[jt]*rate_eval.He4_t__Li7
       +5.00000000000000e-01*rho*2*Y[jhe4]*rate_eval.He4_He4__p_Li7
       +5.00000000000000e-01*rho**2*Y[jn]*2*Y[jhe4]*rate_eval.n_He4_He4__d_Li7
       +2.50000000000000e-01*rho**3*Y[jn]**2*2*Y[jhe4]*rate_eval.n_n_He4_He4__t_Li7
       +5.00000000000000e-01*rho**3*Y[jn]*Y[jp]*2*Y[jhe4]*rate_eval.n_p_He4_He4__He3_Li7
       )

    jac[jli7, jli6] = (
       +rho*Y[jn]*rate_eval.n_Li6__Li7
       +rho*Y[jd]*rate_eval.d_Li6__p_Li7
       )

    jac[jli7, jli7] = (
       -rate_eval.Li7__n_Li6
       -rate_eval.Li7__He4_t
       -rho*Y[jp]*rate_eval.p_Li7__n_Be7
       -rho*Y[jp]*rate_eval.p_Li7__d_Li6
       -rho*Y[jp]*rate_eval.p_Li7__He4_He4
       -rho*Y[jt]*rate_eval.t_Li7__n_Be9
       -rho*Y[jd]*rate_eval.d_Li7__n_He4_He4
       -rho*Y[jt]*rate_eval.t_Li7__n_n_He4_He4
       -rho*Y[jhe3]*rate_eval.He3_Li7__n_p_He4_He4
       )

    jac[jli7, jbe7] = (
       +rho*ye(Y)*rate_eval.Be7__Li7__weak__electron_capture
       +rho*Y[jn]*rate_eval.n_Be7__p_Li7
       )

    jac[jli7, jbe9] = (
       +rho*Y[jn]*rate_eval.n_Be9__t_Li7
       )

    jac[jbe7, jn] = (
       -rho*Y[jbe7]*rate_eval.n_Be7__p_Li7
       -rho*Y[jbe7]*rate_eval.n_Be7__d_Li6
       -rho*Y[jbe7]*rate_eval.n_Be7__He4_He4
       +5.00000000000000e-01*rho**3*Y[jp]*Y[jhe4]**2*rate_eval.n_p_He4_He4__t_Be7
       )

    jac[jbe7, jp] = (
       -rho*Y[jbe7]*rate_eval.p_Be7__B8
       +rho*Y[jli6]*rate_eval.p_Li6__Be7
       +rho*Y[jli7]*rate_eval.p_Li7__n_Be7
       +5.00000000000000e-01*rho**2*Y[jhe4]**2*rate_eval.p_He4_He4__d_Be7
       +5.00000000000000e-01*rho**3*Y[jn]*Y[jhe4]**2*rate_eval.n_p_He4_He4__t_Be7
       +2.50000000000000e-01*rho**3*2*Y[jp]*Y[jhe4]**2*rate_eval.p_p_He4_He4__He3_Be7
       )

    jac[jbe7, jd] = (
       -rho*Y[jbe7]*rate_eval.d_Be7__p_He4_He4
       +rho*Y[jli6]*rate_eval.d_Li6__n_Be7
       )

    jac[jbe7, jt] = (
       -rho*Y[jbe7]*rate_eval.t_Be7__n_p_He4_He4
       )

    jac[jbe7, jhe3] = (
       -rho*Y[jbe7]*rate_eval.He3_Be7__p_p_He4_He4
       +rho*Y[jhe4]*rate_eval.He4_He3__Be7
       )

    jac[jbe7, jhe4] = (
       -rho*Y[jbe7]*rate_eval.He4_Be7__C11
       +rho*Y[jhe3]*rate_eval.He4_He3__Be7
       +5.00000000000000e-01*rho*2*Y[jhe4]*rate_eval.He4_He4__n_Be7
       +5.00000000000000e-01*rho**2*Y[jp]*2*Y[jhe4]*rate_eval.p_He4_He4__d_Be7
       +5.00000000000000e-01*rho**3*Y[jn]*Y[jp]*2*Y[jhe4]*rate_eval.n_p_He4_He4__t_Be7
       +2.50000000000000e-01*rho**3*Y[jp]**2*2*Y[jhe4]*rate_eval.p_p_He4_He4__He3_Be7
       )

    jac[jbe7, jli6] = (
       +rho*Y[jp]*rate_eval.p_Li6__Be7
       +rho*Y[jd]*rate_eval.d_Li6__n_Be7
       )

    jac[jbe7, jli7] = (
       +rho*Y[jp]*rate_eval.p_Li7__n_Be7
       )

    jac[jbe7, jbe7] = (
       -rho*ye(Y)*rate_eval.Be7__Li7__weak__electron_capture
       -rate_eval.Be7__p_Li6
       -rate_eval.Be7__He4_He3
       -rho*Y[jp]*rate_eval.p_Be7__B8
       -rho*Y[jhe4]*rate_eval.He4_Be7__C11
       -rho*Y[jn]*rate_eval.n_Be7__p_Li7
       -rho*Y[jn]*rate_eval.n_Be7__d_Li6
       -rho*Y[jn]*rate_eval.n_Be7__He4_He4
       -rho*Y[jd]*rate_eval.d_Be7__p_He4_He4
       -rho*Y[jt]*rate_eval.t_Be7__n_p_He4_He4
       -rho*Y[jhe3]*rate_eval.He3_Be7__p_p_He4_He4
       )

    jac[jbe7, jb8] = (
       +rate_eval.B8__p_Be7
       )

    jac[jbe7, jc11] = (
       +rate_eval.C11__He4_Be7
       )

    jac[jbe9, jn] = (
       -rho*Y[jbe9]*rate_eval.n_Be9__t_Li7
       +rho*Y[jc12]*rate_eval.n_C12__He4_Be9
       +5.00000000000000e-01*rho**2*Y[jhe4]**2*rate_eval.n_He4_He4__Be9
       +5.00000000000000e-01*rho**3*Y[jp]*Y[jhe4]**2*rate_eval.n_p_He4_He4__p_Be9
       )

    jac[jbe9, jp] = (
       -rho*Y[jbe9]*rate_eval.p_Be9__He4_Li6
       -rho*Y[jbe9]*rate_eval.p_Be9__d_He4_He4
       -rho*Y[jbe9]*rate_eval.p_Be9__n_p_He4_He4
       +5.00000000000000e-01*rho**3*Y[jn]*Y[jhe4]**2*rate_eval.n_p_He4_He4__p_Be9
       )

    jac[jbe9, jd] = (
       +5.00000000000000e-01*rho**2*Y[jhe4]**2*rate_eval.d_He4_He4__p_Be9
       )

    jac[jbe9, jt] = (
       +rho*Y[jli7]*rate_eval.t_Li7__n_Be9
       )

    jac[jbe9, jhe4] = (
       -rho*Y[jbe9]*rate_eval.He4_Be9__n_C12
       +rho*Y[jli6]*rate_eval.He4_Li6__p_Be9
       +5.00000000000000e-01*rho**2*Y[jn]*2*Y[jhe4]*rate_eval.n_He4_He4__Be9
       +5.00000000000000e-01*rho**2*Y[jd]*2*Y[jhe4]*rate_eval.d_He4_He4__p_Be9
       +5.00000000000000e-01*rho**3*Y[jn]*Y[jp]*2*Y[jhe4]*rate_eval.n_p_He4_He4__p_Be9
       )

    jac[jbe9, jli6] = (
       +rho*Y[jhe4]*rate_eval.He4_Li6__p_Be9
       )

    jac[jbe9, jli7] = (
       +rho*Y[jt]*rate_eval.t_Li7__n_Be9
       )

    jac[jbe9, jbe9] = (
       -rate_eval.Be9__n_He4_He4
       -rho*Y[jn]*rate_eval.n_Be9__t_Li7
       -rho*Y[jp]*rate_eval.p_Be9__He4_Li6
       -rho*Y[jhe4]*rate_eval.He4_Be9__n_C12
       -rho*Y[jp]*rate_eval.p_Be9__d_He4_He4
       -rho*Y[jp]*rate_eval.p_Be9__n_p_He4_He4
       )

    jac[jbe9, jc12] = (
       +rho*Y[jn]*rate_eval.n_C12__He4_Be9
       )

    jac[jb8, jn] = (
       -rho*Y[jb8]*rate_eval.n_B8__p_He4_He4
       )

    jac[jb8, jp] = (
       +rho*Y[jbe7]*rate_eval.p_Be7__B8
       +rho*Y[jc11]*rate_eval.p_C11__He4_B8
       +5.00000000000000e-01*rho**2*Y[jhe4]**2*rate_eval.p_He4_He4__n_B8
       )

    jac[jb8, jhe4] = (
       -rho*Y[jb8]*rate_eval.He4_B8__p_C11
       +5.00000000000000e-01*rho**2*Y[jp]*2*Y[jhe4]*rate_eval.p_He4_He4__n_B8
       )

    jac[jb8, jbe7] = (
       +rho*Y[jp]*rate_eval.p_Be7__B8
       )

    jac[jb8, jb8] = (
       -rate_eval.B8__p_Be7
       -rate_eval.B8__He4_He4__weak__wc12
       -rho*Y[jhe4]*rate_eval.He4_B8__p_C11
       -rho*Y[jn]*rate_eval.n_B8__p_He4_He4
       )

    jac[jb8, jc11] = (
       +rho*Y[jp]*rate_eval.p_C11__He4_B8
       )

    jac[jc11, jn] = (
       -rho*Y[jc11]*rate_eval.n_C11__C12
       -rho*Y[jc11]*rate_eval.n_C11__He4_He4_He4
       +rho*Y[jo14]*rate_eval.n_O14__He4_C11
       )

    jac[jc11, jp] = (
       -rho*Y[jc11]*rate_eval.p_C11__He4_B8
       +rho*Y[jn14]*rate_eval.p_N14__He4_C11
       )

    jac[jc11, jhe4] = (
       -rho*Y[jc11]*rate_eval.He4_C11__n_O14
       -rho*Y[jc11]*rate_eval.He4_C11__p_N14
       +rho*Y[jbe7]*rate_eval.He4_Be7__C11
       +rho*Y[jb8]*rate_eval.He4_B8__p_C11
       +1.66666666666667e-01*rho**2*3*Y[jhe4]**2*rate_eval.He4_He4_He4__n_C11
       )

    jac[jc11, jbe7] = (
       +rho*Y[jhe4]*rate_eval.He4_Be7__C11
       )

    jac[jc11, jb8] = (
       +rho*Y[jhe4]*rate_eval.He4_B8__p_C11
       )

    jac[jc11, jc11] = (
       -rate_eval.C11__He4_Be7
       -rho*Y[jn]*rate_eval.n_C11__C12
       -rho*Y[jp]*rate_eval.p_C11__He4_B8
       -rho*Y[jhe4]*rate_eval.He4_C11__n_O14
       -rho*Y[jhe4]*rate_eval.He4_C11__p_N14
       -rho*Y[jn]*rate_eval.n_C11__He4_He4_He4
       )

    jac[jc11, jc12] = (
       +rate_eval.C12__n_C11
       )

    jac[jc11, jn14] = (
       +rho*Y[jp]*rate_eval.p_N14__He4_C11
       )

    jac[jc11, jo14] = (
       +rho*Y[jn]*rate_eval.n_O14__He4_C11
       )

    jac[jc12, jn] = (
       -rho*Y[jc12]*rate_eval.n_C12__C13
       -rho*Y[jc12]*rate_eval.n_C12__He4_Be9
       +rho*Y[jc11]*rate_eval.n_C11__C12
       +rho*Y[jo15]*rate_eval.n_O15__He4_C12
       +2*rho*Y[jmg23]*rate_eval.n_Mg23__C12_C12
       +rho*Y[js31]*rate_eval.n_S31__C12_Ne20
       )

    jac[jc12, jp] = (
       -rho*Y[jc12]*rate_eval.p_C12__N13
       +rho*Y[jn15]*rate_eval.p_N15__He4_C12
       +2*rho*Y[jna23]*rate_eval.p_Na23__C12_C12
       +rho*Y[jal27]*rate_eval.p_Al27__C12_O16
       +rho*Y[jp31]*rate_eval.p_P31__C12_Ne20
       )

    jac[jc12, jhe4] = (
       -rho*Y[jc12]*rate_eval.He4_C12__O16
       -rho*Y[jc12]*rate_eval.He4_C12__n_O15
       -rho*Y[jc12]*rate_eval.He4_C12__p_N15
       +rho*Y[jbe9]*rate_eval.He4_Be9__n_C12
       +2*rho*Y[jne20]*rate_eval.He4_Ne20__C12_C12
       +rho*Y[jmg24]*rate_eval.He4_Mg24__C12_O16
       +rho*Y[jsi28]*rate_eval.He4_Si28__C12_Ne20
       +1.66666666666667e-01*rho**2*3*Y[jhe4]**2*rate_eval.He4_He4_He4__C12
       )

    jac[jc12, jbe9] = (
       +rho*Y[jhe4]*rate_eval.He4_Be9__n_C12
       )

    jac[jc12, jc11] = (
       +rho*Y[jn]*rate_eval.n_C11__C12
       )

    jac[jc12, jc12] = (
       -rate_eval.C12__n_C11
       -rate_eval.C12__He4_He4_He4
       -rho*Y[jn]*rate_eval.n_C12__C13
       -rho*Y[jp]*rate_eval.p_C12__N13
       -rho*Y[jhe4]*rate_eval.He4_C12__O16
       -rho*Y[jn]*rate_eval.n_C12__He4_Be9
       -rho*Y[jhe4]*rate_eval.He4_C12__n_O15
       -rho*Y[jhe4]*rate_eval.He4_C12__p_N15
       -2*5.00000000000000e-01*rho*2*Y[jc12]*rate_eval.C12_C12__n_Mg23
       -2*5.00000000000000e-01*rho*2*Y[jc12]*rate_eval.C12_C12__p_Na23
       -2*5.00000000000000e-01*rho*2*Y[jc12]*rate_eval.C12_C12__He4_Ne20
       -rho*Y[jo16]*rate_eval.C12_O16__p_Al27
       -rho*Y[jo16]*rate_eval.C12_O16__He4_Mg24
       -rho*Y[jne20]*rate_eval.C12_Ne20__n_S31
       -rho*Y[jne20]*rate_eval.C12_Ne20__p_P31
       -rho*Y[jne20]*rate_eval.C12_Ne20__He4_Si28
       )

    jac[jc12, jc13] = (
       +rate_eval.C13__n_C12
       )

    jac[jc12, jn13] = (
       +rate_eval.N13__p_C12
       )

    jac[jc12, jn15] = (
       +rho*Y[jp]*rate_eval.p_N15__He4_C12
       )

    jac[jc12, jo15] = (
       +rho*Y[jn]*rate_eval.n_O15__He4_C12
       )

    jac[jc12, jo16] = (
       -rho*Y[jc12]*rate_eval.C12_O16__p_Al27
       -rho*Y[jc12]*rate_eval.C12_O16__He4_Mg24
       +rate_eval.O16__He4_C12
       )

    jac[jc12, jne20] = (
       -rho*Y[jc12]*rate_eval.C12_Ne20__n_S31
       -rho*Y[jc12]*rate_eval.C12_Ne20__p_P31
       -rho*Y[jc12]*rate_eval.C12_Ne20__He4_Si28
       +2*rho*Y[jhe4]*rate_eval.He4_Ne20__C12_C12
       )

    jac[jc12, jna23] = (
       +2*rho*Y[jp]*rate_eval.p_Na23__C12_C12
       )

    jac[jc12, jmg23] = (
       +2*rho*Y[jn]*rate_eval.n_Mg23__C12_C12
       )

    jac[jc12, jmg24] = (
       +rho*Y[jhe4]*rate_eval.He4_Mg24__C12_O16
       )

    jac[jc12, jal27] = (
       +rho*Y[jp]*rate_eval.p_Al27__C12_O16
       )

    jac[jc12, jsi28] = (
       +rho*Y[jhe4]*rate_eval.He4_Si28__C12_Ne20
       )

    jac[jc12, jp31] = (
       +rho*Y[jp]*rate_eval.p_P31__C12_Ne20
       )

    jac[jc12, js31] = (
       +rho*Y[jn]*rate_eval.n_S31__C12_Ne20
       )

    jac[jc13, jn] = (
       -rho*Y[jc13]*rate_eval.n_C13__C14
       +rho*Y[jc12]*rate_eval.n_C12__C13
       +rho*Y[jn13]*rate_eval.n_N13__p_C13
       +rho*Y[jn14]*rate_eval.n_N14__d_C13
       +rho*Y[jo16]*rate_eval.n_O16__He4_C13
       )

    jac[jc13, jp] = (
       -rho*Y[jc13]*rate_eval.p_C13__N14
       -rho*Y[jc13]*rate_eval.p_C13__n_N13
       )

    jac[jc13, jd] = (
       -rho*Y[jc13]*rate_eval.d_C13__n_N14
       )

    jac[jc13, jhe4] = (
       -rho*Y[jc13]*rate_eval.He4_C13__n_O16
       )

    jac[jc13, jc12] = (
       +rho*Y[jn]*rate_eval.n_C12__C13
       )

    jac[jc13, jc13] = (
       -rate_eval.C13__n_C12
       -rho*Y[jn]*rate_eval.n_C13__C14
       -rho*Y[jp]*rate_eval.p_C13__N14
       -rho*Y[jp]*rate_eval.p_C13__n_N13
       -rho*Y[jd]*rate_eval.d_C13__n_N14
       -rho*Y[jhe4]*rate_eval.He4_C13__n_O16
       )

    jac[jc13, jc14] = (
       +rate_eval.C14__n_C13
       )

    jac[jc13, jn13] = (
       +rate_eval.N13__C13__weak__wc12
       +rho*Y[jn]*rate_eval.n_N13__p_C13
       )

    jac[jc13, jn14] = (
       +rate_eval.N14__p_C13
       +rho*Y[jn]*rate_eval.n_N14__d_C13
       )

    jac[jc13, jo16] = (
       +rho*Y[jn]*rate_eval.n_O16__He4_C13
       )

    jac[jc14, jn] = (
       +rho*Y[jc13]*rate_eval.n_C13__C14
       +rho*Y[jn14]*rate_eval.n_N14__p_C14
       +rho*Y[jn15]*rate_eval.n_N15__d_C14
       +rho*Y[jo17]*rate_eval.n_O17__He4_C14
       )

    jac[jc14, jp] = (
       -rho*Y[jc14]*rate_eval.p_C14__N15
       -rho*Y[jc14]*rate_eval.p_C14__n_N14
       )

    jac[jc14, jd] = (
       -rho*Y[jc14]*rate_eval.d_C14__n_N15
       )

    jac[jc14, jhe4] = (
       -rho*Y[jc14]*rate_eval.He4_C14__O18
       -rho*Y[jc14]*rate_eval.He4_C14__n_O17
       )

    jac[jc14, jc13] = (
       +rho*Y[jn]*rate_eval.n_C13__C14
       )

    jac[jc14, jc14] = (
       -rate_eval.C14__N14__weak__wc12
       -rate_eval.C14__n_C13
       -rho*Y[jp]*rate_eval.p_C14__N15
       -rho*Y[jhe4]*rate_eval.He4_C14__O18
       -rho*Y[jp]*rate_eval.p_C14__n_N14
       -rho*Y[jd]*rate_eval.d_C14__n_N15
       -rho*Y[jhe4]*rate_eval.He4_C14__n_O17
       )

    jac[jc14, jn14] = (
       +rho*Y[jn]*rate_eval.n_N14__p_C14
       )

    jac[jc14, jn15] = (
       +rate_eval.N15__p_C14
       +rho*Y[jn]*rate_eval.n_N15__d_C14
       )

    jac[jc14, jo17] = (
       +rho*Y[jn]*rate_eval.n_O17__He4_C14
       )

    jac[jc14, jo18] = (
       +rate_eval.O18__He4_C14
       )

    jac[jn13, jn] = (
       -rho*Y[jn13]*rate_eval.n_N13__N14
       -rho*Y[jn13]*rate_eval.n_N13__p_C13
       +rho*Y[jf16]*rate_eval.n_F16__He4_N13
       )

    jac[jn13, jp] = (
       -rho*Y[jn13]*rate_eval.p_N13__O14
       +rho*Y[jc12]*rate_eval.p_C12__N13
       +rho*Y[jc13]*rate_eval.p_C13__n_N13
       +rho*Y[jo16]*rate_eval.p_O16__He4_N13
       )

    jac[jn13, jhe4] = (
       -rho*Y[jn13]*rate_eval.He4_N13__n_F16
       -rho*Y[jn13]*rate_eval.He4_N13__p_O16
       )

    jac[jn13, jc12] = (
       +rho*Y[jp]*rate_eval.p_C12__N13
       )

    jac[jn13, jc13] = (
       +rho*Y[jp]*rate_eval.p_C13__n_N13
       )

    jac[jn13, jn13] = (
       -rate_eval.N13__C13__weak__wc12
       -rate_eval.N13__p_C12
       -rho*Y[jn]*rate_eval.n_N13__N14
       -rho*Y[jp]*rate_eval.p_N13__O14
       -rho*Y[jn]*rate_eval.n_N13__p_C13
       -rho*Y[jhe4]*rate_eval.He4_N13__n_F16
       -rho*Y[jhe4]*rate_eval.He4_N13__p_O16
       )

    jac[jn13, jn14] = (
       +rate_eval.N14__n_N13
       )

    jac[jn13, jo14] = (
       +rate_eval.O14__p_N13
       )

    jac[jn13, jo16] = (
       +rho*Y[jp]*rate_eval.p_O16__He4_N13
       )

    jac[jn13, jf16] = (
       +rho*Y[jn]*rate_eval.n_F16__He4_N13
       )

    jac[jn14, jn] = (
       -rho*Y[jn14]*rate_eval.n_N14__N15
       -rho*Y[jn14]*rate_eval.n_N14__p_C14
       -rho*Y[jn14]*rate_eval.n_N14__d_C13
       +rho*Y[jn13]*rate_eval.n_N13__N14
       +rho*Y[jo14]*rate_eval.n_O14__p_N14
       +rho*Y[jf17]*rate_eval.n_F17__He4_N14
       )

    jac[jn14, jp] = (
       -rho*Y[jn14]*rate_eval.p_N14__O15
       -rho*Y[jn14]*rate_eval.p_N14__n_O14
       -rho*Y[jn14]*rate_eval.p_N14__He4_C11
       +rho*Y[jc13]*rate_eval.p_C13__N14
       +rho*Y[jc14]*rate_eval.p_C14__n_N14
       +rho*Y[jo17]*rate_eval.p_O17__He4_N14
       )

    jac[jn14, jd] = (
       +rho*Y[jc13]*rate_eval.d_C13__n_N14
       )

    jac[jn14, jhe4] = (
       -rho*Y[jn14]*rate_eval.He4_N14__F18
       -rho*Y[jn14]*rate_eval.He4_N14__n_F17
       -rho*Y[jn14]*rate_eval.He4_N14__p_O17
       +rho*Y[jc11]*rate_eval.He4_C11__p_N14
       )

    jac[jn14, jc11] = (
       +rho*Y[jhe4]*rate_eval.He4_C11__p_N14
       )

    jac[jn14, jc13] = (
       +rho*Y[jp]*rate_eval.p_C13__N14
       +rho*Y[jd]*rate_eval.d_C13__n_N14
       )

    jac[jn14, jc14] = (
       +rate_eval.C14__N14__weak__wc12
       +rho*Y[jp]*rate_eval.p_C14__n_N14
       )

    jac[jn14, jn13] = (
       +rho*Y[jn]*rate_eval.n_N13__N14
       )

    jac[jn14, jn14] = (
       -rate_eval.N14__n_N13
       -rate_eval.N14__p_C13
       -rho*Y[jn]*rate_eval.n_N14__N15
       -rho*Y[jp]*rate_eval.p_N14__O15
       -rho*Y[jhe4]*rate_eval.He4_N14__F18
       -rho*Y[jn]*rate_eval.n_N14__p_C14
       -rho*Y[jn]*rate_eval.n_N14__d_C13
       -rho*Y[jp]*rate_eval.p_N14__n_O14
       -rho*Y[jp]*rate_eval.p_N14__He4_C11
       -rho*Y[jhe4]*rate_eval.He4_N14__n_F17
       -rho*Y[jhe4]*rate_eval.He4_N14__p_O17
       )

    jac[jn14, jn15] = (
       +rate_eval.N15__n_N14
       )

    jac[jn14, jo14] = (
       +rate_eval.O14__N14__weak__wc12
       +rho*Y[jn]*rate_eval.n_O14__p_N14
       )

    jac[jn14, jo15] = (
       +rate_eval.O15__p_N14
       )

    jac[jn14, jo17] = (
       +rho*Y[jp]*rate_eval.p_O17__He4_N14
       )

    jac[jn14, jf17] = (
       +rho*Y[jn]*rate_eval.n_F17__He4_N14
       )

    jac[jn14, jf18] = (
       +rate_eval.F18__He4_N14
       )

    jac[jn15, jn] = (
       -rho*Y[jn15]*rate_eval.n_N15__d_C14
       +rho*Y[jn14]*rate_eval.n_N14__N15
       +rho*Y[jo15]*rate_eval.n_O15__p_N15
       +rho*Y[jf18]*rate_eval.n_F18__He4_N15
       )

    jac[jn15, jp] = (
       -rho*Y[jn15]*rate_eval.p_N15__O16
       -rho*Y[jn15]*rate_eval.p_N15__n_O15
       -rho*Y[jn15]*rate_eval.p_N15__He4_C12
       +rho*Y[jc14]*rate_eval.p_C14__N15
       +rho*Y[jo18]*rate_eval.p_O18__He4_N15
       )

    jac[jn15, jd] = (
       +rho*Y[jc14]*rate_eval.d_C14__n_N15
       )

    jac[jn15, jhe4] = (
       -rho*Y[jn15]*rate_eval.He4_N15__F19
       -rho*Y[jn15]*rate_eval.He4_N15__n_F18
       -rho*Y[jn15]*rate_eval.He4_N15__p_O18
       +rho*Y[jc12]*rate_eval.He4_C12__p_N15
       )

    jac[jn15, jc12] = (
       +rho*Y[jhe4]*rate_eval.He4_C12__p_N15
       )

    jac[jn15, jc14] = (
       +rho*Y[jp]*rate_eval.p_C14__N15
       +rho*Y[jd]*rate_eval.d_C14__n_N15
       )

    jac[jn15, jn14] = (
       +rho*Y[jn]*rate_eval.n_N14__N15
       )

    jac[jn15, jn15] = (
       -rate_eval.N15__n_N14
       -rate_eval.N15__p_C14
       -rho*Y[jp]*rate_eval.p_N15__O16
       -rho*Y[jhe4]*rate_eval.He4_N15__F19
       -rho*Y[jn]*rate_eval.n_N15__d_C14
       -rho*Y[jp]*rate_eval.p_N15__n_O15
       -rho*Y[jp]*rate_eval.p_N15__He4_C12
       -rho*Y[jhe4]*rate_eval.He4_N15__n_F18
       -rho*Y[jhe4]*rate_eval.He4_N15__p_O18
       )

    jac[jn15, jo15] = (
       +rate_eval.O15__N15__weak__wc12
       +rho*Y[jn]*rate_eval.n_O15__p_N15
       )

    jac[jn15, jo16] = (
       +rate_eval.O16__p_N15
       )

    jac[jn15, jo18] = (
       +rho*Y[jp]*rate_eval.p_O18__He4_N15
       )

    jac[jn15, jf18] = (
       +rho*Y[jn]*rate_eval.n_F18__He4_N15
       )

    jac[jn15, jf19] = (
       +rate_eval.F19__He4_N15
       )

    jac[jo14, jn] = (
       -rho*Y[jo14]*rate_eval.n_O14__O15
       -rho*Y[jo14]*rate_eval.n_O14__p_N14
       -rho*Y[jo14]*rate_eval.n_O14__He4_C11
       )

    jac[jo14, jp] = (
       -rho*Y[jo14]*rate_eval.p_O14__F15
       +rho*Y[jn13]*rate_eval.p_N13__O14
       +rho*Y[jn14]*rate_eval.p_N14__n_O14
       +rho*Y[jf17]*rate_eval.p_F17__He4_O14
       )

    jac[jo14, jhe4] = (
       -rho*Y[jo14]*rate_eval.He4_O14__Ne18
       -rho*Y[jo14]*rate_eval.He4_O14__p_F17
       +rho*Y[jc11]*rate_eval.He4_C11__n_O14
       )

    jac[jo14, jc11] = (
       +rho*Y[jhe4]*rate_eval.He4_C11__n_O14
       )

    jac[jo14, jn13] = (
       +rho*Y[jp]*rate_eval.p_N13__O14
       )

    jac[jo14, jn14] = (
       +rho*Y[jp]*rate_eval.p_N14__n_O14
       )

    jac[jo14, jo14] = (
       -rate_eval.O14__N14__weak__wc12
       -rate_eval.O14__p_N13
       -rho*Y[jn]*rate_eval.n_O14__O15
       -rho*Y[jp]*rate_eval.p_O14__F15
       -rho*Y[jhe4]*rate_eval.He4_O14__Ne18
       -rho*Y[jn]*rate_eval.n_O14__p_N14
       -rho*Y[jn]*rate_eval.n_O14__He4_C11
       -rho*Y[jhe4]*rate_eval.He4_O14__p_F17
       )

    jac[jo14, jo15] = (
       +rate_eval.O15__n_O14
       )

    jac[jo14, jf15] = (
       +rate_eval.F15__p_O14
       )

    jac[jo14, jf17] = (
       +rho*Y[jp]*rate_eval.p_F17__He4_O14
       )

    jac[jo14, jne18] = (
       +rate_eval.Ne18__He4_O14
       )

    jac[jo15, jn] = (
       -rho*Y[jo15]*rate_eval.n_O15__O16
       -rho*Y[jo15]*rate_eval.n_O15__p_N15
       -rho*Y[jo15]*rate_eval.n_O15__He4_C12
       +rho*Y[jo14]*rate_eval.n_O14__O15
       +rho*Y[jf15]*rate_eval.n_F15__p_O15
       +rho*Y[jne18]*rate_eval.n_Ne18__He4_O15
       )

    jac[jo15, jp] = (
       -rho*Y[jo15]*rate_eval.p_O15__F16
       -rho*Y[jo15]*rate_eval.p_O15__n_F15
       +rho*Y[jn14]*rate_eval.p_N14__O15
       +rho*Y[jn15]*rate_eval.p_N15__n_O15
       +rho*Y[jf18]*rate_eval.p_F18__He4_O15
       )

    jac[jo15, jhe4] = (
       -rho*Y[jo15]*rate_eval.He4_O15__Ne19
       -rho*Y[jo15]*rate_eval.He4_O15__n_Ne18
       -rho*Y[jo15]*rate_eval.He4_O15__p_F18
       +rho*Y[jc12]*rate_eval.He4_C12__n_O15
       )

    jac[jo15, jc12] = (
       +rho*Y[jhe4]*rate_eval.He4_C12__n_O15
       )

    jac[jo15, jn14] = (
       +rho*Y[jp]*rate_eval.p_N14__O15
       )

    jac[jo15, jn15] = (
       +rho*Y[jp]*rate_eval.p_N15__n_O15
       )

    jac[jo15, jo14] = (
       +rho*Y[jn]*rate_eval.n_O14__O15
       )

    jac[jo15, jo15] = (
       -rate_eval.O15__N15__weak__wc12
       -rate_eval.O15__n_O14
       -rate_eval.O15__p_N14
       -rho*Y[jn]*rate_eval.n_O15__O16
       -rho*Y[jp]*rate_eval.p_O15__F16
       -rho*Y[jhe4]*rate_eval.He4_O15__Ne19
       -rho*Y[jn]*rate_eval.n_O15__p_N15
       -rho*Y[jn]*rate_eval.n_O15__He4_C12
       -rho*Y[jp]*rate_eval.p_O15__n_F15
       -rho*Y[jhe4]*rate_eval.He4_O15__n_Ne18
       -rho*Y[jhe4]*rate_eval.He4_O15__p_F18
       )

    jac[jo15, jo16] = (
       +rate_eval.O16__n_O15
       )

    jac[jo15, jf15] = (
       +rho*Y[jn]*rate_eval.n_F15__p_O15
       )

    jac[jo15, jf16] = (
       +rate_eval.F16__p_O15
       )

    jac[jo15, jf18] = (
       +rho*Y[jp]*rate_eval.p_F18__He4_O15
       )

    jac[jo15, jne18] = (
       +rho*Y[jn]*rate_eval.n_Ne18__He4_O15
       )

    jac[jo15, jne19] = (
       +rate_eval.Ne19__He4_O15
       )

    jac[jo16, jn] = (
       -rho*Y[jo16]*rate_eval.n_O16__O17
       -rho*Y[jo16]*rate_eval.n_O16__He4_C13
       +rho*Y[jo15]*rate_eval.n_O15__O16
       +rho*Y[jf16]*rate_eval.n_F16__p_O16
       +rho*Y[jne19]*rate_eval.n_Ne19__He4_O16
       +2*rho*Y[js31]*rate_eval.n_S31__O16_O16
       )

    jac[jo16, jp] = (
       -rho*Y[jo16]*rate_eval.p_O16__F17
       -rho*Y[jo16]*rate_eval.p_O16__n_F16
       -rho*Y[jo16]*rate_eval.p_O16__He4_N13
       +rho*Y[jn15]*rate_eval.p_N15__O16
       +rho*Y[jf19]*rate_eval.p_F19__He4_O16
       +rho*Y[jal27]*rate_eval.p_Al27__C12_O16
       +2*rho*Y[jp31]*rate_eval.p_P31__O16_O16
       )

    jac[jo16, jhe4] = (
       -rho*Y[jo16]*rate_eval.He4_O16__Ne20
       -rho*Y[jo16]*rate_eval.He4_O16__n_Ne19
       -rho*Y[jo16]*rate_eval.He4_O16__p_F19
       +rho*Y[jc12]*rate_eval.He4_C12__O16
       +rho*Y[jc13]*rate_eval.He4_C13__n_O16
       +rho*Y[jn13]*rate_eval.He4_N13__p_O16
       +rho*Y[jmg24]*rate_eval.He4_Mg24__C12_O16
       +2*rho*Y[jsi28]*rate_eval.He4_Si28__O16_O16
       )

    jac[jo16, jc12] = (
       -rho*Y[jo16]*rate_eval.C12_O16__p_Al27
       -rho*Y[jo16]*rate_eval.C12_O16__He4_Mg24
       +rho*Y[jhe4]*rate_eval.He4_C12__O16
       )

    jac[jo16, jc13] = (
       +rho*Y[jhe4]*rate_eval.He4_C13__n_O16
       )

    jac[jo16, jn13] = (
       +rho*Y[jhe4]*rate_eval.He4_N13__p_O16
       )

    jac[jo16, jn15] = (
       +rho*Y[jp]*rate_eval.p_N15__O16
       )

    jac[jo16, jo15] = (
       +rho*Y[jn]*rate_eval.n_O15__O16
       )

    jac[jo16, jo16] = (
       -rate_eval.O16__n_O15
       -rate_eval.O16__p_N15
       -rate_eval.O16__He4_C12
       -rho*Y[jn]*rate_eval.n_O16__O17
       -rho*Y[jp]*rate_eval.p_O16__F17
       -rho*Y[jhe4]*rate_eval.He4_O16__Ne20
       -rho*Y[jn]*rate_eval.n_O16__He4_C13
       -rho*Y[jp]*rate_eval.p_O16__n_F16
       -rho*Y[jp]*rate_eval.p_O16__He4_N13
       -rho*Y[jhe4]*rate_eval.He4_O16__n_Ne19
       -rho*Y[jhe4]*rate_eval.He4_O16__p_F19
       -rho*Y[jc12]*rate_eval.C12_O16__p_Al27
       -rho*Y[jc12]*rate_eval.C12_O16__He4_Mg24
       -2*5.00000000000000e-01*rho*2*Y[jo16]*rate_eval.O16_O16__n_S31
       -2*5.00000000000000e-01*rho*2*Y[jo16]*rate_eval.O16_O16__p_P31
       -2*5.00000000000000e-01*rho*2*Y[jo16]*rate_eval.O16_O16__He4_Si28
       )

    jac[jo16, jo17] = (
       +rate_eval.O17__n_O16
       )

    jac[jo16, jf16] = (
       +rho*Y[jn]*rate_eval.n_F16__p_O16
       )

    jac[jo16, jf17] = (
       +rate_eval.F17__p_O16
       )

    jac[jo16, jf19] = (
       +rho*Y[jp]*rate_eval.p_F19__He4_O16
       )

    jac[jo16, jne19] = (
       +rho*Y[jn]*rate_eval.n_Ne19__He4_O16
       )

    jac[jo16, jne20] = (
       +rate_eval.Ne20__He4_O16
       )

    jac[jo16, jmg24] = (
       +rho*Y[jhe4]*rate_eval.He4_Mg24__C12_O16
       )

    jac[jo16, jal27] = (
       +rho*Y[jp]*rate_eval.p_Al27__C12_O16
       )

    jac[jo16, jsi28] = (
       +2*rho*Y[jhe4]*rate_eval.He4_Si28__O16_O16
       )

    jac[jo16, jp31] = (
       +2*rho*Y[jp]*rate_eval.p_P31__O16_O16
       )

    jac[jo16, js31] = (
       +2*rho*Y[jn]*rate_eval.n_S31__O16_O16
       )

    jac[jo17, jn] = (
       -rho*Y[jo17]*rate_eval.n_O17__O18
       -rho*Y[jo17]*rate_eval.n_O17__He4_C14
       +rho*Y[jo16]*rate_eval.n_O16__O17
       +rho*Y[jf17]*rate_eval.n_F17__p_O17
       +rho*Y[jne20]*rate_eval.n_Ne20__He4_O17
       )

    jac[jo17, jp] = (
       -rho*Y[jo17]*rate_eval.p_O17__F18
       -rho*Y[jo17]*rate_eval.p_O17__n_F17
       -rho*Y[jo17]*rate_eval.p_O17__He4_N14
       )

    jac[jo17, jhe4] = (
       -rho*Y[jo17]*rate_eval.He4_O17__Ne21
       -rho*Y[jo17]*rate_eval.He4_O17__n_Ne20
       +rho*Y[jc14]*rate_eval.He4_C14__n_O17
       +rho*Y[jn14]*rate_eval.He4_N14__p_O17
       )

    jac[jo17, jc14] = (
       +rho*Y[jhe4]*rate_eval.He4_C14__n_O17
       )

    jac[jo17, jn14] = (
       +rho*Y[jhe4]*rate_eval.He4_N14__p_O17
       )

    jac[jo17, jo16] = (
       +rho*Y[jn]*rate_eval.n_O16__O17
       )

    jac[jo17, jo17] = (
       -rate_eval.O17__n_O16
       -rho*Y[jn]*rate_eval.n_O17__O18
       -rho*Y[jp]*rate_eval.p_O17__F18
       -rho*Y[jhe4]*rate_eval.He4_O17__Ne21
       -rho*Y[jn]*rate_eval.n_O17__He4_C14
       -rho*Y[jp]*rate_eval.p_O17__n_F17
       -rho*Y[jp]*rate_eval.p_O17__He4_N14
       -rho*Y[jhe4]*rate_eval.He4_O17__n_Ne20
       )

    jac[jo17, jo18] = (
       +rate_eval.O18__n_O17
       )

    jac[jo17, jf17] = (
       +rate_eval.F17__O17__weak__wc12
       +rho*Y[jn]*rate_eval.n_F17__p_O17
       )

    jac[jo17, jf18] = (
       +rate_eval.F18__p_O17
       )

    jac[jo17, jne20] = (
       +rho*Y[jn]*rate_eval.n_Ne20__He4_O17
       )

    jac[jo17, jne21] = (
       +rate_eval.Ne21__He4_O17
       )

    jac[jo18, jn] = (
       +rho*Y[jo17]*rate_eval.n_O17__O18
       +rho*Y[jf18]*rate_eval.n_F18__p_O18
       +rho*Y[jne21]*rate_eval.n_Ne21__He4_O18
       )

    jac[jo18, jp] = (
       -rho*Y[jo18]*rate_eval.p_O18__F19
       -rho*Y[jo18]*rate_eval.p_O18__n_F18
       -rho*Y[jo18]*rate_eval.p_O18__He4_N15
       )

    jac[jo18, jhe4] = (
       -rho*Y[jo18]*rate_eval.He4_O18__n_Ne21
       +rho*Y[jc14]*rate_eval.He4_C14__O18
       +rho*Y[jn15]*rate_eval.He4_N15__p_O18
       )

    jac[jo18, jc14] = (
       +rho*Y[jhe4]*rate_eval.He4_C14__O18
       )

    jac[jo18, jn15] = (
       +rho*Y[jhe4]*rate_eval.He4_N15__p_O18
       )

    jac[jo18, jo17] = (
       +rho*Y[jn]*rate_eval.n_O17__O18
       )

    jac[jo18, jo18] = (
       -rate_eval.O18__n_O17
       -rate_eval.O18__He4_C14
       -rho*Y[jp]*rate_eval.p_O18__F19
       -rho*Y[jp]*rate_eval.p_O18__n_F18
       -rho*Y[jp]*rate_eval.p_O18__He4_N15
       -rho*Y[jhe4]*rate_eval.He4_O18__n_Ne21
       )

    jac[jo18, jf18] = (
       +rate_eval.F18__O18__weak__wc12
       +rho*Y[jn]*rate_eval.n_F18__p_O18
       )

    jac[jo18, jf19] = (
       +rate_eval.F19__p_O18
       )

    jac[jo18, jne21] = (
       +rho*Y[jn]*rate_eval.n_Ne21__He4_O18
       )

    jac[jf15, jn] = (
       -rho*Y[jf15]*rate_eval.n_F15__F16
       -rho*Y[jf15]*rate_eval.n_F15__p_O15
       )

    jac[jf15, jp] = (
       +rho*Y[jo14]*rate_eval.p_O14__F15
       +rho*Y[jo15]*rate_eval.p_O15__n_F15
       +rho*Y[jne18]*rate_eval.p_Ne18__He4_F15
       )

    jac[jf15, jhe4] = (
       -rho*Y[jf15]*rate_eval.He4_F15__p_Ne18
       )

    jac[jf15, jo14] = (
       +rho*Y[jp]*rate_eval.p_O14__F15
       )

    jac[jf15, jo15] = (
       +rho*Y[jp]*rate_eval.p_O15__n_F15
       )

    jac[jf15, jf15] = (
       -rate_eval.F15__p_O14
       -rho*Y[jn]*rate_eval.n_F15__F16
       -rho*Y[jn]*rate_eval.n_F15__p_O15
       -rho*Y[jhe4]*rate_eval.He4_F15__p_Ne18
       )

    jac[jf15, jf16] = (
       +rate_eval.F16__n_F15
       )

    jac[jf15, jne18] = (
       +rho*Y[jp]*rate_eval.p_Ne18__He4_F15
       )

    jac[jf16, jn] = (
       -rho*Y[jf16]*rate_eval.n_F16__F17
       -rho*Y[jf16]*rate_eval.n_F16__p_O16
       -rho*Y[jf16]*rate_eval.n_F16__He4_N13
       +rho*Y[jf15]*rate_eval.n_F15__F16
       )

    jac[jf16, jp] = (
       +rho*Y[jo15]*rate_eval.p_O15__F16
       +rho*Y[jo16]*rate_eval.p_O16__n_F16
       +rho*Y[jne19]*rate_eval.p_Ne19__He4_F16
       )

    jac[jf16, jhe4] = (
       -rho*Y[jf16]*rate_eval.He4_F16__p_Ne19
       +rho*Y[jn13]*rate_eval.He4_N13__n_F16
       )

    jac[jf16, jn13] = (
       +rho*Y[jhe4]*rate_eval.He4_N13__n_F16
       )

    jac[jf16, jo15] = (
       +rho*Y[jp]*rate_eval.p_O15__F16
       )

    jac[jf16, jo16] = (
       +rho*Y[jp]*rate_eval.p_O16__n_F16
       )

    jac[jf16, jf15] = (
       +rho*Y[jn]*rate_eval.n_F15__F16
       )

    jac[jf16, jf16] = (
       -rate_eval.F16__n_F15
       -rate_eval.F16__p_O15
       -rho*Y[jn]*rate_eval.n_F16__F17
       -rho*Y[jn]*rate_eval.n_F16__p_O16
       -rho*Y[jn]*rate_eval.n_F16__He4_N13
       -rho*Y[jhe4]*rate_eval.He4_F16__p_Ne19
       )

    jac[jf16, jf17] = (
       +rate_eval.F17__n_F16
       )

    jac[jf16, jne19] = (
       +rho*Y[jp]*rate_eval.p_Ne19__He4_F16
       )

    jac[jf17, jn] = (
       -rho*Y[jf17]*rate_eval.n_F17__F18
       -rho*Y[jf17]*rate_eval.n_F17__p_O17
       -rho*Y[jf17]*rate_eval.n_F17__He4_N14
       +rho*Y[jf16]*rate_eval.n_F16__F17
       )

    jac[jf17, jp] = (
       -rho*Y[jf17]*rate_eval.p_F17__Ne18
       -rho*Y[jf17]*rate_eval.p_F17__He4_O14
       +rho*Y[jo16]*rate_eval.p_O16__F17
       +rho*Y[jo17]*rate_eval.p_O17__n_F17
       +rho*Y[jne20]*rate_eval.p_Ne20__He4_F17
       )

    jac[jf17, jhe4] = (
       -rho*Y[jf17]*rate_eval.He4_F17__Na21
       -rho*Y[jf17]*rate_eval.He4_F17__p_Ne20
       +rho*Y[jn14]*rate_eval.He4_N14__n_F17
       +rho*Y[jo14]*rate_eval.He4_O14__p_F17
       )

    jac[jf17, jn14] = (
       +rho*Y[jhe4]*rate_eval.He4_N14__n_F17
       )

    jac[jf17, jo14] = (
       +rho*Y[jhe4]*rate_eval.He4_O14__p_F17
       )

    jac[jf17, jo16] = (
       +rho*Y[jp]*rate_eval.p_O16__F17
       )

    jac[jf17, jo17] = (
       +rho*Y[jp]*rate_eval.p_O17__n_F17
       )

    jac[jf17, jf16] = (
       +rho*Y[jn]*rate_eval.n_F16__F17
       )

    jac[jf17, jf17] = (
       -rate_eval.F17__O17__weak__wc12
       -rate_eval.F17__n_F16
       -rate_eval.F17__p_O16
       -rho*Y[jn]*rate_eval.n_F17__F18
       -rho*Y[jp]*rate_eval.p_F17__Ne18
       -rho*Y[jhe4]*rate_eval.He4_F17__Na21
       -rho*Y[jn]*rate_eval.n_F17__p_O17
       -rho*Y[jn]*rate_eval.n_F17__He4_N14
       -rho*Y[jp]*rate_eval.p_F17__He4_O14
       -rho*Y[jhe4]*rate_eval.He4_F17__p_Ne20
       )

    jac[jf17, jf18] = (
       +rate_eval.F18__n_F17
       )

    jac[jf17, jne18] = (
       +rate_eval.Ne18__p_F17
       )

    jac[jf17, jne20] = (
       +rho*Y[jp]*rate_eval.p_Ne20__He4_F17
       )

    jac[jf17, jna21] = (
       +rate_eval.Na21__He4_F17
       )

    jac[jf18, jn] = (
       -rho*Y[jf18]*rate_eval.n_F18__F19
       -rho*Y[jf18]*rate_eval.n_F18__p_O18
       -rho*Y[jf18]*rate_eval.n_F18__He4_N15
       +rho*Y[jf17]*rate_eval.n_F17__F18
       +rho*Y[jne18]*rate_eval.n_Ne18__p_F18
       +rho*Y[jna21]*rate_eval.n_Na21__He4_F18
       )

    jac[jf18, jp] = (
       -rho*Y[jf18]*rate_eval.p_F18__Ne19
       -rho*Y[jf18]*rate_eval.p_F18__n_Ne18
       -rho*Y[jf18]*rate_eval.p_F18__He4_O15
       +rho*Y[jo17]*rate_eval.p_O17__F18
       +rho*Y[jo18]*rate_eval.p_O18__n_F18
       +rho*Y[jne21]*rate_eval.p_Ne21__He4_F18
       )

    jac[jf18, jhe4] = (
       -rho*Y[jf18]*rate_eval.He4_F18__Na22
       -rho*Y[jf18]*rate_eval.He4_F18__n_Na21
       -rho*Y[jf18]*rate_eval.He4_F18__p_Ne21
       +rho*Y[jn14]*rate_eval.He4_N14__F18
       +rho*Y[jn15]*rate_eval.He4_N15__n_F18
       +rho*Y[jo15]*rate_eval.He4_O15__p_F18
       )

    jac[jf18, jn14] = (
       +rho*Y[jhe4]*rate_eval.He4_N14__F18
       )

    jac[jf18, jn15] = (
       +rho*Y[jhe4]*rate_eval.He4_N15__n_F18
       )

    jac[jf18, jo15] = (
       +rho*Y[jhe4]*rate_eval.He4_O15__p_F18
       )

    jac[jf18, jo17] = (
       +rho*Y[jp]*rate_eval.p_O17__F18
       )

    jac[jf18, jo18] = (
       +rho*Y[jp]*rate_eval.p_O18__n_F18
       )

    jac[jf18, jf17] = (
       +rho*Y[jn]*rate_eval.n_F17__F18
       )

    jac[jf18, jf18] = (
       -rate_eval.F18__O18__weak__wc12
       -rate_eval.F18__n_F17
       -rate_eval.F18__p_O17
       -rate_eval.F18__He4_N14
       -rho*Y[jn]*rate_eval.n_F18__F19
       -rho*Y[jp]*rate_eval.p_F18__Ne19
       -rho*Y[jhe4]*rate_eval.He4_F18__Na22
       -rho*Y[jn]*rate_eval.n_F18__p_O18
       -rho*Y[jn]*rate_eval.n_F18__He4_N15
       -rho*Y[jp]*rate_eval.p_F18__n_Ne18
       -rho*Y[jp]*rate_eval.p_F18__He4_O15
       -rho*Y[jhe4]*rate_eval.He4_F18__n_Na21
       -rho*Y[jhe4]*rate_eval.He4_F18__p_Ne21
       )

    jac[jf18, jf19] = (
       +rate_eval.F19__n_F18
       )

    jac[jf18, jne18] = (
       +rate_eval.Ne18__F18__weak__wc12
       +rho*Y[jn]*rate_eval.n_Ne18__p_F18
       )

    jac[jf18, jne19] = (
       +rate_eval.Ne19__p_F18
       )

    jac[jf18, jne21] = (
       +rho*Y[jp]*rate_eval.p_Ne21__He4_F18
       )

    jac[jf18, jna21] = (
       +rho*Y[jn]*rate_eval.n_Na21__He4_F18
       )

    jac[jf18, jna22] = (
       +rate_eval.Na22__He4_F18
       )

    jac[jf19, jn] = (
       +rho*Y[jf18]*rate_eval.n_F18__F19
       +rho*Y[jne19]*rate_eval.n_Ne19__p_F19
       +rho*Y[jna22]*rate_eval.n_Na22__He4_F19
       )

    jac[jf19, jp] = (
       -rho*Y[jf19]*rate_eval.p_F19__Ne20
       -rho*Y[jf19]*rate_eval.p_F19__n_Ne19
       -rho*Y[jf19]*rate_eval.p_F19__He4_O16
       +rho*Y[jo18]*rate_eval.p_O18__F19
       )

    jac[jf19, jhe4] = (
       -rho*Y[jf19]*rate_eval.He4_F19__Na23
       -rho*Y[jf19]*rate_eval.He4_F19__n_Na22
       +rho*Y[jn15]*rate_eval.He4_N15__F19
       +rho*Y[jo16]*rate_eval.He4_O16__p_F19
       )

    jac[jf19, jn15] = (
       +rho*Y[jhe4]*rate_eval.He4_N15__F19
       )

    jac[jf19, jo16] = (
       +rho*Y[jhe4]*rate_eval.He4_O16__p_F19
       )

    jac[jf19, jo18] = (
       +rho*Y[jp]*rate_eval.p_O18__F19
       )

    jac[jf19, jf18] = (
       +rho*Y[jn]*rate_eval.n_F18__F19
       )

    jac[jf19, jf19] = (
       -rate_eval.F19__n_F18
       -rate_eval.F19__p_O18
       -rate_eval.F19__He4_N15
       -rho*Y[jp]*rate_eval.p_F19__Ne20
       -rho*Y[jhe4]*rate_eval.He4_F19__Na23
       -rho*Y[jp]*rate_eval.p_F19__n_Ne19
       -rho*Y[jp]*rate_eval.p_F19__He4_O16
       -rho*Y[jhe4]*rate_eval.He4_F19__n_Na22
       )

    jac[jf19, jne19] = (
       +rate_eval.Ne19__F19__weak__wc12
       +rho*Y[jn]*rate_eval.n_Ne19__p_F19
       )

    jac[jf19, jne20] = (
       +rate_eval.Ne20__p_F19
       )

    jac[jf19, jna22] = (
       +rho*Y[jn]*rate_eval.n_Na22__He4_F19
       )

    jac[jf19, jna23] = (
       +rate_eval.Na23__He4_F19
       )

    jac[jne18, jn] = (
       -rho*Y[jne18]*rate_eval.n_Ne18__Ne19
       -rho*Y[jne18]*rate_eval.n_Ne18__p_F18
       -rho*Y[jne18]*rate_eval.n_Ne18__He4_O15
       )

    jac[jne18, jp] = (
       -rho*Y[jne18]*rate_eval.p_Ne18__He4_F15
       +rho*Y[jf17]*rate_eval.p_F17__Ne18
       +rho*Y[jf18]*rate_eval.p_F18__n_Ne18
       +rho*Y[jna21]*rate_eval.p_Na21__He4_Ne18
       )

    jac[jne18, jhe4] = (
       -rho*Y[jne18]*rate_eval.He4_Ne18__Mg22
       -rho*Y[jne18]*rate_eval.He4_Ne18__p_Na21
       +rho*Y[jo14]*rate_eval.He4_O14__Ne18
       +rho*Y[jo15]*rate_eval.He4_O15__n_Ne18
       +rho*Y[jf15]*rate_eval.He4_F15__p_Ne18
       )

    jac[jne18, jo14] = (
       +rho*Y[jhe4]*rate_eval.He4_O14__Ne18
       )

    jac[jne18, jo15] = (
       +rho*Y[jhe4]*rate_eval.He4_O15__n_Ne18
       )

    jac[jne18, jf15] = (
       +rho*Y[jhe4]*rate_eval.He4_F15__p_Ne18
       )

    jac[jne18, jf17] = (
       +rho*Y[jp]*rate_eval.p_F17__Ne18
       )

    jac[jne18, jf18] = (
       +rho*Y[jp]*rate_eval.p_F18__n_Ne18
       )

    jac[jne18, jne18] = (
       -rate_eval.Ne18__F18__weak__wc12
       -rate_eval.Ne18__p_F17
       -rate_eval.Ne18__He4_O14
       -rho*Y[jn]*rate_eval.n_Ne18__Ne19
       -rho*Y[jhe4]*rate_eval.He4_Ne18__Mg22
       -rho*Y[jn]*rate_eval.n_Ne18__p_F18
       -rho*Y[jn]*rate_eval.n_Ne18__He4_O15
       -rho*Y[jp]*rate_eval.p_Ne18__He4_F15
       -rho*Y[jhe4]*rate_eval.He4_Ne18__p_Na21
       )

    jac[jne18, jne19] = (
       +rate_eval.Ne19__n_Ne18
       )

    jac[jne18, jna21] = (
       +rho*Y[jp]*rate_eval.p_Na21__He4_Ne18
       )

    jac[jne18, jmg22] = (
       +rate_eval.Mg22__He4_Ne18
       )

    jac[jne19, jn] = (
       -rho*Y[jne19]*rate_eval.n_Ne19__Ne20
       -rho*Y[jne19]*rate_eval.n_Ne19__p_F19
       -rho*Y[jne19]*rate_eval.n_Ne19__He4_O16
       +rho*Y[jne18]*rate_eval.n_Ne18__Ne19
       +rho*Y[jmg22]*rate_eval.n_Mg22__He4_Ne19
       )

    jac[jne19, jp] = (
       -rho*Y[jne19]*rate_eval.p_Ne19__He4_F16
       +rho*Y[jf18]*rate_eval.p_F18__Ne19
       +rho*Y[jf19]*rate_eval.p_F19__n_Ne19
       +rho*Y[jna22]*rate_eval.p_Na22__He4_Ne19
       )

    jac[jne19, jhe4] = (
       -rho*Y[jne19]*rate_eval.He4_Ne19__Mg23
       -rho*Y[jne19]*rate_eval.He4_Ne19__n_Mg22
       -rho*Y[jne19]*rate_eval.He4_Ne19__p_Na22
       +rho*Y[jo15]*rate_eval.He4_O15__Ne19
       +rho*Y[jo16]*rate_eval.He4_O16__n_Ne19
       +rho*Y[jf16]*rate_eval.He4_F16__p_Ne19
       )

    jac[jne19, jo15] = (
       +rho*Y[jhe4]*rate_eval.He4_O15__Ne19
       )

    jac[jne19, jo16] = (
       +rho*Y[jhe4]*rate_eval.He4_O16__n_Ne19
       )

    jac[jne19, jf16] = (
       +rho*Y[jhe4]*rate_eval.He4_F16__p_Ne19
       )

    jac[jne19, jf18] = (
       +rho*Y[jp]*rate_eval.p_F18__Ne19
       )

    jac[jne19, jf19] = (
       +rho*Y[jp]*rate_eval.p_F19__n_Ne19
       )

    jac[jne19, jne18] = (
       +rho*Y[jn]*rate_eval.n_Ne18__Ne19
       )

    jac[jne19, jne19] = (
       -rate_eval.Ne19__F19__weak__wc12
       -rate_eval.Ne19__n_Ne18
       -rate_eval.Ne19__p_F18
       -rate_eval.Ne19__He4_O15
       -rho*Y[jn]*rate_eval.n_Ne19__Ne20
       -rho*Y[jhe4]*rate_eval.He4_Ne19__Mg23
       -rho*Y[jn]*rate_eval.n_Ne19__p_F19
       -rho*Y[jn]*rate_eval.n_Ne19__He4_O16
       -rho*Y[jp]*rate_eval.p_Ne19__He4_F16
       -rho*Y[jhe4]*rate_eval.He4_Ne19__n_Mg22
       -rho*Y[jhe4]*rate_eval.He4_Ne19__p_Na22
       )

    jac[jne19, jne20] = (
       +rate_eval.Ne20__n_Ne19
       )

    jac[jne19, jna22] = (
       +rho*Y[jp]*rate_eval.p_Na22__He4_Ne19
       )

    jac[jne19, jmg22] = (
       +rho*Y[jn]*rate_eval.n_Mg22__He4_Ne19
       )

    jac[jne19, jmg23] = (
       +rate_eval.Mg23__He4_Ne19
       )

    jac[jne20, jn] = (
       -rho*Y[jne20]*rate_eval.n_Ne20__Ne21
       -rho*Y[jne20]*rate_eval.n_Ne20__He4_O17
       +rho*Y[jne19]*rate_eval.n_Ne19__Ne20
       +rho*Y[jmg23]*rate_eval.n_Mg23__He4_Ne20
       +rho*Y[js31]*rate_eval.n_S31__C12_Ne20
       )

    jac[jne20, jp] = (
       -rho*Y[jne20]*rate_eval.p_Ne20__Na21
       -rho*Y[jne20]*rate_eval.p_Ne20__He4_F17
       +rho*Y[jf19]*rate_eval.p_F19__Ne20
       +rho*Y[jna23]*rate_eval.p_Na23__He4_Ne20
       +rho*Y[jp31]*rate_eval.p_P31__C12_Ne20
       )

    jac[jne20, jhe4] = (
       -rho*Y[jne20]*rate_eval.He4_Ne20__Mg24
       -rho*Y[jne20]*rate_eval.He4_Ne20__n_Mg23
       -rho*Y[jne20]*rate_eval.He4_Ne20__p_Na23
       -rho*Y[jne20]*rate_eval.He4_Ne20__C12_C12
       +rho*Y[jo16]*rate_eval.He4_O16__Ne20
       +rho*Y[jo17]*rate_eval.He4_O17__n_Ne20
       +rho*Y[jf17]*rate_eval.He4_F17__p_Ne20
       +rho*Y[jsi28]*rate_eval.He4_Si28__C12_Ne20
       )

    jac[jne20, jc12] = (
       -rho*Y[jne20]*rate_eval.C12_Ne20__n_S31
       -rho*Y[jne20]*rate_eval.C12_Ne20__p_P31
       -rho*Y[jne20]*rate_eval.C12_Ne20__He4_Si28
       +5.00000000000000e-01*rho*2*Y[jc12]*rate_eval.C12_C12__He4_Ne20
       )

    jac[jne20, jo16] = (
       +rho*Y[jhe4]*rate_eval.He4_O16__Ne20
       )

    jac[jne20, jo17] = (
       +rho*Y[jhe4]*rate_eval.He4_O17__n_Ne20
       )

    jac[jne20, jf17] = (
       +rho*Y[jhe4]*rate_eval.He4_F17__p_Ne20
       )

    jac[jne20, jf19] = (
       +rho*Y[jp]*rate_eval.p_F19__Ne20
       )

    jac[jne20, jne19] = (
       +rho*Y[jn]*rate_eval.n_Ne19__Ne20
       )

    jac[jne20, jne20] = (
       -rate_eval.Ne20__n_Ne19
       -rate_eval.Ne20__p_F19
       -rate_eval.Ne20__He4_O16
       -rho*Y[jn]*rate_eval.n_Ne20__Ne21
       -rho*Y[jp]*rate_eval.p_Ne20__Na21
       -rho*Y[jhe4]*rate_eval.He4_Ne20__Mg24
       -rho*Y[jn]*rate_eval.n_Ne20__He4_O17
       -rho*Y[jp]*rate_eval.p_Ne20__He4_F17
       -rho*Y[jhe4]*rate_eval.He4_Ne20__n_Mg23
       -rho*Y[jhe4]*rate_eval.He4_Ne20__p_Na23
       -rho*Y[jhe4]*rate_eval.He4_Ne20__C12_C12
       -rho*Y[jc12]*rate_eval.C12_Ne20__n_S31
       -rho*Y[jc12]*rate_eval.C12_Ne20__p_P31
       -rho*Y[jc12]*rate_eval.C12_Ne20__He4_Si28
       )

    jac[jne20, jne21] = (
       +rate_eval.Ne21__n_Ne20
       )

    jac[jne20, jna21] = (
       +rate_eval.Na21__p_Ne20
       )

    jac[jne20, jna23] = (
       +rho*Y[jp]*rate_eval.p_Na23__He4_Ne20
       )

    jac[jne20, jmg23] = (
       +rho*Y[jn]*rate_eval.n_Mg23__He4_Ne20
       )

    jac[jne20, jmg24] = (
       +rate_eval.Mg24__He4_Ne20
       )

    jac[jne20, jsi28] = (
       +rho*Y[jhe4]*rate_eval.He4_Si28__C12_Ne20
       )

    jac[jne20, jp31] = (
       +rho*Y[jp]*rate_eval.p_P31__C12_Ne20
       )

    jac[jne20, js31] = (
       +rho*Y[jn]*rate_eval.n_S31__C12_Ne20
       )

    jac[jne21, jn] = (
       -rho*Y[jne21]*rate_eval.n_Ne21__He4_O18
       +rho*Y[jne20]*rate_eval.n_Ne20__Ne21
       +rho*Y[jna21]*rate_eval.n_Na21__p_Ne21
       +rho*Y[jmg24]*rate_eval.n_Mg24__He4_Ne21
       )

    jac[jne21, jp] = (
       -rho*Y[jne21]*rate_eval.p_Ne21__Na22
       -rho*Y[jne21]*rate_eval.p_Ne21__n_Na21
       -rho*Y[jne21]*rate_eval.p_Ne21__He4_F18
       +rho*Y[jna24]*rate_eval.p_Na24__He4_Ne21
       )

    jac[jne21, jhe4] = (
       -rho*Y[jne21]*rate_eval.He4_Ne21__Mg25
       -rho*Y[jne21]*rate_eval.He4_Ne21__n_Mg24
       -rho*Y[jne21]*rate_eval.He4_Ne21__p_Na24
       +rho*Y[jo17]*rate_eval.He4_O17__Ne21
       +rho*Y[jo18]*rate_eval.He4_O18__n_Ne21
       +rho*Y[jf18]*rate_eval.He4_F18__p_Ne21
       )

    jac[jne21, jo17] = (
       +rho*Y[jhe4]*rate_eval.He4_O17__Ne21
       )

    jac[jne21, jo18] = (
       +rho*Y[jhe4]*rate_eval.He4_O18__n_Ne21
       )

    jac[jne21, jf18] = (
       +rho*Y[jhe4]*rate_eval.He4_F18__p_Ne21
       )

    jac[jne21, jne20] = (
       +rho*Y[jn]*rate_eval.n_Ne20__Ne21
       )

    jac[jne21, jne21] = (
       -rate_eval.Ne21__n_Ne20
       -rate_eval.Ne21__He4_O17
       -rho*Y[jp]*rate_eval.p_Ne21__Na22
       -rho*Y[jhe4]*rate_eval.He4_Ne21__Mg25
       -rho*Y[jn]*rate_eval.n_Ne21__He4_O18
       -rho*Y[jp]*rate_eval.p_Ne21__n_Na21
       -rho*Y[jp]*rate_eval.p_Ne21__He4_F18
       -rho*Y[jhe4]*rate_eval.He4_Ne21__n_Mg24
       -rho*Y[jhe4]*rate_eval.He4_Ne21__p_Na24
       )

    jac[jne21, jna21] = (
       +rate_eval.Na21__Ne21__weak__wc12
       +rho*Y[jn]*rate_eval.n_Na21__p_Ne21
       )

    jac[jne21, jna22] = (
       +rate_eval.Na22__p_Ne21
       )

    jac[jne21, jna24] = (
       +rho*Y[jp]*rate_eval.p_Na24__He4_Ne21
       )

    jac[jne21, jmg24] = (
       +rho*Y[jn]*rate_eval.n_Mg24__He4_Ne21
       )

    jac[jne21, jmg25] = (
       +rate_eval.Mg25__He4_Ne21
       )

    jac[jna21, jn] = (
       -rho*Y[jna21]*rate_eval.n_Na21__Na22
       -rho*Y[jna21]*rate_eval.n_Na21__p_Ne21
       -rho*Y[jna21]*rate_eval.n_Na21__He4_F18
       )

    jac[jna21, jp] = (
       -rho*Y[jna21]*rate_eval.p_Na21__Mg22
       -rho*Y[jna21]*rate_eval.p_Na21__He4_Ne18
       +rho*Y[jne20]*rate_eval.p_Ne20__Na21
       +rho*Y[jne21]*rate_eval.p_Ne21__n_Na21
       +rho*Y[jmg24]*rate_eval.p_Mg24__He4_Na21
       )

    jac[jna21, jhe4] = (
       -rho*Y[jna21]*rate_eval.He4_Na21__Al25
       -rho*Y[jna21]*rate_eval.He4_Na21__p_Mg24
       +rho*Y[jf17]*rate_eval.He4_F17__Na21
       +rho*Y[jf18]*rate_eval.He4_F18__n_Na21
       +rho*Y[jne18]*rate_eval.He4_Ne18__p_Na21
       )

    jac[jna21, jf17] = (
       +rho*Y[jhe4]*rate_eval.He4_F17__Na21
       )

    jac[jna21, jf18] = (
       +rho*Y[jhe4]*rate_eval.He4_F18__n_Na21
       )

    jac[jna21, jne18] = (
       +rho*Y[jhe4]*rate_eval.He4_Ne18__p_Na21
       )

    jac[jna21, jne20] = (
       +rho*Y[jp]*rate_eval.p_Ne20__Na21
       )

    jac[jna21, jne21] = (
       +rho*Y[jp]*rate_eval.p_Ne21__n_Na21
       )

    jac[jna21, jna21] = (
       -rate_eval.Na21__Ne21__weak__wc12
       -rate_eval.Na21__p_Ne20
       -rate_eval.Na21__He4_F17
       -rho*Y[jn]*rate_eval.n_Na21__Na22
       -rho*Y[jp]*rate_eval.p_Na21__Mg22
       -rho*Y[jhe4]*rate_eval.He4_Na21__Al25
       -rho*Y[jn]*rate_eval.n_Na21__p_Ne21
       -rho*Y[jn]*rate_eval.n_Na21__He4_F18
       -rho*Y[jp]*rate_eval.p_Na21__He4_Ne18
       -rho*Y[jhe4]*rate_eval.He4_Na21__p_Mg24
       )

    jac[jna21, jna22] = (
       +rate_eval.Na22__n_Na21
       )

    jac[jna21, jmg22] = (
       +rate_eval.Mg22__p_Na21
       )

    jac[jna21, jmg24] = (
       +rho*Y[jp]*rate_eval.p_Mg24__He4_Na21
       )

    jac[jna21, jal25] = (
       +rate_eval.Al25__He4_Na21
       )

    jac[jna22, jn] = (
       -rho*Y[jna22]*rate_eval.n_Na22__Na23
       -rho*Y[jna22]*rate_eval.n_Na22__He4_F19
       +rho*Y[jna21]*rate_eval.n_Na21__Na22
       +rho*Y[jmg22]*rate_eval.n_Mg22__p_Na22
       +rho*Y[jal25]*rate_eval.n_Al25__He4_Na22
       )

    jac[jna22, jp] = (
       -rho*Y[jna22]*rate_eval.p_Na22__Mg23
       -rho*Y[jna22]*rate_eval.p_Na22__n_Mg22
       -rho*Y[jna22]*rate_eval.p_Na22__He4_Ne19
       +rho*Y[jne21]*rate_eval.p_Ne21__Na22
       +rho*Y[jmg25]*rate_eval.p_Mg25__He4_Na22
       )

    jac[jna22, jhe4] = (
       -rho*Y[jna22]*rate_eval.He4_Na22__Al26
       -rho*Y[jna22]*rate_eval.He4_Na22__n_Al25
       -rho*Y[jna22]*rate_eval.He4_Na22__p_Mg25
       +rho*Y[jf18]*rate_eval.He4_F18__Na22
       +rho*Y[jf19]*rate_eval.He4_F19__n_Na22
       +rho*Y[jne19]*rate_eval.He4_Ne19__p_Na22
       )

    jac[jna22, jf18] = (
       +rho*Y[jhe4]*rate_eval.He4_F18__Na22
       )

    jac[jna22, jf19] = (
       +rho*Y[jhe4]*rate_eval.He4_F19__n_Na22
       )

    jac[jna22, jne19] = (
       +rho*Y[jhe4]*rate_eval.He4_Ne19__p_Na22
       )

    jac[jna22, jne21] = (
       +rho*Y[jp]*rate_eval.p_Ne21__Na22
       )

    jac[jna22, jna21] = (
       +rho*Y[jn]*rate_eval.n_Na21__Na22
       )

    jac[jna22, jna22] = (
       -rate_eval.Na22__n_Na21
       -rate_eval.Na22__p_Ne21
       -rate_eval.Na22__He4_F18
       -rho*Y[jn]*rate_eval.n_Na22__Na23
       -rho*Y[jp]*rate_eval.p_Na22__Mg23
       -rho*Y[jhe4]*rate_eval.He4_Na22__Al26
       -rho*Y[jn]*rate_eval.n_Na22__He4_F19
       -rho*Y[jp]*rate_eval.p_Na22__n_Mg22
       -rho*Y[jp]*rate_eval.p_Na22__He4_Ne19
       -rho*Y[jhe4]*rate_eval.He4_Na22__n_Al25
       -rho*Y[jhe4]*rate_eval.He4_Na22__p_Mg25
       )

    jac[jna22, jna23] = (
       +rate_eval.Na23__n_Na22
       )

    jac[jna22, jmg22] = (
       +rate_eval.Mg22__Na22__weak__wc12
       +rho*Y[jn]*rate_eval.n_Mg22__p_Na22
       )

    jac[jna22, jmg23] = (
       +rate_eval.Mg23__p_Na22
       )

    jac[jna22, jmg25] = (
       +rho*Y[jp]*rate_eval.p_Mg25__He4_Na22
       )

    jac[jna22, jal25] = (
       +rho*Y[jn]*rate_eval.n_Al25__He4_Na22
       )

    jac[jna22, jal26] = (
       +rate_eval.Al26__He4_Na22
       )

    jac[jna23, jn] = (
       -rho*Y[jna23]*rate_eval.n_Na23__Na24
       +rho*Y[jna22]*rate_eval.n_Na22__Na23
       +rho*Y[jmg23]*rate_eval.n_Mg23__p_Na23
       +rho*Y[jal26]*rate_eval.n_Al26__He4_Na23
       )

    jac[jna23, jp] = (
       -rho*Y[jna23]*rate_eval.p_Na23__Mg24
       -rho*Y[jna23]*rate_eval.p_Na23__n_Mg23
       -rho*Y[jna23]*rate_eval.p_Na23__He4_Ne20
       -rho*Y[jna23]*rate_eval.p_Na23__C12_C12
       +rho*Y[jmg26]*rate_eval.p_Mg26__He4_Na23
       )

    jac[jna23, jhe4] = (
       -rho*Y[jna23]*rate_eval.He4_Na23__Al27
       -rho*Y[jna23]*rate_eval.He4_Na23__n_Al26
       -rho*Y[jna23]*rate_eval.He4_Na23__p_Mg26
       +rho*Y[jf19]*rate_eval.He4_F19__Na23
       +rho*Y[jne20]*rate_eval.He4_Ne20__p_Na23
       )

    jac[jna23, jc12] = (
       +5.00000000000000e-01*rho*2*Y[jc12]*rate_eval.C12_C12__p_Na23
       )

    jac[jna23, jf19] = (
       +rho*Y[jhe4]*rate_eval.He4_F19__Na23
       )

    jac[jna23, jne20] = (
       +rho*Y[jhe4]*rate_eval.He4_Ne20__p_Na23
       )

    jac[jna23, jna22] = (
       +rho*Y[jn]*rate_eval.n_Na22__Na23
       )

    jac[jna23, jna23] = (
       -rate_eval.Na23__n_Na22
       -rate_eval.Na23__He4_F19
       -rho*Y[jn]*rate_eval.n_Na23__Na24
       -rho*Y[jp]*rate_eval.p_Na23__Mg24
       -rho*Y[jhe4]*rate_eval.He4_Na23__Al27
       -rho*Y[jp]*rate_eval.p_Na23__n_Mg23
       -rho*Y[jp]*rate_eval.p_Na23__He4_Ne20
       -rho*Y[jp]*rate_eval.p_Na23__C12_C12
       -rho*Y[jhe4]*rate_eval.He4_Na23__n_Al26
       -rho*Y[jhe4]*rate_eval.He4_Na23__p_Mg26
       )

    jac[jna23, jna24] = (
       +rate_eval.Na24__n_Na23
       )

    jac[jna23, jmg23] = (
       +rate_eval.Mg23__Na23__weak__wc12
       +rho*Y[jn]*rate_eval.n_Mg23__p_Na23
       )

    jac[jna23, jmg24] = (
       +rate_eval.Mg24__p_Na23
       )

    jac[jna23, jmg26] = (
       +rho*Y[jp]*rate_eval.p_Mg26__He4_Na23
       )

    jac[jna23, jal26] = (
       +rho*Y[jn]*rate_eval.n_Al26__He4_Na23
       )

    jac[jna23, jal27] = (
       +rate_eval.Al27__He4_Na23
       )

    jac[jna24, jn] = (
       +rho*Y[jna23]*rate_eval.n_Na23__Na24
       +rho*Y[jmg24]*rate_eval.n_Mg24__p_Na24
       +rho*Y[jal27]*rate_eval.n_Al27__He4_Na24
       )

    jac[jna24, jp] = (
       -rho*Y[jna24]*rate_eval.p_Na24__Mg25
       -rho*Y[jna24]*rate_eval.p_Na24__n_Mg24
       -rho*Y[jna24]*rate_eval.p_Na24__He4_Ne21
       )

    jac[jna24, jhe4] = (
       -rho*Y[jna24]*rate_eval.He4_Na24__Al28
       -rho*Y[jna24]*rate_eval.He4_Na24__n_Al27
       +rho*Y[jne21]*rate_eval.He4_Ne21__p_Na24
       )

    jac[jna24, jne21] = (
       +rho*Y[jhe4]*rate_eval.He4_Ne21__p_Na24
       )

    jac[jna24, jna23] = (
       +rho*Y[jn]*rate_eval.n_Na23__Na24
       )

    jac[jna24, jna24] = (
       -rate_eval.Na24__Mg24__weak__wc12
       -rate_eval.Na24__n_Na23
       -rho*Y[jp]*rate_eval.p_Na24__Mg25
       -rho*Y[jhe4]*rate_eval.He4_Na24__Al28
       -rho*Y[jp]*rate_eval.p_Na24__n_Mg24
       -rho*Y[jp]*rate_eval.p_Na24__He4_Ne21
       -rho*Y[jhe4]*rate_eval.He4_Na24__n_Al27
       )

    jac[jna24, jmg24] = (
       +rho*Y[jn]*rate_eval.n_Mg24__p_Na24
       )

    jac[jna24, jmg25] = (
       +rate_eval.Mg25__p_Na24
       )

    jac[jna24, jal27] = (
       +rho*Y[jn]*rate_eval.n_Al27__He4_Na24
       )

    jac[jna24, jal28] = (
       +rate_eval.Al28__He4_Na24
       )

    jac[jmg22, jn] = (
       -rho*Y[jmg22]*rate_eval.n_Mg22__Mg23
       -rho*Y[jmg22]*rate_eval.n_Mg22__p_Na22
       -rho*Y[jmg22]*rate_eval.n_Mg22__He4_Ne19
       )

    jac[jmg22, jp] = (
       +rho*Y[jna21]*rate_eval.p_Na21__Mg22
       +rho*Y[jna22]*rate_eval.p_Na22__n_Mg22
       +rho*Y[jal25]*rate_eval.p_Al25__He4_Mg22
       )

    jac[jmg22, jhe4] = (
       -rho*Y[jmg22]*rate_eval.He4_Mg22__Si26
       -rho*Y[jmg22]*rate_eval.He4_Mg22__p_Al25
       +rho*Y[jne18]*rate_eval.He4_Ne18__Mg22
       +rho*Y[jne19]*rate_eval.He4_Ne19__n_Mg22
       )

    jac[jmg22, jne18] = (
       +rho*Y[jhe4]*rate_eval.He4_Ne18__Mg22
       )

    jac[jmg22, jne19] = (
       +rho*Y[jhe4]*rate_eval.He4_Ne19__n_Mg22
       )

    jac[jmg22, jna21] = (
       +rho*Y[jp]*rate_eval.p_Na21__Mg22
       )

    jac[jmg22, jna22] = (
       +rho*Y[jp]*rate_eval.p_Na22__n_Mg22
       )

    jac[jmg22, jmg22] = (
       -rate_eval.Mg22__Na22__weak__wc12
       -rate_eval.Mg22__p_Na21
       -rate_eval.Mg22__He4_Ne18
       -rho*Y[jn]*rate_eval.n_Mg22__Mg23
       -rho*Y[jhe4]*rate_eval.He4_Mg22__Si26
       -rho*Y[jn]*rate_eval.n_Mg22__p_Na22
       -rho*Y[jn]*rate_eval.n_Mg22__He4_Ne19
       -rho*Y[jhe4]*rate_eval.He4_Mg22__p_Al25
       )

    jac[jmg22, jmg23] = (
       +rate_eval.Mg23__n_Mg22
       )

    jac[jmg22, jal25] = (
       +rho*Y[jp]*rate_eval.p_Al25__He4_Mg22
       )

    jac[jmg22, jsi26] = (
       +rate_eval.Si26__He4_Mg22
       )

    jac[jmg23, jn] = (
       -rho*Y[jmg23]*rate_eval.n_Mg23__Mg24
       -rho*Y[jmg23]*rate_eval.n_Mg23__p_Na23
       -rho*Y[jmg23]*rate_eval.n_Mg23__He4_Ne20
       -rho*Y[jmg23]*rate_eval.n_Mg23__C12_C12
       +rho*Y[jmg22]*rate_eval.n_Mg22__Mg23
       +rho*Y[jsi26]*rate_eval.n_Si26__He4_Mg23
       )

    jac[jmg23, jp] = (
       +rho*Y[jna22]*rate_eval.p_Na22__Mg23
       +rho*Y[jna23]*rate_eval.p_Na23__n_Mg23
       +rho*Y[jal26]*rate_eval.p_Al26__He4_Mg23
       )

    jac[jmg23, jhe4] = (
       -rho*Y[jmg23]*rate_eval.He4_Mg23__n_Si26
       -rho*Y[jmg23]*rate_eval.He4_Mg23__p_Al26
       +rho*Y[jne19]*rate_eval.He4_Ne19__Mg23
       +rho*Y[jne20]*rate_eval.He4_Ne20__n_Mg23
       )

    jac[jmg23, jc12] = (
       +5.00000000000000e-01*rho*2*Y[jc12]*rate_eval.C12_C12__n_Mg23
       )

    jac[jmg23, jne19] = (
       +rho*Y[jhe4]*rate_eval.He4_Ne19__Mg23
       )

    jac[jmg23, jne20] = (
       +rho*Y[jhe4]*rate_eval.He4_Ne20__n_Mg23
       )

    jac[jmg23, jna22] = (
       +rho*Y[jp]*rate_eval.p_Na22__Mg23
       )

    jac[jmg23, jna23] = (
       +rho*Y[jp]*rate_eval.p_Na23__n_Mg23
       )

    jac[jmg23, jmg22] = (
       +rho*Y[jn]*rate_eval.n_Mg22__Mg23
       )

    jac[jmg23, jmg23] = (
       -rate_eval.Mg23__Na23__weak__wc12
       -rate_eval.Mg23__n_Mg22
       -rate_eval.Mg23__p_Na22
       -rate_eval.Mg23__He4_Ne19
       -rho*Y[jn]*rate_eval.n_Mg23__Mg24
       -rho*Y[jn]*rate_eval.n_Mg23__p_Na23
       -rho*Y[jn]*rate_eval.n_Mg23__He4_Ne20
       -rho*Y[jn]*rate_eval.n_Mg23__C12_C12
       -rho*Y[jhe4]*rate_eval.He4_Mg23__n_Si26
       -rho*Y[jhe4]*rate_eval.He4_Mg23__p_Al26
       )

    jac[jmg23, jmg24] = (
       +rate_eval.Mg24__n_Mg23
       )

    jac[jmg23, jal26] = (
       +rho*Y[jp]*rate_eval.p_Al26__He4_Mg23
       )

    jac[jmg23, jsi26] = (
       +rho*Y[jn]*rate_eval.n_Si26__He4_Mg23
       )

    jac[jmg24, jn] = (
       -rho*Y[jmg24]*rate_eval.n_Mg24__Mg25
       -rho*Y[jmg24]*rate_eval.n_Mg24__p_Na24
       -rho*Y[jmg24]*rate_eval.n_Mg24__He4_Ne21
       +rho*Y[jmg23]*rate_eval.n_Mg23__Mg24
       )

    jac[jmg24, jp] = (
       -rho*Y[jmg24]*rate_eval.p_Mg24__Al25
       -rho*Y[jmg24]*rate_eval.p_Mg24__He4_Na21
       +rho*Y[jna23]*rate_eval.p_Na23__Mg24
       +rho*Y[jna24]*rate_eval.p_Na24__n_Mg24
       +rho*Y[jal27]*rate_eval.p_Al27__He4_Mg24
       )

    jac[jmg24, jhe4] = (
       -rho*Y[jmg24]*rate_eval.He4_Mg24__Si28
       -rho*Y[jmg24]*rate_eval.He4_Mg24__p_Al27
       -rho*Y[jmg24]*rate_eval.He4_Mg24__C12_O16
       +rho*Y[jne20]*rate_eval.He4_Ne20__Mg24
       +rho*Y[jne21]*rate_eval.He4_Ne21__n_Mg24
       +rho*Y[jna21]*rate_eval.He4_Na21__p_Mg24
       )

    jac[jmg24, jc12] = (
       +rho*Y[jo16]*rate_eval.C12_O16__He4_Mg24
       )

    jac[jmg24, jo16] = (
       +rho*Y[jc12]*rate_eval.C12_O16__He4_Mg24
       )

    jac[jmg24, jne20] = (
       +rho*Y[jhe4]*rate_eval.He4_Ne20__Mg24
       )

    jac[jmg24, jne21] = (
       +rho*Y[jhe4]*rate_eval.He4_Ne21__n_Mg24
       )

    jac[jmg24, jna21] = (
       +rho*Y[jhe4]*rate_eval.He4_Na21__p_Mg24
       )

    jac[jmg24, jna23] = (
       +rho*Y[jp]*rate_eval.p_Na23__Mg24
       )

    jac[jmg24, jna24] = (
       +rate_eval.Na24__Mg24__weak__wc12
       +rho*Y[jp]*rate_eval.p_Na24__n_Mg24
       )

    jac[jmg24, jmg23] = (
       +rho*Y[jn]*rate_eval.n_Mg23__Mg24
       )

    jac[jmg24, jmg24] = (
       -rate_eval.Mg24__n_Mg23
       -rate_eval.Mg24__p_Na23
       -rate_eval.Mg24__He4_Ne20
       -rho*Y[jn]*rate_eval.n_Mg24__Mg25
       -rho*Y[jp]*rate_eval.p_Mg24__Al25
       -rho*Y[jhe4]*rate_eval.He4_Mg24__Si28
       -rho*Y[jn]*rate_eval.n_Mg24__p_Na24
       -rho*Y[jn]*rate_eval.n_Mg24__He4_Ne21
       -rho*Y[jp]*rate_eval.p_Mg24__He4_Na21
       -rho*Y[jhe4]*rate_eval.He4_Mg24__p_Al27
       -rho*Y[jhe4]*rate_eval.He4_Mg24__C12_O16
       )

    jac[jmg24, jmg25] = (
       +rate_eval.Mg25__n_Mg24
       )

    jac[jmg24, jal25] = (
       +rate_eval.Al25__p_Mg24
       )

    jac[jmg24, jal27] = (
       +rho*Y[jp]*rate_eval.p_Al27__He4_Mg24
       )

    jac[jmg24, jsi28] = (
       +rate_eval.Si28__He4_Mg24
       )

    jac[jmg25, jn] = (
       -rho*Y[jmg25]*rate_eval.n_Mg25__Mg26
       +rho*Y[jmg24]*rate_eval.n_Mg24__Mg25
       +rho*Y[jal25]*rate_eval.n_Al25__p_Mg25
       +rho*Y[jsi28]*rate_eval.n_Si28__He4_Mg25
       )

    jac[jmg25, jp] = (
       -rho*Y[jmg25]*rate_eval.p_Mg25__Al26
       -rho*Y[jmg25]*rate_eval.p_Mg25__n_Al25
       -rho*Y[jmg25]*rate_eval.p_Mg25__He4_Na22
       +rho*Y[jna24]*rate_eval.p_Na24__Mg25
       +rho*Y[jal28]*rate_eval.p_Al28__He4_Mg25
       )

    jac[jmg25, jhe4] = (
       -rho*Y[jmg25]*rate_eval.He4_Mg25__Si29
       -rho*Y[jmg25]*rate_eval.He4_Mg25__n_Si28
       -rho*Y[jmg25]*rate_eval.He4_Mg25__p_Al28
       +rho*Y[jne21]*rate_eval.He4_Ne21__Mg25
       +rho*Y[jna22]*rate_eval.He4_Na22__p_Mg25
       )

    jac[jmg25, jne21] = (
       +rho*Y[jhe4]*rate_eval.He4_Ne21__Mg25
       )

    jac[jmg25, jna22] = (
       +rho*Y[jhe4]*rate_eval.He4_Na22__p_Mg25
       )

    jac[jmg25, jna24] = (
       +rho*Y[jp]*rate_eval.p_Na24__Mg25
       )

    jac[jmg25, jmg24] = (
       +rho*Y[jn]*rate_eval.n_Mg24__Mg25
       )

    jac[jmg25, jmg25] = (
       -rate_eval.Mg25__n_Mg24
       -rate_eval.Mg25__p_Na24
       -rate_eval.Mg25__He4_Ne21
       -rho*Y[jn]*rate_eval.n_Mg25__Mg26
       -rho*Y[jp]*rate_eval.p_Mg25__Al26
       -rho*Y[jhe4]*rate_eval.He4_Mg25__Si29
       -rho*Y[jp]*rate_eval.p_Mg25__n_Al25
       -rho*Y[jp]*rate_eval.p_Mg25__He4_Na22
       -rho*Y[jhe4]*rate_eval.He4_Mg25__n_Si28
       -rho*Y[jhe4]*rate_eval.He4_Mg25__p_Al28
       )

    jac[jmg25, jmg26] = (
       +rate_eval.Mg26__n_Mg25
       )

    jac[jmg25, jal25] = (
       +rate_eval.Al25__Mg25__weak__wc12
       +rho*Y[jn]*rate_eval.n_Al25__p_Mg25
       )

    jac[jmg25, jal26] = (
       +rate_eval.Al26__p_Mg25
       )

    jac[jmg25, jal28] = (
       +rho*Y[jp]*rate_eval.p_Al28__He4_Mg25
       )

    jac[jmg25, jsi28] = (
       +rho*Y[jn]*rate_eval.n_Si28__He4_Mg25
       )

    jac[jmg25, jsi29] = (
       +rate_eval.Si29__He4_Mg25
       )

    jac[jmg26, jn] = (
       +rho*Y[jmg25]*rate_eval.n_Mg25__Mg26
       +rho*Y[jal26]*rate_eval.n_Al26__p_Mg26
       +rho*Y[jsi29]*rate_eval.n_Si29__He4_Mg26
       )

    jac[jmg26, jp] = (
       -rho*Y[jmg26]*rate_eval.p_Mg26__Al27
       -rho*Y[jmg26]*rate_eval.p_Mg26__n_Al26
       -rho*Y[jmg26]*rate_eval.p_Mg26__He4_Na23
       )

    jac[jmg26, jhe4] = (
       -rho*Y[jmg26]*rate_eval.He4_Mg26__Si30
       -rho*Y[jmg26]*rate_eval.He4_Mg26__n_Si29
       +rho*Y[jna23]*rate_eval.He4_Na23__p_Mg26
       )

    jac[jmg26, jna23] = (
       +rho*Y[jhe4]*rate_eval.He4_Na23__p_Mg26
       )

    jac[jmg26, jmg25] = (
       +rho*Y[jn]*rate_eval.n_Mg25__Mg26
       )

    jac[jmg26, jmg26] = (
       -rate_eval.Mg26__n_Mg25
       -rho*Y[jp]*rate_eval.p_Mg26__Al27
       -rho*Y[jhe4]*rate_eval.He4_Mg26__Si30
       -rho*Y[jp]*rate_eval.p_Mg26__n_Al26
       -rho*Y[jp]*rate_eval.p_Mg26__He4_Na23
       -rho*Y[jhe4]*rate_eval.He4_Mg26__n_Si29
       )

    jac[jmg26, jal26] = (
       +rate_eval.Al26__Mg26__weak__wc12
       +rho*Y[jn]*rate_eval.n_Al26__p_Mg26
       )

    jac[jmg26, jal27] = (
       +rate_eval.Al27__p_Mg26
       )

    jac[jmg26, jsi29] = (
       +rho*Y[jn]*rate_eval.n_Si29__He4_Mg26
       )

    jac[jmg26, jsi30] = (
       +rate_eval.Si30__He4_Mg26
       )

    jac[jal25, jn] = (
       -rho*Y[jal25]*rate_eval.n_Al25__Al26
       -rho*Y[jal25]*rate_eval.n_Al25__p_Mg25
       -rho*Y[jal25]*rate_eval.n_Al25__He4_Na22
       )

    jac[jal25, jp] = (
       -rho*Y[jal25]*rate_eval.p_Al25__Si26
       -rho*Y[jal25]*rate_eval.p_Al25__He4_Mg22
       +rho*Y[jmg24]*rate_eval.p_Mg24__Al25
       +rho*Y[jmg25]*rate_eval.p_Mg25__n_Al25
       +rho*Y[jsi28]*rate_eval.p_Si28__He4_Al25
       )

    jac[jal25, jhe4] = (
       -rho*Y[jal25]*rate_eval.He4_Al25__P29
       -rho*Y[jal25]*rate_eval.He4_Al25__p_Si28
       +rho*Y[jna21]*rate_eval.He4_Na21__Al25
       +rho*Y[jna22]*rate_eval.He4_Na22__n_Al25
       +rho*Y[jmg22]*rate_eval.He4_Mg22__p_Al25
       )

    jac[jal25, jna21] = (
       +rho*Y[jhe4]*rate_eval.He4_Na21__Al25
       )

    jac[jal25, jna22] = (
       +rho*Y[jhe4]*rate_eval.He4_Na22__n_Al25
       )

    jac[jal25, jmg22] = (
       +rho*Y[jhe4]*rate_eval.He4_Mg22__p_Al25
       )

    jac[jal25, jmg24] = (
       +rho*Y[jp]*rate_eval.p_Mg24__Al25
       )

    jac[jal25, jmg25] = (
       +rho*Y[jp]*rate_eval.p_Mg25__n_Al25
       )

    jac[jal25, jal25] = (
       -rate_eval.Al25__Mg25__weak__wc12
       -rate_eval.Al25__p_Mg24
       -rate_eval.Al25__He4_Na21
       -rho*Y[jn]*rate_eval.n_Al25__Al26
       -rho*Y[jp]*rate_eval.p_Al25__Si26
       -rho*Y[jhe4]*rate_eval.He4_Al25__P29
       -rho*Y[jn]*rate_eval.n_Al25__p_Mg25
       -rho*Y[jn]*rate_eval.n_Al25__He4_Na22
       -rho*Y[jp]*rate_eval.p_Al25__He4_Mg22
       -rho*Y[jhe4]*rate_eval.He4_Al25__p_Si28
       )

    jac[jal25, jal26] = (
       +rate_eval.Al26__n_Al25
       )

    jac[jal25, jsi26] = (
       +rate_eval.Si26__p_Al25
       )

    jac[jal25, jsi28] = (
       +rho*Y[jp]*rate_eval.p_Si28__He4_Al25
       )

    jac[jal25, jp29] = (
       +rate_eval.P29__He4_Al25
       )

    jac[jal26, jn] = (
       -rho*Y[jal26]*rate_eval.n_Al26__Al27
       -rho*Y[jal26]*rate_eval.n_Al26__p_Mg26
       -rho*Y[jal26]*rate_eval.n_Al26__He4_Na23
       +rho*Y[jal25]*rate_eval.n_Al25__Al26
       +rho*Y[jsi26]*rate_eval.n_Si26__p_Al26
       +rho*Y[jp29]*rate_eval.n_P29__He4_Al26
       )

    jac[jal26, jp] = (
       -rho*Y[jal26]*rate_eval.p_Al26__n_Si26
       -rho*Y[jal26]*rate_eval.p_Al26__He4_Mg23
       +rho*Y[jmg25]*rate_eval.p_Mg25__Al26
       +rho*Y[jmg26]*rate_eval.p_Mg26__n_Al26
       +rho*Y[jsi29]*rate_eval.p_Si29__He4_Al26
       )

    jac[jal26, jhe4] = (
       -rho*Y[jal26]*rate_eval.He4_Al26__P30
       -rho*Y[jal26]*rate_eval.He4_Al26__n_P29
       -rho*Y[jal26]*rate_eval.He4_Al26__p_Si29
       +rho*Y[jna22]*rate_eval.He4_Na22__Al26
       +rho*Y[jna23]*rate_eval.He4_Na23__n_Al26
       +rho*Y[jmg23]*rate_eval.He4_Mg23__p_Al26
       )

    jac[jal26, jna22] = (
       +rho*Y[jhe4]*rate_eval.He4_Na22__Al26
       )

    jac[jal26, jna23] = (
       +rho*Y[jhe4]*rate_eval.He4_Na23__n_Al26
       )

    jac[jal26, jmg23] = (
       +rho*Y[jhe4]*rate_eval.He4_Mg23__p_Al26
       )

    jac[jal26, jmg25] = (
       +rho*Y[jp]*rate_eval.p_Mg25__Al26
       )

    jac[jal26, jmg26] = (
       +rho*Y[jp]*rate_eval.p_Mg26__n_Al26
       )

    jac[jal26, jal25] = (
       +rho*Y[jn]*rate_eval.n_Al25__Al26
       )

    jac[jal26, jal26] = (
       -rate_eval.Al26__Mg26__weak__wc12
       -rate_eval.Al26__n_Al25
       -rate_eval.Al26__p_Mg25
       -rate_eval.Al26__He4_Na22
       -rho*Y[jn]*rate_eval.n_Al26__Al27
       -rho*Y[jhe4]*rate_eval.He4_Al26__P30
       -rho*Y[jn]*rate_eval.n_Al26__p_Mg26
       -rho*Y[jn]*rate_eval.n_Al26__He4_Na23
       -rho*Y[jp]*rate_eval.p_Al26__n_Si26
       -rho*Y[jp]*rate_eval.p_Al26__He4_Mg23
       -rho*Y[jhe4]*rate_eval.He4_Al26__n_P29
       -rho*Y[jhe4]*rate_eval.He4_Al26__p_Si29
       )

    jac[jal26, jal27] = (
       +rate_eval.Al27__n_Al26
       )

    jac[jal26, jsi26] = (
       +rate_eval.Si26__Al26__weak__wc12
       +rho*Y[jn]*rate_eval.n_Si26__p_Al26
       )

    jac[jal26, jsi29] = (
       +rho*Y[jp]*rate_eval.p_Si29__He4_Al26
       )

    jac[jal26, jp29] = (
       +rho*Y[jn]*rate_eval.n_P29__He4_Al26
       )

    jac[jal26, jp30] = (
       +rate_eval.P30__He4_Al26
       )

    jac[jal27, jn] = (
       -rho*Y[jal27]*rate_eval.n_Al27__Al28
       -rho*Y[jal27]*rate_eval.n_Al27__He4_Na24
       +rho*Y[jal26]*rate_eval.n_Al26__Al27
       +rho*Y[jp30]*rate_eval.n_P30__He4_Al27
       )

    jac[jal27, jp] = (
       -rho*Y[jal27]*rate_eval.p_Al27__Si28
       -rho*Y[jal27]*rate_eval.p_Al27__He4_Mg24
       -rho*Y[jal27]*rate_eval.p_Al27__C12_O16
       +rho*Y[jmg26]*rate_eval.p_Mg26__Al27
       +rho*Y[jsi30]*rate_eval.p_Si30__He4_Al27
       )

    jac[jal27, jhe4] = (
       -rho*Y[jal27]*rate_eval.He4_Al27__P31
       -rho*Y[jal27]*rate_eval.He4_Al27__n_P30
       -rho*Y[jal27]*rate_eval.He4_Al27__p_Si30
       +rho*Y[jna23]*rate_eval.He4_Na23__Al27
       +rho*Y[jna24]*rate_eval.He4_Na24__n_Al27
       +rho*Y[jmg24]*rate_eval.He4_Mg24__p_Al27
       )

    jac[jal27, jc12] = (
       +rho*Y[jo16]*rate_eval.C12_O16__p_Al27
       )

    jac[jal27, jo16] = (
       +rho*Y[jc12]*rate_eval.C12_O16__p_Al27
       )

    jac[jal27, jna23] = (
       +rho*Y[jhe4]*rate_eval.He4_Na23__Al27
       )

    jac[jal27, jna24] = (
       +rho*Y[jhe4]*rate_eval.He4_Na24__n_Al27
       )

    jac[jal27, jmg24] = (
       +rho*Y[jhe4]*rate_eval.He4_Mg24__p_Al27
       )

    jac[jal27, jmg26] = (
       +rho*Y[jp]*rate_eval.p_Mg26__Al27
       )

    jac[jal27, jal26] = (
       +rho*Y[jn]*rate_eval.n_Al26__Al27
       )

    jac[jal27, jal27] = (
       -rate_eval.Al27__n_Al26
       -rate_eval.Al27__p_Mg26
       -rate_eval.Al27__He4_Na23
       -rho*Y[jn]*rate_eval.n_Al27__Al28
       -rho*Y[jp]*rate_eval.p_Al27__Si28
       -rho*Y[jhe4]*rate_eval.He4_Al27__P31
       -rho*Y[jn]*rate_eval.n_Al27__He4_Na24
       -rho*Y[jp]*rate_eval.p_Al27__He4_Mg24
       -rho*Y[jp]*rate_eval.p_Al27__C12_O16
       -rho*Y[jhe4]*rate_eval.He4_Al27__n_P30
       -rho*Y[jhe4]*rate_eval.He4_Al27__p_Si30
       )

    jac[jal27, jal28] = (
       +rate_eval.Al28__n_Al27
       )

    jac[jal27, jsi28] = (
       +rate_eval.Si28__p_Al27
       )

    jac[jal27, jsi30] = (
       +rho*Y[jp]*rate_eval.p_Si30__He4_Al27
       )

    jac[jal27, jp30] = (
       +rho*Y[jn]*rate_eval.n_P30__He4_Al27
       )

    jac[jal27, jp31] = (
       +rate_eval.P31__He4_Al27
       )

    jac[jal28, jn] = (
       +rho*Y[jal27]*rate_eval.n_Al27__Al28
       +rho*Y[jsi28]*rate_eval.n_Si28__p_Al28
       +rho*Y[jp31]*rate_eval.n_P31__He4_Al28
       )

    jac[jal28, jp] = (
       -rho*Y[jal28]*rate_eval.p_Al28__Si29
       -rho*Y[jal28]*rate_eval.p_Al28__n_Si28
       -rho*Y[jal28]*rate_eval.p_Al28__He4_Mg25
       )

    jac[jal28, jhe4] = (
       -rho*Y[jal28]*rate_eval.He4_Al28__P32
       -rho*Y[jal28]*rate_eval.He4_Al28__n_P31
       +rho*Y[jna24]*rate_eval.He4_Na24__Al28
       +rho*Y[jmg25]*rate_eval.He4_Mg25__p_Al28
       )

    jac[jal28, jna24] = (
       +rho*Y[jhe4]*rate_eval.He4_Na24__Al28
       )

    jac[jal28, jmg25] = (
       +rho*Y[jhe4]*rate_eval.He4_Mg25__p_Al28
       )

    jac[jal28, jal27] = (
       +rho*Y[jn]*rate_eval.n_Al27__Al28
       )

    jac[jal28, jal28] = (
       -rate_eval.Al28__Si28__weak__wc12
       -rate_eval.Al28__n_Al27
       -rate_eval.Al28__He4_Na24
       -rho*Y[jp]*rate_eval.p_Al28__Si29
       -rho*Y[jhe4]*rate_eval.He4_Al28__P32
       -rho*Y[jp]*rate_eval.p_Al28__n_Si28
       -rho*Y[jp]*rate_eval.p_Al28__He4_Mg25
       -rho*Y[jhe4]*rate_eval.He4_Al28__n_P31
       )

    jac[jal28, jsi28] = (
       +rho*Y[jn]*rate_eval.n_Si28__p_Al28
       )

    jac[jal28, jsi29] = (
       +rate_eval.Si29__p_Al28
       )

    jac[jal28, jp31] = (
       +rho*Y[jn]*rate_eval.n_P31__He4_Al28
       )

    jac[jal28, jp32] = (
       +rate_eval.P32__He4_Al28
       )

    jac[jsi26, jn] = (
       -rho*Y[jsi26]*rate_eval.n_Si26__p_Al26
       -rho*Y[jsi26]*rate_eval.n_Si26__He4_Mg23
       )

    jac[jsi26, jp] = (
       +rho*Y[jal25]*rate_eval.p_Al25__Si26
       +rho*Y[jal26]*rate_eval.p_Al26__n_Si26
       +rho*Y[jp29]*rate_eval.p_P29__He4_Si26
       )

    jac[jsi26, jhe4] = (
       -rho*Y[jsi26]*rate_eval.He4_Si26__p_P29
       +rho*Y[jmg22]*rate_eval.He4_Mg22__Si26
       +rho*Y[jmg23]*rate_eval.He4_Mg23__n_Si26
       )

    jac[jsi26, jmg22] = (
       +rho*Y[jhe4]*rate_eval.He4_Mg22__Si26
       )

    jac[jsi26, jmg23] = (
       +rho*Y[jhe4]*rate_eval.He4_Mg23__n_Si26
       )

    jac[jsi26, jal25] = (
       +rho*Y[jp]*rate_eval.p_Al25__Si26
       )

    jac[jsi26, jal26] = (
       +rho*Y[jp]*rate_eval.p_Al26__n_Si26
       )

    jac[jsi26, jsi26] = (
       -rate_eval.Si26__Al26__weak__wc12
       -rate_eval.Si26__p_Al25
       -rate_eval.Si26__He4_Mg22
       -rho*Y[jn]*rate_eval.n_Si26__p_Al26
       -rho*Y[jn]*rate_eval.n_Si26__He4_Mg23
       -rho*Y[jhe4]*rate_eval.He4_Si26__p_P29
       )

    jac[jsi26, jp29] = (
       +rho*Y[jp]*rate_eval.p_P29__He4_Si26
       )

    jac[jsi28, jn] = (
       -rho*Y[jsi28]*rate_eval.n_Si28__Si29
       -rho*Y[jsi28]*rate_eval.n_Si28__p_Al28
       -rho*Y[jsi28]*rate_eval.n_Si28__He4_Mg25
       +rho*Y[js31]*rate_eval.n_S31__He4_Si28
       )

    jac[jsi28, jp] = (
       -rho*Y[jsi28]*rate_eval.p_Si28__P29
       -rho*Y[jsi28]*rate_eval.p_Si28__He4_Al25
       +rho*Y[jal27]*rate_eval.p_Al27__Si28
       +rho*Y[jal28]*rate_eval.p_Al28__n_Si28
       +rho*Y[jp31]*rate_eval.p_P31__He4_Si28
       )

    jac[jsi28, jhe4] = (
       -rho*Y[jsi28]*rate_eval.He4_Si28__S32
       -rho*Y[jsi28]*rate_eval.He4_Si28__n_S31
       -rho*Y[jsi28]*rate_eval.He4_Si28__p_P31
       -rho*Y[jsi28]*rate_eval.He4_Si28__C12_Ne20
       -rho*Y[jsi28]*rate_eval.He4_Si28__O16_O16
       +rho*Y[jmg24]*rate_eval.He4_Mg24__Si28
       +rho*Y[jmg25]*rate_eval.He4_Mg25__n_Si28
       +rho*Y[jal25]*rate_eval.He4_Al25__p_Si28
       )

    jac[jsi28, jc12] = (
       +rho*Y[jne20]*rate_eval.C12_Ne20__He4_Si28
       )

    jac[jsi28, jo16] = (
       +5.00000000000000e-01*rho*2*Y[jo16]*rate_eval.O16_O16__He4_Si28
       )

    jac[jsi28, jne20] = (
       +rho*Y[jc12]*rate_eval.C12_Ne20__He4_Si28
       )

    jac[jsi28, jmg24] = (
       +rho*Y[jhe4]*rate_eval.He4_Mg24__Si28
       )

    jac[jsi28, jmg25] = (
       +rho*Y[jhe4]*rate_eval.He4_Mg25__n_Si28
       )

    jac[jsi28, jal25] = (
       +rho*Y[jhe4]*rate_eval.He4_Al25__p_Si28
       )

    jac[jsi28, jal27] = (
       +rho*Y[jp]*rate_eval.p_Al27__Si28
       )

    jac[jsi28, jal28] = (
       +rate_eval.Al28__Si28__weak__wc12
       +rho*Y[jp]*rate_eval.p_Al28__n_Si28
       )

    jac[jsi28, jsi28] = (
       -rate_eval.Si28__p_Al27
       -rate_eval.Si28__He4_Mg24
       -rho*Y[jn]*rate_eval.n_Si28__Si29
       -rho*Y[jp]*rate_eval.p_Si28__P29
       -rho*Y[jhe4]*rate_eval.He4_Si28__S32
       -rho*Y[jn]*rate_eval.n_Si28__p_Al28
       -rho*Y[jn]*rate_eval.n_Si28__He4_Mg25
       -rho*Y[jp]*rate_eval.p_Si28__He4_Al25
       -rho*Y[jhe4]*rate_eval.He4_Si28__n_S31
       -rho*Y[jhe4]*rate_eval.He4_Si28__p_P31
       -rho*Y[jhe4]*rate_eval.He4_Si28__C12_Ne20
       -rho*Y[jhe4]*rate_eval.He4_Si28__O16_O16
       )

    jac[jsi28, jsi29] = (
       +rate_eval.Si29__n_Si28
       )

    jac[jsi28, jp29] = (
       +rate_eval.P29__p_Si28
       )

    jac[jsi28, jp31] = (
       +rho*Y[jp]*rate_eval.p_P31__He4_Si28
       )

    jac[jsi28, js31] = (
       +rho*Y[jn]*rate_eval.n_S31__He4_Si28
       )

    jac[jsi28, js32] = (
       +rate_eval.S32__He4_Si28
       )

    jac[jsi28, jcl32] = (
       +rate_eval.Cl32__He4_Si28__weak__wc12
       )

    jac[jsi29, jn] = (
       -rho*Y[jsi29]*rate_eval.n_Si29__Si30
       -rho*Y[jsi29]*rate_eval.n_Si29__He4_Mg26
       +rho*Y[jsi28]*rate_eval.n_Si28__Si29
       +rho*Y[jp29]*rate_eval.n_P29__p_Si29
       +rho*Y[js32]*rate_eval.n_S32__He4_Si29
       )

    jac[jsi29, jp] = (
       -rho*Y[jsi29]*rate_eval.p_Si29__P30
       -rho*Y[jsi29]*rate_eval.p_Si29__n_P29
       -rho*Y[jsi29]*rate_eval.p_Si29__He4_Al26
       +rho*Y[jal28]*rate_eval.p_Al28__Si29
       +rho*Y[jp32]*rate_eval.p_P32__He4_Si29
       )

    jac[jsi29, jhe4] = (
       -rho*Y[jsi29]*rate_eval.He4_Si29__n_S32
       -rho*Y[jsi29]*rate_eval.He4_Si29__p_P32
       +rho*Y[jmg25]*rate_eval.He4_Mg25__Si29
       +rho*Y[jmg26]*rate_eval.He4_Mg26__n_Si29
       +rho*Y[jal26]*rate_eval.He4_Al26__p_Si29
       )

    jac[jsi29, jmg25] = (
       +rho*Y[jhe4]*rate_eval.He4_Mg25__Si29
       )

    jac[jsi29, jmg26] = (
       +rho*Y[jhe4]*rate_eval.He4_Mg26__n_Si29
       )

    jac[jsi29, jal26] = (
       +rho*Y[jhe4]*rate_eval.He4_Al26__p_Si29
       )

    jac[jsi29, jal28] = (
       +rho*Y[jp]*rate_eval.p_Al28__Si29
       )

    jac[jsi29, jsi28] = (
       +rho*Y[jn]*rate_eval.n_Si28__Si29
       )

    jac[jsi29, jsi29] = (
       -rate_eval.Si29__n_Si28
       -rate_eval.Si29__p_Al28
       -rate_eval.Si29__He4_Mg25
       -rho*Y[jn]*rate_eval.n_Si29__Si30
       -rho*Y[jp]*rate_eval.p_Si29__P30
       -rho*Y[jn]*rate_eval.n_Si29__He4_Mg26
       -rho*Y[jp]*rate_eval.p_Si29__n_P29
       -rho*Y[jp]*rate_eval.p_Si29__He4_Al26
       -rho*Y[jhe4]*rate_eval.He4_Si29__n_S32
       -rho*Y[jhe4]*rate_eval.He4_Si29__p_P32
       )

    jac[jsi29, jsi30] = (
       +rate_eval.Si30__n_Si29
       )

    jac[jsi29, jp29] = (
       +rate_eval.P29__Si29__weak__wc12
       +rho*Y[jn]*rate_eval.n_P29__p_Si29
       )

    jac[jsi29, jp30] = (
       +rate_eval.P30__p_Si29
       )

    jac[jsi29, jp32] = (
       +rho*Y[jp]*rate_eval.p_P32__He4_Si29
       )

    jac[jsi29, js32] = (
       +rho*Y[jn]*rate_eval.n_S32__He4_Si29
       )

    jac[jsi30, jn] = (
       +rho*Y[jsi29]*rate_eval.n_Si29__Si30
       +rho*Y[jp30]*rate_eval.n_P30__p_Si30
       )

    jac[jsi30, jp] = (
       -rho*Y[jsi30]*rate_eval.p_Si30__P31
       -rho*Y[jsi30]*rate_eval.p_Si30__n_P30
       -rho*Y[jsi30]*rate_eval.p_Si30__He4_Al27
       )

    jac[jsi30, jhe4] = (
       -rho*Y[jsi30]*rate_eval.He4_Si30__S34
       +rho*Y[jmg26]*rate_eval.He4_Mg26__Si30
       +rho*Y[jal27]*rate_eval.He4_Al27__p_Si30
       )

    jac[jsi30, jmg26] = (
       +rho*Y[jhe4]*rate_eval.He4_Mg26__Si30
       )

    jac[jsi30, jal27] = (
       +rho*Y[jhe4]*rate_eval.He4_Al27__p_Si30
       )

    jac[jsi30, jsi29] = (
       +rho*Y[jn]*rate_eval.n_Si29__Si30
       )

    jac[jsi30, jsi30] = (
       -rate_eval.Si30__n_Si29
       -rate_eval.Si30__He4_Mg26
       -rho*Y[jp]*rate_eval.p_Si30__P31
       -rho*Y[jhe4]*rate_eval.He4_Si30__S34
       -rho*Y[jp]*rate_eval.p_Si30__n_P30
       -rho*Y[jp]*rate_eval.p_Si30__He4_Al27
       )

    jac[jsi30, jp30] = (
       +rate_eval.P30__Si30__weak__wc12
       +rho*Y[jn]*rate_eval.n_P30__p_Si30
       )

    jac[jsi30, jp31] = (
       +rate_eval.P31__p_Si30
       )

    jac[jsi30, js34] = (
       +rate_eval.S34__He4_Si30
       )

    jac[jp29, jn] = (
       -rho*Y[jp29]*rate_eval.n_P29__P30
       -rho*Y[jp29]*rate_eval.n_P29__p_Si29
       -rho*Y[jp29]*rate_eval.n_P29__He4_Al26
       +rho*Y[jcl32]*rate_eval.n_Cl32__He4_P29
       )

    jac[jp29, jp] = (
       -rho*Y[jp29]*rate_eval.p_P29__He4_Si26
       +rho*Y[jsi28]*rate_eval.p_Si28__P29
       +rho*Y[jsi29]*rate_eval.p_Si29__n_P29
       +rho*Y[js32]*rate_eval.p_S32__He4_P29
       )

    jac[jp29, jhe4] = (
       -rho*Y[jp29]*rate_eval.He4_P29__n_Cl32
       -rho*Y[jp29]*rate_eval.He4_P29__p_S32
       +rho*Y[jal25]*rate_eval.He4_Al25__P29
       +rho*Y[jal26]*rate_eval.He4_Al26__n_P29
       +rho*Y[jsi26]*rate_eval.He4_Si26__p_P29
       )

    jac[jp29, jal25] = (
       +rho*Y[jhe4]*rate_eval.He4_Al25__P29
       )

    jac[jp29, jal26] = (
       +rho*Y[jhe4]*rate_eval.He4_Al26__n_P29
       )

    jac[jp29, jsi26] = (
       +rho*Y[jhe4]*rate_eval.He4_Si26__p_P29
       )

    jac[jp29, jsi28] = (
       +rho*Y[jp]*rate_eval.p_Si28__P29
       )

    jac[jp29, jsi29] = (
       +rho*Y[jp]*rate_eval.p_Si29__n_P29
       )

    jac[jp29, jp29] = (
       -rate_eval.P29__Si29__weak__wc12
       -rate_eval.P29__p_Si28
       -rate_eval.P29__He4_Al25
       -rho*Y[jn]*rate_eval.n_P29__P30
       -rho*Y[jn]*rate_eval.n_P29__p_Si29
       -rho*Y[jn]*rate_eval.n_P29__He4_Al26
       -rho*Y[jp]*rate_eval.p_P29__He4_Si26
       -rho*Y[jhe4]*rate_eval.He4_P29__n_Cl32
       -rho*Y[jhe4]*rate_eval.He4_P29__p_S32
       )

    jac[jp29, jp30] = (
       +rate_eval.P30__n_P29
       )

    jac[jp29, js32] = (
       +rho*Y[jp]*rate_eval.p_S32__He4_P29
       )

    jac[jp29, jcl32] = (
       +rho*Y[jn]*rate_eval.n_Cl32__He4_P29
       )

    jac[jp30, jn] = (
       -rho*Y[jp30]*rate_eval.n_P30__P31
       -rho*Y[jp30]*rate_eval.n_P30__p_Si30
       -rho*Y[jp30]*rate_eval.n_P30__He4_Al27
       +rho*Y[jp29]*rate_eval.n_P29__P30
       )

    jac[jp30, jp] = (
       -rho*Y[jp30]*rate_eval.p_P30__S31
       +rho*Y[jsi29]*rate_eval.p_Si29__P30
       +rho*Y[jsi30]*rate_eval.p_Si30__n_P30
       )

    jac[jp30, jhe4] = (
       -rho*Y[jp30]*rate_eval.He4_P30__Cl34
       +rho*Y[jal26]*rate_eval.He4_Al26__P30
       +rho*Y[jal27]*rate_eval.He4_Al27__n_P30
       )

    jac[jp30, jal26] = (
       +rho*Y[jhe4]*rate_eval.He4_Al26__P30
       )

    jac[jp30, jal27] = (
       +rho*Y[jhe4]*rate_eval.He4_Al27__n_P30
       )

    jac[jp30, jsi29] = (
       +rho*Y[jp]*rate_eval.p_Si29__P30
       )

    jac[jp30, jsi30] = (
       +rho*Y[jp]*rate_eval.p_Si30__n_P30
       )

    jac[jp30, jp29] = (
       +rho*Y[jn]*rate_eval.n_P29__P30
       )

    jac[jp30, jp30] = (
       -rate_eval.P30__Si30__weak__wc12
       -rate_eval.P30__n_P29
       -rate_eval.P30__p_Si29
       -rate_eval.P30__He4_Al26
       -rho*Y[jn]*rate_eval.n_P30__P31
       -rho*Y[jp]*rate_eval.p_P30__S31
       -rho*Y[jhe4]*rate_eval.He4_P30__Cl34
       -rho*Y[jn]*rate_eval.n_P30__p_Si30
       -rho*Y[jn]*rate_eval.n_P30__He4_Al27
       )

    jac[jp30, jp31] = (
       +rate_eval.P31__n_P30
       )

    jac[jp30, js31] = (
       +rate_eval.S31__p_P30
       )

    jac[jp30, jcl34] = (
       +rate_eval.Cl34__He4_P30
       )

    jac[jp31, jn] = (
       -rho*Y[jp31]*rate_eval.n_P31__P32
       -rho*Y[jp31]*rate_eval.n_P31__He4_Al28
       +rho*Y[jp30]*rate_eval.n_P30__P31
       +rho*Y[js31]*rate_eval.n_S31__p_P31
       +rho*Y[jcl34]*rate_eval.n_Cl34__He4_P31
       )

    jac[jp31, jp] = (
       -rho*Y[jp31]*rate_eval.p_P31__S32
       -rho*Y[jp31]*rate_eval.p_P31__n_S31
       -rho*Y[jp31]*rate_eval.p_P31__He4_Si28
       -rho*Y[jp31]*rate_eval.p_P31__C12_Ne20
       -rho*Y[jp31]*rate_eval.p_P31__O16_O16
       +rho*Y[jsi30]*rate_eval.p_Si30__P31
       +rho*Y[js34]*rate_eval.p_S34__He4_P31
       )

    jac[jp31, jhe4] = (
       -rho*Y[jp31]*rate_eval.He4_P31__Cl35
       -rho*Y[jp31]*rate_eval.He4_P31__n_Cl34
       -rho*Y[jp31]*rate_eval.He4_P31__p_S34
       +rho*Y[jal27]*rate_eval.He4_Al27__P31
       +rho*Y[jal28]*rate_eval.He4_Al28__n_P31
       +rho*Y[jsi28]*rate_eval.He4_Si28__p_P31
       )

    jac[jp31, jc12] = (
       +rho*Y[jne20]*rate_eval.C12_Ne20__p_P31
       )

    jac[jp31, jo16] = (
       +5.00000000000000e-01*rho*2*Y[jo16]*rate_eval.O16_O16__p_P31
       )

    jac[jp31, jne20] = (
       +rho*Y[jc12]*rate_eval.C12_Ne20__p_P31
       )

    jac[jp31, jal27] = (
       +rho*Y[jhe4]*rate_eval.He4_Al27__P31
       )

    jac[jp31, jal28] = (
       +rho*Y[jhe4]*rate_eval.He4_Al28__n_P31
       )

    jac[jp31, jsi28] = (
       +rho*Y[jhe4]*rate_eval.He4_Si28__p_P31
       )

    jac[jp31, jsi30] = (
       +rho*Y[jp]*rate_eval.p_Si30__P31
       )

    jac[jp31, jp30] = (
       +rho*Y[jn]*rate_eval.n_P30__P31
       )

    jac[jp31, jp31] = (
       -rate_eval.P31__n_P30
       -rate_eval.P31__p_Si30
       -rate_eval.P31__He4_Al27
       -rho*Y[jn]*rate_eval.n_P31__P32
       -rho*Y[jp]*rate_eval.p_P31__S32
       -rho*Y[jhe4]*rate_eval.He4_P31__Cl35
       -rho*Y[jn]*rate_eval.n_P31__He4_Al28
       -rho*Y[jp]*rate_eval.p_P31__n_S31
       -rho*Y[jp]*rate_eval.p_P31__He4_Si28
       -rho*Y[jp]*rate_eval.p_P31__C12_Ne20
       -rho*Y[jp]*rate_eval.p_P31__O16_O16
       -rho*Y[jhe4]*rate_eval.He4_P31__n_Cl34
       -rho*Y[jhe4]*rate_eval.He4_P31__p_S34
       )

    jac[jp31, jp32] = (
       +rate_eval.P32__n_P31
       )

    jac[jp31, js31] = (
       +rate_eval.S31__P31__weak__wc12
       +rho*Y[jn]*rate_eval.n_S31__p_P31
       )

    jac[jp31, js32] = (
       +rate_eval.S32__p_P31
       )

    jac[jp31, js34] = (
       +rho*Y[jp]*rate_eval.p_S34__He4_P31
       )

    jac[jp31, jcl32] = (
       +rate_eval.Cl32__p_P31__weak__wc12
       )

    jac[jp31, jcl34] = (
       +rho*Y[jn]*rate_eval.n_Cl34__He4_P31
       )

    jac[jp31, jcl35] = (
       +rate_eval.Cl35__He4_P31
       )

    jac[jp32, jn] = (
       +rho*Y[jp31]*rate_eval.n_P31__P32
       +rho*Y[js32]*rate_eval.n_S32__p_P32
       +rho*Y[jcl35]*rate_eval.n_Cl35__He4_P32
       )

    jac[jp32, jp] = (
       -rho*Y[jp32]*rate_eval.p_P32__n_S32
       -rho*Y[jp32]*rate_eval.p_P32__He4_Si29
       )

    jac[jp32, jhe4] = (
       -rho*Y[jp32]*rate_eval.He4_P32__n_Cl35
       +rho*Y[jal28]*rate_eval.He4_Al28__P32
       +rho*Y[jsi29]*rate_eval.He4_Si29__p_P32
       )

    jac[jp32, jal28] = (
       +rho*Y[jhe4]*rate_eval.He4_Al28__P32
       )

    jac[jp32, jsi29] = (
       +rho*Y[jhe4]*rate_eval.He4_Si29__p_P32
       )

    jac[jp32, jp31] = (
       +rho*Y[jn]*rate_eval.n_P31__P32
       )

    jac[jp32, jp32] = (
       -rate_eval.P32__S32__weak__wc12
       -rate_eval.P32__n_P31
       -rate_eval.P32__He4_Al28
       -rho*Y[jp]*rate_eval.p_P32__n_S32
       -rho*Y[jp]*rate_eval.p_P32__He4_Si29
       -rho*Y[jhe4]*rate_eval.He4_P32__n_Cl35
       )

    jac[jp32, js32] = (
       +rho*Y[jn]*rate_eval.n_S32__p_P32
       )

    jac[jp32, jcl35] = (
       +rho*Y[jn]*rate_eval.n_Cl35__He4_P32
       )

    jac[js31, jn] = (
       -rho*Y[js31]*rate_eval.n_S31__S32
       -rho*Y[js31]*rate_eval.n_S31__p_P31
       -rho*Y[js31]*rate_eval.n_S31__He4_Si28
       -rho*Y[js31]*rate_eval.n_S31__C12_Ne20
       -rho*Y[js31]*rate_eval.n_S31__O16_O16
       )

    jac[js31, jp] = (
       -rho*Y[js31]*rate_eval.p_S31__Cl32
       +rho*Y[jp30]*rate_eval.p_P30__S31
       +rho*Y[jp31]*rate_eval.p_P31__n_S31
       +rho*Y[jcl34]*rate_eval.p_Cl34__He4_S31
       )

    jac[js31, jhe4] = (
       -rho*Y[js31]*rate_eval.He4_S31__Ar35
       -rho*Y[js31]*rate_eval.He4_S31__p_Cl34
       +rho*Y[jsi28]*rate_eval.He4_Si28__n_S31
       )

    jac[js31, jc12] = (
       +rho*Y[jne20]*rate_eval.C12_Ne20__n_S31
       )

    jac[js31, jo16] = (
       +5.00000000000000e-01*rho*2*Y[jo16]*rate_eval.O16_O16__n_S31
       )

    jac[js31, jne20] = (
       +rho*Y[jc12]*rate_eval.C12_Ne20__n_S31
       )

    jac[js31, jsi28] = (
       +rho*Y[jhe4]*rate_eval.He4_Si28__n_S31
       )

    jac[js31, jp30] = (
       +rho*Y[jp]*rate_eval.p_P30__S31
       )

    jac[js31, jp31] = (
       +rho*Y[jp]*rate_eval.p_P31__n_S31
       )

    jac[js31, js31] = (
       -rate_eval.S31__P31__weak__wc12
       -rate_eval.S31__p_P30
       -rho*Y[jn]*rate_eval.n_S31__S32
       -rho*Y[jp]*rate_eval.p_S31__Cl32
       -rho*Y[jhe4]*rate_eval.He4_S31__Ar35
       -rho*Y[jn]*rate_eval.n_S31__p_P31
       -rho*Y[jn]*rate_eval.n_S31__He4_Si28
       -rho*Y[jn]*rate_eval.n_S31__C12_Ne20
       -rho*Y[jn]*rate_eval.n_S31__O16_O16
       -rho*Y[jhe4]*rate_eval.He4_S31__p_Cl34
       )

    jac[js31, js32] = (
       +rate_eval.S32__n_S31
       )

    jac[js31, jcl32] = (
       +rate_eval.Cl32__p_S31
       )

    jac[js31, jcl34] = (
       +rho*Y[jp]*rate_eval.p_Cl34__He4_S31
       )

    jac[js31, jar35] = (
       +rate_eval.Ar35__He4_S31
       )

    jac[js32, jn] = (
       -rho*Y[js32]*rate_eval.n_S32__p_P32
       -rho*Y[js32]*rate_eval.n_S32__He4_Si29
       +rho*Y[js31]*rate_eval.n_S31__S32
       +rho*Y[jcl32]*rate_eval.n_Cl32__p_S32
       +rho*Y[jar35]*rate_eval.n_Ar35__He4_S32
       )

    jac[js32, jp] = (
       -rho*Y[js32]*rate_eval.p_S32__n_Cl32
       -rho*Y[js32]*rate_eval.p_S32__He4_P29
       +rho*Y[jp31]*rate_eval.p_P31__S32
       +rho*Y[jp32]*rate_eval.p_P32__n_S32
       +rho*Y[jcl35]*rate_eval.p_Cl35__He4_S32
       )

    jac[js32, jhe4] = (
       -rho*Y[js32]*rate_eval.He4_S32__n_Ar35
       -rho*Y[js32]*rate_eval.He4_S32__p_Cl35
       +rho*Y[jsi28]*rate_eval.He4_Si28__S32
       +rho*Y[jsi29]*rate_eval.He4_Si29__n_S32
       +rho*Y[jp29]*rate_eval.He4_P29__p_S32
       )

    jac[js32, jsi28] = (
       +rho*Y[jhe4]*rate_eval.He4_Si28__S32
       )

    jac[js32, jsi29] = (
       +rho*Y[jhe4]*rate_eval.He4_Si29__n_S32
       )

    jac[js32, jp29] = (
       +rho*Y[jhe4]*rate_eval.He4_P29__p_S32
       )

    jac[js32, jp31] = (
       +rho*Y[jp]*rate_eval.p_P31__S32
       )

    jac[js32, jp32] = (
       +rate_eval.P32__S32__weak__wc12
       +rho*Y[jp]*rate_eval.p_P32__n_S32
       )

    jac[js32, js31] = (
       +rho*Y[jn]*rate_eval.n_S31__S32
       )

    jac[js32, js32] = (
       -rate_eval.S32__n_S31
       -rate_eval.S32__p_P31
       -rate_eval.S32__He4_Si28
       -rho*Y[jn]*rate_eval.n_S32__p_P32
       -rho*Y[jn]*rate_eval.n_S32__He4_Si29
       -rho*Y[jp]*rate_eval.p_S32__n_Cl32
       -rho*Y[jp]*rate_eval.p_S32__He4_P29
       -rho*Y[jhe4]*rate_eval.He4_S32__n_Ar35
       -rho*Y[jhe4]*rate_eval.He4_S32__p_Cl35
       )

    jac[js32, jcl32] = (
       +rate_eval.Cl32__S32__weak__wc12
       +rho*Y[jn]*rate_eval.n_Cl32__p_S32
       )

    jac[js32, jcl35] = (
       +rho*Y[jp]*rate_eval.p_Cl35__He4_S32
       )

    jac[js32, jar35] = (
       +rho*Y[jn]*rate_eval.n_Ar35__He4_S32
       )

    jac[js34, jn] = (
       +rho*Y[jcl34]*rate_eval.n_Cl34__p_S34
       )

    jac[js34, jp] = (
       -rho*Y[js34]*rate_eval.p_S34__Cl35
       -rho*Y[js34]*rate_eval.p_S34__n_Cl34
       -rho*Y[js34]*rate_eval.p_S34__He4_P31
       )

    jac[js34, jhe4] = (
       +rho*Y[jsi30]*rate_eval.He4_Si30__S34
       +rho*Y[jp31]*rate_eval.He4_P31__p_S34
       )

    jac[js34, jsi30] = (
       +rho*Y[jhe4]*rate_eval.He4_Si30__S34
       )

    jac[js34, jp31] = (
       +rho*Y[jhe4]*rate_eval.He4_P31__p_S34
       )

    jac[js34, js34] = (
       -rate_eval.S34__He4_Si30
       -rho*Y[jp]*rate_eval.p_S34__Cl35
       -rho*Y[jp]*rate_eval.p_S34__n_Cl34
       -rho*Y[jp]*rate_eval.p_S34__He4_P31
       )

    jac[js34, jcl34] = (
       +rate_eval.Cl34__S34__weak__wc12
       +rho*Y[jn]*rate_eval.n_Cl34__p_S34
       )

    jac[js34, jcl35] = (
       +rate_eval.Cl35__p_S34
       )

    jac[jcl32, jn] = (
       -rho*Y[jcl32]*rate_eval.n_Cl32__p_S32
       -rho*Y[jcl32]*rate_eval.n_Cl32__He4_P29
       )

    jac[jcl32, jp] = (
       +rho*Y[js31]*rate_eval.p_S31__Cl32
       +rho*Y[js32]*rate_eval.p_S32__n_Cl32
       +rho*Y[jar35]*rate_eval.p_Ar35__He4_Cl32
       )

    jac[jcl32, jhe4] = (
       -rho*Y[jcl32]*rate_eval.He4_Cl32__p_Ar35
       +rho*Y[jp29]*rate_eval.He4_P29__n_Cl32
       )

    jac[jcl32, jp29] = (
       +rho*Y[jhe4]*rate_eval.He4_P29__n_Cl32
       )

    jac[jcl32, js31] = (
       +rho*Y[jp]*rate_eval.p_S31__Cl32
       )

    jac[jcl32, js32] = (
       +rho*Y[jp]*rate_eval.p_S32__n_Cl32
       )

    jac[jcl32, jcl32] = (
       -rate_eval.Cl32__S32__weak__wc12
       -rate_eval.Cl32__p_S31
       -rate_eval.Cl32__p_P31__weak__wc12
       -rate_eval.Cl32__He4_Si28__weak__wc12
       -rho*Y[jn]*rate_eval.n_Cl32__p_S32
       -rho*Y[jn]*rate_eval.n_Cl32__He4_P29
       -rho*Y[jhe4]*rate_eval.He4_Cl32__p_Ar35
       )

    jac[jcl32, jar35] = (
       +rho*Y[jp]*rate_eval.p_Ar35__He4_Cl32
       )

    jac[jcl34, jn] = (
       -rho*Y[jcl34]*rate_eval.n_Cl34__Cl35
       -rho*Y[jcl34]*rate_eval.n_Cl34__p_S34
       -rho*Y[jcl34]*rate_eval.n_Cl34__He4_P31
       )

    jac[jcl34, jp] = (
       -rho*Y[jcl34]*rate_eval.p_Cl34__Ar35
       -rho*Y[jcl34]*rate_eval.p_Cl34__He4_S31
       +rho*Y[js34]*rate_eval.p_S34__n_Cl34
       )

    jac[jcl34, jhe4] = (
       +rho*Y[jp30]*rate_eval.He4_P30__Cl34
       +rho*Y[jp31]*rate_eval.He4_P31__n_Cl34
       +rho*Y[js31]*rate_eval.He4_S31__p_Cl34
       )

    jac[jcl34, jp30] = (
       +rho*Y[jhe4]*rate_eval.He4_P30__Cl34
       )

    jac[jcl34, jp31] = (
       +rho*Y[jhe4]*rate_eval.He4_P31__n_Cl34
       )

    jac[jcl34, js31] = (
       +rho*Y[jhe4]*rate_eval.He4_S31__p_Cl34
       )

    jac[jcl34, js34] = (
       +rho*Y[jp]*rate_eval.p_S34__n_Cl34
       )

    jac[jcl34, jcl34] = (
       -rate_eval.Cl34__S34__weak__wc12
       -rate_eval.Cl34__He4_P30
       -rho*Y[jn]*rate_eval.n_Cl34__Cl35
       -rho*Y[jp]*rate_eval.p_Cl34__Ar35
       -rho*Y[jn]*rate_eval.n_Cl34__p_S34
       -rho*Y[jn]*rate_eval.n_Cl34__He4_P31
       -rho*Y[jp]*rate_eval.p_Cl34__He4_S31
       )

    jac[jcl34, jcl35] = (
       +rate_eval.Cl35__n_Cl34
       )

    jac[jcl34, jar35] = (
       +rate_eval.Ar35__p_Cl34
       )

    jac[jcl35, jn] = (
       -rho*Y[jcl35]*rate_eval.n_Cl35__He4_P32
       +rho*Y[jcl34]*rate_eval.n_Cl34__Cl35
       +rho*Y[jar35]*rate_eval.n_Ar35__p_Cl35
       )

    jac[jcl35, jp] = (
       -rho*Y[jcl35]*rate_eval.p_Cl35__n_Ar35
       -rho*Y[jcl35]*rate_eval.p_Cl35__He4_S32
       +rho*Y[js34]*rate_eval.p_S34__Cl35
       )

    jac[jcl35, jhe4] = (
       +rho*Y[jp31]*rate_eval.He4_P31__Cl35
       +rho*Y[jp32]*rate_eval.He4_P32__n_Cl35
       +rho*Y[js32]*rate_eval.He4_S32__p_Cl35
       )

    jac[jcl35, jp31] = (
       +rho*Y[jhe4]*rate_eval.He4_P31__Cl35
       )

    jac[jcl35, jp32] = (
       +rho*Y[jhe4]*rate_eval.He4_P32__n_Cl35
       )

    jac[jcl35, js32] = (
       +rho*Y[jhe4]*rate_eval.He4_S32__p_Cl35
       )

    jac[jcl35, js34] = (
       +rho*Y[jp]*rate_eval.p_S34__Cl35
       )

    jac[jcl35, jcl34] = (
       +rho*Y[jn]*rate_eval.n_Cl34__Cl35
       )

    jac[jcl35, jcl35] = (
       -rate_eval.Cl35__n_Cl34
       -rate_eval.Cl35__p_S34
       -rate_eval.Cl35__He4_P31
       -rho*Y[jn]*rate_eval.n_Cl35__He4_P32
       -rho*Y[jp]*rate_eval.p_Cl35__n_Ar35
       -rho*Y[jp]*rate_eval.p_Cl35__He4_S32
       )

    jac[jcl35, jar35] = (
       +rate_eval.Ar35__Cl35__weak__wc12
       +rho*Y[jn]*rate_eval.n_Ar35__p_Cl35
       )

    jac[jar35, jn] = (
       -rho*Y[jar35]*rate_eval.n_Ar35__p_Cl35
       -rho*Y[jar35]*rate_eval.n_Ar35__He4_S32
       )

    jac[jar35, jp] = (
       -rho*Y[jar35]*rate_eval.p_Ar35__He4_Cl32
       +rho*Y[jcl34]*rate_eval.p_Cl34__Ar35
       +rho*Y[jcl35]*rate_eval.p_Cl35__n_Ar35
       )

    jac[jar35, jhe4] = (
       +rho*Y[js31]*rate_eval.He4_S31__Ar35
       +rho*Y[js32]*rate_eval.He4_S32__n_Ar35
       +rho*Y[jcl32]*rate_eval.He4_Cl32__p_Ar35
       )

    jac[jar35, js31] = (
       +rho*Y[jhe4]*rate_eval.He4_S31__Ar35
       )

    jac[jar35, js32] = (
       +rho*Y[jhe4]*rate_eval.He4_S32__n_Ar35
       )

    jac[jar35, jcl32] = (
       +rho*Y[jhe4]*rate_eval.He4_Cl32__p_Ar35
       )

    jac[jar35, jcl34] = (
       +rho*Y[jp]*rate_eval.p_Cl34__Ar35
       )

    jac[jar35, jcl35] = (
       +rho*Y[jp]*rate_eval.p_Cl35__n_Ar35
       )

    jac[jar35, jar35] = (
       -rate_eval.Ar35__Cl35__weak__wc12
       -rate_eval.Ar35__p_Cl34
       -rate_eval.Ar35__He4_S31
       -rho*Y[jn]*rate_eval.n_Ar35__p_Cl35
       -rho*Y[jn]*rate_eval.n_Ar35__He4_S32
       -rho*Y[jp]*rate_eval.p_Ar35__He4_Cl32
       )

    return jac
