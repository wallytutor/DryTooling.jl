# -*- coding: utf-8 -*-
module Elements

""" Represents a chemical element. """
struct ElementData
    symbol::String
    name::String
    atomicmass::Float64
    standardenthalpy::Float64
    standardentropy::Float64

    function ElementData(
            symbol::String,
            name::String,
            atomicmass::Float64;
            h0::Float64 = 0.0,
            s0::Float64 = 0.0
        )
        return new(symbol, name, atomicmass, h0, s0)        
    end
end

""" Retrieve an element by name. """
element(s::String) = getfield(Elements, Symbol(s))

""" Retrieve atomic mass of element. """
mass(e::ElementData) = e.atomicmass
mass(s::String) = mass(element(s))

const H  = ElementData("H",  "hydrogen",        1.008)
const He = ElementData("He", "helium",          4.002602)
const Li = ElementData("Li", "lithium",         6.94)
const Be = ElementData("Be", "beryllium",       9.0121831)
const B  = ElementData("B",  "boron",          10.81)
const C  = ElementData("C",  "carbon",         12.011)
const N  = ElementData("N",  "nitrogen",       14.007)
const O  = ElementData("O",  "oxygen",         15.999)
const F  = ElementData("F",  "fluorine",       18.998403163)
const Ne = ElementData("Ne", "neon",           20.1797)
const Na = ElementData("Na", "sodium",         22.98976928)
const Mg = ElementData("Mg", "magnesium",      24.305)
const Al = ElementData("Al", "aluminum",       26.9815384)
const Si = ElementData("Si", "silicon",        28.085)
const P  = ElementData("P",  "phosphorus",     30.973761998)
const S  = ElementData("S",  "sulfur",         32.06)
const Cl = ElementData("Cl", "chlorine",       35.45)
const Ar = ElementData("Ar", "argon",          39.95)
const K  = ElementData("K",  "potassium",      39.0983)
const Ca = ElementData("Ca", "calcium",        40.078)
const Sc = ElementData("Sc", "scandium",       44.955908)
const Ti = ElementData("Ti", "titanium",       47.867)
const V  = ElementData("V",  "vanadium",       50.9415)
const Cr = ElementData("Cr", "chromium",       51.9961)
const Mn = ElementData("Mn", "manganese",      54.938043)
const Fe = ElementData("Fe", "iron",           55.845)
const Co = ElementData("Co", "cobalt",         58.933194)
const Ni = ElementData("Ni", "nickel",         58.6934)
const Cu = ElementData("Cu", "copper",         63.546)
const Zn = ElementData("Zn", "zinc",           65.38)
const Ga = ElementData("Ga", "gallium",        69.723)
const Ge = ElementData("Ge", "germanium",      72.630)
const As = ElementData("As", "arsenic",        74.921595)
const Se = ElementData("Se", "selenium",       78.971)
const Br = ElementData("Br", "bromine",        79.904)
const Kr = ElementData("Kr", "krypton",        83.798)
const Rb = ElementData("Rb", "rubidium",       85.4678)
const Sr = ElementData("Sr", "strontium",      87.62)
const Y  = ElementData("Y",  "yttrium",        88.90584)
const Zr = ElementData("Zr", "zirconium",      91.224)
const Nb = ElementData("Nb", "nobelium",       92.90637)
const Mo = ElementData("Mo", "molybdenum",     95.95)
const Ru = ElementData("Ru", "ruthenium",     101.07)
const Rh = ElementData("Rh", "rhodium",       102.90549)
const Pd = ElementData("Pd", "palladium",     106.42)
const Ag = ElementData("Ag", "silver",        107.8682)
const Cd = ElementData("Cd", "cadmium",       112.414)
const In = ElementData("In", "indium",        114.818)
const Sn = ElementData("Sn", "tin",           118.710)
const Sb = ElementData("Sb", "antimony",      121.760)
const Te = ElementData("Te", "tellurium",     127.60 )
const I  = ElementData("I",  "iodine",        126.90447)
const Xe = ElementData("Xe", "xenon",         131.293)
const Cs = ElementData("Cs", "cesium",        132.90545196)
const Ba = ElementData("Ba", "barium",        137.327)
const La = ElementData("La", "lanthanum",     138.90547)
const Ce = ElementData("Ce", "cerium",        140.116)
const Pr = ElementData("Pr", "praseodymium",  140.90766)
const Nd = ElementData("Nd", "neodymium",     144.242)
const Sm = ElementData("Sm", "samarium",      150.36)
const Eu = ElementData("Eu", "europium",      151.964)
const Gd = ElementData("Gd", "gadolinium",    157.25)
const Tb = ElementData("Tb", "terbium",       158.925354)
const Dy = ElementData("Dy", "dysprosium",    162.500)
const Ho = ElementData("Ho", "holmium",       164.930328)
const Er = ElementData("Er", "erbium",        167.259)
const Tm = ElementData("Tm", "thulium",       168.934218)
const Yb = ElementData("Yb", "ytterbium",     173.045)
const Lu = ElementData("Lu", "lutetium",      174.9668)
const Hf = ElementData("Hf", "hafnium",       178.49)
const Ta = ElementData("Ta", "tantalum",      180.94788)
const W  = ElementData("W",  "tungsten",      183.84)
const Re = ElementData("Re", "rhenium",       186.207)
const Os = ElementData("Os", "osmium",        190.23 )
const Ir = ElementData("Ir", "iridium",       192.217)
const Pt = ElementData("Pt", "platinum",      195.084)
const Au = ElementData("Au", "gold",          196.966570)
const Hg = ElementData("Hg", "mercury",       200.592)
const Tl = ElementData("Tl", "thallium",      204.38)
const Pb = ElementData("Pb", "lead",          207.2 )
const Bi = ElementData("Bi", "bismuth",       208.98040)
const Th = ElementData("Th", "thorium",       232.0377)
const Pa = ElementData("Pa", "protactinium",  231.03588)
const U  = ElementData("U",  "uranium",       238.02891)

end # (module Elements)