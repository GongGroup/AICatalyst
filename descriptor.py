import re

PeriodicTable = [
    'H', 'He',
    'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
    'Na', 'Mg', 'Al', 'P', 'S', 'Cl', 'Ar', 'Si',
    'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
    'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
    'Cs', 'Ba', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
    'Fr', 'Ra', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og',
    'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',
    'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr',
]


class Descriptor(object):
    """
    Base Descriptor, <get, set, del> method of each param is unlimited
    """

    def __new__(cls, *args, **kwargs):
        if cls is Descriptor:
            raise TypeError(f"<{cls.__name__} class> may not be instantiated")
        return super().__new__(cls)

    def __init__(self, name):
        self.name = name

    def __get__(self, instance, cls):
        if instance is None:
            return self
        elif self.name in instance.__dict__.keys():
            return instance.__dict__[self.name]
        else:
            return cls.__getattr__(instance, self.name)

    def __set__(self, instance, value):
        instance.__dict__[self.name] = value

    def __delete__(self, instance):
        del instance.__dict__[self.name]


class ValueDescriptor(Descriptor):
    """
    Value Descriptor, limit param's value
    """

    def __init__(self, name, value):
        super(ValueDescriptor, self).__init__(name=name)
        self.value = value


class FormulaDescriptor(Descriptor):
    """
    FormulaDescriptor, check name is valid as chemical formula
    """

    def __set__(self, instance, value: str):

        # check digital
        if value.isdigit():
            raise ValueError(f"`{value}` is invalid name")

        # check special symbols
        symbols = [",", " ", "[", "]", "WSC"]
        for symbol in symbols:
            if symbol in value:
                raise ValueError(f"`{value}` is invalid name")

        # check lower
        if value.islower():
            raise ValueError(f"`{value}` is invalid name")

        # check element
        check_value = value
        check_value = re.sub(r'[0-9()*η+-]+', '', check_value)
        for element in PeriodicTable[::-1]:
            check_value = re.sub(element, '', check_value)
        if len(check_value):
            raise ValueError(f"`{value}` is invalid name")

        super().__set__(instance, value)


class CatalystDescriptor(ValueDescriptor):
    """
    CatalystDescriptor, check name is valid as catalyst
    """

    # TODO: need add other conditions, e.g. Amino not but True, N,N-diethylisonicotinamide not but true
    def __set__(self, instance, value):
        for item in self.value:
            if item in value:
                break
        else:
            raise ValueError(f"`{value}` is invalid name")
        super().__set__(instance, value)
