import re

from AICatalysis.common.constant import PeriodicTable, Abbr


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
        symbols = [",", "WSC"]
        for symbol in symbols:
            if symbol in value:
                raise ValueError(f"`{value}` is invalid name")

        # check lower
        if value.islower():
            raise ValueError(f"`{value}` is invalid name")

        check_value = value

        # check abbr
        for element in Abbr[::-1]:
            check_value = re.sub(element, '', check_value)

        # check element
        for element in PeriodicTable[::-1]:
            check_value = re.sub(element, '', check_value)

        check_value = re.sub(r'[\d()\[\]*.·=/Xxnη+ -{}]+', '', check_value)
        if len(check_value):
            raise ValueError(f"`{value}` is invalid name")

        super().__set__(instance, value)


def not_have_exclude(value, excludes):
    for item in excludes:
        if item in value:
            return False
    else:
        return True


class ReactantDescriptor(Descriptor):
    """
    ReactantDescriptor, check name is valid as time
    """

    def __set__(self, instance, value):
        if re.search('reactant', value) is None:  # match with no-groups
            raise ValueError(f"`{value}` is invalid name")
        super().__set__(instance, value)


class ProductDescriptor(Descriptor):
    """
    ProductDescriptor, check name is valid as time
    """

    def __set__(self, instance, value):
        if re.search('^ product', value) is None:  # match with no-groups
            raise ValueError(f"`{value}` is invalid name")
        super().__set__(instance, value)


global_exclude = ["Reaction", "Conditions", "General", "1a"]


class ReagentDescriptor(ValueDescriptor):
    """
    MetalDescriptor, check name has metal
    """
    excludes = ["Ligand", "Yield", "Table", "Run", "Ref"] + global_exclude

    def __set__(self, instance, value):
        for item in self.value:
            if item in value and not_have_exclude(value, self.excludes):
                break
        else:
            raise ValueError(f"`{value}` is invalid `{instance.__class__.__name__}` name")
        super().__set__(instance, value)


class SolDescriptor(ReagentDescriptor):
    """
    SolDescriptor, check name is valid as solvent
    """
    excludes = ["1a", "2a", "TBD", "3a", "iodobenzene", "CO (1", "bromobenzene", 'base', "1 bar", "PhI(OAc)2",
                'MePh2SiCO2H', 'NaBPh4']

    def __set__(self, instance, value):
        super(SolDescriptor, self).__set__(instance, value)


class GasDescriptor(ReagentDescriptor):
    """
    SolDescriptor, check name is valid as solvent
    """
    excludes = ["HCOOH", 'carbonyl', 'MePh2SiCO2H', 'Rh4(CO)12', 'HCO2H', 'K2CO3']

    def __set__(self, instance, value):
        super(GasDescriptor, self).__set__(instance, value)


class AdditiveDescriptor(ReagentDescriptor):
    """
    SolDescriptor, check name is valid as solvent
    """
    excludes = ["1a", "2a", "3a", 'carbonyl', "Yield"] + global_exclude

    def __set__(self, instance, value):
        super(AdditiveDescriptor, self).__set__(instance, value)


class OxidantDescriptor(ReagentDescriptor):
    """
    SolDescriptor, check name is valid as solvent
    """
    excludes = ["1a", "Reaction"]

    def __set__(self, instance, value):
        super(OxidantDescriptor, self).__set__(instance, value)


class LigandDescriptor(ReagentDescriptor):
    """
    LigandDescriptor, check name is valid as solvent
    """
    excludes = ["bromobenzene", "1a", 'solvent', "CO (1 atm)", "base", "Et3SiH", "1 bar", "PhI(OAc)2"]

    def __set__(self, instance, value):
        super(LigandDescriptor, self).__set__(instance, value)


class TimeDescriptor(Descriptor):
    """
    TimeDescriptor, check name is valid as time
    """

    def __set__(self, instance, value):
        if re.search('[0-9]*–?[0-9]+\s*(?:h|min)', value) is None:  # match with no-groups
            raise ValueError(f"`{value}` is invalid name")
        super().__set__(instance, value)
