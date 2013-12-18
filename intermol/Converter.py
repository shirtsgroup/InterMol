import intermol.unit as units

def convert_units(arg, unit):
    """Checks compatibility and converts units using simtk.units package

    Args:
        arg (Quantity): quantity to be converted
        unit (Unit): Unit to be converted to

    Returns:
        arg (Quantity): Quantity scaled to the new unit
    """
    conversionFactor = (arg.unit).conversion_factor_to(unit)
    arg = arg * conversionFactor
    return arg._value * unit
