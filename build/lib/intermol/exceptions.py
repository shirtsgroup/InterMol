class InterMolError(Exception):
    """"""


class MultipleValidationErrors(InterMolError):
    """"""
    def __str__(self):
        return '\n\n{0}\n\n'.format('\n'.join(self.args))


class ConversionError(InterMolError):
    """"""
    def __init__(self, could_not_convert, engine):
        Exception.__init__(self)
        self.could_not_convert = could_not_convert
        self.engine = engine


class UnsupportedFunctional(ConversionError):
    """Force functional that is not supported in a specific engine."""
    def __str__(self):
        return "The {} interaction type is not supported in {}.".format(
            self.could_not_convert.__class__.__name__, self.engine.upper())


class UnimplementedFunctional(ConversionError):
    """Force functional that is not yet implemented in intermol for a specific engine."""
    def __str__(self):
        return ('{} conversion has not yet been implemented in InterMol '
                'for {}.'.format(self.could_not_convert.__class__.__name__, self.engine.upper()))


class UnsupportedSetting(ConversionError):
    """Any setting that is not supported in a specific engine."""
    def __str__(self):
        return "{} is not supported in {}.".format(
            self.could_not_convert, self.engine.upper())


class UnimplementedSetting(ConversionError):
    """Any setting that is not yet implemented in intermol for a specific engine."""
    def __str__(self):
        return ('{} has not yet been implemented in InterMol for {}.'.format(
            self.could_not_convert, self.engine.upper()))


class ParsingError(InterMolError):
    """"""


class GromacsError(ParsingError):
    """"""


class AmberError(ParsingError):
    """"""


class DesmondError(ParsingError):
    """"""


class LammpsError(ParsingError):
    """"""
