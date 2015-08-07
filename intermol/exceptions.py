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
    """"""
    def __str__(self):
        return "{}'s are not supported in {}.".format(
            self.could_not_convert.__class__.__name__, self.engine)


class UnimplementedFunctional(ConversionError):
    """"""
    def __str__(self):
        return ('{} conversion has not yet been implemented in InterMol '
                'for {}.'.format(self.could_not_convert.__class__.__name__, self.engine))


class UnsupportedSetting(ConversionError):
    """"""
    def __str__(self):
        return "{} is not supported in {}.".format(
            self.could_not_convert, self.engine)


class UnimplementedSetting(ConversionError):
    """"""
    def __str__(self):
        return ('{} has not yet been implemented in InterMol for {}.'.format(
            self.could_not_convert, self.engine))


class ParsingError(InterMolError):
    """"""


class GromacsError(ParsingError):
    """"""


class DesmondError(ParsingError):
    """"""


class LammpsError(ParsingError):
    """"""
