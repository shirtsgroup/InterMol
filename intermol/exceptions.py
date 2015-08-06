class InterMolError(Exception):
    """"""


class MultipleValidationErrors(InterMolError):
    """"""
    def __str__(self):
        return '\n\n{0}\n\n'.format('\n'.join(self.args))


class ConversionError(InterMolError):
    """"""
    def __init__(self, functional, engine):
        Exception.__init__(self)
        self.functional = functional
        self.engine = engine


class UnsupportedConversion(ConversionError):
    """"""
    def __str__(self):
        return "{}'s are not supported in {}.".format(
            self.functional.__class__.__name__, self.engine)


class UnimplementedConversion(ConversionError):
    """"""
    def __str__(self):
        return ('{} conversion has not yet been implemented in InterMol '
                'for {}.'.format(self.functional.__class__.__name__, self.engine))


class ParsingError(InterMolError):
    """"""


class GromacsError(ParsingError):
    """"""


class DesmondError(ParsingError):
    """"""


class LammpsError(ParsingError):
    """"""
