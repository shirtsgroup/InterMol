class MultipleValidationErrors(Exception):
    """"""
    def __str__(self):
        return '\n\n{0}\n\n'.format('\n'.join(self.args))


class ConversionError(Exception):
    """"""
    def __init__(self, functional, engine):
        Exception.__init__(self)
        self.functional = functional
        self.engine = engine


class UnsupportedConversion(ConversionError):
    """"""
    def __str__(self):
        return '{} are not supported in {}.'.format(
            self.functional.__name__, self.engine)


class UnimplementedConversion(Exception):
    """"""
    def __str__(self):
        return ('{} conversion has not yet been implemented in InterMol '
                'for {}.'.format(self.functional.__name__, self.engine))


class ParsingError(Exception):
    """"""


class GromacsError(ParsingError):
    """"""


class DesmondError(ParsingError):
    """"""


class LammpsError(ParsingError):
    """"""
