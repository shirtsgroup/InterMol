class MultipleValidationErrors(Exception):
    """"""
    def __str__(self):
        return '\n\n{0}\n\n'.format('\n'.join(self.args))


class UnsupportedConversion(Exception):
    """"""


class UnimplementedConversion(Exception):
    """"""


class ParsingError(Exception):
    """"""


class GromacsError(ParsingError):
    """"""


class DesmondError(ParsingError):
    """"""


class LammpsError(ParsingError):
    """"""
