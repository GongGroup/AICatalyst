class FileFormatError(IOError):
    pass


class StructureError(TypeError):
    pass


class ParseError(FileFormatError):
    pass
