"""Exception types raised by corset-tools."""


class CorsetToolsError(Exception):
    """Base class for all errors raised by corset-tools."""


class DuplicateTranscriptError(CorsetToolsError):
    """Raised when the same transcript name occurs in more than one input set."""


class InputFormatError(CorsetToolsError):
    """Raised when an input file cannot be parsed in the expected format."""
