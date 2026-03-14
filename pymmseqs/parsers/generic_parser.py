# pymmseqs/parsers/generic_parser.py


class GenericParser:
    """Generic parser for auto-generated MMseqs2 command wrappers.

    Finds the output path from the config's YAML defaults by locating
    the first path parameter with should_exist=False.
    """

    def __init__(self, config):
        self._output_path = None
        for param_name, param_info in config._defaults.items():
            if param_info['type'] == 'path' and not param_info['should_exist']:
                self._output_path = str(getattr(config, param_name))
                break

    def to_path(self) -> str:
        """Return the output path prefix.

        Returns
        -------
        str
            The resolved output path for the command results.
        """
        return self._output_path
