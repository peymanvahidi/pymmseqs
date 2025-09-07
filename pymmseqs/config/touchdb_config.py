from pathlib import Path
from typing import Union

from .base import BaseConfig
from ..defaults import loader
from ..utils import (
    get_caller_dir,
    run_mmseqs_command
)

DEFAULTS = loader.load("touchdb")


class TouchDBConfig(BaseConfig):
    """
    Load an MMseqs2 database index into memory using the touchdb module.

    This speeds up subsequent searches by ensuring the index is resident in RAM
    (via the OS page cache / mmap touching).
    """

    def __init__(
        self,
        # Required parameters
        sequence_db: Union[str, Path],
    ):
        super().__init__()

        # Required parameter
        self.sequence_db = Path(sequence_db)

        self._defaults = DEFAULTS
        self._path_params = [param for param, info in DEFAULTS.items() if info['type'] == 'path']
        self._caller_dir = get_caller_dir()

    def _validate(self) -> None:
        self._check_required_files()
        self._validate_choices()

    def run(self) -> None:
        self._resolve_all_path(self._caller_dir)

        self._validate()

        args = self._get_command_args("touchdb")
        mmseqs_output = run_mmseqs_command(args)

        self._handle_command_output(
            mmseqs_output=mmseqs_output,
            output_identifier="TouchDB",
            output_path=str(self.sequence_db)
        )


