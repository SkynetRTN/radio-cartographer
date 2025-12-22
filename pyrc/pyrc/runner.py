import subprocess
import logging
from typing import Optional, List
from .parameters import RadioCartographerConfig

class RadioCartographerRunner:
    def __init__(self, executable_path: str = "radio-cartographer"):
        """
        Initialize the runner with the path to the executable.
        :param executable_path: Path to the radio-cartographer binary. 
                                Default assumes it's in the system PATH.
        """
        self.executable_path = executable_path
        self.logger = logging.getLogger(__name__)

    def build_command(self, config: RadioCartographerConfig, input_file: str) -> List[str]:
        """
        Constructs the full subprocess command list.
        :param config: The configuration object.
        :param input_file: Path to the input FITS file (argv[1]).
        :return: List of strings representing the command.
        """
        args = [self.executable_path, input_file]
        args.extend(config.to_args())
        return args

    def run(self, config: RadioCartographerConfig, input_file: str, 
            cwd: Optional[str] = None, check: bool = True,
            capture_output: bool = False) -> subprocess.CompletedProcess:
        """
        Runs the radio-cartographer executable with the given configuration.
        
        :param config: The configuration object.
        :param input_file: Path to the input FITS file.
        :param cwd: Working directory for the subprocess.
        :param check: If True, raise CalledProcessError on non-zero exit.
        :param capture_output: If True, capture stdout/stderr.
        :return: subprocess.CompletedProcess object.
        """
        command = self.build_command(config, input_file)
        self.logger.info(f"Executing: {' '.join(command)}")
        
        return subprocess.run(command, cwd=cwd, check=check, 
                              text=True, capture_output=capture_output)
