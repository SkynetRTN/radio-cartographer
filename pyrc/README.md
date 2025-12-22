# Radio Cartographer Python Wrapper (PyRC)

This directory contains `PyRC` (Python Radio Cartographer, pronounced "perk"), the Python wrapper for the Radio Cartographer C++ library.

## Prerequisites

Install the required Python packages:

```bash
pip install -r requirements.txt
```

## Usage

The `radio_cartographer` package allows you to configure and run the cartographer from Python.

### Example

See `examples/run_example.py` for a complete example of how to configure and run the cartographer.

To run the example:

```bash
python3 examples/run_example.py
```

You can also specify the path to the binary:

```bash
python3 examples/run_example.py --bin /path/to/radio-cartographer
```

### Basic Import structure

```python
from radio_cartographer import (
    RadioCartographerConfig,
    RadioCartographerRunner,
    Channel,
    CalibrationMethod,
    CoordinateSystem,
    # ... other enums
)

config = RadioCartographerConfig(...)
runner = RadioCartographerRunner("/path/to/executable")
runner.run(config, "input_file.fits")
```

## Docker Environment

If you are developing or running within the provided Docker container:

### Default Build Directory

The default build directory where the C++ project is compiled is:

```
/skynet/radio-cartographer
```

The executable is usually found at `/skynet/radio-cartographer/radio-cartographer`.

### Python Path

The Docker container is configured with the correct `PYTHONPATH`:

```bash
export PYTHONPATH="${PYTHONPATH}:/skynet/radio-cartographer/pyrc"
```

This ensures that `import radio_cartographer` works correctly without additional configuration.
