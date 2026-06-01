# input/

Drop your Skynet SDFITS files here (`*.fits`). `run.py` (in the repo
root) processes every `.fits` file in this directory and writes the
results to `output/`.

```
python run.py                       # process every .fits in input/
python run.py input/my_target.fits  # process a specific file
```

This directory is tracked by git (so the directory exists for new
clones) but the `*.fits` files inside it are not — see `.gitignore`.
Each lab member's data stays local.
