# Wheel Generation
This directory contains scripts to generate precompiled python wheels for GridFire

# Notes
- MacOS wheels can only be generated on macos
- aarch64 wheels can only be generated on aarch64 machines
- x86_64 wheels can only be generated on x86_64 machines
- linux wheels can be generated on any linux machine, but the target architecture must match the machine architecture
- Running each script will take **a very long time** (could be upwards of half of a day depending on your system) and will require roughly 2GB of disk space
- When generating MacOS wheels, you must have all the correct versions of python installed with `pyenv`. Run the script `utils/wheels/installPyEnvVersions.sh` to install the correct versions of python.

# Usage
Once you know you are on the correct machine, run the script for your desired architecture and operating system. For example, to generate a macos x86_64 wheel, run:

```bash
./build-wheels-macos-aarch64.sh https://github.com/4D-STAR/GridFire
```

Once you have all the wheels generated (which will likely require multiple systems), copy all the wheels into a single
directory (lets assume its called `wheels` and in the root of the directory) and then run (from the root of the repository):

```bash
python -m pip install --upgrade build
python -m build --sdist --outdir wheels
twine upload wheels/*
```

Thie will also take a while (it needs to upload all the wheels to PyPI) but will result in all the wheels being uploaded to PyPI.