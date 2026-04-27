# Makefile Run Target

## Overview

The `make run` target provides a convenient way to execute the full spatial-tk pipeline using a single TOML configuration file.

## Usage

```bash
make run ROOT=/path/to/project/directory
```

The `ROOT` directory must contain a `config.toml` file with all pipeline parameters.

## Requirements

- The `ROOT` directory must exist
- The `ROOT` directory must contain `config.toml`
- The config file must have sections for all commands: `[concat]`, `[normalize]`, `[cluster]`, `[annotate]`, and `[differential]`

## What It Does

The `run` target executes all five pipeline steps sequentially:

1. **Concatenate samples** - Reads `[concat]` section from config
2. **Normalize data** - Reads `[normalize]` section from config
3. **Cluster cells** - Reads `[cluster]` section from config
4. **Annotate cell types** - Reads `[annotate]` section from config
5. **Differential expression** - Reads `[differential]` section from config

Each step runs from the `ROOT` directory, so paths in the config file should be relative to that directory.

## Example

```bash
# Project structure:
# projects/PDAC_HIV/
#   ├── config.toml
#   ├── samples.csv
#   └── markers.csv

# Run full pipeline
make run ROOT=projects/PDAC_HIV
```

## Error Handling

The target stops immediately if:
- ROOT is not specified
- ROOT directory doesn't exist
- config.toml is missing
- Any pipeline step fails

## Benefits

- **Reproducibility**: All parameters stored in version-controlled config file
- **Simplicity**: Single command runs entire pipeline
- **Consistency**: Same config used for all steps
- **Convenience**: No need to remember all command-line arguments

## See Also

- `notes/TOML_CONFIG_GUIDE.md` - Complete guide to TOML config files
- `example_config.toml` - Example configuration file
- `README.md` - General usage documentation

