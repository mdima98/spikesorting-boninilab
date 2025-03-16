# Deuteron to Kilosort4 Binary Converter

## Overview

This document describes the parameters and file format for the deuteron2kilosort-bin.py script, which converts Deuteron .DTx files to binary format compatible with Kilosort4.

## TOML Configuration File

The script requires a TOML (Tom's Obvious, Minimal Language) configuration file that specifies input data sources, processing parameters, and output settings.

### Basic Structure

```toml
# Input parameters for the binary conversion script
[input_data]
path = '/path/to/source/data'
session_name = "SessionName"
tasks_names = ["TaskName1", "TaskName2"]

[deuteron]
file_type = "DT6"
start_times_s = [10.5, 20.0]
stop_times_s = [100.8, 120.5]

[output]
path = '/path/to/output/directory'
file_name = "recording.bin"
data_type = "uint16"

# Optional performance parameters
process_count = 4
buffer_size = 16777216
```

### Parameters Explanation

#### [input_data] Section

| Parameter | Type | Description | Required |
|-----------|------|-------------|----------|
| `path` | String | Root path where recording sessions are stored | Yes |
| `session_name` | String | Name of the recording session | Yes |
| `tasks_names` | String or Array | Names of recording tasks to process | Yes |

#### [deuteron] Section

| Parameter | Type | Description | Required |
|-----------|------|-------------|----------|
| `file_type` | String | Deuteron file type (DT2, DT4, DT6, DT8, DAT) | Yes |
| `start_times_s` | Number or Array | Start time in seconds for each task | Yes |
| `stop_times_s` | Number or Array | Stop time in seconds for each task | Yes |

#### [output] Section

| Parameter | Type | Description | Required |
|-----------|------|-------------|----------|
| `path` | String | Directory where binary file will be saved | Yes |
| `file_name` | String | Name of output binary file | Yes |
| `data_type` | String | Data type (uint16, int16, int32, float32) | No, defaults to uint16 |

#### Optional Root Parameters

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `process_count` | Integer | Number of parallel processes | CPU count |
| `buffer_size` | Integer | Buffer size for file operations in bytes | 16777216 (16MB) |

## File Structure Requirements

The script expects the following directory structure for input files:

```plaintext
/path/to/source/data/
    └── SessionName/
        ├── TaskName1/
        │   ├── file1.DT6
        │   ├── file2.DT6
        │   └── ...
        └── TaskName2/
            ├── file1.DT6
            ├── file2.DT6
            └── ...
```

## Deuteron File Types

| Type | Channels | Resolution | Sampling Rate | Description |
|------|----------|------------|--------------|-------------|
| DT2 | 32 | 16-bit | 32 kHz | Standard 32-channel logger |
| DT4 | 64 | 16-bit | 32 kHz | 64-channel logger |
| DT6 | 128 | 16-bit | 32 kHz | 128-channel logger |
| DT8 | 8 | 15-bit | 4 kHz | 8-channel low sampling rate logger |
| DAT | 16 | 12-bit | 31.25 kHz | Legacy 16-channel logger |

## Example Configurations

### Single Task Recording

```toml
[input_data]
path = '/data/recordings'
session_name = "Mouse01"
tasks_names = "OpenField"

[deuteron]
file_type = "DT6"
start_times_s = 10.5
stop_times_s = 600.0

[output]
path = '/data/kilosort_files'
file_name = "Mouse01_OpenField.bin"
data_type = "uint16"
```

### Multiple Tasks Combined

```toml
[input_data]
path = '/data/recordings'
session_name = "Rat02"
tasks_names = ["OpenField", "SleepRecording", "ObjectExploration"]

[deuteron]
file_type = "DT4"
start_times_s = [5.0, 10.0, 8.0]
stop_times_s = [305.0, 3600.0, 900.0]

[output]
path = '/data/kilosort_files'
file_name = "Rat02_combined.bin"
data_type = "int16"

# Use fewer processes to save memory
process_count = 2
```

## Troubleshooting

If you receive errors about parameter validation:

1. Check that all required parameters are present
2. Ensure file paths exist and are accessible
3. Verify that the number of tasks, start times, and stop times match
4. Check that start times are less than stop times
5. Make sure the specified file type is one of the supported types

## Performance Tips

- The `process_count` parameter controls parallelization. Set it lower on machines with limited RAM.
- `buffer_size` affects file I/O efficiency. The default (16MB) works well for most systems.
- For very large datasets, ensure you have at least 10% more free disk space than the estimated output size.

## Notes on Path Formatting

When specifying paths in the TOML file:

- Use forward slashes (/) even on Windows systems
- Enclose paths in single quotes to avoid escape character issues
- Absolute paths are recommended to avoid ambiguity

## Output Format

The output binary file will contain neural data in the format:

- Arranged as channels × samples matrix
- Each value stored in the specified data type (default: uint16)
- No header information (raw binary)

This format is directly compatible with Kilosort4 and similar spike sorting tools.

## Memory Considerations

The script uses memory-mapped file operations to efficiently handle large datasets without loading everything into RAM. However, processing very large recordings with many channels may still require significant memory resources. In these cases:

- Reduce `process_count` to limit concurrent memory usage
- Process tasks separately instead of combining multiple tasks
- Ensure your system has adequate swap/virtual memory configure
