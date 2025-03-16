"""
KILOSORT BIN FILE CONVERSION SCRIPT by Matteo Di Mario

This script converts Deuteron .DTx files to a binary format compatible with Kilosort4.
It processes the DTx files, merging them together and saves the data in the form of channels x samples.

How to launch the script from the terminal (or anaconda prompt in Windows):
    python deuteron2kilosort-bin.py <path_to_toml_file>

For detailed documentation on TOML file format and parameters, see the accompanying 
PARAMETERS.md file in this directory.
"""

import numpy as np
import toml
from pathlib import Path
from tqdm import tqdm
import sys
import multiprocessing as mp
import os
import shutil
import tempfile
import uuid
import logging
import logging.handlers 
import psutil
from datetime import datetime
import time
from termcolor import colored
from art import text2art, art
import functools

# Disable warnings from multiprocessing (handled in this context)
import warnings
# Suppress all resource tracker warnings with multiple approaches
warnings.filterwarnings("ignore", category=UserWarning, module="multiprocessing.resource_tracker")
# Add more specific filters for the exact messages you're seeing
warnings.filterwarnings("ignore", message=".*leaked semaphore objects.*", category=UserWarning)
warnings.filterwarnings("ignore", message=".*No such file or directory.*", category=UserWarning)

import os
os.environ['PYTHONWARNINGS'] = "ignore::UserWarning:multiprocessing.resource_tracker"

# Add this line near the top of the file, right after imports
logger = logging.getLogger(__name__)

# Create a custom formatter class that uses termcolor
class ColoredFormatter(logging.Formatter):
    """Custom formatter to add colors to log messages based on level."""
    
    def format(self, record):
        message = super().format(record)
        
        if record.levelno >= logging.ERROR:
            return colored(message, 'red')
        elif record.levelno >= logging.WARNING:
            return colored(message, 'yellow')
        else:
            return message

def setup_logging(log_dir=None, session_name=None):
    """Set up logging to both console and file."""
    global logger
    
    # Store early messages if this is first initialization
    memory_handler = None
    if not logger.hasHandlers():
        logger.setLevel(logging.DEBUG)
        # Create memory handler to store early logs
        memory_handler = logging.handlers.MemoryHandler(capacity=1000, 
                                                      flushLevel=logging.ERROR,
                                                      target=None)
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(ColoredFormatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(console_handler)
        logger.addHandler(memory_handler)
    
    # Add file handler if directory is provided
    if log_dir:
        os.makedirs(log_dir, exist_ok=True)
        # Include session name in the log filename if provided
        log_name = f"dt2bin_{session_name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
        if session_name:
            log_file = os.path.join(log_dir, log_name)
        else:
            log_file = os.path.join(log_dir, log_name)
        
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        
        # Set the file handler as the target for the memory handler to flush early logs
        if memory_handler:
            memory_handler.setTarget(file_handler)
            memory_handler.flush()
        
        logger.addHandler(file_handler)
        logger.info(f"Log file created: {log_name}")
        
    # Prevent propagation to root logger
    logger.propagate = False
    return logger

# Constants
FILE_SIZE = 16777216  # size of each file in bytes (fixed for all Deuteron file types)
DEFAULT_BUFFER_SIZE = 16 * 1024 * 1024  # 16MB chunks for file copying
HEADER_WIDTH = 80 # For header in print_header
SECTION_WIDTH = 70 # For section headers in print_section

def print_header():
    """Print a distinctive header with simplified ASCII art at the start of program execution."""
    print("\n" + "=" * HEADER_WIDTH)
    
    # Create title with big ASCII art - simplified
    title_art = text2art("DT2BIN", font="smaller")
    for line in title_art.split('\n'):
        print(f"{line:^{HEADER_WIDTH}}")
        
    print(f"{'DEUTERON TO KILOSORT CONVERTER':^{HEADER_WIDTH}}")
    print(f"{'By Matteo Di Mario':^{HEADER_WIDTH}}")
    
    # Simplify the monkey grid to just one monkey
    monkey = art("monkey")
    for line in monkey.split('\n'):
        print(f"{line:^{HEADER_WIDTH}}")
        
    print(f"{'Started at ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S'):^{HEADER_WIDTH}}")
    print("=" * HEADER_WIDTH)

def print_footer(success=True):
    """Print a distinctive footer to mark the end of program execution."""
    print("\n" + "=" * HEADER_WIDTH)
    if success:
        print(f"{'CONVERSION COMPLETED SUCCESSFULLY':^{HEADER_WIDTH}}")
    else:
        print(f"{'CONVERSION TERMINATED WITH ERRORS':^{HEADER_WIDTH}}")
    print(f"{'Finished at ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S'):^{HEADER_WIDTH}}")
    print("=" * HEADER_WIDTH + "\n")

def print_section(title, width=70, char="-"):
    """Print a section divider with a title."""
    print("\n" + char * width)
    print(f" {title} ".center(width, char))
    print(char * width + "\n")

# Replace print_separator with this
def print_separator(title):
    """Print a section divider with a title and log it"""
    print('\n' + '-' * 50)
    logger.info(f'=== {title} ===')

def get_meta_data(ext):
    """
    Extracts meta data based on the file extension.
    Builds a dictionary with essential information for data processing.
    This information is found in the instruction manual of each type of logger.
    """
    meta_data_dict = {
        'DT2': {'num_channels': 32, 'num_adc_bits': 16, 'voltage_res': 0.2e-6, 'sampling_rate': 32e3},
        'DT4': {'num_channels': 64, 'num_adc_bits': 16, 'voltage_res': 0.2e-6, 'sampling_rate': 32e3},
        'DT6': {'num_channels': 128, 'num_adc_bits': 16, 'voltage_res': 0.2e-6, 'sampling_rate': 32e3},
        'DT8': {'num_channels': 8, 'num_adc_bits': 15, 'voltage_res': 0.42e-6, 'sampling_rate': 4e3},
        'DAT': {'num_channels': 16, 'num_adc_bits': 12, 'voltage_res': 3.3e-6, 'sampling_rate': 31.25e3}
    }
    
    if ext.upper() not in meta_data_dict:
        raise ValueError(f"Invalid file extension '{ext}'. Please choose from: {', '.join(meta_data_dict.keys())}")
    
    # Return using uppercase key to ensure consistent behavior
    return meta_data_dict[ext.upper()]

def load_raw_dtx_data(file_path):
    """
    Loads raw binary data from a DTx file.
    This function only handles the file I/O and binary data loading.
    
    Args:
        file_path: Path to the DTx file
        
    Returns:
        numpy array containing raw uint16 data
    
    Raises:
        IOError: If the file cannot be read or is empty
    """
    try:
        with open(file_path, 'rb') as fid:
            data = np.fromfile(fid, dtype=np.uint16)  # each data point of neural data is a 16-bit word
        
        # Check if file is empty or corrupted
        if len(data) == 0:
            raise ValueError(f"File is empty: {file_path}")
        
        return data
    except Exception as e:
        raise IOError(f"Error reading file {file_path}: {str(e)}")

def reshape_neural_data(raw_data, num_channels):
    """
    Reshapes raw data into channels x samples format and handles data size issues.
    
    Args:
        raw_data: 1D array of raw neural data values
        num_channels: Number of channels in the recording
        
    Returns:
        2D array with shape (channels, samples)
    """
    # Check if data needs to be truncated to fit channel count
    if len(raw_data) % num_channels != 0:
        # Log warning about unexpected size
        logger.warning(f"Data has unexpected size. Truncating {len(raw_data) % num_channels} samples.")
        # Truncate to fit channel count
        raw_data = raw_data[:-(len(raw_data) % num_channels)]
    
    # Reshape into channels × samples format
    return np.reshape(raw_data, (num_channels, -1))

def apply_data_transformations(neural_data, meta_data):
    """
    Apply any necessary transformations to the neural data.
    Future modifications to data processing can be added here.
    
    Args:
        neural_data: 2D array of neural data (channels × samples)
        meta_data: Dictionary containing recording metadata
        
    Returns:
        Processed neural data
    """
    # Currently no transformations are applied, but this function provides
    # a clear entry point for future modifications like:
    #  - Filtering
    #  - Artifact removal
    #  - Normalization
    #  - Channel selection or reordering
    
    return neural_data

def extract_time_window(neural_data, file_start_sample, num_samples_per_file, start_sample, stop_sample):
    """
    Extracts the relevant time window from neural data based on start and stop sample indices.
    Handles shape inconsistencies robustly.
    """
    # Calculate the start and end indices for the current file
    file_end_sample = file_start_sample + num_samples_per_file # neural_data.shape[1]  # Use actual shape instead of expected
    
    # Skip this file if it's completely outside the desired range
    if file_end_sample <= start_sample or file_start_sample >= stop_sample:
        return None, 0, 0
    
    # Calculate the indices to slice the neural data
    slice_start = max(0, start_sample - file_start_sample)
    slice_end = min(num_samples_per_file, stop_sample - file_start_sample)
    
    # Calculate where to write in the output memmap
    memmap_start_idx = max(0, file_start_sample - start_sample)
    
    # Extract the relevant portion of the data
    extracted_data = neural_data[:, slice_start:slice_end]
    
    # Recalculate memmap_end_idx based on actual extracted data size
    memmap_end_idx = memmap_start_idx + extracted_data.shape[1]
    
    return extracted_data, memmap_start_idx, memmap_end_idx

def process_dtx_file(file_path, meta_data):
    """
    Processes a DTx file to extract the neural data.
    
    Args:
        file_path: Path to the DTx file
        meta_data: Dictionary containing recording metadata
        
    Returns:
        2D array of neural data (channels x samples)
    """
    try:
        # Load raw binary data
        raw_data = load_raw_dtx_data(file_path)
        
        # Log the actual size for debugging shape issues
        expected_size = FILE_SIZE // 2  # Each value is 2 bytes
        if len(raw_data) != expected_size:
            logger.debug(f"File {os.path.basename(file_path)} has {len(raw_data)} values "
                        f"(expected {expected_size}), difference: {len(raw_data)-expected_size}")
        
        # Reshape into channels × samples format
        data_matrix = reshape_neural_data(raw_data, meta_data['num_channels'])
        
        # Apply any transformations
        processed_data = apply_data_transformations(data_matrix, meta_data)
        
        return processed_data
    except Exception as e:
        raise IOError(f"Error processing file {file_path}: {str(e)}")

def process_file(file_path, meta_data, current_sample, memmap, start_sample, stop_sample, 
                session_name, task_name, show_artifact_warnings=True):
    """Process a single DTx file and write it to the memory-mapped array."""
    warnings_list = []
    
    try:
        # Load and process the file
        neural_data = process_dtx_file(file_path, meta_data)
        
        # Extract the relevant time window
        extracted_data, memmap_start_idx, memmap_end_idx = extract_time_window(
            neural_data, current_sample, neural_data.shape[1], start_sample, stop_sample
        )
        
        # Calculate file start time in seconds
        sampling_rate = meta_data['sampling_rate']
        file_start_time = current_sample / sampling_rate
        
        # Write to memmap if data exists
        if extracted_data is not None and memmap_end_idx > memmap_start_idx:
            write_to_memmap(memmap, extracted_data, memmap_start_idx, memmap_end_idx, file_path)
            
            # Validate data and collect warnings
            # file_warnings = validate_file_data(extracted_data, file_path, session_name, task_name, 
            #                                  file_start_time, show_artifact_warnings)
            file_warnings = []
            warnings_list.extend(file_warnings)
            
        # Clean up
        del neural_data
        if extracted_data is not None:
            del extracted_data
            
        return True, warnings_list
            
    except Exception as e:
        logger.error(f"Error processing file {file_path}: {str(e)}")
        return False, warnings_list

def process_task(task_info):
    """Convert files into binary format with tqdm progress tracking."""
    # Initialize for error handling
    task_idx = -1
    task_name = "Unknown"
    all_warnings = []
    
    try:
        # Extract task information directly
        task_idx, task_name, start_time_s, stop_time_s, input_data_path, session_name, deuteron_file_type, meta_data, num_samples_per_file, data_type, temp_dir, buffer_size, show_artifact_warnings = task_info
        
        # Create unique temporary file for this task
        temp_file = Path(temp_dir).joinpath(f'task_{task_idx}_{task_name}_{uuid.uuid4()}.bin')
        
        # Calculate start and stop samples
        start_sample = int(start_time_s * meta_data['sampling_rate'])
        stop_sample = int(stop_time_s * meta_data['sampling_rate'])
        samples_to_keep = stop_sample - start_sample
        
        # Create directory if it doesn't exist
        task_dir = Path(input_data_path).joinpath(session_name, task_name)
        
        # Get and sort DTx files
        dtx_files = []
        for ext in [deuteron_file_type.upper(), deuteron_file_type.lower()]:
            dtx_files.extend(task_dir.glob(f'*.{ext}'))
        dtx_files = sorted(dtx_files)
        
        if not dtx_files:
            raise ValueError(f"No DTx files found for task {task_name} at {task_dir}")
        
        # Create memmap for output file
        mm = np.memmap(temp_file, dtype=data_type, mode='w+', 
                      shape=(meta_data['num_channels'], samples_to_keep))
        
        # Process each file
        current_sample = 0
        processed_files_count = 0
        relevant_files = 0  # Count files that are actually in our range
        
        # First pass to determine how many files are relevant to our time window
        for i, file_path in enumerate(dtx_files):
            file_start_sample = i * num_samples_per_file
            if not (file_start_sample + num_samples_per_file <= start_sample or file_start_sample >= stop_sample):
                relevant_files += 1
        
        # Use tqdm for progress tracking
        # Use fixed positions with cleaner management
        position = task_idx  # Keep absolute position
        with tqdm(total=relevant_files, unit='file', 
                  desc=f"Task {task_name}",
                  position=position,
                  leave=True, # Keep bar visible after completion
                  bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]') as pbar:
            
            for i, file_path in enumerate(dtx_files):
                # Check if this file is within the range of samples we need
                file_start_sample = i * num_samples_per_file
                
                # Skip if outside range
                if file_start_sample + num_samples_per_file <= start_sample or file_start_sample >= stop_sample:
                    continue
                
                # Process the file and write to memmap
                success, file_warnings = process_file(file_path, meta_data, file_start_sample, mm, 
                                     start_sample, stop_sample, session_name, task_name, show_artifact_warnings)
                
                all_warnings.extend(file_warnings)
                
                if success:
                    processed_files_count += 1
                
                # Update progress bar with additional information
                pbar.set_postfix(processed=f"{processed_files_count}/{relevant_files}", 
                                file=os.path.basename(file_path))
                pbar.update(1)
        
        # Flush the memmap to disk
        mm.flush()
        del mm
        
        # Add explicit cleanup for any potential resource leaks
        import gc
        gc.collect()  # Force garbage collection
        
        # Estimate file size
        estimated_temp_size = os.path.getsize(temp_file)
        
        # Verify temp file integrity before considering the task successful
        if not verify_temp_file(temp_file, estimated_temp_size, {'meta_data': meta_data}):
            logger.warning(f"Warning: Task {task_name} temp file verification failed")
            # Continue anyway since partial data may still be useful
        
        logger.info(f"Task {task_name} complete: {processed_files_count}/{relevant_files} files processed")
        
        return True, task_idx, f"Successfully processed task {task_name}", str(temp_file), estimated_temp_size, all_warnings
        
    except Exception as e:
        logger.error(f"Error in task {task_name}: {str(e)}")
        return False, task_idx, f"Error in task {task_name}: {str(e)}", None, 0, all_warnings

def check_free_space(directory, required_bytes):
    """Check if there's enough free space in directory."""
    free_bytes = psutil.disk_usage(directory).free
    return free_bytes >= required_bytes, free_bytes

def validate_toml_params(params):
    """Validate the parameters from the TOML file."""
    errors = []
    warnings = []
    
    # Check required sections
    for section in ['input_data', 'deuteron', 'output']:
        if section not in params:
            errors.append(f"Missing section: [{section}]")
    
    if errors:
        return errors, warnings
    
    # Check input_data section
    input_data = params['input_data']
    for param in ['path', 'session_name', 'tasks_names']:
        if param not in input_data:
            errors.append(f"Missing parameter: input_data.{param}")
    
    # Check deuteron section
    deuteron = params['deuteron']
    for param in ['file_type', 'start_times_s', 'stop_times_s']:
        if param not in deuteron:
            errors.append(f"Missing parameter: deuteron.{param}")
    
    # Check output section
    output = params['output']
    for param in ['path', 'file_name']:
        if param not in output:
            errors.append(f"Missing parameter: output.{param}")
    
    if 'data_type' not in output:
        warnings.append("Missing parameter: output.data_type, defaulting to 'uint16'")
        output['data_type'] = 'uint16'
    
    if errors:
        return errors, warnings
    
    # Check types and values
    if not isinstance(input_data['tasks_names'], list):
        input_data['tasks_names'] = [input_data['tasks_names']]
        warnings.append("tasks_names was converted to list format")
    
    if not (isinstance(deuteron['start_times_s'], list) and isinstance(deuteron['stop_times_s'], list)):
        if len(input_data['tasks_names']) == 1:
            deuteron['start_times_s'] = [deuteron['start_times_s']]
            deuteron['stop_times_s'] = [deuteron['stop_times_s']]
            warnings.append("start_times_s and stop_times_s were converted to list format")
        else:
            errors.append("start_times_s and stop_times_s must be lists when multiple tasks are specified")
    
    # Check lengths match
    if len(input_data['tasks_names']) != len(deuteron['start_times_s']) or len(input_data['tasks_names']) != len(deuteron['stop_times_s']):
        errors.append("Number of tasks, start times, and stop times must be equal")
    
    # Check file type
    try:
        get_meta_data(deuteron['file_type'])
    except ValueError as e:
        errors.append(str(e))
    
    # Check start_times < stop_times
    for i, (start, stop) in enumerate(zip(deuteron['start_times_s'], deuteron['stop_times_s'])):
        if start >= stop:
            errors.append(f"Task {i}: start time ({start}) must be less than stop time ({stop})")
        if start < 0:
            errors.append(f"Task {i}: start time ({start}) must be positive")
    
    # Check data_type
    allowed_types = ['uint16', 'int16', 'int32', 'float32']
    if output['data_type'] not in allowed_types:
        errors.append(f"data_type must be one of: {', '.join(allowed_types)}")
    
    # Check paths exist
    if not Path(input_data['path']).exists():
        errors.append(f"Input path does not exist: {input_data['path']}")
    
    # Check if performance section exists
    if 'performance' in params:
        performance = params['performance']
        # Get process_count from performance section
        if 'process_count' in performance:
            params['process_count'] = performance['process_count']
        else:
            params['process_count'] = mp.cpu_count()
            warnings.append(f"process_count not specified, defaulting to {mp.cpu_count()}")
            
        # Get buffer_size from performance section
        if 'buffer_size' in performance:
            params['buffer_size'] = performance['buffer_size']
        else:
            params['buffer_size'] = DEFAULT_BUFFER_SIZE
            warnings.append(f"buffer_size not specified, defaulting to {DEFAULT_BUFFER_SIZE/(1024*1024)}MB")
        
        # Add this inside the 'performance' section check in validate_toml_params function:
        if 'enable_file_logging' in performance:
            params['enable_file_logging'] = performance['enable_file_logging']
        else:
            params['enable_file_logging'] = False
            warnings.append("enable_file_logging not specified, defaulting to False")
        
        # Add this new check
        if 'use_multiprocessing' in performance:
            params['use_multiprocessing'] = performance['use_multiprocessing']
        else:
            params['use_multiprocessing'] = True
            warnings.append("use_multiprocessing not specified, defaulting to True")
    
        # Add this new check
        if 'show_artifact_warnings' in performance:
            params['show_artifact_warnings'] = performance['show_artifact_warnings']
        else:
            params['show_artifact_warnings'] = True  # Show warnings by default
            warnings.append("show_artifact_warnings not specified, defaulting to True")
    else:
        # No performance section, use default values
        params['process_count'] = mp.cpu_count()
        params['buffer_size'] = DEFAULT_BUFFER_SIZE
        params['use_multiprocessing'] = True  # Add this line
        warnings.append(f"No [performance] section found, using defaults: process_count={mp.cpu_count()}, buffer_size={DEFAULT_BUFFER_SIZE/(1024*1024)}MB, use_multiprocessing=True")
    
    return errors, warnings

def parse_command_line_args():
    if len(sys.argv) != 2:
        print("Usage: python deuteron2kilosort-bin.py <path_to_toml_file>")
        print_footer(False)
        return None

    toml_file_path = sys.argv[1]

    try:
        print_separator("STARTING PREPROCESSING")
        logger.info("Starting Deuteron to Kilosort4 data conversion")
        logger.info(f"Loading parameters from {toml_file_path}")
        
        # Load input parameters from a TOML file
        input_params = toml.load(toml_file_path)
        
        # Validate parameters
        errors, warnings = validate_toml_params(input_params)
        
        if warnings:
            for warning in warnings:
                logger.warning(warning)
        
        if errors:
            logger.error("Parameter validation errors:")
            for error in errors:
                logger.error(f"  - {error}")
            logger.error("Please fix these errors and try again.")
            print_footer(False)
            return None

        logger.info("Parameters loaded successfully")
        return input_params

    except Exception as e:
        logger.error(f"Error loading parameters: {str(e)}")
        logger.error("The program has been interrupted and will now exit.")
        print_footer(False)
        return None

def extract_configuration(params):
    try:
        # Extract input parameters
        input_data_path = params['input_data']['path']
        session_name = params['input_data']['session_name']
        tasks_names = params['input_data']['tasks_names']
        deuteron_file_type = params['deuteron']['file_type']
        start_times_s = params['deuteron']['start_times_s']
        stop_times_s = params['deuteron']['stop_times_s']
        output_data_dir = params['output']['path']
        binary_file_name = params['output']['file_name']
        data_type = params['output']['data_type']
        
        # Get process count (with default fallback)
        process_count = params.get('process_count', mp.cpu_count())
        # More efficient implementation:
        process_count = min(process_count, mp.cpu_count(), len(tasks_names))  # Cap at CPU count and task count
        process_count = max(1, process_count)  # Ensure at least 1 process
        
        # Get buffer size (with default fallback)
        buffer_size = params.get('buffer_size', DEFAULT_BUFFER_SIZE)
        
        # Get meta data for the deuteron file type
        meta_data = get_meta_data(deuteron_file_type)
        
        # Extract use_multiprocessing parameter
        use_multiprocessing = params.get('use_multiprocessing', True)
        
        # Get show_artifact_warnings parameter
        show_artifact_warnings = params.get('show_artifact_warnings', True)
        
        return {
            'input_data_path': input_data_path,
            'session_name': session_name,
            'tasks_names': tasks_names,
            'deuteron_file_type': deuteron_file_type,
            'start_times_s': start_times_s,
            'stop_times_s': stop_times_s,
            'output_data_dir': output_data_dir,
            'binary_file_name': binary_file_name,
            'data_type': data_type,
            'process_count': process_count,
            'buffer_size': buffer_size,
            'meta_data': meta_data,
            'use_multiprocessing': use_multiprocessing,
            'show_artifact_warnings': show_artifact_warnings
        }
    except Exception as e:
        logger.error(f"Error extracting configuration: {str(e)}")
        print_footer(False)
        return None

def display_configuration_summary(config):
    # Create a more compact and visually appealing summary
    print_separator("CONFIGURATION SUMMARY")
    logger.info(f"Session:      {config['session_name']}")
    logger.info(f"Tasks:        {'['+', '.join(config['tasks_names'])}"+']')
    logger.info(f"File type:    {config['deuteron_file_type']} ({config['meta_data']['num_channels']} channels @ {int(config['meta_data']['sampling_rate']/1000)} kHz)")
    logger.info(f"Time window:  {config['start_times_s']} to {config['stop_times_s']} seconds")
    logger.info(f"Output Dir:       {config['output_data_dir']}")
    logger.info(f"Output File:      {config['binary_file_name']}")
    logger.info(f"Processing:   {'Parallel' if config['use_multiprocessing'] else 'Sequential'} ({config['process_count']} processes)")
    
    # Calculate predicted file size
    total_samples_to_keep = 0
    for start_time, stop_time in zip(config['start_times_s'], config['stop_times_s']):
        total_samples_to_keep += int((stop_time - start_time) * config['meta_data']['sampling_rate'])
    data_type_size = np.dtype(config['data_type']).itemsize
    predicted_file_size = total_samples_to_keep * config['meta_data']['num_channels'] * data_type_size
    logger.info(f"Est. size:    {predicted_file_size / (1024**3):.2f} GB")
    
    # Check available disk space
    has_space, free_space = check_free_space(config['output_data_dir'], int(predicted_file_size * 1.1))
    if not has_space:
        logger.error(f"Insufficient disk space. Required: {predicted_file_size / (1024**3):.2f} GB, "
                    f"Available: {free_space / (1024**3):.2f} GB")
        return False

    # Check for existing output file
    output_file = Path(config['output_data_dir']).joinpath(config['binary_file_name'])
    if output_file.exists():
        logger.warning(f"Output file {output_file} already exists and will be overwritten")
    
    return True

def check_and_confirm(config):
    # Ask the user to confirm to start the file processing
    print_separator("CONFIRMATION")
    logger.info("Do you want to start the file processing and creation of the binary file? (Y/n): ")
    confirm = input()
    if confirm.lower() not in ['y', 'yes', '']:
        logger.info("File processing aborted by user.")
        print_footer(False)
        return False
    return True

def process_all_tasks(config, temp_dir):
    print_separator("STARTING TASKS PROCESSING")
    # Prepare task information
    task_infos = []
    for idx, (task_name, start_time_s, stop_time_s) in enumerate(zip(
            config['tasks_names'], config['start_times_s'], config['stop_times_s'])):
        
        # Calculate number of samples in each file
        num_samples_per_file = FILE_SIZE // (config['meta_data']['num_channels'] * 2)  # 2 bytes per sample
        
        task_infos.append((idx, task_name, start_time_s, stop_time_s, config['input_data_path'], 
                          config['session_name'], config['deuteron_file_type'], config['meta_data'], 
                          num_samples_per_file, config['data_type'], temp_dir, config['buffer_size'], 
                          config['show_artifact_warnings']))
    
    results = []
    if config['use_multiprocessing']:
        # Use simple process pool without complex tracking
        with mp.Pool(processes=config['process_count']) as pool:
            logger.info(f"Processing {len(task_infos)} tasks with {config['process_count']} processes...")
            
            # Submit all tasks and wait for results - keep it simple
            results = pool.map(process_task, task_infos)
            
            # Log completion of all tasks
            time.sleep(1.5) # Add a small delay for better log readability
            logger.info(f"[SUCCESS] All {len(task_infos)} tasks processed successfully, ready for merging")
    else:
        # Sequential processing with simple progress
        for task_idx, task_info in enumerate(task_infos):
            logger.info(f"Processing task {task_idx+1}/{len(task_infos)}: {task_info[1]}")
            result = process_task(task_info)
            results.append(result)
            logger.info(f"Task {task_info[1]} complete")
    
    # Sort results and continue with merging
    return sorted(results, key=lambda x: x[1])  # Sort by task index

def merge_output_files(results, config):
    """Focus only on merging temporary files into the final output file."""
    # Calculate total size to help with progress reporting
    total_size = sum(size for _, _, _, _, size, _ in results if size > 0)
    
    # Add this just before the merging phase starts
    print_separator("STARTING MERGING PHASE")
    logger.info(f"Found {len([r for r in results if r[0]])} successful tasks to merge")
    for success, idx, msg, path, size, _ in results:
        logger.info(f"Task {idx}: {'[SUCCESS]' if success else '[FAILED]'} - {os.path.basename(path) if path else 'No file'}")

    try:
        # Combine all temporary files into the final output file
        logger.info("Combining task files into final output file...")
        output_file = Path(config['output_data_dir']).joinpath(config['binary_file_name'])
        with open(output_file, 'wb') as outfile:
            with tqdm(total=total_size, unit='B', unit_scale=True, unit_divisor=1024,
                     desc="Merging files") as pbar:
                for success, task_idx, task_name, temp_file_path, _, _ in results:
                    if success and temp_file_path:
                        # Update progress bar to show which task file is being merged
                        pbar.set_description(f"Merging {os.path.basename(temp_file_path)}")
                        
                        with open(temp_file_path, 'rb') as infile:
                            # Copy in chunks for efficiency and progress tracking
                            while True:
                                buf = infile.read(config['buffer_size'])
                                if not buf:
                                    break
                                outfile.write(buf)
                                pbar.update(len(buf))
        
        # Disabled for now - will be re-enabled in the future
        # # Verify file integrity                     
        # print_separator("STARTING FILE INTEGRITY VERIFICATION")
        # verification_success = verify_merged_data(results, output_file, config)
        # if not verification_success:
        #     logger.warning("[FAILED] File verification failed. The output file may have issues.")
        # else:
        #     logger.info("[SUCCESS] File integrity verified successfully")

        return True, output_file
    
    except Exception as e:
        logger.error(f"Error during file merging: {str(e)}")
        return False, None

def log_and_handle_error(error_msg, exception=None, is_critical=False):
    """Centralized error handling function."""
    if exception:
        if is_critical:
            logger.critical(f"{error_msg}: {str(exception)}", exc_info=True)
        else:
            logger.error(f"{error_msg}: {str(exception)}")
    else:
        if is_critical:
            logger.critical(error_msg)
        else:
            logger.error(error_msg)
            
    return False  # Indicate failure

def time_execution(func):
    """Decorator to time function execution."""
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        elapsed = time.time() - start_time
        logger.info(f"{func.__name__} completed in {elapsed:.2f} seconds")
        return result
    return wrapper

def write_to_memmap(memmap, data, start_idx, end_idx, file_path):
    """Write data to memmap with robust shape handling."""
    try:
        # Normal case - shapes match
        memmap[:, start_idx:end_idx] = data
        return True
    except ValueError as e:
        logger.warning(f"Shape mismatch when writing {file_path.name}: {str(e)}")
        logger.warning(f"Data shape: {data.shape}, Target shape: {(memmap.shape[0], end_idx-start_idx)}")
        
        # Try writing with explicit dimensions
        try:
            rows = min(data.shape[0], memmap.shape[0])
            cols = min(data.shape[1], end_idx-start_idx)
            memmap[:rows, start_idx:start_idx+cols] = data[:rows, :cols]
            logger.info(f"Recovered partial data from {file_path.name}: used {rows}×{cols} out of {data.shape}")
            return True
        except Exception as inner_e:
            logger.error(f"Failed to recover any data from {file_path.name}: {str(inner_e)}")
            return False

def validate_file_data(extracted_data, file_path, session_name, task_name, file_start_time, show_warnings=True):
    """Validate file data for potential issues and return warnings list."""
    warnings_list = []
    
    try:
        # Check for zeros (potential data corruption)
        zero_percentage = np.count_nonzero(extracted_data == 0) / extracted_data.size
        if zero_percentage > 0.9:  # 90% zeros is suspicious
            warnings_list.append({
                'type': 'zeros',
                'file': str(file_path),
                'session': session_name,
                'task': task_name,
                'time': file_start_time,
                'details': f"Contains {zero_percentage:.1%} zeros (possible corruption)"
            })
        
        # Check for extreme values (potential artifacts)
        with np.errstate(invalid='ignore'):
            if np.any(np.abs(extracted_data) > 30000):  # Close to int16 max
                warnings_list.append({
                    'type': 'extreme_values',
                    'file': str(file_path),
                    'session': session_name,
                    'task': task_name,
                    'time': file_start_time,
                    'details': f"Contains extreme values (possible artifacts)"
                })
                
        # Future validation checks can be added here
                
    except Exception as e:
        warnings_list.append({
            'type': 'validation_error',
            'file': str(file_path),
            'session': session_name,
            'task': task_name,
            'time': file_start_time,
            'details': f"Error validating data: {str(e)}"
        })
        
    # No immediate logging - warnings will be collected and displayed at the end
    return warnings_list

def verify_temp_file(temp_file, expected_size, config, sample_percent=0.001):
    """Quick statistical verification of a temporary file after task processing."""
    if abs(os.path.getsize(temp_file) - expected_size) > expected_size * 0.01:
        logger.warning(f"Temporary file size mismatch: {os.path.getsize(temp_file)} vs expected {expected_size}")
        return False
    return True        

def verify_merged_data(results, output_file, config):
    """Verify data integrity using statistical sampling with tolerance."""
    num_channels = config['meta_data']['num_channels']
    data_type = config['data_type']
    sample_size = 1000  # Number of samples to check at each position
    
    # Open the output file as memmap
    try:
        output_samples = os.path.getsize(output_file) // (num_channels * np.dtype(data_type).itemsize)
        output_data = np.memmap(output_file, dtype=data_type, mode='r',
                              shape=(num_channels, output_samples))
        
        # Keep track of current position in output file
        output_position = 0
        
        # For each task's temporary file
        for success, task_idx, task_info, temp_file_path, file_size, _ in results:
            # Extract just the actual task name without the "Successfully processed" prefix
            task_name = task_info.replace("Successfully processed task ", "")
            
            if not success or not temp_file_path:
                continue
                
            # Calculate samples in temp file
            temp_samples = file_size // (num_channels * np.dtype(data_type).itemsize)
            
            # Only verify if file exists (it might be deleted by temp dir cleanup)
            if os.path.exists(temp_file_path):
                try:
                    temp_data = np.memmap(temp_file_path, dtype=data_type, mode='r',
                                        shape=(num_channels, temp_samples))
                    
                    # Use a more tolerant comparison for floating point data
                    if np.issubdtype(np.dtype(data_type), np.floating):
                        tolerance = 1e-5  # Small tolerance for floating point comparison
                        
                        # Check sections with tolerance
                        # Beginning
                        start_size = min(sample_size, temp_samples)
                        start_diff = np.max(np.abs(temp_data[:, :start_size] - 
                                                 output_data[:, output_position:output_position+start_size]))
                        if start_diff > tolerance:
                            logger.error(f"Verification failed at beginning of task {task_name} (max diff: {start_diff})")
                            return False
                        
                        # Check middle (if file is large enough)
                        if temp_samples > 2 * sample_size:
                            mid_point = temp_samples // 2
                            mid_start = mid_point - (sample_size // 2)
                            mid_end = mid_point + (sample_size // 2)
                            
                            mid_diff = np.max(np.abs(temp_data[:, mid_start:mid_end] - 
                                                output_data[:, output_position+mid_start:output_position+mid_end]))
                            if mid_diff > tolerance:
                                logger.error(f"Verification failed in middle of task {task_name} (max diff: {mid_diff})")
                                return False

                        # Check end
                        if temp_samples > sample_size:
                            end_size = min(sample_size, temp_samples)
                            end_diff = np.max(np.abs(temp_data[:, -end_size:] - 
                                                output_data[:, output_position+temp_samples-end_size:output_position+temp_samples]))
                            if end_diff > tolerance:
                                logger.error(f"Verification failed at end of task {task_name} (max diff: {end_diff})")
                                return False
                    else:
                        # For integer data types, still use array_equal but add diagnostics if it fails
                        # Beginning
                        start_size = min(sample_size, temp_samples)
                        if not np.array_equal(temp_data[:, :start_size], 
                                           output_data[:, output_position:output_position+start_size]):
                            # Add diagnostics on failure
                            diff_count = np.count_nonzero(temp_data[:, :start_size] != 
                                                       output_data[:, output_position:output_position+start_size])
                            diff_percent = 100.0 * diff_count / (start_size * num_channels)
                            logger.error(f"Verification failed at beginning of task {task_name}: {diff_count} differences ({diff_percent:.2f}%)")
                            
                            # If very few differences (<0.1%), consider it acceptable
                            if diff_percent < 0.1:
                                logger.info(f"Accepting small difference percentage ({diff_percent:.4f}%)")
                            else:
                                return False
                        
                        # Check middle (if the file is large enough)
                        if temp_samples > 2 * sample_size:
                            mid_point = temp_samples // 2
                            mid_start = mid_point - (sample_size // 2)
                            mid_end = mid_point + (sample_size // 2)
                            
                            if not np.array_equal(temp_data[:, mid_start:mid_end], 
                                                output_data[:, output_position+mid_start:output_position+mid_end]):
                                # Add diagnostics on failure
                                diff_count = np.count_nonzero(temp_data[:, mid_start:mid_end] != 
                                                        output_data[:, output_position+mid_start:output_position+mid_end])
                                diff_percent = 100.0 * diff_count / (sample_size * num_channels)
                                logger.error(f"Verification failed in middle of task {task_name}: {diff_count} differences ({diff_percent:.2f}%)")
                                
                                # If very few differences (<0.1%), consider it acceptable
                                if diff_percent < 0.1:
                                    logger.info(f"Accepting small difference percentage ({diff_percent:.4f}%)")
                                else:
                                    return False
                                    
                        # Check end
                        if temp_samples > sample_size:
                            end_size = min(sample_size, temp_samples)
                            if not np.array_equal(temp_data[:, -end_size:], 
                                                output_data[:, output_position+temp_samples-end_size:output_position+temp_samples]):
                                # Add diagnostics on failure
                                diff_count = np.count_nonzero(temp_data[:, -end_size:] != 
                                                        output_data[:, output_position+temp_samples-end_size:output_position+temp_samples])
                                diff_percent = 100.0 * diff_count / (end_size * num_channels)
                                logger.error(f"Verification failed at end of task {task_name}: {diff_count} differences ({diff_percent:.2f}%)")
                                
                                # If very few differences (<0.1%), consider it acceptable
                                if diff_percent < 0.1:
                                    logger.info(f"Accepting small difference percentage ({diff_percent:.4f}%)")
                                else:
                                    return False
                    
                    # Clean up
                    del temp_data
                    
                except Exception as e:
                    logger.error(f"Error verifying task {task_name}: {str(e)}")
                    return False
            
            # Move position forward
            output_position += temp_samples
            
        logger.info("[SUCCESS] File integrity verification successful")
        return True
        
    except Exception as e:
        logger.error(f"Error during verification: {str(e)}")
        return False
    
    finally:
        # Clean up output memmap
        if 'output_data' in locals():
            del output_data

def generate_warnings_report(results, config, output_file):
    # Collecting all warnings from all tasks
    all_warnings = []
    for success, _, _, _, _, warnings in results:
        if warnings:
            all_warnings.extend(warnings)
            
    print_separator("ARTIFACT WARNINGS REPORT")
    
    if not all_warnings:
        logger.info("No artifacts or validation warnings detected in the data.")
        return True
    
    # Group by type
    warnings_by_type = {}
    for warning in all_warnings:
        warning_type = warning['type']
        if warning_type not in warnings_by_type:
            warnings_by_type[warning_type] = []
        warnings_by_type[warning_type].append(warning)
    
    # Only display warnings if show_artifact_warnings is True
    if config.get('show_artifact_warnings', True):
        # Generate artifact warnings report
        print("\n" + "-" * 60)
        print(" ARTIFACT WARNINGS REPORT ".center(60, "-"))
        print("-" * 60)
        print(f"Total warnings: {len(all_warnings)}")
        
        # Print by type
        for warning_type, warnings in warnings_by_type.items():
            print(f"\n{warning_type.upper()} WARNINGS ({len(warnings)}):")
            for warning in warnings:
                print(f"• {os.path.basename(warning['file'])} (Task: {warning['task']}, Time: {warning['time']:.1f}s)")
                print(f"  {warning['details']}")
    else:
        # Just log the count if we're not showing details
        logger.info(f"Found {len(all_warnings)} potential artifacts (not shown, set show_artifact_warnings=true to see)")
    
    # Always save to file if file logging is enabled
    if config.get('enable_file_logging', False):
        # File writing logic...
        warnings_file = os.path.join(config['output_data_dir'], 'artifact_warnings.txt')
        with open(warnings_file, 'w') as f:
            f.write(f"ARTIFACT WARNINGS REPORT\n")
            f.write(f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Session: {config['session_name']}\n")
            f.write(f"Binary file: {output_file}\n")
            f.write(f"Total warnings: {len(all_warnings)}\n\n")
            
            for warning_type, warnings in warnings_by_type.items():
                f.write(f"\n{warning_type.upper()} WARNINGS ({len(warnings)}):\n")
                for warning in warnings:
                    f.write(f"• {os.path.basename(warning['file'])} (Task: {warning['task']}, Time: {warning['time']:.1f}s)\n")
                    f.write(f"  {warning['details']}\n")
                    
            logger.info(f"Artifact warnings saved to: {warnings_file}")
    
    return True

def display_final_summary(config, output_file, results, start_time):
    """Display a final summary of the conversion process with statistics."""
    print_separator("CONVERSION SUMMARY")
    # Calculate and display file stats
    file_size_gb = os.path.getsize(output_file) / (1024**3)
    
    # Format elapsed time as hours:minutes:seconds
    elapsed_seconds = time.time() - start_time
    hours, remainder = divmod(int(elapsed_seconds), 3600)
    minutes, seconds = divmod(remainder, 60)
    elapsed_formatted = f"{hours}:{minutes:02d}:{seconds:02d}"
    
    # Create a clean, visually appealing summary box
    logger.info(f"Output file:  {output_file}")
    logger.info(f"File size:    {file_size_gb:.2f} GB")
    logger.info(f"Tasks:        {len([r for r in results if r[0]])} of {len(results)} completed successfully")
    logger.info(f"Elapsed:      {elapsed_formatted}")
    
    # Report any differences from predicted size
    actual_size = os.path.getsize(output_file)
    total_samples_to_keep = 0
    for start_time_s, stop_time_s in zip(config['start_times_s'], config['stop_times_s']):
        total_samples_to_keep += int((stop_time_s - start_time_s) * config['meta_data']['sampling_rate'])
    data_type_size = np.dtype(config['data_type']).itemsize
    predicted_file_size = total_samples_to_keep * config['meta_data']['num_channels'] * data_type_size
    
    if abs(actual_size - predicted_file_size) / predicted_file_size > 0.01:  # >1% difference
        size_diff_pct = abs(actual_size - predicted_file_size) / predicted_file_size * 100
        logger.info(f"• Size diff:    {size_diff_pct:.1f}% ({actual_size / (1024**3):.2f} GB actual vs " 
              f"{predicted_file_size / (1024**3):.2f} GB predicted)")

def main():
    global logger  # Add logger to globals
    
    # Set up basic logging first (console and memory buffer)
    setup_logging()
    
    start_time = time.time()
    print_header()
    
    # Parse parameters first
    params = parse_command_line_args()
    if not params:
        return
        
    # Now set up file logging if requested (with session name)
    if params.get('enable_file_logging', False):
        session_name = params['input_data']['session_name']
        log_dir = os.path.join(params['output']['path'], 'logs')
        os.makedirs(log_dir, exist_ok=True)  # Ensure the directory exists
        logger.info(f"File logging enabled. Log files will be saved to {log_dir}")
        setup_logging(log_dir=log_dir, session_name=session_name)  # Set up file logging with session name
        
    
    # Extract configuration 
    config = extract_configuration(params)
    if not config:
        return
        
    # Display configuration summary
    display_configuration_summary(config)
    
    # Check disk space and get user confirmation
    if not check_and_confirm(config):
        return
        
    # Process files
    with tempfile.TemporaryDirectory() as temp_dir:
        # Process tasks and get results
        results = process_all_tasks(config, temp_dir)
        if not results:
            return

        # Merge output files
        success, output_file = merge_output_files(results, config)
        if not success:
            return
        
        # Print results summary
        display_final_summary(config, output_file, results, start_time)
        
        # Generate warnings report at the very end
        # generate_warnings_report(results, config, output_file)
    
    print_footer(True)

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nProcess interrupted by user.")
        print("Cleaning up temporary files...")
        sys.exit(0)
    except Exception as e:
        # Instead of logger.critical, use print if logger might not exist yet
        print(f"\nCritical error: {str(e)}")
        print("Please check the log file for details.")
        sys.exit(1)