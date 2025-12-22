# ForCommunication

MATLAB project for communication-based control system analysis.

## Structure

- `+cfg/`: Configuration constants
- `basic_src/`: Basic simulation and visualization functions
  - `+datasim/`: Data simulation functions
  - `+visualize/`: Visualization functions
- `com_src/`: Communication-related source code
  - `+attack/`: Attack-related functions
  - `+original_thesis/`: Original thesis implementation
  - `+takakiSenpai/`: Alternative implementation
- `Experiment/`: Experiment scripts and data
- `scripts/`: Demo scripts
- `resources/`: Project resources

## Setup

Run `startup.m` to add all necessary paths to MATLAB.

```matlab
startup
```

## Usage

See `scripts/demo_sdp.m` for a basic example.

## Requirements

- MATLAB R2025a or later
- YALMIP
- SeDuMi or MOSEK solver

