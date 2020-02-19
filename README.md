# HCP Suite

## Description

This is a set of scripts to perform common analyses on the Human Connectome Project's neuroimaging data. WIP.

## Usage

### Task analysis

```
Usage: task.sh <action> [action-specific arguments] [common options]
  Actions:
    prepare
      -s <arg>: specify file with subject IDs or a space-separated list of subjects
      -t <arg>: specify a task
      -p <arg>: specify a parcellation
      You might want to consider using -e to prevent exiting when one subject fails.
    design
      -s <arg>: specify file with subject IDs or a space-separated list of subjects
      -o <arg>: specify an output file (e.g. design.mat)
      -V <arg>: specify variables (e.g. "Age_in_Yrs Gender BMI"; repeatable)
      -f <arg>: specify CSV file containing the specified variables (repeatable; e.g.
                "-f restricted.csv -f unrestricted.csv"\)
    eb
      -s <arg>: specify file with subject IDs or a space-separated list of subjects
      -o <arg>: specify an output file (e.g. design.mat)
      -f <arg>: specify CSV file containing the specified variables (repeatable; e.g.
                "-f restricted.csv -f unrestricted.csv"\)
    analyse
      -s <arg>: specify file with subject IDs or a space-separated list of subjects
      -t <arg>: specify a task
      -c <arg>: specify a COPE number
      -p <arg>: specify a parcellation
      -d <arg>: specify target directory
      [-n <arg>]: number of permutations (default: 5000)
      [-o <arg>]: additional PALM arguments
      Specify design, contrast and EB via additional PALM arguments or place them in
      <target directory>/../ as design_<task>.mat, design.con and eb_<task>.csv

  Common options:
    [-l <log file>]: specify log file instead of default /tmp/${0}-${RANDOM}.log
    [-e]: turn off error handling (do not exit on error)
```

### Task connectivity analysis

Code exists, needs to be properly integrated

### Resting-state connectome analysis

Code exists, needs to be properly integrated

### Seed-based correlation analysis

Code exists, needs to be properly integrated

### Dual regression

Code exists, needs to be properly integrated

## TODO

- Add, clean, and integrate more code
- Better in-code documentation
