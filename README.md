# Provably efficient machine learning for quantum many-body problems

This is an open source implementation for the paper Provably efficient machine learning for quantum many-body problems.

We require `g++` and `python` version 3.

### Quick Start

The following is for using a command line interface to run the C++ implementation. The most important code in this repo is the creation of measurement schemes, see `data_acquisition_shadow.cpp` or `data_acquisition_shadow.py`.

```shell
# Compile the codes
> g++ -std=c++0x -O3 data_acquisition_shadow.cpp -o data_acquisition_shadow
> g++ -std=c++0x -O3 prediction_shadow.cpp -o prediction_shadow

# Generate observables you want to predict
> g++ -O3 -std=c++0x generate_observables.cpp -o generate_observables
> ./generate_observables

# Create measurement scheme (stored in scheme.txt) using derandomized version of classical shadows
> ./data_acquisition_shadow -d 100 generated_observables.txt 1> scheme.txt

# Do the physical experiments
# Store the data in measurement.txt

# Predicting many few-body observables
> ./prediction_shadow -o measurement.txt observables.txt
# Predicting many subsystem entanglement entropy
> ./prediction_shadow -e measurement.txt subsystems.txt
```

Since many people are using Python, we have implemented `data_acquisition_shadow.py` which is the Python version of `data_acquisition_shadow.cpp` and `prediction_shadow.py` which is the Python version of `prediction_shadow.cpp`. The purpose of the two codes is only to facilitate understanding of the procedure and it could be orders of magnitude slower than the C++ implementation. It can be used through the command line interface
```shell
> python data_acquisition_shadow.py -d 10 observables.txt
> python prediction_shadow.py -o measurement.txt observables.txt
```
or by importing into a Python code
```python
import data_acquisition_shadow
import prediction_shadow

# randomized classical shadow consisting of 100 parallel measurements in a 20-qubit system
measurement_procedure = data_acquisition_shadow.randomized_classical_shadow(100, 20)

# measurement_procedure = [a list of 100 parallel measurements, each being [a list of 20 single-qubit Pauli bases]]
print(measurement_procedure)

# full_measurement = [a list of [list of (XYZ basis, +-1 outcome) for each qubit]]
# one_observable = [a list of (Pauli-XYZ, index for the qubit)]
estimate_exp(full_measurement, one_observable)
```
Currently, `prediction_shadow.py` only support the `-o` option for prediction expectation value of observables.

### Step 1: Compile the code
In your terminal, perform the following to compile the C++ codes to executable files:
```shell
> g++ -std=c++0x -O3 data_acquisition_shadow.cpp -o data_acquisition_shadow
> g++ -std=c++0x -O3 prediction_shadow.cpp -o prediction_shadow
```

### Step 2: Prepare the measurements
The executable `data_acquisition_shadow` could be used to produce an efficient measurement scheme for predicting many local properties from very few measurements. There are two ways to use this program:

#### 1. Randomized measurements:
```shell
> ./data_acquisition_shadow -r [number of measurements] [system size]
```
This generates random Pauli measurements. There would be `[number of measurements]` repetitions on a system with `[system size]` qubits.
You may then use this set of randomized measurements to perform the experiment.
