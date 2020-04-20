# Microgrid Dispatch Simulator

# Overview
This project provides tools to simulate energy management and various dispatch algorithms in community microgrids with distributed energy resources (DERs). The primary features are:
* A quasi-static simulation of steady-state DER frequency response and active power sharing using tie-line bias control
* A bottom-up model of loads that includes a demand-response model for electricity users to optimize energy and power usage subject to constraints imposed by the microgrid operator
* Receding horizon control loops for energy management, load control, and power dipatch
* Implementations of optimal predictive load control algorithms using mixed-integer quadratic programming with a Gurobi solver

We recommend the paper below for a more comprehensive discussion of the modeling.

# Terms of use
The code is available under the MIT license (see license file). In addition, we request that any publications using this code directly or following from the program structure or algorithms within cite the following paper:
> Lee, Jonathan T., Anderson, Sean, Vergara, Claudio, and Callaway, Duncan. "Non-Intrusive Load Management Under Forecast Uncertainty in Energy Constrained Microgrids." *Electric Power Systems Research*, (Under Review).

# Contributing
The project has been developed primarily to support the research of the contributors, and will continue ad-hoc until broader interest is expressed. Please contact Jonathan Lee through GitHub if you are interested in contributing or have a potential use case that could be supported by extensions to the project.

# Guide
## Setup
The project was developed in MATLAB 2018A, and requires the optimization toolbox. To use, clone the repository into a local folder. Either add this folder to the MATLAB path or use the folder as MATLAB's working directory. Add the `scripts` subfolder to the MATLAB path to run the example scripts.

Dependencies:
* MATLAB 2018A with the Optimization Toolbox
* Gurobi version 9 is required for the predictive load control algorithms (installed, licensed, and added to the MATLAB path). See [Gurobi MATLAB installation instructions](https://www.gurobi.com/documentation/9.0/quickstart_mac/matlab_setting_up_grb_for_.html)
* [CVX](http://cvxr.com/cvx/doc/install.html) can be used as a modelling language with Gurobi as a solver. It is not currently required because the current version uses the Gurobi API directly, but previous implementations of the control algorithms using CVX are included in the repository.

## Program structure
The package `MicrogridDispatchController` consists of the following subpackages
* `DataParsing`: Functions for reading configuration and time series data from the file system, and creating models
* `DispatchControllers`: Optimization functions to compute control actions. These are called by the `MicrogridController` object.
* `Models`: Classes to represent objects within the microgrid. Most of these are implemented as handle classes. Many of the objects (excluding the ones with the `Type` suffix) have a state and are meant to be used as subsystems with state evolving over time. They define an `update(u,t,deltaT)` method that takes a struct of inputs called `u` by convention and evolves the model's state from absolute time `t` by an amount of time `deltaT`. By convention time is specified in seconds.
  * `Microgrid`: The physical microgrid. It consists of buses (currently modeled as connected via a copper plate), which have loads and DERs attached to them.
  * `MicrogridController`: A controller that sets load limits and power injection setpoints.
  * `User`: An end user of electricity. Users are of a certain user type, and can have DERs, loads, and a collection of activities. Users adjust their activities in response to signals from the microgrid to maximize their utility of electricity use.
  * `Activity`: Activities belong to users and are of a certain activity type corresponding to specific loads. Activities evolve as a finite state machine that can be queued, in progress, completed, interrupted, or cancelled. The state transitions through user action, by time passing as activities move to completion, or via interruptions in the power availability. Activities are associated with value from completion and costs from interruption that determine the users' utility. Transitions in activity state correspond to transition in the power state of loads.
  * `Load`: Loads correspond to activity types and have a power demand. They have multiple states relating to whether or not they are 'on', specifically whether they are `Connected`, which is controlled by the user, and whether they are connected to a bus that has a voltage greater than zero. At present, only voltage of 1 or 0 is used at the bus - this can be extened to a model with a continous voltage and loads that are dependent on this.
  * `Bus`: A bus serves to model the physical association of loads to the microgrid. The bus has a voltage state `V` that is controlled by the microgrid, and can return the downstream connected load (power demand given the current load state) as a dependent property. Buses also have DERs attached, which includes stored energy as a state.
  * `MeterRelay`: This models a meter with a relay that associate users with the microgrid buses. A users loads are connected to its `LoadBus` and the `SourceBus` is one of the microgrid buses. The relay can be controlled directly, and the object also accepts power limits (where the relay will be opened if the limit is exceeded) and energy limits (where the relay will open if the energy limit is exceeded over the specified time). When the relay is opened, the load bus is set to 0, and the loads (and hence the activities) respond accordingly to the loss of power.
  * `UserType`, `ActivityType`, and `LoadType` describe a typology of these models as their name implies. These values are defined in the `data` folder. The user type defines a probabilistic schedule for which activity types users will engage in and when. The activity types are one-to-one with loads and defined for each user type. These specify the value and costs associated with activities, the power demand of the loads, and the range of times activities take to complete.
  * `SimulationOutputs`: Time series outputs from the simulation in a structured format
  * `Simulation`: Contains two function `simRHC` and `simLoad` to simulate the microgrid with DERs and receding horizon control and users' load only, respectively.
* `Utilities`: Helper functions that are used in multiple places
* `Visualization`: Plotting tools to visualize simulation outputs

The folder `data` contains solar time series data, user type data, and DER parameters. See the `DataParsing` subpackage to see how these data are parsed.
## Example usage
See the `scripts` folder for example usage scripts.
