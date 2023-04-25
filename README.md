
# Thompson Problem in Julia

This repository contains the implementation of the Thompson Problem optimization algorithm in Julia, using various optimization techniques such as gradient descent, BFGS, and backtracking line search with GD. Each of these techniques is implemented in a separate notebook file.

## Problem Statement

The Thompson Problem is an optimization problem where we aim to find the minimum enclosing circle of a set of points in a plane. The problem is stated as follows:

Given a set of n points in a plane, find the smallest circle that encloses all the points.

## Implementation

The implementation of the Thompson Problem in Julia is done using various optimization techniques. The implementation includes the following:

- Gradient Descent (GD)
- Broyden-Fletcher-Goldfarb-Shanno (BFGS)
- Backtracking Line Search with GD

Each of these techniques is implemented in a separate notebook file. The code is well-documented, and examples are provided to demonstrate how to use the functions.

## Web-Scraping

Additionally, a Python script is attached to this repository, which can be used to web-scrape known solutions of the Thompson Problem from the Wikipedia page. This script can be used to compare the solutions obtained from the implementation with the known solutions. The obtained CSV with the data is also attached.

## Usage

To use the implementation, clone this repository and run the notebook files in a Julia environment. The Python script can be run independently and requires the `beautifulsoup4` and `requests` libraries to be installed.
