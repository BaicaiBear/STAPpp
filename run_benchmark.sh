#!/bin/bash

# Script to build and run the STAP++ project with different solvers
# and compare performance

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}Building STAP++ project...${NC}"
mkdir -p build
cd build
cmake ../src
make -j4

if [ $? -ne 0 ]; then
    echo -e "${RED}Build failed!${NC}"
    exit 1
fi

echo -e "${GREEN}Build successful!${NC}"

# Function to run the program with a specific solver and input file
run_test() {
    local input_file=$1
    local solver_type=$2
    local solver_name=$3
    
    echo -e "${BLUE}Running with $solver_name solver...${NC}"
    echo "Command: ./stap++ $input_file $solver_type"
    
    # Run the program and capture the output
    output=$(./stap++ $input_file $solver_type 2>&1)
    
    # Extract the timing information
    input_time=$(echo "$output" | grep "TIME FOR INPUT PHASE" | awk -F'=' '{print $2}' | tr -d ' ')
    matrix_time=$(echo "$output" | grep "TIME FOR CALCULATION OF STIFFNESS MATRIX" | awk -F'=' '{print $2}' | tr -d ' ')
    solution_time=$(echo "$output" | grep "TIME FOR FACTORIZATION AND LOAD CASE SOLUTIONS" | awk -F'=' '{print $2}' | tr -d ' ')
    total_time=$(echo "$output" | grep "T O T A L   S O L U T I O N   T I M E" | awk -F'=' '{print $2}' | tr -d ' ')
    
    echo -e "${GREEN}Results for $solver_name:${NC}"
    echo "  Input time: $input_time seconds"
    echo "  Matrix assembly time: $matrix_time seconds"
    echo "  Solution time: $solution_time seconds"
    echo "  Total time: $total_time seconds"
    echo ""
    
    # Return the solution time for comparison
    echo $solution_time
}

# Check if an input file was provided
if [ $# -eq 0 ]; then
    echo -e "${RED}No input file provided!${NC}"
    echo "Usage: $0 <input_file_without_extension>"
    exit 1
fi

input_file=$1

# Run tests with different solvers
skyline_time=$(run_test $input_file 0 "Skyline LDLT (Original)")
eigen_direct_time=$(run_test $input_file 1 "Eigen Direct (SimplicialLDLT)")
eigen_cg_time=$(run_test $input_file 2 "Eigen Iterative (Conjugate Gradient)")

# Compare performance
echo -e "${BLUE}Performance Comparison:${NC}"
echo "Skyline LDLT solution time: $skyline_time seconds"
echo "Eigen Direct solution time: $eigen_direct_time seconds"
echo "Eigen CG solution time: $eigen_cg_time seconds"

# Calculate speedup
skyline_time_val=$(echo $skyline_time | sed 's/[^0-9.]*//g')
eigen_direct_time_val=$(echo $eigen_direct_time | sed 's/[^0-9.]*//g')
eigen_cg_time_val=$(echo $eigen_cg_time | sed 's/[^0-9.]*//g')

if [ $(echo "$skyline_time_val > 0" | bc -l) -eq 1 ]; then
    eigen_direct_speedup=$(echo "scale=2; $skyline_time_val / $eigen_direct_time_val" | bc -l)
    eigen_cg_speedup=$(echo "scale=2; $skyline_time_val / $eigen_cg_time_val" | bc -l)
    
    echo -e "${GREEN}Speedup with Eigen Direct: ${eigen_direct_speedup}x${NC}"
    echo -e "${GREEN}Speedup with Eigen CG: ${eigen_cg_speedup}x${NC}"
fi

echo -e "${BLUE}Benchmark complete!${NC}"
