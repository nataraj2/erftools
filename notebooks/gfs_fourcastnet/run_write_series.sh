#!/bin/bash

# Specify the root directory for traversal
root_directory="./"

# Create a temporary file to store the list of directories
temp_file=$(mktemp)
find "$root_directory" -type d -print0 > "$temp_file"

# Initialize an empty string to store file information
files_list=""

# file counter to be used as time entry
count=0

ls -d */ | sort -z > temporary_file

# Loop through all .vtk files in the current directory
for file in *.vtk; do
    # Skip if no .vtk files are found
    [ -e "$file" ] || continue

    # Extract base name without extension
    base_name=$(basename "$file" .vtk)

    echo "$base_name"

    # Create file information
    files_list+="$(printf "{ \"name\": \"%s\", \"time\": %d }," "$file" "$count")"
    files_list+=$'\n'

    ((count++))
done

# Remove trailing comma from the last entry
files_list="${files_list%,}"

# Create the final JSON structure
# Header line
header_line="{ \"file-series-version\": \"1.0\", \"files\": ["
# Write the files list
all_files="$(printf '%s\n' "$files_list") ] }"

file_series_data="$header_line"
file_series_data+=$'\n'
file_series_data+="$all_files"

# Write the generated JSON structure to a file named plot_files.series
echo "$file_series_data" > plot_files.series

# Remove the temporary file
rm "$temp_file"

echo "JSON structure has been written to plot_files.series"
