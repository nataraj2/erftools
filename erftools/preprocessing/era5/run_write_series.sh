#!/bin/bash

# Specify the root directory for traversal
root_directory="./"

# Create a temporary file to store the list of directories
temp_file=$(mktemp)
find "$root_directory" -type d -print0 > "$temp_file"

# Initialize an empty string to store file information
files_list=""

# File counter to be used as time entry
count=0

# Sort matching files and save to a temporary file
ls ERA5*.vtk | sort -V -o temporary_file

files_list=""

# Read from the temporary file
while IFS= read -r file; do
    filename=$(basename "$file")   # Keep extension here

    # Build the file list entry with extension included
    files_list+="$(printf "{ \"name\": \"%s\", \"time\": %d}," "$filename" "$count")"
    files_list+=$'\n'

    ((count++))
done < temporary_file

# Optionally remove the trailing comma and print
files_list=$(echo "$files_list" | sed '$ s/,$//')
echo "$files_list"

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
