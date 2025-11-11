import os

# Define the file path
file_path = "/workspace/sequitur/data/test.txt"

# Ensure the directory exists
os.makedirs(os.path.dirname(file_path), exist_ok=True)

# Open the file in write mode
with open(file_path, "w") as file:
    # Write some text to the file
    file.write("Hello, world!")