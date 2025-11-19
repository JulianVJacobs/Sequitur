#!/bin/bash

if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
    echo "Usage: $0 <script_name> <input_file>"
    exit 1
fi

SCRIPT_NAME=$1
INPUT_FILE=$2

if [ ! -d ".venv" ]; then
    echo "Virtual environment not found. Creating .venv and installing requirements..."
    python3 -m venv .venv
    source .venv/bin/activate
    pip install -r requirements.txt
else
    source .venv/bin/activate
fi

python "$SCRIPT_NAME" "$INPUT_FILE"