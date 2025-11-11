#!/bin/bash

if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
    echo "Usage: $0 <script_name> [inputs]"
    exit 1
fi

SCRIPT_NAME=$1
INPUTS=$2

if [ ! -d ".venv" ]; then
    echo "Virtual environment not found. Creating .venv and installing requirements..."
    python3 -m venv .venv
    source .venv/bin/activate
    pip install -r requirements.txt
else
    source .venv/bin/activate
fi

if [ -z "$INPUTS" ]; then
    echo "Running $SCRIPT_NAME without inputs..."
    python "$SCRIPT_NAME"
else
    echo "Running $SCRIPT_NAME with inputs: $INPUTS"
    python "$SCRIPT_NAME" "$INPUTS"
fi