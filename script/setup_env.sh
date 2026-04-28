#!/usr/bin/env bash

if [ ! -d "./.venv/" ]; then 
    python3 -m venv ./.venv

    source ./.venv/bin/activate

    pip install --upgrade pip
    pip --version

    pip install requests

    pip install -U pandas numpy seaborn matplotlib uv

else

    source ./.venv/bin/activate

fi