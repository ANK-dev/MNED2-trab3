#!/bin/bash
gcc main.c -o bin/main.out -Wall -Wextra -Werror -pedantic -ansi -O2 && \
./bin/main.out && \
python3 plot_graph.py