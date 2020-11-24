@echo off
gcc.exe main.c -o bin\main.exe -Wall -Wextra -Werror -pedantic -ansi -O2 && ^
.\bin\main.exe && ^
python3 plot_graph.py