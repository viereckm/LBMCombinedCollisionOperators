#!/bin/bash

. .venv/bin/activate
cd build

./benchmark/obstacle/squareObstacle ./benchmark/obstacle/BreuerGuo/Re10/Nx500Ny80/config.prm
./benchmark/obstacle/squareObstacle ./benchmark/obstacle/BreuerGuo/Re10/Nx1000Ny160/config.prm

./benchmark/obstacle/squareObstacle ./benchmark/obstacle/BreuerGuo/Re30/Nx500Ny80/config.prm
./benchmark/obstacle/squareObstacle ./benchmark/obstacle/BreuerGuo/Re30/Nx1000Ny160/config.prm

./benchmark/obstacle/squareObstacle ./benchmark/obstacle/BreuerGuo/Re60/Nx500Ny80/config.prm
./benchmark/obstacle/squareObstacle ./benchmark/obstacle/BreuerGuo/Re60/Nx1000Ny160/config.prm

./benchmark/obstacle/squareObstacle ./benchmark/obstacle/BreuerGuo/Re100/Nx500Ny80/config.prm
./benchmark/obstacle/squareObstacle ./benchmark/obstacle/BreuerGuo/Re100/Nx1000Ny160/config.prm

./benchmark/obstacle/squareObstacle ./benchmark/obstacle/BreuerGuo/Re133/Nx500Ny80/config.prm
./benchmark/obstacle/squareObstacle ./benchmark/obstacle/BreuerGuo/Re133/Nx1000Ny160/config.prm

./benchmark/obstacle/squareObstacle ./benchmark/obstacle/BreuerGuo/Re150/Nx500Ny80/config.prm
./benchmark/obstacle/squareObstacle ./benchmark/obstacle/BreuerGuo/Re150/Nx1000Ny160/config.prm

python3 fixedGrid.py
