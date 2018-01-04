g++ -x c++ -flto -march=native -O3 -fopenmp -ffast-math -funroll-loops -funroll-all-loops -fforce-addr -fprefetch-loop-arrays -finline-limit=6000 main.cpp -I "SFML-2.4.2/include/" -L "SFML-2.4.2/lib/" -o Program.exe -lsfml-system -lsfml-window -lsfml-graphics

pause

Program.exe
