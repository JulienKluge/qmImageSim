clang -flto -O3 -march=native -fopenmp -ffast-math main.cpp -I "SFML-2.4.2/include/" -L "SFML-2.4.2/lib/" -o Program.exe -lsfml-system -lsfml-window -lsfml-graphics
