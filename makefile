seamcarver: *.cpp
	g++ *.cpp -lsfml-graphics -lsfml-window -lsfml-system -o seamcarver $(shell python-config --libs --cflags)
clean:
	rm seamcarver
