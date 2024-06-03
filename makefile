APPS    = ./apps
BIN     = ./bin
INCLUDE = ./include
OBJ     = ./obj
LIBS     = ./libs

#If compiled with gfortran
include $(INCLUDE)/gfortran.def
 
all:
	g++ -std=c++17 -c $(LIBS)/library.cpp -I $(INCLUDE) -o $(OBJ)/library.o
	g++ -std=c++17 -c $(LIBS)/Plot2D.cpp -I $(INCLUDE) -o $(OBJ)/Plot2D.o
	g++ -std=c++17 -c $(LIBS)/MakeHistogram.cpp -I $(INCLUDE) -o $(OBJ)/MakeHistogram.o
	g++ -std=c++17 -c $(LIBS)/MakeScatter.cpp -I $(INCLUDE) -o $(OBJ)/MakeScatter.o
	g++ -std=c++17 -c $(LIBS)/secbeam.cpp -I $(INCLUDE) -o $(OBJ)/secbeam.o
	g++ -std=c++17 -c $(LIBS)/bipa.cpp -I $(INCLUDE) -o $(OBJ)/bipa.o
	g++ $(APPS)/SecRIBRAS.cpp $(OBJ)/*.o -I $(INCLUDE) -o $(BIN)/SecRIBRAS
	g++ $(APPS)/main.cpp $(OBJ)/*.o -I $(INCLUDE) -o $(BIN)/Apps/main
	gfortran $(APPS)/cinema.for -o $(BIN)/Apps/cinema
	gfortran $(APPS)/twsp.f90 -o $(BIN)/Apps/twsp
	sudo chmod +x $(BIN)/Apps/stopx2
	sudo chmod +x $(BIN)/Apps/stopx
	sudo chmod +x $(BIN)/Apps/kineq

clean:
	rm $(OBJ)/*.o
