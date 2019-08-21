EXE = a5
DLIBS = -L/Users/anzuhakone/glui-2.36/src/lib 
INC = -I/Users/anzuhakone/glui-2.36/src/include -framework OpenGL -framework GLUT -lglui

$(EXE) : 
	g++ -g -O2 $(DLIBS) *.cpp -o $@ $(INC)

clean:
	rm -rf *.o *~ $(EXE)

