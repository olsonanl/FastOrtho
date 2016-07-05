SRC=FastOrtho.cpp \
  Executer.cpp \
  OptionPass.cpp \
  OrthoTop.cpp \
  Classify/Normal.cpp \
  Classify/Ortho.cpp \
  FileIO/LinePull.cpp \
  LinkParser/LinkParser.cpp \
  LinkParser/PickAll.cpp \
  LinkParser/PickE.cpp \
  LinkParser/PickPiE.cpp \
  List/BlockMarch.cpp \
  List/Blocks.cpp \
  List/Links.cpp \
  List/LinkFinder.cpp \
  List/StringStore.cpp \
  Parser/BpoAllParser.cpp \
  Parser/CharSplitParser.cpp \
  Parser/GgParser.cpp \
  Parser/LineParser.cpp
    
OBJ=$(SRC:.cpp=.o)
EXE=FastOrtho

CC=g++
CFLAGS=-Wall -O3
LDFLAGS=-enable-auto-import
RM=rm
	
%.o: %.cpp 
	$(CC) $(CFLAGS) -o $@ -c $<

.PHONY : all
all: $(EXE)

$(EXE): $(OBJ)
	$(CC) $(OBJ) $(LDFLAGS) -o $@

.PHONY : clean
clean:
	-$(RM) $(OBJ) core
