SHELL = /bin/sh
CC = g++
LIBS = -lm  -lboost_random -lPocoFoundation -lpthread -std=c++0x -L/usr/local/qwt-6.0.2/lib -lqwt -lQtGui -lQtCore -lQt3Support -lQtOpenGL
CPPFLAGS = -g -O2 -Wall

VPATH=%.h ./include
VPATH=%.o ./obj

OBJDIR = ./obj
INCLUDE_DIR = ./include /usr/include/Qt /usr/include/Qt3Support
INCLUDES  := $(addprefix -I,$(INCLUDE_DIR))

objects = $(addprefix $(OBJDIR)/, main.o individual.o costfunc.o population.o random.o worker.o window.o)

MY_APPS = test

$(MY_APPS) : $(objects)
	${CC} -o app1 ${objects} ${CPPFLAGS} ${LIBS}

$(OBJDIR)/%.o: %.cpp
	$(CC) -c $(CPPFLAGS) ${LIBS} ${INCLUDES} $< -o $@


$(OBJDIR):
	mkdir $(OBJDIR)

.PHONY : clean
clean:
	rm -f ${MY_APPS}
	rm -f ${objects}
