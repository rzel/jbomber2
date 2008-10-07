SRCDIR   :=./src/
DEPSDIR  :=./.deps
CLASSDIR :=./build/classes
JARDIR   :=./build/jar
BUILDDIR :=./build

JAVAFILES  :=$(shell ls $(SRCDIR)/*.java)
CLASSFILES :=$(JAVAFILES:$(SRCDIR)/%.java=$(CLASSDIR)/%.class)
MAKEFILES  :=Makefile Makefile.inc Makefile.config $(shell cat Makefile.inc 2> /dev/null | sed "s:-include ::")

ifeq ($(CLASSPATH),)
	CLASSPATH:=.:..
endif

CLASSPATH:=$(CLASSPATH):/usr/share/bcel/lib/bcel.jar:$(shell pwd)/lib/jogl/linux-i586/lib/gluegen-rt.jar:$(shell pwd)/lib/jna.jar:$(shell pwd)/lib/jogl/linux-i586/lib/jogl.jar:${CLASSDIR}
SRCPATH  :=$(SRCDIR)

LD_HACK :=LD_LIBRARY_PATH=../../lib/:../../lib/jogl/linux-i586/lib/
JAVA_XCB_HACK := $(shell ls LIBXCB_ALLOW_SLOPPY_LOCK > /dev/null 2>&1 /dev/null && echo "LIBXCB_ALLOW_SLOPPY_LOCK=1")
AOSS:=$(shell which aoss 2> /dev/null)
#JAVA_BIN := $(LD_HACK) $(JAVA_XCB_HACK) $(AOSS) java -cp "${CLASSPATH}" -Dsun.java2d.opengl=True -Dsun.java3d.opengl=True
JAVA_BIN := $(LD_HACK) $(JAVA_XCB_HACK) $(AOSS) java -cp "${CLASSPATH}"


-include Makefile.config
-include Makefile.inc

# compiler choose targets
help:
	@echo '**********************************************'
	@echo Use one of the following to choose a compiler:
	@echo - make javac
	@echo - make jikes
	@echo - make gcj
	@echo Use the following to display this message:
	@echo - make help
	@echo '**********************************************'

javac: clean
	@echo "-include Makefile.javac" > Makefile.inc

jikes: clean
	@echo "-include Makefile.jikes" > Makefile.inc

gcj: clean
	@echo "-include Makefile.gcj" > Makefile.inc

# common targets

run: all
	cd $(CLASSDIR) && $(JAVA_BIN) Smoke


jarrun: $(JARDIR)/smoke.jar
	cd $(JARDIR) && $(JAVA_BIN) -jar smoke.jar

$(JARDIR)/smoke.jar: $(JAVAFILES) $(CLASSFILES)
	@echo "[packing] smoke.jar"
	@echo "Manifest-Version: 1.2" > $(CLASSDIR)/MANIFEST.MF.in
	@echo "Class-Path: / /lib/jna.jar /lib/gluegen-rt.jar /lib/jogl.jar" >> $(CLASSDIR)/MANIFEST.MF.in
	@echo "Main-Class: Smoke" >> $(CLASSDIR)/MANIFEST.MF.in
	@echo "Bundle-NativeCode: /lib/libfftw.so; /lib/librfftw.so; osname=Linux; processor=i586" >> $(CLASSDIR)/MANIFEST.MF.in
	@mkdir -p $(CLASSDIR)/lib
	@cp $(shell pwd)/dist/lib/* $(CLASSDIR)/lib
	@cp $(shell pwd)/lib/lib* $(shell pwd)/lib/*.dll $(CLASSDIR)/lib
	@cp $(shell pwd)/lib/jogl/linux-i586/lib/* $(CLASSDIR)/lib
	@cd $(CLASSDIR) && jar cmf MANIFEST.MF.in smoke.jar *.class lib/*
	@mkdir -p $(JARDIR)
	@mv $(CLASSDIR)/smoke.jar $(JARDIR)/smoke.jar

jar: $(JARDIR)/smoke.jar

jarzip: $(CLASSFILES)
	@echo "[packing] smoke.jar (zip) - (Warning: needs 1.6 JVM)"
	@mkdir -p $(CLASSDIR)/META-INF
	@echo "Manifest-Version: 1.2\r\n" > $(CLASSDIR)/META-INF/MANIFEST.MF
	@echo "Created-By: 1.5.0_12 (Sun Microsystems Inc.)\r\n" >> $(CLASSDIR)/META-INF/MANIFEST.MF
	@echo "Main-Class: swinggui" >> $(CLASSDIR)/META-INF/MANIFEST.MF
	@echo "" >> $(CLASSDIR)/META-INF/MANIFEST.MF
	@cd $(CLASSDIR) && zip -r -9 smoke.zip META-INF smoke
	@mkdir -p $(JARDIR)
	@mv $(CLASSDIR)/smoke.zip $(JARDIR)/smoke.jar

clean:
	rm -f $(CLASSDIR)/*.class
	rm -f $(JARDIR)/smoke.jar
	rm -f $(CLASSDIR)/smoke.zip
	rm -f $(DEPSDIR)/*.d

cleaner: clean
	rm -f $(CLASSDIR)/*.class
	rm -f Makefile.inc
	rm -f $(DEPSDIR)/*

$(CLASSFILES): $(MAKEFILES)
