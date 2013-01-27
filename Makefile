all:
	$(MAKE) -C bwa
	$(MAKE) -C pixz
	cp bwa/bwa resources/usr/bin/bwa
	cp pixz/pixz resources/usr/bin/pixz
