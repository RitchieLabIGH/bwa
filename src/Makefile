all: pixz install_pixz

PIXZ_TAG=bd56f57
pixz:
	wget -O vasi-pixz.tar.gz "https://github.com/vasi/pixz/tarball/$(PIXZ_TAG)"
	tar xvzf vasi-pixz.tar.gz
	$(MAKE) -C vasi-pixz-$(PIXZ_TAG)

install_pixz:
	cp -f vasi-pixz-$(PIXZ_TAG)/pixz ../resources/usr/bin/

.PHONY: pixz install_pixz
