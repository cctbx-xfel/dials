#!/bin/sh
build/bin/libtbx.python modules/cctbx_project/libtbx/auto_build/create_installer.py \
	--dest            tmp/dials-installer-dev \
	--install_script  modules/dials/installer/dials_installer.py \
	--version         dev \
	--binary
