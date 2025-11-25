#!/bin/bash

## This file just puts a running version of fs into your PATH. You can do this manually if you know what you are doing.

## Detect shell
test=`echo $SHELL | grep -F '/bash' | wc -l`
if [ "$test" -eq 1 ] ; then
    shell="bash"
else
    shell="other"
fi


## Identify which fs you want:
fsversions="fs_linux_glibc2.3 fs_mac fs_mac_clang"
for fs in $fsversions; do
    ./$fs -V &> test.log
    test=`cat test.log | grep "fs build date" | wc -l`
    if [ "$test" -eq 1 ];  then
	echo "Using fs version called ${fs}. This is now named 'fs' in this directory."
	cp $fs fs

	## Make a personal bin folder:

	echo "Copying fs into $HOME/bin (which will be created if it did not already exist)"
	mkdir -p $HOME/bin
	cp fs $HOME/bin/

	export PATH="$PATH:$HOME/bin"

	if [ $shell == "other" ] ; then
	    echo "You are not using bash so we have not modified your PATH. You should set 'PATH=$PATH:$HOME/bin' manually for fs to be found."
	    rm -f test.log
	    exit 0
	fi

	## Test whether we already have bin in our path
	test=`grep '$HOME/bin' $HOME/.bashrc  | grep PATH  | grep -v ' §#' | wc -l`
	if [ "$test" -gt 0 ] ;then
	    echo "You appear to have $HOME/bin in your path and therefore fs should work"
	else
	    echo "Attempting to update $HOME/.bashrc to include $HOME/bin in the PATH."
	    echo 'export PATH="$PATH:$HOME/bin"' >> $HOME/.bashrc
	fi
       
	## Test fs in the path is the newest version
	fs -V &> test2.log
	test=`diff test.log test2.log | wc -l`
	if [ "$test" -eq 0 ] ; then
	    echo "fs appears to have been successfully installed into your path as the newest version. Run with 'fs -h' to get help. 
Make sure to look at the examples in this folder!"
	    echo "NOTE: run 'source ~/.bashrc' to update your path in this terminal. It should happen automatically when you open a new terminal."
	else
	    echo "fs was installed into $HOME/bin but was not in your PATH. Set this variable manually, or call fs using its full path. For hpc mode you will have to specify this path manually, using fs -exec <path/to>/fs"
	fi

	## Tidy up
	rm -f test.log test2.log
	exit 0
    fi
done

echo "Unable to find a working version of fs. Check them manually."
rm -f test.log
