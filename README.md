# DYNAC V6.0 (WINDOWS and LINUX/MAC versions)             12-Nov-2015


## GETTING STARTED
After unzipping you should have the following directories in the
current directory:
datafiles  help  plot  source

The compiling is based on CMake. To build the executables, run
```
mkdir build
cd build
cmake ..
make
```

On UNIX systems you may run `make install` to get the binaries in your path. On a Windows system you may set the PATH variable. If you do this then you do not have to specify the full path to the binaries as shown in the commands below.

## EXAMPLES
Go to the datafiles directory.

You are now ready to run DYNAC. An example file is given for 2
different types of accelerators:  
sns_mebt_dtl1.in     -> THE MEBT and first DTL tank for the SNS  
egun_example2.in     -> Electron gun (you will also need egun_field.txt)

Run it on linux/MAC by typing:
```
../build/source/dynac sns/mebt_dtl1.in
```
or
```
../build/source/dynac egun/egun.in
```

Run it on Windows by typing:
```
..\build\source\dynac sns\mebt_dtl1.in  
```
or
```
..\build\source\dynac egun\egun.in
```

You may view the plots by typing `plotit`. There is a minor problem
with plotit for the WINDOWS version: the first time you use it it may
complain about not finding 3 files. Ignore this message. Details on
`plotit` may be found in chapter 7 of the DYNAC help file.


Input files contain a sequence of keywords or type code entries. The
help file containing the user instructions for these type code entries
is in the help directory.

## CONTACT
Please feel free to send any suggestions, comments, modifications or
new routines to Eugene Tanke <dynac.support@cern.ch>.

