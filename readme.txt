DYNAC V6.0 (WINDOWS and LINUX/MAC versions)             11-Aug-2014

Software License Agreement

NOTICE TO USERS: 
This End User License Agreement (the "Agreement") is a legally binding 
agreement between you (either an individual or an entity, the "User"), and the
DYNAC development team (currently Saby Valero and Eugene Tanke) regarding the
DYNAC software (the "Software" or the" Program"), and b) all successor
upgrades, revisions, patches, fixes, modifications, copies, additions or
maintenance releases of the Software, if any, licensed to you by the DYNAC
development team (collectively, the "Updates"), and c) related user
documentation and explanatory materials or files provided in written, "online"
or electronic form (the "Documentation" and together with the Software and
Updates, the "PRODUCT" or the "Distribution Package"). 

CAREFULLY READ THE FOLLOWING LEGAL AGREEMENT. USE OF THE PRODUCT PROVIDED WITH
THIS AGREEMENT CONSTITUTES YOUR ACCEPTANCE OF THESE TERMS. IF YOU DO NOT AGREE
TO THE TERMS OF THIS AGREEMENT, DO NOT INSTALL AND/OR USE THE PRODUCT. YOUR 
USE OF THE PRODUCT IS CONDITIONED UPON COMPLIANCE WITH THE TERMS OF THIS 
AGREEMENT.

1. Intellectual property rights  

The DYNAC Program and associated software ("PRODUCT") is owned by the DYNAC 
development team. Your possession, installation or use of the PRODUCT does
not transfer to you any title to the intellectual property in the PRODUCT, and
you will not acquire any rights in the PRODUCT except as expressly set forth 
in this Agreement. 

The DYNAC Program and associated software is freeware and as such the authors
cannot be held accountable for any problems and or costs arising from its
usage. As freeware, DYNAC and associated software may be used for commercial 
goals, with the exception of commercializing DYNAC PRODUCT, parts of the DYNAC
Program and/or associated software itself.

2. Scope of the License 

You are granted a non-exclusive license to use the PRODUCT as set forth 
herein. You can use the PRODUCT as set forth in the Agreement for 
non-commercial purposes in non-business, non-commercial environment. As 
freeware, DYNAC and associated software may be used for commercial goals, with
the exception of commercializing DYNAC, parts of DYNAC and/or associated 
software itself. 

The Software may be freely distributed, but no person or company may charge a
fee for the distribution of the PRODUCT without written permission from the 
DYNAC development team. 

3. Limited warranties 

The DYNAC development team DOES NOT WARRANT THAT THE SOFTWARE IS FIT FOR ANY 
PARTICULAR PURPOSE. The DYNAC development team DISCLAIMS ALL OTHER WARRANTIES
WITH RESPECT TO THE SOFTWARE, EITHER EXPRESS OR IMPLIED. SOME JURISDICTIONS 
DO NOT ALLOW THE EXCLUSION OF IMPLIED WARRANTIES OR LIMITATIONS ON HOW LONG 
AN IMPLIED WARRANTY MAY LAST, SO THE ABOVE LIMITATIONS OR EXCLUSIONS MAY NOT 
APPLY TO YOU. 

4. Legality statement

The program that is licensed to you is absolutely legal and you can use it 
provided that you are the legal owner of all files or data you are going to
use with our software or have permission from the legitimate owner to perform
these acts. Any illegal use of our software will be solely your 
responsibility. Accordingly, you affirm that you have the legal right to 
access all data, information and files. 

You further attest that the data and/or files will not be used for any illegal
purpose. 
 
5. Final provisions 

All rights not expressly granted here are reserved by the DYNAC development
team. 


GETTING STARTED:
After unzipping you should have the following directories in the
current directory:
bin  datafiles  help  plot  source

In order to create a DYNAC executable, go to the source directory and
type comv6_0 (in the case of WINDOWS, it is assumed here that you have
already installed the gfortran compiler and set the path for it).

In order to create an executable for the plotting post-processor, which
is based on GNU plot, go to the plot directory and type complt
Then go to the datafiles directory.

You are now ready to run DYNAC. An example file is given for 2
different types of accelerators:
sns_mebt_dtl1.in     -> THE MEBT and first DTL tank for the SNS
egun_example2.in     -> Electron gun (you will also need egun_field.txt)

Run it on linux/MAC by typing:
../bin/dynacV6_0 sns_mebt_dtl1.in
or
../bin/dynacV6_0 egun_example2.in

Run it on Windows by typing:
..\bin\dynacv6_0 sns_mebt_dtl1.in
or
..\bin\dynacv6_0 egun_example2.in

Under Windows you can set the Windows PATH variable such as to point
to the bin directory that contains dynacv6_0.exe (i.e. you would then
only need to type dynacv6_0 as opposed to ..\bin\dynacv6_0)


You may view the plots by typing " plotit ". There is a minor problem
with plotit for the WINDOWS version: the first time you use it it may
complain about not finding 3 files. Ignore this message. Details on
"plotit" may be found in chapter 7 of the DYNAC help file.


Input files contain a sequence of keywords or type code entries. The
help file containing the user instructions for these type code entries
is in the help directory.


Please feel free to send any suggestions, comments, modifications or
new routines to Eugene Tanke ( dynac.support@cern.ch ).
