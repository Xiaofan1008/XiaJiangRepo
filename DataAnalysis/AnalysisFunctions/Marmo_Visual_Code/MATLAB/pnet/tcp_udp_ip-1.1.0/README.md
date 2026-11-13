TCP/UDP/IP Toolbox 2.0.6  
Example of remote tcp-ip connection:  
 con=pnet('tcpconnect','remote-server.xxx.com',1677);  
 pnet(con,'printf','Hello world!\n');  
 pnet(con,'close');

This toolbox can be used to set up TCP/IP connections or  
send/receive UDP/IP packets in MATLAB. It can transmit data over the Intranet/Internet between MATLAB processes or other applications. It is possible to act as server and/or client and transmit textstrings, arrays of any datatype, files or MATLAB variables.

NOTE 2015-05-08:  
 This toolbox has not rely been updated towards latest versions of MATLAB, I know that the MEX API has some minor changes that may need updates in this pnet tool. Maybe some adjustments i c code and recompile towards latest API (mex.h). I do not have any license of latest matlab versions. Maybe some one can support me with that, and I will upload it as updated package. /Peter Rydesäter peter@rydesater.se

New in 2.0.6 2008-03-06

-Included a precompiled pnet.mexw32 for support to in matlab R2007 (R2008?)  
 No need to compile pnet.c to run at the R2007 platform.

NEW in 2.0.5

More "proper" license specification makes it possible to be useded as (dynamical) linked extention to non-free stand alones of your MATLAB commercial applications.   
The precompiled pnet.dll is possible to run in matlab 5.x for Windows.

NEW in 2.0.4  
\* THE bug when running at Windows platform is removed.  
\* pnet\_remote can call arrays of slaves with just one call.

NEW in 2.0.3  
\* Support for the single datatype.

NEW in 2.0.2  
\* Important Bugfix for UDP support.  
\* Calls to drawnow makes blocking operations breakable also in windows.

NEW in 2.0.1  
\* Bugfix for windows version, Sleep() delay is div. with 1000 to correspond to unix sleep().

NEW in 2.0:  
\* UDP packets.  
\* FAST Datatype support for int\*, uint\*, double, char.  
\* FAST INTERNAL byte swapping.  
\* UDP/IP-demos.  
\* TCP/IP-demos.  
\* WEB-SERVER demo.  
\* New nonblocking remote "eval" of matlab expression on remote host.  
  
Tested to work with MATLAB 5.3/6.x for Linux, Solaris, Windows Compiled mex files for tested platforms included. C-source should compile under any Unix/Windows system. Also basic support for remote evaluation and paralelization of MATLAB is included.

This zip package includes the zip package containing version 1.2.4 with mex files (\*.dll) for backward compatibility.

Toolbox by:

Peter Rydesäter   
Peter@Rydesater.se

With support from Mario Bergeron & Mike Medeiros

Any donations from commercial use are welcome to sponsor this project.

TCP/IP connectivity is also included in The MathWorks Instrument Control Toolbox:

[http://www.mathworks.com/products/instrument/](http://www.mathworks.com/products/instrument/)