# Data Collection
[Term] - Terminal on UNIX;
[Rpi] - Raspberry Pi command line;
[QGC] - QGroundControl;
[RC] - Radio Controller;

Machine running [Term] must be on the same network as [RPi]. Connection made between Pixhawk and [QGC] via telemetry (connect via USB).

## COLLECTION OF DATA
 - [RPi] Open Putty, type raspberrypi.local into host name and open
		username = pi
		password = raspberry
 - [QGC] Run Widget > Analyse:
       select  `SYS_STATUS.battery_current`
               `SYS_STATUS.battery_voltage`
 - [RC] Check MANUAL KILL SWITH ENGAGED (G)
 - [RPi] Run STATIC_TEST.py with fan stationary, wait for calibration msg.
     ```
     cd drone
     python STATIC_TEST.py
     ```
 - [QGC] Begin fan by typing from the Mavlink console inserting PWM value between 1000 and 2000.
     ```
     pwm test -p [1000 - 2000] -c 3
     ```
 - [RC] Start fan by disengaging MANUAL KILL SWITCH
 - [QGC] Click `START LOGGING` on QGC and save as 'power' in ./Logs
 - [RPi] Run as default on STATIC_TEST.py by typing `[ENTER]`
 - [RC] Run for desired period of time, engage kill switch to stop fan
 
 ## STOP RECORDING
 - [QGC] Click `STOP LOGGING`. Click `Yes` and `OK`.
 - [RPi] `Ctrl-C` on STATIC_TEST.py
 
 ## GET DATA TO PC
- [Term] Transfer output using scp.
     ```
     scp pi@raspberrypi.local;~/drone/STATIC_TEST.txt ./STATIC_TEST.txt
     ```
 ## PROCESS DATA
 - [Matlab] Run EXP_COLLECT.m (Will run EXP_UPDATE.m and EXP_PLOT.m as
 well)
