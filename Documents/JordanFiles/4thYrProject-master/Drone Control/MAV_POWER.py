import time
from pymavlink import mavutil

mavutil.set_dialect("ardupilotmega")

autopilot = mavutil.mavlink_connection('/dev/ttyAMA0', baud=57600)

msg = None

# wait for autopilot connection
while msg is None:
        msg = autopilot.recv_msg()

print msg

# The values of these heartbeat fields is not really important here
# I just used the same numbers that QGC uses
# It is standard practice for any system communicating via mavlink emit the HEARTBEAT message at 1Hz! Your autopilot may not behave the way you want otherwise!
autopilot.mav.heartbeat_send(
6, # type
8, # autopilot
192, # base_mode
0, # custom_mode
4, # system_status
3  # mavlink_version
)

autopilot.mav.command_long_send(
1, # autopilot system id
1, # autopilot component id
400, # command id, ARM/DISARM
0, # confirmation
1, # arm!
0,0,0,0,0,0 # unused parameters for this command
)

time.sleep(2)

autopilot.mav.command_long_send(
1, # autopilot system id
1, # autopilot component id
mavutil.mavlink.system_status.current_battery, # command id, ARM/DISARM
0, # confirmation
0, # disarm!
0,0,0,0,0,0 # unused parameters for this command
)

[a,b,c,d,e,f,g] = autopilot.mav.system_status.recv_msg()

autopilot.mav.param_request_read(
    1,
    1,
    20,
)