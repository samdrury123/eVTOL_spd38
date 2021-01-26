# Position hold using sonar sensors
from ADCPi import ADCPi
import time
import numpy as np


# Set desired position
# Determine current position
# Send required local coordinate position

class drone:

    def __init__(self, ADCchannels=[8,7]):

        self.xchan = ADCchannels[0]
        self.ychan = ADCchannels[1]

        # I2C addresses set by jumpers on ADC board
        # Last argument determines ADC bit number (12, 14, 16, 18)
        self.adc = ADCPi(0x68, 0x69, 12)

        # Setup MAVPROXY here too
        #
        #
        #

    def xdist(self):
        # Return distance in meters or 0 if fail
        go = True
        fail = 0

        while go == True:

            try:
                # Get distance to wall in x direction (meters)
                raw = self.adc.read_voltage(self.xchan)

                # Handle ADC == 0
                if raw == 0:
                    if fail == 0:
                        raise ValueError
                    else:
                        raise RuntimeError
                else:
                    return raw/0.976


            except ValueError:
                fail += 1
                time.sleep(0.001)

            except RuntimeError:
                # Ignore new position sound error and enable landing
                go = False
                return 0


    def ydist(self):
        # Return distance in meters or 0 if fail
        go = True
        fail = 0

        while go == True:

            try:
                # Get distance to wall in x direction (meters)
                raw = self.adc.read_voltage(self.ychan)

                # Handle ADC == 0
                if raw == 0:
                    if fail == 0:
                        raise ValueError
                    else:
                        raise RuntimeError
                else:
                    return raw/0.976


            except ValueError:
                fail += 1
                time.sleep(0.001)

            except RuntimeError:
                # Ignore new position sound error and enable landing
                go = False
                return 0


    def set_dist(self, target=[1.0,1.0]):
        # Send command to drone to move to new local position
        err = 0
        try:
            # Get distance, handle if 0 (ie fail)
            xdist = self.xdist()
            ydist = self.ydist()

            # Check for error
            if xdist == 0 and ydist != 0:
                err = 1
                yoffset = ydist - target[1]
                raise ValueError

            elif ydist == 0 and xdist != 0:
                err = 2
                xoffset = xdist - target[0]
                raise ValueError

            elif ydist == 0 and xdist == 0:
                err = 3
                raise ValueError

            else:
                xoffset = xdist - target[0]
                yoffset = ydist - target[1]

            # MAVPROXY send new coordinates in the local frame
            #
            #
            #

            print([xoffset, yoffset])

        except ValueError:
            # Hold position, send error tone/land, do not send location changes
            #
            #
            #

            if err == 1:
                print(['ERROR', yoffset])

            elif err == 2:
                print([xoffset, 'ERROR'])

            elif err == 3:
                print(['ERROR', 'ERROR'])


kd = drone()
yes = True

while yes == True:

    try:
        kd.set_dist()
        time.sleep(2)

    except KeyboardInterrupt:
        yes = False
