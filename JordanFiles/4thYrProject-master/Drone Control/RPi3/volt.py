# Voltmeter using ADC
from ADCPi import ADCPi
import time
import numpy as np


# Get voltage

class voltmeter:

    def __init__(self, ADCchannels=[1]):

        self.input = ADCchannels[0]

        # I2C addresses set by jumpers on ADC board
        # Last argument determines ADC bit number (12, 14, 16, 18)
        self.adc = ADCPi(0x68, 0x69, 12)

    def raw(self):
        # Return raw ADC ouput
        go = True
        fail = 0

        while go == True:

            try:
                # Get distance to wall in x direction (meters)
                raw = self.adc.read_voltage(self.input)

                # Handle ADC == 0
                if raw == 0:
                    if fail == 0:
                        raise ValueError
                    else:
                        raise RuntimeError
                else:
                    return raw


            except ValueError:
                fail += 1
                time.sleep(0.001)

            except RuntimeError:
                # Ignore new position sound error and enable landing
                go = False
                return 0


    def volt(self):
        # Return distance in meters or 0 if fail
        go = True
        fail = 0

        while go == True:

            try:
                # Get distance to wall in x direction (meters)
                raw = self.adc.read_voltage(self.input)

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

v = voltmeter()
yes = True

var = bool(input('Choose raw (0) or volts (1)'))

if type(var) != bool:
    var = bool(input('INCORRECT VALUE. Choose raw (0) or volts (1)'))
else:
    if var == 1:
        previous = v.volt()
    else:
        previous = v.raw()

while yes == True:
    
    try:

        if var == 1:
            current = v.volt()
        else:
            current = v.raw()
         
        if current != previous:
            print(current)
        else:
            pass

        time.sleep(0.5)

    except KeyboardInterrupt:
        yes = False
