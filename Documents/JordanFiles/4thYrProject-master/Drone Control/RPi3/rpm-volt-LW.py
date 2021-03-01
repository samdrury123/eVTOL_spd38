from __future__ import print_function
from ADCPi import ADCPi
import time
import numpy as np


##
##
## SET THRESHOLD
##
##

THRESHOLD = 0.25 # Volts

class voltmeter:

    def __init__(self, ADCchannels=[1]):

        self.input1 = ADCchannels[0]

        # I2C addresses set by jumpers on ADC board
        # Last argument determines ADC bit number (12, 14, 16, 18)
        self.adc1 = ADCPi(0x68, 0x69, 12)

    def raw(self):
        raw1 = self.adc1.read_voltage(self.input1)
        return raw1



    def volt(self):
        raw1 = self.adc1.read_voltage(self.input1)
        return raw1/0.976


# Setup sensor and set q = 0
adc1 = ADCPi(0x68, 0x69, 12)
volt = adc1.read_voltage(1)/0.976

print(volt)

go = True

# Count number of readings taken
n = 0

# Count number of revolutions
r1 = 0

# Average number of revs (in 10)
av1 = 0

# If the time between adjacent readings increases by a significant amount
# take it as one revolution
# Sensor at 0 for section not on tape

# Array to store rpms to take average from
rpmav1 = np.array([0])

go = True

# Use to determine number of readings per second
start = time.time()

# Set low counters to zero
zero1 = 0


# Determine gap between revolutions
rev1 = time.time()

print("Begin scanning \n ")



# Run until KeyboardInterrupt
while go==True:

    try:
        x = adc1.read_voltage(1)/0.976
        line = str(x) + " \n"
        # Save RPMs to file
        with open("rpm.txt", "a") as myfile:
            myfile.write(line)

        # # Get sensor reading from GPIO pin (high or low)
        # x = adc1.read_voltage(1)/0.976
        
        # # time.sleep(2)
        # # print(x)
        # # How many readings are low, more than 70

        # #
        # # MOTOR1
        # #
        # if x < THRESHOLD:
        #     # Count readings in low
        #     zero1 += 1

        # else:
        #     if zero1 > 1:
        #         # Determine how long revolution took
        #         revtime1 = time.time() - rev1
        #         rev1 = time.time()

        #         # Check for 10th reading to get average
        #         if r1 == 0:
        #             r1 = 0
        #             av1 = int(np.mean(rpmav1))
        #             line = str(av1) + " \n"
        #             # Save RPMs to file
        #             with open("rpm.txt", "a") as myfile:
        #                 myfile.write(line)
        #         else:
        #             r1 += 1

        #         # Calculate RPM
        #         rpm1 = 60*(1/revtime1)
        #         rpmav1[r1] = rpm1
        #     else:
        #         pass

        #     # Reset zero counter
        #     zero1 = 0

        # # Count readings
        # n += 1

    except (KeyboardInterrupt, SystemExit):
        stop = time.time()

        print(" ")
        print(n)
        go = False
        
# Determine elapsed time
elapsed = stop - start

# Determine average number of readings per second
avread = n/elapsed

print(avread, ' reads per second')
