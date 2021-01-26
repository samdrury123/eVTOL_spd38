from __future__ import print_function
from ADCPi import ADCPi
import time
import numpy as np


##
##
## SET THRESHOLD
##
##

THRESHOLD = 1 # Volts

class voltmeter:

    def __init__(self, ADCchannels=[1,2,3,4]):

        self.input1 = ADCchannels[0]
        self.input2 = ADCchannels[1]
        self.input3 = ADCchannels[2]
        self.input4 = ADCchannels[3]

        # I2C addresses set by jumpers on ADC board
        # Last argument determines ADC bit number (12, 14, 16, 18)
        self.adc1 = ADCPi(0x68, 0x69, 12)
        self.adc2 = ADCPi(0x68, 0x69, 12)
        self.adc3 = ADCPi(0x68, 0x69, 12)
        self.adc4 = ADCPi(0x68, 0x69, 12)

    def raw(self):
        raw1 = self.adc1.read_voltage(self.input1)
        raw2 = self.adc2.read_voltage(self.input2)
        raw3 = self.adc3.read_voltage(self.input3)
        raw4 = self.adc4.read_voltage(self.input4)
        return [raw1, raw2, raw3, raw4]



    def volt(self):
        raw1 = self.adc1.read_voltage(self.input1)
        raw2 = self.adc2.read_voltage(self.input2)
        raw3 = self.adc3.read_voltage(self.input3)
        raw4 = self.adc4.read_voltage(self.input4)
        return [raw1/0.976, raw2/0.976, raw3/0.976, raw4/0.976]


# Setup sensor and set q = 0
v = voltmeter()

test = v.volt()
print(test)

go = True

# Comparitor to decide if reading has changed ie if a == q
q = [0, 0, 0, 0]

# Count number of readings taken
n = 0

# Count number of revolutions
r1 = 0
r2 = 0
r3 = 0
r4 = 0

# Average number of revs (in 10)
av1 = 0
av2 = 0
av3 = 0
av4 = 0

# If the time between adjacent readings increases by a significant amount
# take it as one revolution
# Sensor at 0 for section not on tape

# Array to store rpms to take average from
rpmav1 = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
rpmav2 = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
rpmav3 = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
rpmav4 = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

go = True

# Use to determine number of readings per second
start = time.time()

# Set low counters to zero
zero1 = 0
zero2 = 0
zero3 = 0
zero4 = 0

# Determine gap between revolutions
rev1 = time.time()
rev2 = time.time()
rev3 = time.time()
rev4 = time.time()

print("Begin scanning \n ")
print("TIME (s) |  M1  |  M2   |    M3    |    M4  ")
print("")

# Run until KeyboardInterrupt
while go==True:

    try:

        print("  ",int(time.time()-start)," | ", av1," | ",av2," | ",av3," | ",av4, end="\r")
        # Get sensor reading from GPIO pin (high or low)
        x = v.volt()
        
        # time.sleep(2)
        # print(x)
        # How many readings are low, more than 70

        #
        # MOTOR1
        #
        if x[0] < THRESHOLD:
            # Count readings in low
            zero1 += 1

        else:
            if zero1 > 70:
                # Determine how long revolution took
                revtime1 = time.time() - rev1
                rev1 = time.time()

                # Check for 10th reading to get average
                if r1 == 9:
                    r1 = 0
                    av1 = int(np.mean(rpmav1))
                else:
                    r1 += 1

                # Calculate RPM
                rpm1 = 60*(1/revtime1)
                rpmav1[r1] = rpm1
            else:
                pass

            # Reset zero counter
            zero1 = 0

        # How many readings are low, more than 70

        ##
        ## MOTOR 2
        ##
        if x[1] < THRESHOLD:
            # Count readings in low
            zero2 += 1

        else:
            if zero2 > 70:
                # Determine how long revolution took
                revtime2 = time.time() - rev2
                rev2 = time.time()

                # Check for 10th reading to get average
                if r2 == 9:
                    r2 = 0
                    av2 = int(np.mean(rpmav2))
                else:
                    r2 += 1

                # Calculate RPM
                rpm2 = 60*(1/revtime2)
                rpmav2[r2] = rpm2
            else:
                pass

            # Reset zero counter
            zero2 = 0

        # How many readings are low, more than 70

        ##
        ## MOTOR 3
        ##
        if x[2] < THRESHOLD:
            # Count readings in low
            zero3 += 1

        else:
            if zero3 > 70:
                # Determine how long revolution took
                revtime3 = time.time() - rev3
                rev3 = time.time()

                # Check for 10th reading to get average
                if r3 == 9:
                    r3 = 0
                    av3 = int(np.mean(rpmav3))
                    
                else:
                    r3 += 1

                # Calculate RPM
                rpm3 = 60*(1/revtime3)
                rpmav3[r3] = rpm3
            else:
                pass

            # Reset zero counter
            zero3 = 0

        # How many readings are low, more than 70

        ##
        ## MOTOR 4
        ##
        if x[3] < THRESHOLD:
            # Count readings in low
            zero4 += 1

        else:
            if zero4 > 70:
                # Determine how long revolution took
                revtime4 = time.time() - rev4
                rev4 = time.time()

                # Check for 10th reading to get average
                if r4 == 9:
                    r4 = 0
                    av4 = int(np.mean(rpmav4))
                    
                else:
                    r4 += 1

                # Calculate RPM
                rpm4 = 60*(1/revtime4)
                rpmav4[r4] = rpm4
            else:
                pass

            # Reset zero counter
            zero4 = 0

        # # # Count readings
        n += 1

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
