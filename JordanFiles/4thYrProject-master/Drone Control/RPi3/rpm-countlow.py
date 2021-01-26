from __future__ import print_function
import RPi.GPIO as GPIO
import time
import numpy as np

class RPM:

    def __init__(self, signalpin):
        self.outpin = signalpin

        GPIO.setmode(GPIO.BCM)
        GPIO.setup(self.outpin, GPIO.IN)

    def get_switch(self):
        return GPIO.input(self.outpin)

def setup_all():
    # Set sensors as RPM class
    motor1 = RPM(26)
    motor2 = RPM(16)
    motor3 = RPM(19)
    motor4 = RPM(20)
    return motor1, motor2, motor3, motor4

def get_switch_all(m1, m2, m3, m4):
    # Get readings from all sensors
    # Motor 1 -> GPIO 4
    # Motor 2 -> GPIO 17
    # Motor 3 -> GPIO 27
    # Motor 4 -> GPIO 22
    ms1 = m1.get_switch()
    ms2 = m2.get_switch()
    ms3 = m3.get_switch()
    ms4 = m4.get_switch()
    return [ms1, ms2, ms3, ms4]

# Setup sensor and set q = 0
m1, m2, m3, m4 = setup_all()
test = get_switch_all(m1, m2, m3, m4)
#print(test)
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


        # Get sensor reading from GPIO pin (high or low)
        x = get_switch_all(m1, m2, m3, m4)
        # time.sleep(2)
        # print(x)
        # How many readings are low, more than 70

        #
        # MOTOR1
        #
        if x[0] == 0:
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
        if x[1] == 0:
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
        if x[2] == 0:
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
                    print("  ",int(time.time()-start),"     ", av1,"   ",av2,"    ",av3,"    ",av4, end="\r")
                    
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
        if x[3] == 0:
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
        # n += 1

    except (KeyboardInterrupt, SystemExit):
        stop = time.time()

        # Clear GPIO allocation
        GPIO.cleanup()

        # print(" ")
        # print(n)
        go = False

# # Determine elapsed time
# elapsed = stop - start
#
# # Determine average number of readings per second
# avread = n/elapsed
#
# print(avread, ' reads per second')
