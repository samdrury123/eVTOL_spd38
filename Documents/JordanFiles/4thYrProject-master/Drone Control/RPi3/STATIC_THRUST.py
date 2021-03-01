from __future__ import print_function
from ADCPi import ADCPi
from hx711 import HX711
import RPi.GPIO as GPIO
import numpy as np
import time
import sys
import os

class loadcell:
    def __init__(self, pins=[24,23], units='mN',calweight=309.5, calfactor=-0.00320306508155):
        self.pins = pins
        self.units = units
        self.calweight = calweight
        self.calfactor = calfactor

        self.lc = HX711(self.pins[0], self.pins[1])

        self.lc.set_reading_format("MSB", "MSB")

        print("Set Zero Weight")
        self.lc.reset()
        self.lc.tare()
        print("Tare done. Add weight now...")

        try:
            self.calin = input("Weight added (grams/ [ ENTER for default ] ): ")
            gravity = 9.81

            print("____")
            print("Input Weight = %f grams" % self.calin)
            # to use both channels, you'll need to tare them both
            #hx.tare_A()
            #hx.tare_B()

            self.calout = self.lc.get_weight(5)
            print("Raw Value = %f" % self.calout)
            # y = Ax + B - Assume B set by tare. Therefore A...

            self.Ax = self.calin*gravity/self.calout

        except SyntaxError: # User hits ENTER, enabling use of previous calibration values
            self.calin = self.calweight
            print("Calibration weight set to default (309.5g)")
            self.Ax = self.calfactor

        print("Calibration factor = %s " % self.Ax)

    def thrust(self, times=1):
        val = self.lc.get_weight(times)
        force = val*self.Ax

        self.lc.power_down()
        self.lc.power_up()

        return force

def cleanAndExit():
    print(" ")
    print(" ")
    print("Cleaning GPIO...")
    GPIO.cleanup()
    print("DONE")
    sys.exit()

# SETUP
thr     = loadcell()

filename = "STATIC_THRUST.txt"
os.system("rm STATIC_THRUST.txt")

currentthr = thr.thrust()

sys.stdout.write('\nForce = %.4e %s    \n' % (currentthr, thr.units))

start = time.time()
now = time.time() - start

while True:
    try:
        now = time.time() - start
        currentthr = thr.thrust()
        
        line = str(now) + "," + str(currentthr) + " \n"
        
        # Save RPMs to file
        with open(filename, "a") as myfile:
            myfile.write(line)


    except KeyboardInterrupt:
        
        elapsed = time.time() - start
        line = str(elapsed) + "," + str(currentthr) + " \n"
        
        # Save RPMs to file
        with open(filename, "a") as myfile:
            myfile.write(line)

        cleanAndExit()
        break