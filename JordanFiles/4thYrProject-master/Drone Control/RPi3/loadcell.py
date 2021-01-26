from __future__ import print_function
from ADCPi import ADCPi
from hx711 import HX711
import RPi.GPIO as GPIO
import numpy as np
import time
import sys
import os

class loadcell:
    def __init__(self, pins=[24,23], units='mN'):
        self.pins = pins
        self.units = units
        # self.calweight = calweight
        # self.calfactor = calfactor

        self.lc = HX711(self.pins[0], self.pins[1])

        self.lc.set_reading_format("MSB", "MSB")

        print("Set Zero Weight")
        self.lc.reset()
        self.lc.tare()
        print("Tare done. Add weight now...")

    def test(self, times=10):
        val = self.lc.get_weight(times)
        # force = val*self.Ax

        self.lc.power_down()
        self.lc.power_up()

        return val

def cleanAndExit():
    print(" ")
    print(" ")
    print("Cleaning GPIO...")
    GPIO.cleanup()
    print("DONE")
    sys.exit()

# SETUP
thr     = loadcell()

filename = "LC_cal.txt"
os.system("rm LC_cal.txt")

currentthr = thr.test()

sys.stdout.write('\nForce = %.4e %s\n' % (currentthr, thr.units))

while True:
    try:
            currentthr = thr.test()
        
        line = str(currentthr) + " \n"
        
        # Save RPMs to file
        with open(filename, "a") as myfile:
            myfile.write(line)

    except KeyboardInterrupt:
        
        line = str(currentthr) + " \n"
        
        # Save RPMs to file
        with open(filename, "a") as myfile:
            myfile.write(line)

        cleanAndExit()
        break