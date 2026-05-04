import time
import minimalmodbus

instrument = minimalmodbus.Instrument('/dev/ttyUSB0', slaveaddress=49, debug=False)
instrument.serial.baudrate = 921600
instrument.serial.baudrate = 115200

instrument.serial.parity = minimalmodbus.serial.PARITY_NONE
instrument.serial.timeout = 2

communication = False
while communication ==  False:
    time.sleep(1)
    try:
        instrument.read_float(43030, 4, 2, 3)
    except:
        print("No communication with sensor, trying again...")
    else:
        communication = True
        print("Communication with sensor established.")

time.sleep(2)

while True:
    val = instrument.read_float(43030, 4, 2, 3)
    mystring = ' '.join([str(x) for x in val])
    myfloat = float(mystring)
    rndfloat = round(myfloat, 2)
    print(rndfloat)
