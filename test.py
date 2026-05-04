
import time
import minimalmodbus


PORT = "/dev/ttyUSB0"
SLAVE_ADDRESS = 45


instrument = minimalmodbus.Instrument(PORT, SLAVE_ADDRESS)
instrument.serial.baudrate = 19200
instrument.serial.bytesize = 8
instrument.serial.parity = minimalmodbus.serial.PARITY_NONE
instrument.serial.stopbits = 1
instrument.serial.timeout = 3
instrument.mode = minimalmodbus.MODE_RTU
instrument.clear_buffers_before_each_transaction = True


while True:
    try:
        value = instrument.read_float(
            registeraddress=40000,
            functioncode=4,
            number_of_registers=2,
            byteorder=3,
        )

        print("Methane:", value)

    except Exception as exc:
        print("Read error:", exc)

    time.sleep(1)
