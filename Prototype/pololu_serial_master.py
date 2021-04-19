import serial
import numpy as np

class servo_serial_master:
    def __init__ (self, port = '/dev/ttyACM0', baud_rate = 9600):
        # Dictionary with polulu's servo master commands
        self.commands = {}  
        self.commands['set_target'] = 132
        self.commands['set_speed'] = 135

        # PWM value, in microseconds, corresponding to each servo home position
        self.servo_home = np.array([1480, 1585, 1569, 1500, 1455, 1650])
        
        # (Delta u)/(Delta alpha)
        self.gamma = 2000/np.pi

        # Open the serial port
        self.serial_port = serial.Serial (port, baud_rate)

       

    def set_target_rad (self, target_vec):      
        for servo_n, angle in np.ndenumerate(target_vec):        
            target =  np.uint16(np.around(4*(self.servo_home[servo_n[0]] + np.power(-1,servo_n[0]+1) * self.gamma*angle))) # 2000?! Acho que deve ser 1000
            print(servo_n[0], angle, target)
            array = bytearray([self.commands['set_target'], servo_n[0], target & 0x7F, (target >> 7) & 0x7F]) # 132 = command flag, 0 = servo number, [112,46] = 8 bit target
            # in quarter microseconds. 1 quarter microseconds = 0,25 microseconds
            # 1500 micro = 6000 quater micro
            # 6000 = 0b1011101110000
            # 112 = 0b0b1110000
            # 46 = 0b101110
            self.serial_port.write(array)

    
    def set_speed_rad_s (self, speed_vec):
        factor = ((10*10**(-3))/0.25)*(self.gamma)
        for servo_n, speed in np.ndenumerate(speed_vec):
            speed = np.uint16(np.around(factor*speed))
            array = bytearray([self.commands['set_speed'], servo_n[0], speed & 0x7F, (speed >> 7) & 0x7F])
            self.serial_port.write(array)