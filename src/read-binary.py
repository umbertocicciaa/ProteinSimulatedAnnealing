import struct

with open('phi_256_to20_k1_alpha1_sd3.ds2', 'rb') as file:
    data = file.read()

float_size = 4 
num_floats = len(data) // float_size

float_values = struct.unpack(f'{num_floats}f', data)

print(float_values)
