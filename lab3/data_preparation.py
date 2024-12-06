import struct
import numpy as np
import intvalpy as ip


def read_bin_file_with_numpy(file_path):
    with open(file_path, 'rb') as f:
        header_data = f.read(256)
        side, mode_, frame_count = struct.unpack('<BBH', header_data[:4])

        frames = []
        point_dtype = np.dtype('<8H')

        for _ in range(frame_count):
            frame_header_data = f.read(16)
            stop_point, timestamp = struct.unpack('<HL', frame_header_data[:6])
            frame_data = np.frombuffer(f.read(1024 * 16), dtype=point_dtype)
            frames.append(frame_data)
        print("Complete load data")
        return np.array(frames)


def scalar_to_interval(x, rad):
    return ip.Interval(x - rad, x + rad)


scalar_to_interval_vec = np.vectorize(scalar_to_interval)


def get_avg(data):
    avg = [[0]*8]*1024
    for i in range(len(data)): # 100
        avg = np.add(avg, data[i])
    return np.divide(avg, len(data))


def GetData():
    x_data = read_bin_file_with_numpy('-0.205_lvl_side_a_fast_data.bin')
    y_data = read_bin_file_with_numpy('0.225_lvl_side_a_fast_data.bin')

    # x_avg = get_avg(x_data)
    # y_avg = get_avg(y_data)

    x_voltage = x_data / 16384.0 - 0.5
    y_voltage = y_data / 16384.0 - 0.5

    # x_avg = x_avg / 16384.0 - 0.5
    # y_avg = y_avg / 16384.0 - 0.5

    rad = 2 ** (-14)

    X = scalar_to_interval_vec(x_voltage, rad).flatten()
    Y = scalar_to_interval_vec(y_voltage, rad).flatten()

    # X_mode = scalar_to_interval_vec(x_avg, rad).flatten()
    # Y_mode = scalar_to_interval_vec(y_avg, rad).flatten()
    print("Convert X and Y")
    return X, Y

