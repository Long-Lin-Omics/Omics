import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import pod5


def pod5_plot(data, reads_ids,outpng=True):
    if isinstance(data, str) and os.path.isdir(data):
        with pod5.DatasetReader(data) as dataset:
            for read_id in reads_ids:
                read = next(dataset.reads([read_id]))
                sample_rate = read.run_info.sample_rate
                signal = read.signal
                time = np.arange(len(signal)) / sample_rate
                
                plt.figure(figsize=(15, 4))
                plt.plot(time, signal)
                plt.title(f"Read: {read_id}")
                plt.xlabel("Time (s) = len(signal)/fixed_sample_rate")
                plt.ylabel("Signal")
                if outpng == True:
                    plt.savefig(read_id + '.plot.png')
                    plt.close()
    elif isinstance(data, str) and os.path.isfile(data):
        with pod5.Reader(data) as reader:
            for read_id in reads_ids:
                read = next(dataset.reads([read_id]))
                sample_rate = read.run_info.sample_rate
                signal = read.signal
                time = np.arange(len(signal)) / sample_rate
                plt.figure(figsize=(15, 4))
                plt.plot(time, signal)
                plt.title(f"Read: {read_id}")
                plt.xlabel("Time (s) = len(signal)/fixed_sample_rate")
                plt.ylabel("Signal")
                if outpng == True:
                    plt.savefig(read_id + '.plot.png')
                    plt.close()


if __name__ == '__main__':
    if len(sys.argv) < 2:
        data = "/ddn/gs1/project/nextgen/niehsinbox/gridion/Hu_Lab/091124_E14WT_3Users/091124_E14WT_3Users/091124_E14WT_3Users/20240911_1752_X1_FAZ60790_46151abd/pod5/"
        read_ids = ['59e8b4e4-e30c-4794-97e9-19fe7e7a1e1f']
    else:
        data = sys.argv[1]
        read_ids = sys.argv[2:]
    print(read_ids)
    pod5_plot(data, read_ids)